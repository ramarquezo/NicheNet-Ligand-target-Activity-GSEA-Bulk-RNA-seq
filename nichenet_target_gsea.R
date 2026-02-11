# ============================================================
# Bulk NicheNet ligand activity (GitHub-friendly, OS-safe)
# ------------------------------------------------------------
# What this script does:
#  1) Reads two DE tables (Excel): sender + receiver
#  2) Builds receiver DE signature + background
#  3) Loads NicheNet resources (package datasets OR Zenodo files)
#  4) Infers ligand activities (sender -> receiver)
#  5) Adds sender DE evidence (logFC, padj) for top ligands
#  6) Builds ligand->target links for top ligands
#  7) Exports: TSV + plots + QC text + optional GO dotplots per ligand
#
# Notes:
#  - No setwd()/choose.dir(); you pass file paths and outdir.
#  - Uses file.path() everywhere (cross-platform).
#  - Uses explicit dplyr:: prefixes to avoid Bioconductor masking.
#  - Handles nichenetr output schema differences (test_ligand vs ligand).
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(ggplot2)
  library(tibble)
  library(dplyr)
  library(nichenetr)
  library(clusterProfiler)
})

# ---------------------------
# Utilities
# ---------------------------

safe_name <- function(x) gsub("[^A-Za-z0-9_\\-]+", "_", x)

make_outdir <- function(outdir = ".", prefix = "Run", stamp = TRUE) {
  prefix_safe <- safe_name(prefix)
  if (isTRUE(stamp)) {
    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    outdir <- file.path(outdir, paste0(prefix_safe, "_", ts))
  } else {
    outdir <- file.path(outdir, prefix_safe)
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  normalizePath(outdir, winslash = "/", mustWork = FALSE)
}

standardize_de <- function(df, symbol_col="symbol", lfc_col="log2FoldChange", padj_col="padj") {
  df <- as.data.frame(df)
  
  # common alternate column names
  if ("log2Foldchange" %in% colnames(df) && !lfc_col %in% colnames(df)) {
    colnames(df)[colnames(df) == "log2Foldchange"] <- lfc_col
  }
  if ("gene" %in% colnames(df) && !symbol_col %in% colnames(df)) {
    colnames(df)[colnames(df) == "gene"] <- symbol_col
  }
  
  missing <- setdiff(c(symbol_col, lfc_col, padj_col), colnames(df))
  if (length(missing)) stop("DE table missing required columns: ", paste(missing, collapse=", "))
  
  df %>%
    dplyr::mutate(
      !!symbol_col := as.character(.data[[symbol_col]]),
      !!lfc_col := as.numeric(.data[[lfc_col]]),
      !!padj_col := as.numeric(.data[[padj_col]])
    )
}

clean_symbols <- function(x) {
  x <- as.character(x)
  x[x == "" | x == "NA"] <- NA
  x
}

get_DE_genes <- function(df, symbol_col="symbol", lfc_col="log2FoldChange", padj_col="padj",
                         padj_cut = 0.05, lfc_cut = 0.25, direction = c("both","up","down"),
                         to_upper = TRUE) {
  direction <- match.arg(direction)
  
  out <- df %>%
    dplyr::mutate(
      .sym = clean_symbols(.data[[symbol_col]]),
      .lfc = as.numeric(.data[[lfc_col]]),
      .padj = as.numeric(.data[[padj_col]])
    ) %>%
    dplyr::filter(!is.na(.sym), is.finite(.lfc), is.finite(.padj))
  
  out <- if (direction == "both") {
    out %>% dplyr::filter(.padj < padj_cut, abs(.lfc) > lfc_cut)
  } else if (direction == "up") {
    out %>% dplyr::filter(.padj < padj_cut, .lfc > lfc_cut)
  } else {
    out %>% dplyr::filter(.padj < padj_cut, .lfc < -lfc_cut)
  }
  
  genes <- out %>% dplyr::pull(.sym) %>% unique()
  if (to_upper) genes <- toupper(genes)
  genes
}

get_background_genes <- function(df, symbol_col="symbol", lfc_col="log2FoldChange", to_upper = TRUE) {
  genes <- df %>%
    dplyr::mutate(
      .sym = clean_symbols(.data[[symbol_col]]),
      .lfc = as.numeric(.data[[lfc_col]])
    ) %>%
    dplyr::filter(!is.na(.sym), is.finite(.lfc)) %>%
    dplyr::pull(.sym) %>%
    unique()
  if (to_upper) genes <- toupper(genes)
  genes
}

# ---------------------------
# NicheNet resources loader
# ---------------------------

load_nichenet_resources <- function(mode = c("package","zenodo"),
                                    nichenet_dir = "nichenet_resources") {
  mode <- match.arg(mode)
  
  if (mode == "package") {
    message("Loading NicheNet resources from nichenetr package datasets...")
    data("ligand_target_matrix", package = "nichenetr", envir = environment())
    data("lr_network",          package = "nichenetr", envir = environment())
    data("weighted_networks",   package = "nichenetr", envir = environment())
    
    return(list(
      ligand_target_matrix = get("ligand_target_matrix", envir = environment()),
      lr_network = get("lr_network", envir = environment()),
      weighted_networks = get("weighted_networks", envir = environment())
    ))
  }
  
  # Zenodo mode: download files if missing
  dir.create(nichenet_dir, showWarnings = FALSE, recursive = TRUE)
  
  ltm_file <- file.path(nichenet_dir, "ligand_target_matrix.rds")
  lr_file  <- file.path(nichenet_dir, "lr_network.rds")
  wn_file  <- file.path(nichenet_dir, "weighted_networks.rds")
  
  ltm_url <- "https://zenodo.org/records/3260758/files/ligand_target_matrix.rds"
  lr_url  <- "https://zenodo.org/records/3260758/files/lr_network.rds"
  wn_url  <- "https://zenodo.org/records/3260758/files/weighted_networks.rds"
  
  download_if_missing <- function(url, dest) {
    if (file.exists(dest) && file.info(dest)$size > 1000) return(invisible(TRUE))
    message("Downloading: ", basename(dest))
    options(timeout = 600)
    download.file(url, destfile = dest, mode = "wb", quiet = FALSE)
    invisible(TRUE)
  }
  
  download_if_missing(lr_url, lr_file)
  download_if_missing(wn_url, wn_file)
  download_if_missing(ltm_url, ltm_file)
  
  list(
    ligand_target_matrix = readRDS(ltm_file),
    lr_network = readRDS(lr_file),
    weighted_networks = readRDS(wn_file)
  )
}

# ---------------------------
# OrgDb resolver for GO dotplots
# ---------------------------

resolve_orgdb <- function(organism_db = c("org.Mm.eg.db", "org.Hs.eg.db")) {
  organism_db <- match.arg(organism_db)
  
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!requireNamespace(organism_db, quietly = TRUE)) {
    BiocManager::install(organism_db, ask = FALSE, update = FALSE)
  }
  
  if (organism_db == "org.Mm.eg.db") {
    return(org.Mm.eg.db::org.Mm.eg.db)
  } else {
    return(org.Hs.eg.db::org.Hs.eg.db)
  }
}

# ---------------------------
# Main runner
# ---------------------------

run_nichenet_bulk <- function(
    sender_file,
    receiver_file,
    sender_label = "Sender",
    receiver_label = "Receiver",
    
    # output
    prefix = "NicheNet",
    outdir = ".",
    timestamp_dir = TRUE,
    
    # column mapping
    symbol_col = "symbol",
    lfc_col = "log2FoldChange",
    padj_col = "padj",
    to_upper = TRUE,
    
    # receiver signature
    receiver_padj_cut = 0.05,
    receiver_lfc_cut  = 0.25,
    receiver_direction = c("both","up","down"),
    
    # sender ligand prior
    sender_ligand_mode = c("DE_up","DE_both","all_expressed"),
    sender_padj_cut = 0.10,
    sender_lfc_cut  = 0.25,
    
    # NicheNet
    resources_mode = c("package","zenodo"),
    top_n_ligands_plot = 20,
    top_n_ligands_targets = 30,
    top_targets_per_ligand = 200,
    
    # GO dotplots per ligand (optional)
    do_go_dotplots = FALSE,
    organism_db = c("org.Mm.eg.db", "org.Hs.eg.db"),
    go_ont = c("BP","MF","CC"),
    go_padj_cutoff = 0.05,
    min_targets_go = 5,
    go_showCategory = 10,
    
    # exports
    export_plots = TRUE,
    export_tables = TRUE,
    save_rds = TRUE
) {
  receiver_direction <- match.arg(receiver_direction)
  sender_ligand_mode <- match.arg(sender_ligand_mode)
  resources_mode <- match.arg(resources_mode)
  organism_db <- match.arg(organism_db)
  go_ont <- match.arg(go_ont)
  
  outdir_final <- make_outdir(outdir = outdir, prefix = prefix, stamp = timestamp_dir)
  message("Outputs will be written to: ", outdir_final)
  
  # Read + standardize
  sender_df <- readxl::read_excel(sender_file) %>% standardize_de(symbol_col, lfc_col, padj_col)
  receiver_df <- readxl::read_excel(receiver_file) %>% standardize_de(symbol_col, lfc_col, padj_col)
  
  if (to_upper) {
    sender_df[[symbol_col]] <- toupper(sender_df[[symbol_col]])
    receiver_df[[symbol_col]] <- toupper(receiver_df[[symbol_col]])
  }
  
  # Receiver DE + bg
  receiver_DE <- get_DE_genes(
    receiver_df, symbol_col, lfc_col, padj_col,
    padj_cut = receiver_padj_cut, lfc_cut = receiver_lfc_cut,
    direction = receiver_direction, to_upper = to_upper
  )
  receiver_bg <- get_background_genes(receiver_df, symbol_col, lfc_col, to_upper = to_upper)
  message("Receiver DE genes: ", length(receiver_DE))
  message("Receiver background genes: ", length(receiver_bg))
  
  # Load NicheNet resources
  res <- load_nichenet_resources(mode = resources_mode)
  ligand_target_matrix <- res$ligand_target_matrix
  lr_network <- res$lr_network
  
  # normalize case
  lr_network <- lr_network %>%
    dplyr::mutate(
      from = if (to_upper) toupper(.data$from) else .data$from,
      to   = if (to_upper) toupper(.data$to) else .data$to
    )
  rownames(ligand_target_matrix) <- if (to_upper) toupper(rownames(ligand_target_matrix)) else rownames(ligand_target_matrix)
  colnames(ligand_target_matrix) <- if (to_upper) toupper(colnames(ligand_target_matrix)) else colnames(ligand_target_matrix)
  
  all_ligands <- unique(lr_network$from)
  
  # Sender candidates
  sender_bg <- get_background_genes(sender_df, symbol_col, lfc_col, to_upper = to_upper)
  
  sender_de_tbl <- sender_df %>%
    dplyr::transmute(
      ligand = .data[[symbol_col]],
      sender_logFC = .data[[lfc_col]],
      sender_padj  = .data[[padj_col]]
    ) %>%
    dplyr::filter(!is.na(ligand), is.finite(sender_logFC), is.finite(sender_padj)) %>%
    dplyr::group_by(ligand) %>%
    dplyr::slice_max(order_by = abs(sender_logFC), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  sender_DE_up <- sender_de_tbl %>%
    dplyr::filter(sender_padj < sender_padj_cut, sender_logFC > sender_lfc_cut) %>%
    dplyr::pull(ligand) %>% unique()
  
  sender_DE_both <- sender_de_tbl %>%
    dplyr::filter(sender_padj < sender_padj_cut, abs(sender_logFC) > sender_lfc_cut) %>%
    dplyr::pull(ligand) %>% unique()
  
  sender_candidates <- switch(
    sender_ligand_mode,
    "all_expressed" = intersect(sender_bg, all_ligands),
    "DE_up"         = intersect(sender_DE_up, all_ligands),
    "DE_both"       = intersect(sender_DE_both, all_ligands),
    stop("Unknown sender_ligand_mode: ", sender_ligand_mode)
  )
  
  message("Sender ligand mode: ", sender_ligand_mode)
  message("Sender candidate ligands used: ", length(sender_candidates))
  
  # Ligand activity
  message("Running NicheNet ligand activity inference...")
  ligand_activities <- nichenetr::predict_ligand_activities(
    geneset = receiver_DE,
    background_expressed_genes = receiver_bg,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = sender_candidates
  ) %>%
    dplyr::arrange(dplyr::desc(aupr_corrected)) %>%
    dplyr::mutate(rank = dplyr::row_number())
  
  # Version-robust ligand column
  lig_col <- dplyr::case_when(
    "test_ligand" %in% colnames(ligand_activities) ~ "test_ligand",
    "ligand"      %in% colnames(ligand_activities) ~ "ligand",
    "from"        %in% colnames(ligand_activities) ~ "from",
    TRUE ~ NA_character_
  )
  if (is.na(lig_col)) {
    stop("Could not find ligand column in ligand_activities. Columns: ",
         paste(colnames(ligand_activities), collapse = ", "))
  }
  
  ligand_activities2 <- ligand_activities %>%
    dplyr::rename(ligand = !!lig_col) %>%
    dplyr::left_join(sender_de_tbl, by = "ligand")
  
  message("Top ligands by AUPR corrected:")
  print(utils::head(ligand_activities2, 10))
  
  # Top ligands lists
  top_ligands_plot <- ligand_activities2 %>%
    dplyr::slice_head(n = top_n_ligands_plot) %>%
    dplyr::pull(ligand)
  
  top_ligs_targets <- ligand_activities2 %>%
    dplyr::slice_head(n = top_n_ligands_targets) %>%
    dplyr::pull(ligand)
  
  # Ligand -> target links
  get_top_targets <- function(ligand, ltm, n = 200) {
    if (!ligand %in% rownames(ltm)) return(tibble::tibble())
    targets <- sort(ltm[ligand, ], decreasing = TRUE)
    tibble::tibble(
      ligand = ligand,
      target = names(targets)[seq_len(min(n, length(targets)))],
      regulatory_potential = unname(targets[seq_len(min(n, length(targets)))])
    )
  }
  
  ligand_target_links <- dplyr::bind_rows(lapply(
    top_ligs_targets,
    get_top_targets,
    ltm = ligand_target_matrix,
    n = top_targets_per_ligand
  )) %>%
    dplyr::mutate(target_in_receiver_DE = .data$target %in% receiver_DE)
  
  # Plots
  p_activity <- ggplot2::ggplot(ligand_activities2, ggplot2::aes(x = rank, y = aupr_corrected)) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = paste0("NicheNet ligand activity (", sender_label, " → ", receiver_label, ")"),
      x = "Ligand rank", y = "AUPR corrected"
    )
  
  p_top20 <- ligand_activities2 %>%
    dplyr::filter(.data$ligand %in% top_ligands_plot) %>%
    dplyr::slice_head(n = top_n_ligands_plot) %>%
    ggplot2::ggplot(ggplot2::aes(x = aupr_corrected, y = stats::reorder(ligand, aupr_corrected))) +
    ggplot2::geom_col(fill = "grey40", width = 0.75) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::labs(
      title = paste0("Top ", top_n_ligands_plot, " predicted active ligands"),
      y = "Ligand", x = "AUPR corrected"
    )
  
  # sender FC evidence bar plot (top ligands)
  p_sender_fc <- ligand_activities2 %>%
    dplyr::filter(.data$ligand %in% top_ligands_plot) %>%
    dplyr::mutate(
      sender_logFC = as.numeric(sender_logFC),
      sender_padj = as.numeric(sender_padj),
      padj_lbl = paste0("p-adj = ", formatC(sender_padj, format = "f", digits = 6)),
      ligand = factor(ligand, levels = rev(top_ligands_plot))
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x = sender_logFC, y = ligand)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, alpha = 0.5) +
    ggplot2::geom_col(fill = "red", width = 0.75) +
    ggplot2::geom_text(
      ggplot2::aes(label = padj_lbl,
                   x = sender_logFC + 0.03 * diff(range(sender_logFC, na.rm = TRUE))),
      hjust = 0, size = 3
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(plot.margin = ggplot2::margin(5.5, 80, 5.5, 5.5)) +
    ggplot2::labs(x = paste0(sender_label, " log2FC"), y = NULL, title = "Sender evidence (top ligands)")
  
  # Optional composite if patchwork exists
  p_summary <- NULL
  if (requireNamespace("patchwork", quietly = TRUE)) {
    p_summary <- (p_top20 | p_sender_fc) + patchwork::plot_layout(widths = c(1, 1))
  }
  
  # GO dotplots per ligand (optional)

  # GO dotplots per ligand (optional)
  dotplots_list <- list()
  if (isTRUE(do_go_dotplots)) {
    message("Running GO enrichment dotplots per top ligands...")
    
    # Require AnnotationDbi for mapIds
    if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install("AnnotationDbi", ask = FALSE, update = FALSE)
    }
    
    OrgDb <- resolve_orgdb(organism_db)
    go_outdir <- file.path(outdir_final, "nichenet_GO_dotplots")
    dir.create(go_outdir, showWarnings = FALSE, recursive = TRUE)
    
    # Helper: convert to TitleCase for mouse if you uppercased everything
    to_title_case_mouse <- function(x) {
      x <- as.character(x)
      x <- ifelse(
        is.na(x), NA_character_,
        paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
      )
      x
    }
    
    # Pull all targets for top ligands (don’t group here)
    lt_sub <- ligand_target_links %>%
      dplyr::filter(.data$ligand %in% top_ligands_plot)
    
    run_go_for_ligand <- function(ligand_name) {
      targets <- lt_sub %>%
        dplyr::filter(.data$ligand == ligand_name) %>%
        dplyr::pull(.data$target) %>%
        unique()
      
      if (length(targets) == 0) return(NULL)
      
      # Keep only targets that are significant in receiver (as in your original)
      receiver_sig <- receiver_df %>%
        dplyr::mutate(.sym = clean_symbols(.data[[symbol_col]])) %>%
        dplyr::filter(!is.na(.sym), is.finite(.data[[padj_col]])) %>%
        dplyr::mutate(.sym = if (to_upper) toupper(.sym) else .sym) %>%
        dplyr::filter(.data[[padj_col]] < go_padj_cutoff) %>%
        dplyr::pull(.sym) %>%
        unique()
      
      targets_sig <- intersect(targets, receiver_sig)
      if (length(targets_sig) < min_targets_go) {
        message("Skipping ", ligand_name, ": too few significant targets (n = ", length(targets_sig), ")")
        return(NULL)
      }
      
      # If mouse + we uppercased earlier, convert to TitleCase for SYMBOL mapping
      targets_for_map <- targets_sig
      if (organism_db == "org.Mm.eg.db") {
        targets_for_map <- to_title_case_mouse(targets_for_map)
      }
      
      # Map SYMBOL -> ENTREZID (robust)
      entrez <- AnnotationDbi::mapIds(
        x = OrgDb,
        keys = targets_for_map,
        column = "ENTREZID",
        keytype = "SYMBOL",
        multiVals = "first"
      )
      
      entrez <- unique(na.omit(unname(entrez)))
      
      if (length(entrez) < min_targets_go) {
        message("Skipping ", ligand_name, ": too few mapped Entrez IDs (n = ", length(entrez), ")")
        return(NULL)
      }
      
      ego <- clusterProfiler::enrichGO(
        gene = entrez,
        OrgDb = OrgDb,
        keyType = "ENTREZID",
        ont = go_ont,
        readable = TRUE
      )
      
      if (is.null(ego) || nrow(ego@result) == 0) {
        message("Skipping ", ligand_name, ": no enriched GO terms")
        return(NULL)
      }
      
      p <- clusterProfiler::dotplot(ego, showCategory = go_showCategory) +
        ggplot2::ggtitle(paste0(
          ligand_name, " predicted targets (NicheNet)\n",
          "Top targets ∩ receiver padj < ", go_padj_cutoff,
          " • mapped=", length(entrez)
        )) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
      
      ggplot2::ggsave(file.path(go_outdir, paste0("dotplot_", ligand_name, ".pdf")), p, width = 7, height = 5)
      ggplot2::ggsave(file.path(go_outdir, paste0("dotplot_", ligand_name, ".png")), p, width = 7, height = 5, dpi = 300)
      p
    }
    
    for (lig in top_ligands_plot) {
      p <- run_go_for_ligand(lig)
      if (!is.null(p)) dotplots_list[[lig]] <- p
    }
    
    if (length(dotplots_list) > 0) {
      grDevices::pdf(file.path(go_outdir, "all_ligands_GO.pdf"), width = 7, height = 5)
      for (p in dotplots_list) print(p)
      grDevices::dev.off()
    }
  }
  
  # Exports
  if (export_plots) {
    ggplot2::ggsave(file.path(outdir_final, "nichenet_ligand_activity_scatter.png"),
                    p_activity, width = 6, height = 4, dpi = 300)
    ggplot2::ggsave(file.path(outdir_final, paste0("nichenet_top", top_n_ligands_plot, "_ligands.png")),
                    p_top20, width = 7, height = 5, dpi = 300)
    ggplot2::ggsave(file.path(outdir_final, paste0("nichenet_sender_fc_top", top_n_ligands_plot, ".png")),
                    p_sender_fc, width = 7, height = 5, dpi = 300)
    
    ggplot2::ggsave(file.path(outdir_final, "nichenet_ligand_activity_scatter.pdf"),
                    p_activity, width = 6, height = 4)
    ggplot2::ggsave(file.path(outdir_final, paste0("nichenet_top", top_n_ligands_plot, "_ligands.pdf")),
                    p_top20, width = 7, height = 5)
    ggplot2::ggsave(file.path(outdir_final, paste0("nichenet_sender_fc_top", top_n_ligands_plot, ".pdf")),
                    p_sender_fc, width = 7, height = 5)
    
    if (!is.null(p_summary)) {
      ggplot2::ggsave(file.path(outdir_final, "p_nichenet_summary.png"),
                      p_summary, width = 14, height = 6, dpi = 600)
      ggplot2::ggsave(file.path(outdir_final, "p_nichenet_summary.pdf"),
                      p_summary, width = 14, height = 6)
    }
  }
  
  if (export_tables) {
    readr::write_tsv(ligand_activities2, file.path(outdir_final, "ligand_activity.tsv"))
    readr::write_tsv(ligand_target_links, file.path(outdir_final, paste0("ligand_target_links_top", top_n_ligands_targets, "x", top_targets_per_ligand, ".tsv")))
    
    qc <- c(
      paste0("Sender label: ", sender_label),
      paste0("Receiver label: ", receiver_label),
      paste0("Receiver DE genes: ", length(receiver_DE)),
      paste0("Receiver background genes: ", length(receiver_bg)),
      paste0("Sender ligand mode: ", sender_ligand_mode),
      paste0("Sender candidate ligands used: ", length(sender_candidates)),
      paste0("Resources mode: ", resources_mode),
      paste0("Top ligands plot: ", top_n_ligands_plot),
      paste0("Top ligands targets: ", top_n_ligands_targets),
      paste0("Top targets per ligand: ", top_targets_per_ligand),
      paste0("GO dotplots: ", do_go_dotplots),
      paste0("Output folder: ", outdir_final)
    )
    writeLines(qc, file.path(outdir_final, "qc_summary.txt"))
  }
  
  if (save_rds) {
    saveRDS(list(
      meta = list(
        sender_file = sender_file,
        receiver_file = receiver_file,
        sender_label = sender_label,
        receiver_label = receiver_label,
        outdir = outdir_final
      ),
      sender_df = sender_df,
      receiver_df = receiver_df,
      receiver_DE = receiver_DE,
      receiver_bg = receiver_bg,
      ligand_activities = ligand_activities2,
      ligand_target_links = ligand_target_links,
      plots = list(
        activity = p_activity,
        top_ligands = p_top20,
        sender_fc = p_sender_fc,
        summary = p_summary
      ),
      go_dotplots = dotplots_list
    ), file = file.path(outdir_final, "run_objects.rds"))
  }
  
  message("==> Done. Outputs written to: ", outdir_final)
  
  invisible(list(
    outdir = outdir_final,
    ligand_activities = ligand_activities2,
    ligand_target_links = ligand_target_links,
    go_dotplots = dotplots_list
  ))
}

# ---------------------------
# Example usage
# ---------------------------
# run_nichenet_bulk(
#   sender_file   = "sender_DE.xlsx",
#   receiver_file = "receiver_DE.xlsx",
#   sender_label = "Astrocytes_RNA",
#   receiver_label = "Tumor_RNA",
#   prefix = "RNA_RNA",
#   outdir = ".",
#   resources_mode = "package",  # or "zenodo"
#   sender_ligand_mode = "DE_up",
#   do_go_dotplots = TRUE,
#   organism_db = "org.Mm.eg.db",
#   go_ont = "BP"
# )
