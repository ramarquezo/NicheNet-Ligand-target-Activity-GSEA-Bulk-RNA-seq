# ============================================================
# NicheNet + Target-GSEA (GitHub-friendly, cross-platform)
# ------------------------------------------------------------
# What this script does:
#  1) Reads two DE tables: sender and receiver (generic labels)
#  2) Runs bulk NicheNet ligand activity (sender -> receiver)
#  3) Plots top ligands + sender logFC evidence
#  4) Builds "ligand target gene sets" (top N targets per ligand)
#  5) Runs GSEA on receiver ranked list using those target sets
#     (i.e., are predicted targets coordinately up/down in receiver?)
#  6) Exports: TSV + PNG/PDF + optional RDS
#
# Notes:
#  - Default resource loading uses nichenetr packaged datasets (portable).
#  - Optional Zenodo download mode is included for robustness.
#  - No choose.dir(); output goes to outdir + timestamped folder.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readxl)
  library(readr)
  library(ggplot2)
  library(clusterProfiler)
  library(nichenetr)
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

  # common alternate casing
  if ("log2Foldchange" %in% colnames(df) && !lfc_col %in% colnames(df)) {
    colnames(df)[colnames(df) == "log2Foldchange"] <- lfc_col
  }
  if ("gene" %in% colnames(df) && !symbol_col %in% colnames(df)) {
    colnames(df)[colnames(df) == "gene"] <- symbol_col
  }

  missing <- setdiff(c(symbol_col, lfc_col, padj_col), colnames(df))
  if (length(missing)) stop("DE table missing required columns: ", paste(missing, collapse=", "))

  df %>%
    mutate(
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

collapse_rank_by_log2fc <- function(df, symbol_col="symbol", lfc_col="log2FoldChange", to_upper=TRUE) {
  out <- df %>%
    mutate(
      .sym = clean_symbols(.data[[symbol_col]]),
      .lfc = as.numeric(.data[[lfc_col]])
    ) %>%
    filter(!is.na(.sym), is.finite(.lfc)) %>%
    mutate(.sym = if (to_upper) toupper(.sym) else .sym) %>%
    group_by(.sym) %>%
    slice_max(order_by = abs(.lfc), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(.lfc)) %>%
    select(.sym, .lfc) %>%
    deframe()

  stopifnot(!anyDuplicated(names(out)))
  out
}

get_DE_genes <- function(df, symbol_col="symbol", lfc_col="log2FoldChange", padj_col="padj",
                         padj_cut=0.05, lfc_cut=0.25, direction=c("both","up","down"),
                         to_upper=TRUE) {
  direction <- match.arg(direction)

  out <- df %>%
    mutate(
      .sym = clean_symbols(.data[[symbol_col]]),
      .lfc = as.numeric(.data[[lfc_col]]),
      .padj = as.numeric(.data[[padj_col]])
    ) %>%
    filter(!is.na(.sym), is.finite(.lfc), is.finite(.padj))

  out <- if (direction == "both") {
    out %>% filter(.padj < padj_cut, abs(.lfc) > lfc_cut)
  } else if (direction == "up") {
    out %>% filter(.padj < padj_cut, .lfc > lfc_cut)
  } else {
    out %>% filter(.padj < padj_cut, .lfc < -lfc_cut)
  }

  genes <- out %>% pull(.sym) %>% unique()
  if (to_upper) genes <- toupper(genes)
  genes
}

get_background_genes <- function(df, symbol_col="symbol", lfc_col="log2FoldChange", to_upper=TRUE) {
  genes <- df %>%
    mutate(
      .sym = clean_symbols(.data[[symbol_col]]),
      .lfc = as.numeric(.data[[lfc_col]])
    ) %>%
    filter(!is.na(.sym), is.finite(.lfc)) %>%
    pull(.sym) %>%
    unique()
  if (to_upper) genes <- toupper(genes)
  genes
}

# ---------------------------
# NicheNet resources loader
# ---------------------------
load_nichenet_resources <- function(mode=c("package","zenodo"),
                                    nichenet_dir="nichenet_resources") {
  mode <- match.arg(mode)

  if (mode == "package") {
    # Most portable: uses datasets shipped with nichenetr
    message("Loading NicheNet resources from the nichenetr package datasets...")
    data("ligand_target_matrix", package = "nichenetr", envir = environment())
    data("lr_network",          package = "nichenetr", envir = environment())
    data("weighted_networks",   package = "nichenetr", envir = environment())

    ligand_target_matrix <- get("ligand_target_matrix", envir = environment())
    lr_network <- get("lr_network", envir = environment())
    weighted_networks <- get("weighted_networks", envir = environment())

    return(list(
      ligand_target_matrix = ligand_target_matrix,
      lr_network = lr_network,
      weighted_networks = weighted_networks
    ))
  }

  # Optional: Zenodo mode (kept for users who prefer explicit files)
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
# Main runner
# ---------------------------
run_nichenet_target_gsea <- function(
    # Inputs
    sender_file = NULL,
    receiver_file = NULL,
    sender_label = "Sender",
    receiver_label = "Receiver",

    # Output
    prefix = "NicheNet_Run",
    outdir = ".",
    timestamp_dir = TRUE,

    # Column mapping
    symbol_col = "symbol",
    lfc_col = "log2FoldChange",
    padj_col = "padj",

    # Receiver DE signature (for NicheNet)
    receiver_padj_cut = 0.05,
    receiver_lfc_cut  = 0.25,
    receiver_direction = c("both","up","down"),

    # Sender ligand prior
    sender_ligand_mode = c("DE_up","DE_both","all_expressed"),
    sender_padj_cut = 0.10,
    sender_lfc_cut  = 0.25,

    # NicheNet settings
    resources_mode = c("package","zenodo"),
    to_upper = TRUE,
    top_n_ligands = 20,
    top_targets_per_ligand = 50,  # for target gene sets

    # Target-GSEA settings (ligand target sets vs receiver rank)
    target_gsea_padj_cutoff = 0.1,
    minGSSize = 10,
    maxGSSize = 500,

    # Exports
    export_plots = TRUE,
    export_tables = TRUE,
    save_rds = TRUE
) {
  receiver_direction <- match.arg(receiver_direction)
  sender_ligand_mode <- match.arg(sender_ligand_mode)
  resources_mode <- match.arg(resources_mode)

  # ---- Resolve inputs ----
  if (is.null(sender_file)) {
    message("Choose SENDER Excel file...")
    sender_file <- file.choose()
  }
  if (is.null(receiver_file)) {
    message("Choose RECEIVER Excel file...")
    receiver_file <- file.choose()
  }

  prefix_safe <- safe_name(prefix)
  sender_label_safe <- safe_name(sender_label)
  receiver_label_safe <- safe_name(receiver_label)

  outdir_final <- make_outdir(outdir = outdir, prefix = prefix_safe, stamp = timestamp_dir)
  message("Outputs will be written to: ", outdir_final)

  # ---- Read + standardize ----
  sender_df <- readxl::read_excel(sender_file) %>%
    standardize_de(symbol_col, lfc_col, padj_col)
  receiver_df <- readxl::read_excel(receiver_file) %>%
    standardize_de(symbol_col, lfc_col, padj_col)

  if (to_upper) {
    sender_df[[symbol_col]] <- toupper(sender_df[[symbol_col]])
    receiver_df[[symbol_col]] <- toupper(receiver_df[[symbol_col]])
  }

  # ---- Receiver signature + background ----
  receiver_DE <- get_DE_genes(
    receiver_df, symbol_col, lfc_col, padj_col,
    padj_cut = receiver_padj_cut, lfc_cut = receiver_lfc_cut,
    direction = receiver_direction, to_upper = to_upper
  )
  receiver_bg <- get_background_genes(receiver_df, symbol_col, lfc_col, to_upper = to_upper)

  message("Receiver DE genes: ", length(receiver_DE))
  message("Receiver background genes: ", length(receiver_bg))

  # ---- Load NicheNet resources ----
  res <- load_nichenet_resources(mode = resources_mode)
  ligand_target_matrix <- res$ligand_target_matrix
  lr_network <- res$lr_network

  # Normalize symbols/case
  lr_network <- lr_network %>%
    mutate(from = if (to_upper) toupper(from) else from,
           to   = if (to_upper) toupper(to) else to)

  rownames(ligand_target_matrix) <- if (to_upper) toupper(rownames(ligand_target_matrix)) else rownames(ligand_target_matrix)
  colnames(ligand_target_matrix) <- if (to_upper) toupper(colnames(ligand_target_matrix)) else colnames(ligand_target_matrix)

  all_ligands <- unique(lr_network$from)

  # ---- Sender ligand candidates ----
  sender_bg <- get_background_genes(sender_df, symbol_col, lfc_col, to_upper = to_upper)

  sender_de_tbl <- sender_df %>%
    transmute(
      ligand = .data[[symbol_col]],
      sender_logFC = .data[[lfc_col]],
      sender_padj  = .data[[padj_col]]
    ) %>%
    filter(!is.na(ligand), is.finite(sender_logFC), is.finite(sender_padj)) %>%
    group_by(ligand) %>%
    slice_max(order_by = abs(sender_logFC), n = 1, with_ties = FALSE) %>%
    ungroup()

  sender_DE_up <- sender_de_tbl %>%
    filter(sender_padj < sender_padj_cut, sender_logFC > sender_lfc_cut) %>%
    pull(ligand) %>% unique()

  sender_DE_both <- sender_de_tbl %>%
    filter(sender_padj < sender_padj_cut, abs(sender_logFC) > sender_lfc_cut) %>%
    pull(ligand) %>% unique()

  sender_candidates <- switch(
    sender_ligand_mode,
    "all_expressed" = intersect(sender_bg, all_ligands),
    "DE_up"         = intersect(sender_DE_up, all_ligands),
    "DE_both"       = intersect(sender_DE_both, all_ligands)
  )

  message("Sender ligand mode: ", sender_ligand_mode)
  message("Sender candidate ligands used: ", length(sender_candidates))

  if (length(sender_candidates) < 5) {
    warning("Very few sender candidate ligands. Consider relaxing sender cutoffs or using sender_ligand_mode='all_expressed'.")
  }

  # ---- NicheNet ligand activity ----
  message("Running NicheNet ligand activity inference...")
  ligand_activities <- predict_ligand_activities(
    geneset = receiver_DE,
    background_expressed_genes = receiver_bg,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = sender_candidates
  ) %>%
    arrange(desc(aupr_corrected)) %>%
    mutate(rank = row_number())

  ligand_activities2 <- ligand_activities %>%
    left_join(sender_de_tbl, by = c("test_ligand" = "ligand")) %>%
    rename(ligand = test_ligand)

  top_ligs <- ligand_activities2 %>%
    slice_head(n = top_n_ligands) %>%
    pull(ligand)

  # ---- Ligand -> targets table ----
  get_top_targets <- function(ligand, ltm, n = 200) {
    if (!ligand %in% rownames(ltm)) return(tibble())
    targets <- sort(ltm[ligand, ], decreasing = TRUE)
    tibble(
      ligand = ligand,
      target = names(targets)[seq_len(min(n, length(targets)))],
      regulatory_potential = unname(targets[seq_len(min(n, length(targets)))])
    )
  }

  ligand_target_links <- bind_rows(lapply(
    top_ligs,
    get_top_targets,
    ltm = ligand_target_matrix,
    n = max(200, top_targets_per_ligand)
  )) %>%
    mutate(target_in_receiver_DE = target %in% receiver_DE)

  # ---- Plots: ligand activity + top ligands + sender FC evidence ----
  p_activity_scatter <- ggplot(ligand_activities2, aes(x = rank, y = aupr_corrected)) +
    geom_point(alpha = 0.9) +
    theme_classic() +
    labs(
      title = paste0("NicheNet ligand activity (", sender_label, " → ", receiver_label, ")"),
      subtitle = paste0("sender_ligand_mode=", sender_ligand_mode,
                        " • receiver_DE=", receiver_direction,
                        " • n_candidates=", length(sender_candidates)),
      x = "Ligand rank",
      y = "AUPR corrected"
    )

  top20_df <- ligand_activities2 %>% slice_head(n = top_n_ligands)

  p_top_ligands <- ggplot(top20_df, aes(x = aupr_corrected, y = reorder(ligand, aupr_corrected))) +
    geom_col(width = 0.75, alpha = 0.9) +
    theme_classic(base_size = 12) +
    labs(
      title = paste0("Top ", top_n_ligands, " predicted ligands"),
      subtitle = paste0(sender_label, " → ", receiver_label),
      x = "AUPR corrected",
      y = "Ligand"
    )

  p_sender_fc <- ggplot(
    top20_df %>%
      mutate(
        rank = as.numeric(rank),
        y_rank = factor(rank, levels = rev(sort(unique(rank)))),
        padj_lbl = paste0("p-adj=", formatC(sender_padj, format = "f", digits = 6))
      ),
    aes(x = sender_logFC, y = y_rank)
  ) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, alpha = 0.5) +
    geom_col(width = 0.75, alpha = 0.9) +
    geom_text(
      aes(label = padj_lbl,
          x = sender_logFC + 0.03 * diff(range(sender_logFC, na.rm = TRUE))),
      hjust = 0, size = 3
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 12) +
    theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin  = margin(5.5, 80, 5.5, 5.5)
    ) +
    labs(x = paste0(sender_label, " log2FC"), y = NULL,
         title = paste0("Sender evidence for top ligands"))

  # Optional composite if patchwork is installed
  p_summary <- NULL
  if (requireNamespace("patchwork", quietly = TRUE)) {
    p_summary <- (p_top_ligands | p_sender_fc) + patchwork::plot_layout(widths = c(1, 1))
  }

  # ---- Target-GSEA: ligand target sets vs receiver rank ----
  message("Running target-GSEA (ligand target sets against receiver rank)...")

  receiver_rank <- collapse_rank_by_log2fc(receiver_df, symbol_col, lfc_col, to_upper)

  # Build TERM2GENE: each ligand is a "geneset" of its top targets (by regulatory potential)
  term2gene_targets <- ligand_target_links %>%
    group_by(ligand) %>%
    slice_max(regulatory_potential, n = top_targets_per_ligand, with_ties = FALSE) %>%
    ungroup() %>%
    select(gs_name = ligand, gene_symbol = target) %>%
    distinct()

  # Run GSEA where "pathways" are ligand target sets
  gsea_targets <- clusterProfiler::GSEA(
    geneList = receiver_rank,
    TERM2GENE = term2gene_targets,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    pvalueCutoff = 1,
    verbose = FALSE
  )

  gsea_targets_df <- as.data.frame(gsea_targets) %>%
    arrange(p.adjust, desc(abs(NES)))

  gsea_targets_sig <- gsea_targets_df %>%
    filter(p.adjust < target_gsea_padj_cutoff)

  # Plot top ligands' target-GSEA (if available)
  gsea_plot_df <- gsea_targets_df %>%
    filter(ID %in% top_ligs) %>%
    mutate(ID = factor(ID, levels = rev(top_ligs)))

  p_target_gsea <- NULL
  if (nrow(gsea_plot_df) > 0) {
    p_target_gsea <- ggplot(gsea_plot_df, aes(x = NES, y = ID)) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, alpha = 0.5) +
      geom_point(aes(size = setSize), alpha = 0.9) +
      theme_classic(base_size = 12) +
      labs(
        title = "Target-GSEA for top ligands",
        subtitle = paste0("Gene sets = top ", top_targets_per_ligand, " predicted targets per ligand • ranked list = receiver log2FC"),
        x = "NES (targets enriched in receiver rank)",
        y = "Ligand"
      )
  }

  # ---- Exports ----
  if (export_plots) {
    ggsave(file.path(outdir_final, paste0("nichenet_activity_scatter_", prefix_safe, ".png")),
           p_activity_scatter, width = 7, height = 5, dpi = 300)
    ggsave(file.path(outdir_final, paste0("nichenet_top_ligands_", prefix_safe, ".png")),
           p_top_ligands, width = 7, height = 5, dpi = 300)
    ggsave(file.path(outdir_final, paste0("nichenet_sender_fc_top_ligands_", prefix_safe, ".png")),
           p_sender_fc, width = 7, height = 5, dpi = 300)

    ggsave(file.path(outdir_final, paste0("nichenet_activity_scatter_", prefix_safe, ".pdf")),
           p_activity_scatter, width = 7, height = 5)
    ggsave(file.path(outdir_final, paste0("nichenet_top_ligands_", prefix_safe, ".pdf")),
           p_top_ligands, width = 7, height = 5)
    ggsave(file.path(outdir_final, paste0("nichenet_sender_fc_top_ligands_", prefix_safe, ".pdf")),
           p_sender_fc, width = 7, height = 5)

    if (!is.null(p_summary)) {
      ggsave(file.path(outdir_final, paste0("nichenet_summary_", prefix_safe, ".png")),
             p_summary, width = 14, height = 6, dpi = 600)
      ggsave(file.path(outdir_final, paste0("nichenet_summary_", prefix_safe, ".pdf")),
             p_summary, width = 14, height = 6)
    }

    if (!is.null(p_target_gsea)) {
      ggsave(file.path(outdir_final, paste0("target_gsea_top_ligands_", prefix_safe, ".png")),
             p_target_gsea, width = 7, height = 5, dpi = 300)
      ggsave(file.path(outdir_final, paste0("target_gsea_top_ligands_", prefix_safe, ".pdf")),
             p_target_gsea, width = 7, height = 5)
    }
  }

  if (export_tables) {
    write_tsv(ligand_activities2, file.path(outdir_final, paste0("ligand_activity_", prefix_safe, ".tsv")))
    write_tsv(ligand_target_links, file.path(outdir_final, paste0("ligand_target_links_", prefix_safe, ".tsv")))
    write_tsv(gsea_targets_df, file.path(outdir_final, paste0("target_gsea_all_", prefix_safe, ".tsv")))
    write_tsv(gsea_targets_sig, file.path(outdir_final, paste0("target_gsea_sig_", prefix_safe, ".tsv")))

    qc <- c(
      paste0("Sender label: ", sender_label),
      paste0("Receiver label: ", receiver_label),
      paste0("Sender ligand mode: ", sender_ligand_mode),
      paste0("Receiver DE direction: ", receiver_direction),
      paste0("Receiver DE genes: ", length(receiver_DE)),
      paste0("Receiver background genes: ", length(receiver_bg)),
      paste0("Sender candidate ligands: ", length(sender_candidates)),
      paste0("Top ligands plotted: ", top_n_ligands),
      paste0("Top targets per ligand (for target-GSEA): ", top_targets_per_ligand),
      paste0("NicheNet resources mode: ", resources_mode),
      paste0("Output folder: ", outdir_final)
    )
    writeLines(qc, file.path(outdir_final, paste0("qc_summary_", prefix_safe, ".txt")))
  }

  if (save_rds) {
    saveRDS(
      list(
        meta = list(
          sender_file = sender_file,
          receiver_file = receiver_file,
          sender_label = sender_label,
          receiver_label = receiver_label,
          prefix = prefix_safe,
          outdir = outdir_final,
          symbol_col = symbol_col,
          lfc_col = lfc_col,
          padj_col = padj_col
        ),
        sender_df = sender_df,
        receiver_df = receiver_df,
        receiver_DE = receiver_DE,
        receiver_bg = receiver_bg,
        ligand_activities = ligand_activities2,
        ligand_target_links = ligand_target_links,
        term2gene_targets = term2gene_targets,
        gsea_targets = gsea_targets,
        gsea_targets_df = gsea_targets_df,
        plots = list(
          activity_scatter = p_activity_scatter,
          top_ligands = p_top_ligands,
          sender_fc = p_sender_fc,
          summary = p_summary,
          target_gsea = p_target_gsea
        )
      ),
      file = file.path(outdir_final, paste0("run_", prefix_safe, "_objects.rds"))
    )
  }

  message("Done. Outputs written to: ", outdir_final)

  invisible(list(
    outdir = outdir_final,
    ligand_activities = ligand_activities2,
    ligand_target_links = ligand_target_links,
    gsea_targets_df = gsea_targets_df,
    gsea_targets_sig = gsea_targets_sig,
    plots = list(
      activity_scatter = p_activity_scatter,
      top_ligands = p_top_ligands,
      sender_fc = p_sender_fc,
      summary = p_summary,
      target_gsea = p_target_gsea
    )
  ))
}

# ---------------------------
# Example usage
# ---------------------------
# source("nichenet_target_gsea.R")
#
# run_nichenet_target_gsea(
#   sender_file   = "sender_DE.xlsx",
#   receiver_file = "receiver_DE.xlsx",
#   sender_label = "Astrocytes_RNA",
#   receiver_label = "Tumor_RNA",
#   prefix = "RNA_RNA",
#   outdir = ".",
#   resources_mode = "package",      # "package" or "zenodo"
#   sender_ligand_mode = "DE_up",
#   top_n_ligands = 20,
#   top_targets_per_ligand = 50,
#   target_gsea_padj_cutoff = 0.1
# )

