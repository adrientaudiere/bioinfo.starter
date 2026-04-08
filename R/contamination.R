#' @export
print.wrapped_pheatmap <- function(x, ...) {
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
}

#' Classify taxa found in negative controls into contamination categories
#'
#' @param physeq A phyloseq object
#' @param neg_control Logical vector identifying negative control samples,
#'   or a character string naming a column in sam_data
#' @param min_reads_artifact Max total reads across all samples to be considered artifact
#' @param max_samples_artifact Max number of samples (total) for artifact classification
#' @param max_ratio_lab_contam Max ratio (non-neg occurrence / neg occurrence) for lab contaminant
#' @param min_neg_samples_lab Min number of neg control samples a taxon must appear in for lab contaminant
#' @return A named list with $classification (data.frame), $summary_plot, $distribution_plot, $heatmap_by_category
classify_neg_control_taxa <- function(
    physeq,
    neg_control,
    min_reads_artifact = 10,
    max_samples_artifact = 2,
    max_ratio_lab_contam = 0.2,
    min_neg_samples_lab = 2) {
  requireNamespace("phyloseq", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("pheatmap", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)

  # Resolve neg_control to logical vector

if (is.character(neg_control) && length(neg_control) == 1) {
    neg_control <- as.logical(phyloseq::sample_data(physeq)[[neg_control]])
  }
  stopifnot(is.logical(neg_control), length(neg_control) == phyloseq::nsamples(physeq))

  otu <- as(phyloseq::otu_table(physeq), "matrix")
  if (!phyloseq::taxa_are_rows(physeq)) otu <- t(otu)
  # otu: taxa x samples

  neg_idx <- which(neg_control)
  non_neg_idx <- which(!neg_control)

  otu_neg <- otu[, neg_idx, drop = FALSE]
  otu_non_neg <- otu[, non_neg_idx, drop = FALSE]

  # Taxa present in at least one negative control
  taxa_in_neg <- rownames(otu_neg)[rowSums(otu_neg > 0) >= 1]

  if (length(taxa_in_neg) == 0) {
    message("No taxa found in negative control samples.")
    return(list(
      classification = data.frame(),
      summary_plot = ggplot2::ggplot(),
      distribution_plot = ggplot2::ggplot(),
      heatmap_by_category = NULL
    ))
  }

  # Build classification data frame
  clf <- data.frame(
    taxon = taxa_in_neg,
    total_reads = rowSums(otu[taxa_in_neg, , drop = FALSE]),
    reads_in_neg = rowSums(otu_neg[taxa_in_neg, , drop = FALSE]),
    reads_in_samples = rowSums(otu_non_neg[taxa_in_neg, , drop = FALSE]),
    n_neg_samples = rowSums(otu_neg[taxa_in_neg, , drop = FALSE] > 0),
    n_non_neg_samples = rowSums(otu_non_neg[taxa_in_neg, , drop = FALSE] > 0),
    stringsAsFactors = FALSE
  )
  clf$ratio_non_neg_to_neg <- clf$n_non_neg_samples / clf$n_neg_samples

  # Classification logic
  clf$category <- NA_character_

  # 1. Artifact
  is_artifact <- clf$total_reads <= min_reads_artifact &
    (clf$n_neg_samples + clf$n_non_neg_samples) <= max_samples_artifact
  clf$category[is_artifact] <- "artifact"

  # 2. Lab contaminant (among non-artifacts)
  remaining <- is.na(clf$category)
  is_lab <- remaining &
    clf$ratio_non_neg_to_neg <= max_ratio_lab_contam &
    clf$n_neg_samples >= min_neg_samples_lab
  clf$category[is_lab] <- "lab_contaminant"

  # 3. Sample contaminant (everything else)
  clf$category[is.na(clf$category)] <- "sample_contaminant"

  clf$category <- factor(clf$category,
    levels = c("artifact", "lab_contaminant", "sample_contaminant")
  )

  # Add taxonomy if available
  tax <- tryCatch(
    as.data.frame(phyloseq::tax_table(physeq)[taxa_in_neg, ]),
    error = function(e) NULL
  )
  if (!is.null(tax)) {
    tax$taxon <- rownames(tax)
    clf <- dplyr::left_join(clf, tax, by = "taxon")
  }

  # Order by category then total_reads descending
  clf <- clf[order(clf$category, -clf$total_reads), ]
  rownames(clf) <- NULL

  # --- Plots ---

  # Summary plot: taxa count + total reads per category
  cat_summary <- clf |>
    dplyr::group_by(category) |>
    dplyr::summarise(
      n_taxa = dplyr::n(),
      total_reads = sum(total_reads),
      .groups = "drop"
    )

  p_count <- ggplot2::ggplot(cat_summary, ggplot2::aes(x = category, y = n_taxa, fill = category)) +
    ggplot2::geom_col() +
    ggplot2::geom_text(ggplot2::aes(label = n_taxa), vjust = -0.3) +
    ggplot2::labs(title = "Number of taxa per category", x = NULL, y = "Number of taxa") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  p_reads <- ggplot2::ggplot(cat_summary, ggplot2::aes(x = category, y = total_reads, fill = category)) +
    ggplot2::geom_col() +
    ggplot2::geom_text(ggplot2::aes(label = total_reads), vjust = -0.3) +
    ggplot2::labs(title = "Total reads per category", x = NULL, y = "Total reads") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  summary_plot <- patchwork::wrap_plots(p_count, p_reads, ncol = 2)

  # Distribution plot: violin of reads per taxon by category + non-contaminant taxa
  taxa_not_in_neg <- setdiff(rownames(otu), taxa_in_neg)
  if (length(taxa_not_in_neg) > 0) {
    non_contam_df <- data.frame(
      taxon = taxa_not_in_neg,
      total_reads = rowSums(otu[taxa_not_in_neg, , drop = FALSE]),
      category = factor("non_contaminant",
        levels = c("artifact", "lab_contaminant", "sample_contaminant", "non_contaminant")
      ),
      stringsAsFactors = FALSE
    )
    distrib_df <- rbind(
      clf[, c("taxon", "total_reads", "category")],
      non_contam_df
    )
    distrib_df$category <- factor(distrib_df$category,
      levels = c("artifact", "lab_contaminant", "sample_contaminant", "non_contaminant")
    )
  } else {
    distrib_df <- clf[, c("taxon", "total_reads", "category")]
  }

  distribution_plot <- ggplot2::ggplot(distrib_df, ggplot2::aes(x = category, y = total_reads + 1, fill = category)) +
    ggplot2::geom_violin(alpha = 0.5, scale = "width") +
    ggplot2::geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_fill_manual(values = c(
      artifact = "#66c2a5", lab_contaminant = "#fc8d62",
      sample_contaminant = "#8da0cb", non_contaminant = "#a6d854"
    )) +
    ggplot2::labs(
      title = "Distribution of total reads per taxon by category",
      x = NULL, y = "Total reads (log scale)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  # Heatmap of log-transformed abundances in neg controls, annotated by category
  mat <- otu_neg[taxa_in_neg, , drop = FALSE]
  mat_log <- log1p(mat)

  # Rename columns to include taxa count per sample
  n_taxa_per_sample <- colSums(mat > 0)
  colnames(mat_log) <- paste0(colnames(mat_log), " (", n_taxa_per_sample, ")")

  annotation_row <- data.frame(
    category = clf$category,
    row.names = clf$taxon
  )
  annotation_row <- annotation_row[rownames(mat_log), , drop = FALSE]

  ann_colors <- list(
    category = c(
      artifact = "#66c2a5",
      lab_contaminant = "#fc8d62",
      sample_contaminant = "#8da0cb"
    )
  )

  heatmap_obj <- pheatmap::pheatmap(
    mat_log,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = nrow(mat_log) <= 60,
    show_colnames = TRUE,
    color = grDevices::colorRampPalette(c("white", "firebrick3"))(100),
    annotation_row = annotation_row,
    annotation_colors = ann_colors,
    main = "Neg. control abundances (log) by contamination category",
    silent = TRUE
  )

  # Wrap pheatmap as a printable object that handles grid properly
  heatmap_wrapped <- structure(
    list(gtable = heatmap_obj$gtable),
    class = "wrapped_pheatmap"
  )

  # Ordination plot highlighting negative controls
  ordination_plot <- tryCatch({
    physeq_clean <- phyloseq::prune_taxa(phyloseq::taxa_sums(physeq) > 0, physeq)
    physeq_clean <- phyloseq::prune_samples(phyloseq::sample_sums(physeq_clean) > 0, physeq_clean)
    physeq_log <- phyloseq::transform_sample_counts(physeq_clean, function(x) log1p(x))
    ord <- phyloseq::ordinate(physeq_log, method = "MDS", distance = "bray")

    samp_df <- data.frame(
      phyloseq::sample_data(physeq_clean),
      check.names = FALSE
    )
    samp_df$neg_control <- neg_control[match(
      phyloseq::sample_names(physeq_clean),
      phyloseq::sample_names(physeq)
    )]
    samp_df$sample_label <- ifelse(samp_df$neg_control,
      phyloseq::sample_names(physeq_clean), NA_character_
    )

    coords <- ord$vectors[, 1:2]
    colnames(coords) <- c("Axis1", "Axis2")
    samp_df <- cbind(samp_df, coords[phyloseq::sample_names(physeq_clean), ])

    eig <- ord$values$Eigenvalues
    var_pct <- round(100 * eig / sum(eig), 1)

    ggplot2::ggplot(samp_df, ggplot2::aes(x = Axis1, y = Axis2)) +
      ggplot2::geom_point(
        ggplot2::aes(color = neg_control, shape = neg_control),
        size = 3, alpha = 0.8
      ) +
      ggplot2::geom_text(
        ggplot2::aes(label = sample_label),
        vjust = -0.8, size = 3, na.rm = TRUE
      ) +
      ggplot2::scale_color_manual(
        values = c("FALSE" = "grey60", "TRUE" = "red"),
        labels = c("Sample", "Negative control")
      ) +
      ggplot2::scale_shape_manual(
        values = c("FALSE" = 16, "TRUE" = 17),
        labels = c("Sample", "Negative control")
      ) +
      ggplot2::labs(
        title = "PCoA ordination (Bray-Curtis) — negative controls highlighted",
        x = paste0("Axis 1 (", var_pct[1], "%)"),
        y = paste0("Axis 2 (", var_pct[2], "%)"),
        color = NULL, shape = NULL
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "bottom")
  }, error = function(e) {
    message("Could not compute ordination plot: ", e$message)
    NULL
  })

  list(
    classification = clf,
    summary_plot = summary_plot,
    distribution_plot = distribution_plot,
    heatmap_by_category = heatmap_wrapped,
    ordination_plot = ordination_plot
  )
}
