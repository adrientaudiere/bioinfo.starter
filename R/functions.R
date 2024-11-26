# Add some functions with documentation here

# fastq_dir <- system.file("extdata", package = "dada2")
# fastqc_agg(fastq_dir, qc.dir = "~/essai_fastqc")
fastqc_agg <- function(fq.dir = getwd(),
                       qc.dir = NULL,
                       threads = 4,
                       fastqc.path = "~/bin/FastQC/fastqc",
                       aggregate = TRUE,
                       progressbar = FALSE,
                       multiqc = FALSE) {
  fastqcr::fastqc(
    fq.dir = fq.dir,
    qc.dir = qc.dir,
    threads = threads,
    fastqc.path = fastqc.path
  )

  if (aggregate) {
    res <- fastqcr::qc_aggregate(qc.dir, progressbar = progressbar)
  }

  if (multiqc) {
    system(paste0("cd ", qc.dir, "; multiqc ."))
  }

  message(paste("Fastqc reports are available in the folder", qc.dir, "."))
  return(res)
}



# fastq_dir <- system.file("extdata", package = "dada2")
# fastqc_agg(fastq_dir, qc.dir = "~/essai_fastqc")
# fastqc_plot("~/essai_fastqc")

fastqc_plot <- function(qc.dir,
                        modules = c(
                          "Basic Statistics",
                          "Per base sequence quality",
                          "Per base sequence content",
                          "Sequence Length Distribution",
                          "Adapter Content"
                        )) {
  p <- fastqcr::qc_plot_collection(
    fastqcr::qc_read_collection(
      list.files(qc.dir, pattern = "*.zip", full.names = TRUE),
      list.files(qc.dir, pattern = "*.zip"),
      modules = modules, verbose = FALSE
    ),
    modules = modules
  )

  return(p)
}



ggvalue_box <- function(
    values = NULL,
    names = NULL,
    colors = NULL,
    icons = NULL,
    icons_png = NULL,
    tile_width = 0.8,
    tile_height = 0.5,
    value_nudge_x = 0.7,
    value_nudge_y = 1.1,
    value_font_size = 25,
    value_font_color = "white",
    lab_nudge_x = 0.7,
    lab_nudge_y = 0.9,
    lab_font_size = 8,
    lab_font_color = "white",
    icon_nudge_x = 0.73,
    icon_nudge_y = 0.9,
    icon_size = 20,
    alpha_icon = 0.5,
    negate_icon_png = FALSE,
    color_icon = "grey20") {
  if (is.null(colors)) {
    colors <- funky_color(length(values))
  }
  df <-
    data.frame(
      "value" = values,
      "lab" = names,
      "color" = colors
    )
  if (!is.null(icons)) {
    df$icon <- icons
  }
  if (!is.null(icons_png)) {
    df$icon <- icons_png
  }
  p <- ggplot(df) +
    geom_tile(aes(lab, 1, fill = lab), width = tile_width, height = tile_height) +
    geom_text(
      aes(
        seq(
          value_nudge_x,
          length(value) + value_nudge_x - 1,
          length.out = length(value)
        ),
        value_nudge_y,
        label = value
      ),
      color = value_font_color,
      size = value_font_size,
      hjust = 0
    ) +
    scale_fill_manual(values = df$color) +
    geom_text(
      aes(
        seq(lab_nudge_x, length(lab) + lab_nudge_x - 1, length.out = length(lab)),
        lab_nudge_y,
        label = lab
      ),
      color = lab_font_color,
      size = lab_font_size,
      hjust = 0
    ) +
    coord_equal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank()
    ) +
    guides(fill = "none") +
    theme_void()

  if (!is.null(icons)) {
    library(emojifont)
    p <-
      p + geom_text(
        size = icon_size,
        aes(
          label = icon,
          family = "fontawesome-webfont",
          seq(
            icon_nudge_x * 1.7,
            length(lab) + icon_nudge_x * 1.7 - 1,
            length.out = length(lab)
          ),
          icon_nudge_y * 1.25
        ),
        alpha = alpha_icon,
        color = color_icon,
      )
  }
  if (!is.null(icons_png)) {
    transparent <- function(img) {
      magick::image_fx(img,
        expression = paste0(alpha_icon, "*a"),
        channel = "alpha"
      )
    }
    transparent_neg <- function(img) {
      magick::image_negate(magick::image_fx(img,
        expression = paste0(alpha_icon, "*a"),
        channel = "alpha"
      ))
    }
    if (negate_icon_png) {
      p <-
        p + ggimage::geom_image(
          aes(
            image = icon,
            seq(
              icon_nudge_x * 1.7,
              length(lab) + icon_nudge_x * 1.7 - 1,
              length.out = length(lab)
            ),
            icon_nudge_y * 1.25
          ),
          size = icon_size,
          image_fun = transparent_neg
        )
    } else {
      p <-
        p + ggimage::geom_image(
          aes(
            image = icon,
            seq(
              icon_nudge_x * 1.7,
              length(lab) + icon_nudge_x * 1.7 - 1,
              length.out = length(lab)
            ),
            icon_nudge_y * 1.25
          ),
          size = icon_size,
          image_fun = transparent
        )
    }
  }

  return(p)
}
