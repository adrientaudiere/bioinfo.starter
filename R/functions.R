# Add some functions with documentation here

#fastq_dir <- system.file("extdata", package = "dada2")
#fastqc_agg(fastq_dir, qc.dir = "~/essai_fastqc")

fastqc_agg <- function(fq.dir = getwd(),
                       qc.dir = NULL,
                       threads = 4,
                       fastqc.path = "~/bin/FastQC/fastqc",
                       aggregate = TRUE) {
  fastqcr::fastqc(
    fq.dir = fq.dir,
    qc.dir = qc.dir,
    threads = threads,
    fastqc.path = fastqc.path
  )

  if (aggregate) {
    res <-  fastqcr::qc_aggregate(qc.dir)
  }

  message(paste("Fastqc reports are available in the folder", qc.dir, "."))
  return(res)
}



#fastq_dir <- system.file("extdata", package = "dada2")
#fastqc_agg(fastq_dir, qc.dir = "~/essai_fastqc")
#fastqc_plot("~/essai_fastqc")

fastqc_plot <- function(qc.dir,
                        modules = c(
                          "Basic Statistics",
                          "Per base sequence quality",
                          "Per base sequence content",
                          "Sequence Length Distribution",
                          "Adapter Content"
                        )) {
  p <- fastqcr::qc_plot_collection(qc_read_collection(
    list.files(qc.dir, pattern = "*.zip", full.names = TRUE),
    list.files(qc.dir, pattern = "*.zip"),
    modules = modules
  ),
  modules = modules)

  return(p)
}




#' ggplot theme: based on [hrbrthemes](https://github.com/hrbrmstr/hrbrthemes/tree/master)
#' `hrbrthemes::theme_ipsum()` by boB Rudis.
#'
#'
#' @param sans_family
#' @param serif_family
#' @param mono_family
#' @param base_size
#' @param plot_title_family
#' @param plot_title_size
#' @param plot_title_face
#' @param plot_title_margin
#' @param subtitle_family
#' @param subtitle_size
#' @param subtitle_face
#' @param subtitle_margin
#' @param strip_text_family
#' @param strip_text_size
#' @param strip_text_face
#' @param caption_family
#' @param caption_size
#' @param caption_face
#' @param caption_margin
#' @param axis_text_size
#' @param axis_text_family
#' @param axis_title_family
#' @param axis_title_size
#' @param axis_title_face
#' @param axis_title_just
#' @param plot_margin
#' @param panel_spacing
#' @param grid_col
#' @param grid
#' @param axis_col
#' @param axis
#' @param ticks
#'
#' @return
#' @export
#'
#' @examples
theme_idest <- function(sans_family = if (.Platform$OS.type == "windows")
  "Roboto Condensed"
  else
    "Roboto Condensed Light",
  serif_family = "Linux Libertine G",
  mono_family = "Fira Code",
  base_size = 11.5,
  plot_title_family = serif_family,
  plot_title_size = 18,
  plot_title_face = "bold",
  plot_title_margin = 10,
  subtitle_family = serif_family,
  subtitle_size = 13,
  subtitle_face = "plain",
  subtitle_margin = 15,
  subtitle_color = "grey30",
  strip_text_family = mono_family,
  strip_text_size = 13,
  strip_text_face = "plain",
  strip_back_grey = FALSE,
  caption_family = sans_family,
  caption_size = 9,
  caption_face = "plain",
  caption_margin = 10,
  axis_text_size = base_size * 0.8,
  axis_text_family = sans_family,
  axis_title_family = mono_family,
  axis_title_size = 12,
  axis_title_face = "plain",
  axis_title_just = "c",
  plot_margin = margin(12, 12, 12, 12),
  panel_spacing = grid::unit(1.2, "lines"),
  grid_col = "#cccccc",
  grid = TRUE,
  axis_col = "#cccccc",
  axis = FALSE,
  ticks = FALSE) {
  ret <- ggplot2::theme_minimal(base_family = sans_family, base_size = base_size)

  ret <- ret + theme(legend.background = element_blank())
  ret <- ret + theme(legend.key = element_blank())

  ret <- ret + theme(plot.margin = plot_margin)
  ret <- ret + theme(panel.spacing = panel_spacing)

  if (inherits(grid, "character") | grid == TRUE) {
    ret <- ret + theme(panel.grid = element_line(color = grid_col, linewidth = 0.2))
    ret <- ret + theme(panel.grid.major = element_line(color = grid_col, linewidth = 0.2))
    ret <- ret + theme(panel.grid.minor = element_line(color = grid_col, linewidth = 0.15))

    if (inherits(grid, "character")) {
      if (regexpr("X", grid)[1] < 0)
        ret <- ret + theme(panel.grid.major.x = element_blank())
      if (regexpr("Y", grid)[1] < 0)
        ret <- ret + theme(panel.grid.major.y = element_blank())
      if (regexpr("x", grid)[1] < 0)
        ret <- ret + theme(panel.grid.minor.x = element_blank())
      if (regexpr("y", grid)[1] < 0)
        ret <- ret + theme(panel.grid.minor.y = element_blank())
    }
  } else {
    ret <- ret + theme(panel.grid = element_blank())
    ret <- ret + theme(panel.grid.major  = element_blank())
    ret <- ret + theme(panel.grid.major.x  = element_blank())
    ret <- ret + theme(panel.grid.major.y  = element_blank())
    ret <- ret + theme(panel.grid.minor  = element_blank())
    ret <- ret + theme(panel.grid.minor.x  = element_blank())
    ret <- ret + theme(panel.grid.minor.y  = element_blank())
  }

  if (inherits(axis, "character") | axis == TRUE) {
    ret <- ret + theme(axis.line = element_line(color = axis_col, linewidth = 0.15))
    if (inherits(axis, "character")) {
      axis <- tolower(axis)
      if (regexpr("x", axis)[1] < 0) {
        ret <- ret + theme(axis.line.x = element_blank())
      } else {
        ret <- ret + theme(axis.line.x = element_line(color = axis_col, linewidth = 0.15))
      }
      if (regexpr("y", axis)[1] < 0) {
        ret <- ret + theme(axis.line.y = element_blank())
      } else {
        ret <- ret + theme(axis.line.y = element_line(color = axis_col, linewidth = 0.15))
      }
    } else {
      ret <- ret + theme(axis.line.x = element_line(color = axis_col, linewidth = 0.15))
      ret <- ret + theme(axis.line.y = element_line(color = axis_col, linewidth = 0.15))
    }
  } else {
    ret <- ret + theme(axis.line = element_blank())
  }

  if (!ticks) {
    ret <- ret + theme(axis.ticks = element_blank())
    ret <- ret + theme(axis.ticks.x = element_blank())
    ret <- ret + theme(axis.ticks.y = element_blank())
  } else {
    ret <- ret + theme(axis.ticks = element_line(linewidth = 0.15))
    ret <- ret + theme(axis.ticks.x = element_line(linewidth = 0.15))
    ret <- ret + theme(axis.ticks.y = element_line(linewidth = 0.15))
    ret <- ret + theme(axis.ticks.length = grid::unit(5, "pt"))
  }

  xj <- switch(
    tolower(substr(axis_title_just, 1, 1)),
    b = 0,
    l = 0,
    m = 0.5,
    c = 0.5,
    r = 1,
    t = 1
  )
  yj <- switch(
    tolower(substr(axis_title_just, 2, 2)),
    b = 0,
    l = 0,
    m = 0.5,
    c = 0.5,
    r = 1,
    t = 1
  )

  ret <- ret + theme(axis.text = element_text(
    size = axis_text_size,
    family = axis_text_family,
    margin = margin(t = 0, r = 0)
  ))
  ret <- ret + theme(axis.text.x = element_text(
    size = axis_text_size,
    family = axis_text_family,
    margin = margin(t = 0)
  ))
  ret <- ret + theme(axis.text.y = element_text(
    size = axis_text_size,
    family = axis_text_family,
    margin = margin(r = 0)
  ))

  ret <- ret + theme(axis.title = element_text(size = axis_title_size, family = axis_title_family))
  ret <- ret + theme(
    axis.title.x = element_text(
      hjust = xj,
      size = axis_title_size,
      family = axis_title_family,
      face = axis_title_face
    )
  )
  ret <- ret + theme(
    axis.title.y = element_text(
      hjust = yj,
      size = axis_title_size,
      family = axis_title_family,
      face = axis_title_face
    )
  )
  ret <- ret + theme(
    axis.title.y.right = element_text(
      hjust = yj,
      size = axis_title_size,
      angle = 90,
      family = axis_title_family,
      face = axis_title_face
    )
  )

  ret <- ret + theme(
    strip.text = element_text(
      hjust = 0,
      size = strip_text_size,
      face = strip_text_face,
      family = strip_text_family
    )
  )

  ret <- ret + theme(
    plot.title = element_text(
      hjust = 0,
      size = plot_title_size,
      margin = margin(b = plot_title_margin),
      family = plot_title_family,
      face = plot_title_face
    )
  )
  ret <- ret + theme(
    plot.subtitle = element_text(
      hjust = 0,
      size = subtitle_size,
      margin = margin(b = subtitle_margin),
      family = subtitle_family,
      face = subtitle_face,
      color = subtitle_color
    )
  )
  ret <- ret + theme(
    plot.caption = element_text(
      hjust = 1,
      size = caption_size,
      margin = margin(t = caption_margin),
      family = caption_family,
      face = caption_face
    )
  )

  if (strip_back_grey) {
    ret <- ret + theme(
      strip.background = element_rect(fill = "grey90", color = NA),
      panel.border = element_rect(color = "grey90", fill = NA)
    )
  }

  ret

}



idest_pal <- list(
  all_color_idest = list(
    c(
      "#dc863b",
      "#faefd1",
      "#2e7891",
      "#8a9da4",
      "#aa4c26",
      "#c8a734",
      "#a6d3e3",
      "#003f5f",
      "#ba63b6",
      "#f20040",
      "#774fa0",
      "#b4dfa7"
    ),
    c(1:11),
    colorblind = FALSE
  ),
  ligth_color_idest =  list(
    c("#faefd1", "#a6d3e3", "#b4dfa7", "#8a9da4"),
    c(1:4),
    colorblind = FALSE
  ),
  dark_color_idest = list(
    c(
      "#dc863b",
      "#2e7891",
      "#aa4c26",
      "#c8a734",
      "#003f5f",
      "#ba63b6",
      "#f20040",
      "#774fa0"
    ),
    c(1:8),
    colorblind = FALSE
  ),
  Picabia = list(
    c(
      "#53362e",
      "#744940",
      "#9f7064",
      "#c99582",
      "#e6bcac",
      "#e2d8d6",
      "#a5a6ae",
      "#858794",
      "#666879",
      "#515260",
      "#3d3d47"
    ),
    c(10, 4, 8, 1, 6, 3, 7, 2, 9, 5, 11),
    colorblind = TRUE
  ),
  # Copy from https://github.com/BlakeRMills/MoMAColors/blob/main/R/Functions.R
  Picasso = list(
    c(
      "#d5968c",
      "#c2676d",
      "#5c363a",
      "#995041",
      "#45939c",
      "#0f6a81"
    ),
    c(6, 3, 4, 2, 1, 5),
    colorblind = TRUE
  ),
  # Copy from https://github.com/BlakeRMills/MoMAColors/blob/main/R/Functions.R
  Levine2  = list(
    c(
      "#E3C1CB",
      "#AD5A6B",
      "#C993A2",
      "#365C83",
      "#384351",
      "#4D8F8B",
      "#CDD6AD"
    ),
    c(7, 1, 5, 3, 6, 2, 4),
    colorblind = TRUE
  ),
  # Copy from https://github.com/BlakeRMills/MoMAColors/blob/main/R/Functions.R
  Rattner = list(
    c(
      "#de8e69",
      "#f1be99",
      "#c1bd38",
      "#7a9132",
      "#4c849a",
      "#184363",
      "#5d5686",
      "#a39fc9"
    ),
    c(1, 5, 6, 2, 3, 7, 8, 4),
    colorblind = TRUE
  ),
  # Copy from https://github.com/BlakeRMills/MoMAColors/blob/main/R/Functions.R
  Sidhu = list(
    c(
      "#af4646",
      "#762b35",
      "#005187",
      "#251c4a",
      "#78adb7",
      "#4c9a77",
      "#1b7975"
    ),
    c(5, 2, 6, 7, 3, 4, 1),
    colorblind = TRUE
  ),# Copy from https://github.com/BlakeRMills/MoMAColors/blob/main/R/Functions.R
  Hokusai2 = list(c("#abc9c8", "#72aeb6", "#4692b0", "#2f70a1", "#134b73", "#0a3351"), c(5, 2, 4, 1, 6, 3), colorblind=TRUE), # copy from https://github.com/BlakeRMills/MetBrewer/blob/main/R/PaletteCode.R
  Hokusai3 = list(c("#d8d97a", "#95c36e", "#74c8c3", "#5a97c1", "#295384", "#0a2e57"), c(4, 2, 5, 3, 1, 6), colorblind=TRUE) # copy from https://github.com/BlakeRMills/MetBrewer/blob/main/R/PaletteCode.R
)

# palette_check(all_color_idest, plot = TRUE)
# palette_check(ligth_color_idest[[1]], plot = TRUE)
# palette_check(dark_color_idest, plot = TRUE)


# Copy from https://github.com/BlakeRMills/MoMAColors/blob/main/R/Functions.R
scale_color_idest_c <- function(palette_name, direction = 1, ...) {
  `%notin%` <- Negate(`%in%`)

  if (direction %notin% c(1, -1)) {
    stop("Direction not valid. Please use 1 for standard palette or -1 for reversed palette.")
  }

  scale_color_gradientn(colors = idest_colors(
    palette_name = palette_name,
    direction = direction,
    override_order = F
  ),
  ...)
}

# Copy from https://github.com/BlakeRMills/MoMAColors/blob/main/R/Functions.R
scale_fill_idest_c <- function(palette_name, direction = 1, ...) {
  `%notin%` <- Negate(`%in%`)

  if (direction %notin% c(1, -1)) {
    stop("Direction not valid. Please use 1 for standard palette or -1 for reversed palette.")
  }

  scale_color_gradientn(colors = idest_colors(
    palette_name = palette_name,
    direction = direction,
    override_order = F
  ),
  ...)
}


# Copy from https://github.com/BlakeRMills/MoMAColors/blob/main/R/Functions.R
scale_color_idest_d <- function(palette_name,
                                direction = 1,
                                override_order = FALSE,
                                ...) {
  discrete_scale(
    aesthetics = "colour",
    scale_name = "moma_d",
    palette = function(n)
      idest_colors(
        palette_name = palette_name,
        n = n,
        direction = direction,
        override_order = override_order
      ),
    ...
  )
}

# Copy from https://github.com/BlakeRMills/MoMAColors/blob/main/R/Functions.R
scale_fill_idest_d <- function(palette_name,
                               direction = 1,
                               override_order = FALSE,
                               ...) {
  discrete_scale(
    aesthetics = "fill",
    scale_name = "moma_d",
    palette = function(n)
      idest_colors(
        palette_name = palette_name,
        n = n,
        direction = direction,
        override_order = override_order
      ),
    ...
  )
}


# Copy from https://github.com/BlakeRMills/MoMAColors/blob/main/R/Functions.R
idest_colors <- function(palette_name,
                         n,
                         type = c("discrete", "continuous"),
                         direction = c(1, -1),
                         override_order = FALSE,
                         return_hex = FALSE) {
  `%notin%` <- Negate(`%in%`)

  palette <- idest_pal[[palette_name]]

  if (is.null(palette) | is.numeric(palette_name)) {
    stop("Palette does not exist.")
  }

  if (missing(n)) {
    n <- length(palette[[1]])
  }

  if (missing(direction)) {
    direction <- 1
  }

  if (direction %notin% c(1, -1)) {
    stop("Direction not valid. Please use 1 for standard palette or -1 for reversed palette.")
  }

  if (missing(type)) {
    if (n > length(palette[[1]])) {
      type <- "continuous"
    }
    else{
      type <- "discrete"
    }
  }

  type <- match.arg(type)


  if (type == "discrete" && n > length(palette[[1]])) {
    stop(
      "Number of requested colors greater than what discrete palette can offer, \n use continuous instead."
    )
  }

  continuous <-  if (direction == 1) {
    grDevices::colorRampPalette(palette[[1]])(n)
  } else{
    grDevices::colorRampPalette(rev(palette[[1]]))(n)
  }

  discrete <- if (direction == 1 & override_order == FALSE) {
    palette[[1]][which(palette[[2]] %in% c(1:n) == TRUE)]
  } else if (direction == -1 & override_order == FALSE) {
    rev(palette[[1]][which(palette[[2]] %in% c(1:n) == TRUE)])
  } else if (direction == 1 & override_order == TRUE) {
    palette[[1]][1:n]
  } else{
    rev(palette[[1]])[1:n]
  }

  out <- switch(type, continuous = continuous, discrete = discrete)
  if (return_hex == T) {
    print(out)
  }
  structure(out, class = "palette", name = palette_name)

}



#
#
# ggplot(data = mpg,
#        mapping = aes(x = displ, y = hwy, color = hwy)) +
#   geom_point(size = 3) +
#   scale_color_idest_c("Sidhu", direction = -1) +
#   facet_wrap(vars(drv)) +
#   labs(
#     x = "Displacement",
#     y = "Highway MPG",
#     color = "Car class",
#     title = "Heavier cars get worse mileage",
#     subtitle = "Except two-seaters?",
#     caption = "Here's a caption"
#   ) +
#   theme_idest()
