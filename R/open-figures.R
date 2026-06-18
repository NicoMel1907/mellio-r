# Internal figure handoff helpers for mellio_open().
#
# These functions keep plots as a web-app handoff concern. They intentionally
# do not expose a package-level figure class or a separate in-R figure API.

mellio_open_figure_object <- function(plot = NULL, image_raw = NULL,
                                      image_path = NULL,
                                      image_format = "png",
                                      style = "apa7", title = NULL,
                                      number = NULL, note = NULL,
                                      source = NULL, width = NULL,
                                      height = NULL, dpi = NULL,
                                      width_px = NULL, height_px = NULL) {
  style <- match.arg(tolower(style), list_styles())

  list(
    plot = plot,
    image_raw = image_raw,
    image_path = image_path,
    image_format = image_format,
    width_in = width,
    height_in = height,
    dpi = if (is.null(dpi)) NULL else as.integer(dpi),
    width_px = if (is.null(width_px)) NULL else as.integer(width_px),
    height_px = if (is.null(height_px)) NULL else as.integer(height_px),
    style = style,
    title = title,
    number = number,
    note = note,
    source = source,
    provenance = ms_provenance_basic()
  )
}

mellio_open_figure_from_gg <- function(x, style = "apa7", title = NULL,
                                       number = NULL, note = NULL,
                                       source = NULL, width = 6,
                                       height = 4, dpi = 300, ...) {
  rlang::check_installed("ggplot2", reason = "to send ggplot objects to Mellio")

  mellio_open_figure_object(
    plot = x,
    image_format = "png",
    style = style,
    title = title,
    number = number,
    note = note,
    source = source,
    width = width,
    height = height,
    dpi = dpi,
    width_px = width * dpi,
    height_px = height * dpi
  )
}

mellio_open_figure_from_png <- function(path, plot = NULL, style = "apa7",
                                        title = NULL, number = NULL,
                                        note = NULL, source = NULL,
                                        width = NULL, height = NULL,
                                        dpi = NULL, width_px = NULL,
                                        height_px = NULL) {
  mellio_open_figure_object(
    plot = plot,
    image_raw = readBin(path, "raw", file.info(path)$size),
    image_path = normalizePath(path, mustWork = FALSE),
    image_format = "png",
    style = style,
    title = title,
    number = number,
    note = note,
    source = source,
    width = width,
    height = height,
    dpi = dpi,
    width_px = width_px,
    height_px = height_px
  )
}

mellio_open_figure_from_printed <- function(x, style = "apa7", title = NULL,
                                            number = NULL, note = NULL,
                                            source = NULL, width = 6,
                                            height = 4, dpi = 300, ...) {
  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)

  grDevices::png(tmp, width = width, height = height, units = "in",
                 res = dpi, bg = "white")
  dev_id <- as.integer(grDevices::dev.cur())
  dev_closed <- FALSE
  on.exit({
    open_devices <- as.integer(grDevices::dev.list() %||% integer())
    if (!dev_closed && dev_id %in% open_devices) {
      grDevices::dev.off(dev_id)
    }
  }, add = TRUE)

  print(x)
  grDevices::dev.off(dev_id)
  dev_closed <- TRUE

  mellio_open_figure_from_png(
    tmp,
    plot = NULL,
    style = style,
    title = title,
    number = number,
    note = note,
    source = source,
    width = width,
    height = height,
    dpi = dpi,
    width_px = width * dpi,
    height_px = height * dpi
  )
}

mellio_open_figure_from_htmlwidget <- function(x, style = "apa7",
                                               title = NULL,
                                               number = NULL,
                                               note = NULL,
                                               source = NULL,
                                               dpi = 150,
                                               vwidth = 1200,
                                               vheight = 800,
                                               selfcontained = TRUE,
                                               ...) {
  rlang::check_installed("htmlwidgets",
    reason = "to save htmlwidgets before sending them to Mellio")
  rlang::check_installed("webshot2",
    reason = "to render htmlwidgets as static images for Mellio")

  tmp_html <- tempfile(fileext = ".html")
  tmp_png <- tempfile(fileext = ".png")
  html_dir <- paste0(tools::file_path_sans_ext(tmp_html), "_files")
  on.exit(unlink(c(tmp_html, tmp_png, html_dir), recursive = TRUE, force = TRUE),
          add = TRUE)

  htmlwidgets::saveWidget(x, file = tmp_html, selfcontained = selfcontained)
  webshot2::webshot(tmp_html, file = tmp_png, vwidth = vwidth,
                    vheight = vheight, ...)

  mellio_open_figure_from_png(
    tmp_png,
    plot = NULL,
    style = style,
    title = title,
    number = number,
    note = note,
    source = source,
    width = vwidth / dpi,
    height = vheight / dpi,
    dpi = dpi,
    width_px = vwidth,
    height_px = vheight
  )
}

mellio_open_figure_from_file <- function(x, style = "apa7", title = NULL,
                                         number = NULL, note = NULL,
                                         source = NULL, ...) {
  if (!file.exists(x)) {
    cli::cli_abort("File not found: {.file {x}}")
  }

  ext <- tolower(tools::file_ext(x))
  supported <- c("png", "jpg", "jpeg", "tiff", "tif", "bmp", "pdf", "svg")
  if (!ext %in% supported) {
    cli::cli_abort(c(
      "Unsupported image format: {.val {ext}}",
      "i" = "Supported formats: {.val {supported}}"
    ))
  }

  fmt <- ext
  if (fmt %in% c("jpg", "jpeg")) fmt <- "jpeg"
  if (fmt %in% c("tiff", "tif")) fmt <- "tiff"

  mellio_open_figure_object(
    image_raw = readBin(x, "raw", file.info(x)$size),
    image_path = normalizePath(x),
    image_format = fmt,
    style = style,
    title = title,
    number = number,
    note = note,
    source = source
  )
}

mellio_open_figure_from_recordedplot <- function(x, style = "apa7",
                                                 title = NULL,
                                                 number = NULL,
                                                 note = NULL,
                                                 source = NULL,
                                                 width = 6,
                                                 height = 4,
                                                 dpi = 300, ...) {
  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)

  grDevices::png(tmp, width = width, height = height, units = "in", res = dpi,
                 bg = "white")
  dev_id <- as.integer(grDevices::dev.cur())
  dev_closed <- FALSE
  on.exit({
    open_devices <- as.integer(grDevices::dev.list() %||% integer())
    if (!dev_closed && dev_id %in% open_devices) {
      grDevices::dev.off(dev_id)
    }
  }, add = TRUE)

  grDevices::replayPlot(x)
  grDevices::dev.off(dev_id)
  dev_closed <- TRUE

  size <- file.info(tmp)$size
  if (!file.exists(tmp) || is.na(size) || size == 0L) {
    cli::cli_abort(c(
      "Could not replay the recorded plot.",
      "i" = "Some non-interactive graphics devices do not keep a display list.",
      "i" = "In scripts, call {.code grDevices::dev.control(\"enable\")} before plotting, or pass a saved image file to {.fn mellio_open}."
    ))
  }

  mellio_open_figure_from_png(
    tmp,
    plot = NULL,
    style = style,
    title = title,
    number = number,
    note = note,
    source = source,
    width = width,
    height = height,
    dpi = dpi,
    width_px = width * dpi,
    height_px = height * dpi
  )
}
