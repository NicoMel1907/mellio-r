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
                                      width_px = NULL, height_px = NULL,
                                      recipe_payload = NULL,
                                      recipe_type = NULL) {
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
    recipe_payload = recipe_payload,
    recipe_type = recipe_type,
    provenance = ms_provenance_basic()
  )
}

mellio_open_figure_from_gg <- function(x, style = "apa7", title = NULL,
                                       number = NULL, note = NULL,
                                       source = NULL, width = 6,
                                       height = 4, dpi = 300, ...,
                                       .call = NULL) {
  rlang::check_installed("ggplot2", reason = "to send ggplot objects to Mellio")

  recipe_payload <- mellio_ggplot_raw_scatter_payload(
    x,
    title = title,
    .call = .call
  )

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
    height_px = height * dpi,
    recipe_payload = recipe_payload,
    recipe_type = if (is.null(recipe_payload)) NULL else "raw_scatter"
  )
}

mellio_ggplot_raw_scatter_payload <- function(plot, title = NULL,
                                              .call = NULL,
                                              max_points = 1000L,
                                              max_groups = 20L) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("rlang", quietly = TRUE)) {
    return(NULL)
  }

  layers <- plot$layers %||% list()
  if (length(layers) != 1L) return(NULL)

  layer <- layers[[1]]
  if (!inherits(layer$geom, "GeomPoint") ||
      !inherits(layer$stat, "StatIdentity")) {
    return(NULL)
  }

  data <- mellio_ggplot_layer_data(plot, layer)
  if (!is.data.frame(data) || nrow(data) < 4L) return(NULL)

  mapping <- mellio_ggplot_layer_mapping(plot, layer)
  if (is.null(mapping$x) || is.null(mapping$y)) return(NULL)

  x_raw <- mellio_ggplot_eval_mapping(mapping$x, data)
  y_raw <- mellio_ggplot_eval_mapping(mapping$y, data)
  if (!is.numeric(x_raw) || !is.numeric(y_raw)) return(NULL)

  x_vals <- as.numeric(x_raw)
  y_vals <- as.numeric(y_raw)
  if (length(x_vals) != nrow(data) || length(y_vals) != nrow(data)) return(NULL)

  built <- tryCatch(ggplot2::ggplot_build(plot), error = function(e) NULL)
  if (is.null(built) || length(built$data) != 1L) return(NULL)
  panel <- built$data[[1]]$PANEL
  if (!is.null(panel) && length(unique(as.character(panel))) > 1L) return(NULL)

  valid <- is.finite(x_vals) & is.finite(y_vals)
  if (sum(valid) < 4L || sum(valid) > max_points) return(NULL)

  group_info <- mellio_ggplot_group_values(plot, mapping, data,
                                           valid = valid,
                                           max_groups = max_groups)
  if (identical(group_info, FALSE)) return(NULL)

  idx <- which(valid)
  points <- lapply(seq_along(idx), function(i) {
    row <- idx[[i]]
    point <- list(
      id = as.integer(i),
      row = as.integer(row),
      x = x_vals[[row]],
      y = y_vals[[row]]
    )
    if (!is.null(group_info)) {
      point$group <- group_info$values[[i]]
    }
    point
  })

  x_var <- mellio_ggplot_mapping_label(mapping$x)
  y_var <- mellio_ggplot_mapping_label(mapping$y)
  x_label <- mellio_ggplot_plot_label(plot, "x", x_var)
  y_label <- mellio_ggplot_plot_label(plot, "y", y_var)

  raw_scatter <- list(
    source = "r_ggplot",
    x = list(variable = x_var, label = x_label),
    y = list(variable = y_var, label = y_label),
    points = points
  )
  fields <- list(table_type = "raw_scatter", x = x_var, y = y_var)

  if (!is.null(group_info)) {
    raw_scatter$group <- list(
      variable = group_info$variable,
      label = group_info$label,
      levels = group_info$levels
    )
    fields$group <- group_info$variable
  }

  payload_title <- title %||%
    mellio_ggplot_label_value(plot, "title") %||%
    paste0("Scatterplot of ", y_label, " by ", x_label)

  list(
    schema_version = "0.1",
    result_id = paste0("r_ggplot_", format(Sys.time(), "%Y%m%d%H%M%OS6",
                                           tz = "UTC")),
    card_kind = "table",
    type = "r_ggplot_raw_scatter",
    type_label = "ggplot scatterplot",
    name = payload_title,
    call = .call %||% "ggplot object",
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    fields = fields,
    figure_data = list(raw_scatter = raw_scatter),
    metadata = list(
      source = "r_ggplot",
      available_figures = list(list(
        type = "raw_scatter",
        label = "Scatterplot",
        default = TRUE
      ))
    ),
    packages = ms_packages_basic("ggplot2")
  )
}

mellio_ggplot_layer_data <- function(plot, layer) {
  data <- layer$data
  if (is.null(data) || inherits(data, "waiver")) {
    data <- plot$data
  }
  if (is.function(data)) return(NULL)
  if (!is.data.frame(data)) return(NULL)
  data
}

mellio_ggplot_layer_mapping <- function(plot, layer) {
  mapping <- if (isTRUE(layer$inherit.aes)) {
    plot$mapping %||% ggplot2::aes()
  } else {
    ggplot2::aes()
  }
  layer_mapping <- layer$mapping %||% ggplot2::aes()
  for (nm in names(layer_mapping)) {
    mapping[[nm]] <- layer_mapping[[nm]]
  }
  mapping
}

mellio_ggplot_eval_mapping <- function(expr, data) {
  value <- tryCatch(
    rlang::eval_tidy(expr, data = data),
    error = function(e) NULL
  )
  if (length(value) == 1L && nrow(data) > 1L) {
    value <- rep(value, nrow(data))
  }
  value
}

mellio_ggplot_group_values <- function(plot, mapping, data, valid,
                                       max_groups = 20L) {
  group_aes <- NULL
  for (candidate in c("colour", "color", "fill", "group")) {
    if (!is.null(mapping[[candidate]])) {
      group_aes <- candidate
      break
    }
  }
  if (is.null(group_aes)) return(NULL)

  raw <- mellio_ggplot_eval_mapping(mapping[[group_aes]], data)
  if (is.null(raw) || length(raw) != nrow(data)) return(NULL)

  text <- as.character(raw)
  text[is.na(text) | !nzchar(text)] <- "Missing"
  kept <- text[valid]
  if (length(unique(kept)) <= 1L) return(NULL)

  if (is.factor(raw)) {
    level_values <- intersect(levels(raw), unique(kept))
    level_values <- c(level_values, setdiff(unique(kept), level_values))
  } else {
    level_values <- unique(kept)
  }
  if (length(level_values) > max_groups) return(FALSE)

  variable <- mellio_ggplot_mapping_label(mapping[[group_aes]])
  label <- mellio_ggplot_plot_label(plot, group_aes, variable)
  list(
    variable = variable,
    label = label,
    values = kept,
    levels = lapply(level_values, function(value) {
      list(value = value, label = value)
    })
  )
}

mellio_ggplot_mapping_label <- function(expr) {
  label <- tryCatch(rlang::as_label(expr), error = function(e) NULL)
  label <- label %||% ""
  label <- sub("^~", "", label)
  if (!nzchar(label)) "value" else label
}

mellio_ggplot_plot_label <- function(plot, aes, fallback) {
  value <- mellio_ggplot_label_value(plot, aes)
  if (is.null(value) && identical(aes, "colour")) {
    value <- mellio_ggplot_label_value(plot, "color")
  }
  if (is.null(value) && identical(aes, "color")) {
    value <- mellio_ggplot_label_value(plot, "colour")
  }
  value %||% fallback
}

mellio_ggplot_label_value <- function(plot, key) {
  labels <- plot$labels %||% list()
  value <- labels[[key]]
  if (is.null(value) || length(value) != 1L) return(NULL)
  value <- tryCatch(as.character(value), error = function(e) NULL)
  if (is.null(value) || is.na(value) || !nzchar(value)) return(NULL)
  value
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
