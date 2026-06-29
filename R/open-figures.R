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

  recipe <- mellio_ggplot_recipe(
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
    recipe_payload = recipe$payload,
    recipe_type = recipe$type
  )
}

mellio_ggplot_recipe <- function(plot, title = NULL, .call = NULL) {
  recipes <- list(
    scatter_plot = mellio_ggplot_lm_scatter_payload,
    adjusted_means = mellio_ggplot_summary_errorbar_payload,
    box_plot = mellio_ggplot_boxplot_payload,
    bar_chart = mellio_ggplot_bar_chart_payload,
    raw_scatter = mellio_ggplot_raw_scatter_payload
  )

  for (type in names(recipes)) {
    payload <- recipes[[type]](plot, title = title, .call = .call)
    if (!is.null(payload)) {
      figure_type <- type
      available <- payload$metadata$available_figures
      if (length(available) && !is.null(available[[1]]$type)) {
        figure_type <- available[[1]]$type
      }
      return(list(type = figure_type, payload = payload))
    }
  }

  list(type = NULL, payload = NULL)
}

mellio_ggplot_raw_scatter_payload <- function(plot, title = NULL,
                                              .call = NULL,
                                              max_points = 1500L,
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

  built <- tryCatch(
    suppressMessages(ggplot2::ggplot_build(plot)),
    error = function(e) NULL
  )
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

mellio_ggplot_bar_chart_payload <- function(plot, title = NULL, .call = NULL,
                                            max_categories = 60L) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("rlang", quietly = TRUE)) {
    return(NULL)
  }

  layers <- plot$layers %||% list()
  if (length(layers) != 1L) return(NULL)
  layer <- layers[[1]]
  if (!inherits(layer$geom, "GeomBar") &&
      !inherits(layer$geom, "GeomCol")) {
    return(NULL)
  }

  data <- mellio_ggplot_layer_data(plot, layer)
  if (!is.data.frame(data) || nrow(data) < 1L) return(NULL)
  mapping <- mellio_ggplot_layer_mapping(plot, layer)
  if (is.null(mapping$x)) return(NULL)
  if (!is.null(mapping$fill) || !is.null(mapping$colour) ||
      !is.null(mapping$color) || !is.null(mapping$group)) {
    return(NULL)
  }

  built <- tryCatch(
    suppressMessages(ggplot2::ggplot_build(plot)),
    error = function(e) NULL
  )
  if (is.null(built) || length(built$data) != 1L) return(NULL)
  panel <- built$data[[1]]$PANEL
  if (!is.null(panel) && length(unique(as.character(panel))) > 1L) return(NULL)

  x_raw <- mellio_ggplot_eval_mapping(mapping$x, data)
  if (is.null(x_raw) || length(x_raw) != nrow(data)) return(NULL)
  x_text <- as.character(x_raw)
  valid_x <- !is.na(x_text) & nzchar(x_text)
  if (!any(valid_x)) return(NULL)

  stat_count <- inherits(layer$stat, "StatCount")
  stat_identity <- inherits(layer$stat, "StatIdentity")
  if (!stat_count && !stat_identity) return(NULL)

  x_var <- mellio_ggplot_mapping_label(mapping$x)
  x_label <- mellio_ggplot_plot_label(plot, "x", x_var)
  y_label <- mellio_ggplot_plot_label(
    plot,
    "y",
    if (stat_count) "Count" else mellio_ggplot_mapping_label(mapping$y)
  )

  if (stat_count) {
    categories <- mellio_ggplot_count_categories(x_raw[valid_x])
    caption <- "Counts per category."
  } else {
    if (is.null(mapping$y)) return(NULL)
    y_raw <- mellio_ggplot_eval_mapping(mapping$y, data)
    if (!is.numeric(y_raw) || length(y_raw) != nrow(data)) return(NULL)
    y_vals <- as.numeric(y_raw)
    valid <- valid_x & is.finite(y_vals)
    if (!any(valid) || any(y_vals[valid] < 0)) return(NULL)
    categories <- mellio_ggplot_value_categories(x_raw[valid], y_vals[valid])
    caption <- "Values by category."
  }

  if (!length(categories) || length(categories) > max_categories) return(NULL)

  payload_title <- title %||%
    mellio_ggplot_label_value(plot, "title") %||%
    if (stat_count) paste0("Frequencies of ", x_label) else paste0(y_label, " by ", x_label)

  metadata <- list(
    source = "r_ggplot",
    available_figures = list(list(
      type = "bar_chart",
      label = "Bar chart",
      default = TRUE
    ))
  )
  if (inherits(plot$coordinates, "CoordFlip")) {
    metadata$figure_settings <- list(barOrientation = "horizontal")
  }

  list(
    schema_version = "0.1",
    result_id = paste0("r_ggplot_", format(Sys.time(), "%Y%m%d%H%M%OS6",
                                           tz = "UTC")),
    card_kind = "table",
    type = if (stat_count) "r_ggplot_bar_count" else "r_ggplot_bar_value",
    type_label = if (stat_count) "ggplot count bar chart" else "ggplot bar chart",
    name = payload_title,
    call = .call %||% "ggplot object",
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    fields = list(table_type = "bar_chart", x = x_var, y = y_label),
    figure_data = list(
      bar_chart = list(
        source = "r_ggplot",
        variable = x_label,
        y_label = y_label,
        title = payload_title,
        caption = caption,
        categories = categories
      )
    ),
    metadata = metadata,
    packages = ms_packages_basic("ggplot2")
  )
}

mellio_ggplot_summary_errorbar_payload <- function(plot, title = NULL,
                                                   .call = NULL,
                                                   max_groups = 60L) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("rlang", quietly = TRUE)) {
    return(NULL)
  }

  layers <- plot$layers %||% list()
  if (inherits(plot$coordinates, "CoordFlip")) return(NULL)

  point_layer <- NULL
  range_layer <- NULL
  if (length(layers) == 1L &&
      inherits(layers[[1]]$geom, "GeomPointrange") &&
      inherits(layers[[1]]$stat, "StatIdentity")) {
    point_layer <- layers[[1]]
    range_layer <- layers[[1]]
  } else if (length(layers) == 2L) {
    point_idx <- which(vapply(layers, function(layer) {
      inherits(layer$geom, "GeomPoint") && inherits(layer$stat, "StatIdentity")
    }, logical(1)))
    range_idx <- which(vapply(layers, function(layer) {
      inherits(layer$geom, "GeomErrorbar") && inherits(layer$stat, "StatIdentity")
    }, logical(1)))
    if (length(point_idx) != 1L || length(range_idx) != 1L) return(NULL)
    point_layer <- layers[[point_idx]]
    range_layer <- layers[[range_idx]]
  } else {
    return(NULL)
  }

  point_data <- mellio_ggplot_layer_data(plot, point_layer)
  range_data <- mellio_ggplot_layer_data(plot, range_layer)
  if (!is.data.frame(point_data) || !is.data.frame(range_data) ||
      nrow(point_data) < 2L || nrow(point_data) != nrow(range_data)) {
    return(NULL)
  }

  point_mapping <- mellio_ggplot_layer_mapping(plot, point_layer)
  range_mapping <- mellio_ggplot_layer_mapping(plot, range_layer)
  if (is.null(point_mapping$x) || is.null(point_mapping$y) ||
      is.null(range_mapping$ymin) || is.null(range_mapping$ymax)) {
    return(NULL)
  }

  x_raw <- mellio_ggplot_eval_mapping(point_mapping$x, point_data)
  y_raw <- mellio_ggplot_eval_mapping(point_mapping$y, point_data)
  ymin_raw <- mellio_ggplot_eval_mapping(range_mapping$ymin, range_data)
  ymax_raw <- mellio_ggplot_eval_mapping(range_mapping$ymax, range_data)
  if (is.null(x_raw) || length(x_raw) != nrow(point_data) ||
      !is.numeric(y_raw) || length(y_raw) != nrow(point_data) ||
      !is.numeric(ymin_raw) || length(ymin_raw) != nrow(range_data) ||
      !is.numeric(ymax_raw) || length(ymax_raw) != nrow(range_data)) {
    return(NULL)
  }

  built <- tryCatch(
    suppressMessages(ggplot2::ggplot_build(plot)),
    error = function(e) NULL
  )
  if (is.null(built) || length(built$data) != length(layers)) return(NULL)
  panels <- unlist(lapply(built$data, function(layer_data) {
    as.character(layer_data$PANEL %||% 1)
  }), use.names = FALSE)
  if (length(unique(panels)) > 1L) return(NULL)

  valid <- !is.na(as.character(x_raw)) &
    nzchar(as.character(x_raw)) &
    is.finite(as.numeric(y_raw)) &
    is.finite(as.numeric(ymin_raw)) &
    is.finite(as.numeric(ymax_raw))
  if (sum(valid) < 2L) return(NULL)

  x_var <- mellio_ggplot_mapping_label(point_mapping$x)
  y_var <- mellio_ggplot_mapping_label(point_mapping$y)
  if (!mellio_ggplot_group_aes_safe(point_mapping, "fill", point_data, valid, x_var) ||
      !mellio_ggplot_group_aes_safe(point_mapping, "colour", point_data, valid, x_var) ||
      !mellio_ggplot_group_aes_safe(point_mapping, "color", point_data, valid, x_var) ||
      !mellio_ggplot_group_aes_safe(point_mapping, "group", point_data, valid, x_var) ||
      !mellio_ggplot_group_aes_safe(range_mapping, "fill", range_data, valid, x_var) ||
      !mellio_ggplot_group_aes_safe(range_mapping, "colour", range_data, valid, x_var) ||
      !mellio_ggplot_group_aes_safe(range_mapping, "color", range_data, valid, x_var) ||
      !mellio_ggplot_group_aes_safe(range_mapping, "group", range_data, valid, x_var)) {
    return(NULL)
  }

  x_text <- as.character(x_raw)
  levels <- mellio_ggplot_category_order(x_raw[valid])
  if (length(levels) < 2L || length(levels) > max_groups) return(NULL)
  if (anyDuplicated(x_text[valid])) return(NULL)

  y_vals <- as.numeric(y_raw)
  ymin_vals <- as.numeric(ymin_raw)
  ymax_vals <- as.numeric(ymax_raw)
  x_label <- mellio_ggplot_plot_label(plot, "x", x_var)
  y_label <- mellio_ggplot_plot_label(plot, "y", y_var)
  group_rows <- lapply(levels, function(level) {
    rows <- which(valid & x_text == level)
    if (length(rows) != 1L) return(NULL)
    row <- rows[[1]]
    lo <- min(ymin_vals[[row]], ymax_vals[[row]])
    hi <- max(ymin_vals[[row]], ymax_vals[[row]])
    list(
      level = level,
      label = level,
      mean = y_vals[[row]],
      ci_lower = lo,
      ci_upper = hi,
      x = match(level, levels)
    )
  })
  group_rows <- Filter(Negate(is.null), group_rows)
  if (length(group_rows) != length(levels)) return(NULL)

  payload_title <- title %||%
    mellio_ggplot_label_value(plot, "title") %||%
    paste0(y_label, " by ", x_label)

  list(
    schema_version = "0.1",
    result_id = paste0("r_ggplot_", format(Sys.time(), "%Y%m%d%H%M%OS6",
                                           tz = "UTC")),
    card_kind = "table",
    type = "r_ggplot_summary_errorbar",
    type_label = "ggplot summary error-bar plot",
    name = payload_title,
    call = .call %||% "ggplot object",
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    fields = list(table_type = "summary_errorbar", factor = x_var, y = y_var),
    figure_data = list(
      adjusted_means = list(
        mean_kind = "observed",
        source = "r_ggplot_summary",
        factor = list(
          variable = x_var,
          term = x_var,
          label = x_label,
          levels = lapply(levels, function(level) {
            list(value = level, label = level)
          })
        ),
        groups = group_rows,
        y_label = y_label,
        outcome = y_label,
        ci_level = 0.95,
        ci_method = "source_bounds"
      )
    ),
    metadata = list(
      source = "r_ggplot",
      available_figures = list(list(
        type = "adjusted_means",
        label = "Error-bar plot",
        default = TRUE
      ))
    ),
    packages = ms_packages_basic("ggplot2")
  )
}

mellio_ggplot_boxplot_payload <- function(plot, title = NULL, .call = NULL,
                                          max_points = 1500L,
                                          max_groups = 40L) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("rlang", quietly = TRUE)) {
    return(NULL)
  }

  layers <- plot$layers %||% list()
  if (length(layers) != 1L) return(NULL)
  if (inherits(plot$coordinates, "CoordFlip")) return(NULL)

  layer <- layers[[1]]
  if (!inherits(layer$geom, "GeomBoxplot") ||
      !inherits(layer$stat, "StatBoxplot")) {
    return(NULL)
  }

  data <- mellio_ggplot_layer_data(plot, layer)
  if (!is.data.frame(data) || nrow(data) < 2L) return(NULL)

  mapping <- mellio_ggplot_layer_mapping(plot, layer)
  if (is.null(mapping$y)) return(NULL)

  y_raw <- mellio_ggplot_eval_mapping(mapping$y, data)
  if (!is.numeric(y_raw) || length(y_raw) != nrow(data)) return(NULL)
  y_vals <- as.numeric(y_raw)

  built <- tryCatch(
    suppressMessages(ggplot2::ggplot_build(plot)),
    error = function(e) NULL
  )
  if (is.null(built) || length(built$data) != 1L) return(NULL)
  panel <- built$data[[1]]$PANEL
  if (!is.null(panel) && length(unique(as.character(panel))) > 1L) return(NULL)

  valid_y <- is.finite(y_vals)
  if (sum(valid_y) < 2L || sum(valid_y) > max_points) return(NULL)

  x_var <- NULL
  x_label <- NULL
  x_raw <- NULL
  x_text <- NULL
  if (!is.null(mapping$x)) {
    x_raw <- mellio_ggplot_eval_mapping(mapping$x, data)
    if (is.null(x_raw) || length(x_raw) != nrow(data)) return(NULL)
    x_text <- as.character(x_raw)
    x_text[is.na(x_text) | !nzchar(x_text)] <- NA_character_
    x_var <- mellio_ggplot_mapping_label(mapping$x)
    x_label <- mellio_ggplot_plot_label(plot, "x", x_var)
  }

  if (!mellio_ggplot_group_aes_safe(mapping, "fill", data, valid_y, x_var) ||
      !mellio_ggplot_group_aes_safe(mapping, "colour", data, valid_y, x_var) ||
      !mellio_ggplot_group_aes_safe(mapping, "color", data, valid_y, x_var) ||
      !mellio_ggplot_group_aes_safe(mapping, "group", data, valid_y, x_var)) {
    return(NULL)
  }

  y_var <- mellio_ggplot_mapping_label(mapping$y)
  y_label <- mellio_ggplot_plot_label(plot, "y", y_var)
  payload_title <- title %||%
    mellio_ggplot_label_value(plot, "title") %||%
    if (!is.null(x_label)) {
      paste0("Distribution of ", y_label, " by ", x_label)
    } else {
      paste0("Distribution of ", y_label)
    }

  if (is.null(x_text)) {
    values <- y_vals[valid_y]
    box <- mellio_ggplot_box_stats(values)
    if (is.null(box)) return(NULL)
    return(list(
      schema_version = "0.1",
      result_id = paste0("r_ggplot_", format(Sys.time(), "%Y%m%d%H%M%OS6",
                                             tz = "UTC")),
      card_kind = "table",
      type = "r_ggplot_boxplot",
      type_label = "ggplot box plot",
      name = payload_title,
      call = .call %||% "ggplot object",
      created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
      fields = list(table_type = "distribution", y = y_var),
      figure_data = list(
        distribution = list(
          source = "r_ggplot_boxplot",
          variable = y_label,
          n = length(values),
          values = values,
          box = box
        )
      ),
      metadata = list(
        source = "r_ggplot",
        available_figures = list(list(
          type = "box_plot",
          label = "Box plot",
          default = TRUE
        ))
      ),
      packages = ms_packages_basic("ggplot2")
    ))
  }

  valid <- valid_y & !is.na(x_text)
  if (sum(valid) < 2L) return(NULL)
  levels <- mellio_ggplot_category_order(x_raw[valid])
  if (length(levels) < 2L || length(levels) > max_groups) return(NULL)

  groups <- lapply(levels, function(level) {
    rows <- which(valid & x_text == level)
    values <- y_vals[rows]
    box <- mellio_ggplot_box_stats(values)
    if (is.null(box)) return(NULL)
    list(
      level = level,
      label = level,
      n = length(values),
      median = box$median,
      q1 = box$q1,
      q3 = box$q3,
      iqr = box$iqr,
      min = box$min,
      max = box$max,
      whisker_lo = box$whisker_lo,
      whisker_hi = box$whisker_hi,
      mean = mean(values)
    )
  })
  groups <- Filter(Negate(is.null), groups)
  if (length(groups) < 2L) return(NULL)

  idx <- which(valid)
  observations <- lapply(seq_along(idx), function(i) {
    row <- idx[[i]]
    list(
      id = as.integer(i),
      row = as.integer(row),
      group = x_text[[row]],
      value = y_vals[[row]]
    )
  })

  list(
    schema_version = "0.1",
    result_id = paste0("r_ggplot_", format(Sys.time(), "%Y%m%d%H%M%OS6",
                                           tz = "UTC")),
    card_kind = "table",
    type = "r_ggplot_grouped_boxplot",
    type_label = "ggplot grouped box plot",
    name = payload_title,
    call = .call %||% "ggplot object",
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    fields = list(table_type = "grouped_distribution", factor = x_var, y = y_var),
    figure_data = list(
      nonparametric_group_plot = list(
        source = "r_ggplot_boxplot",
        factor = list(
          variable = x_var,
          term = x_var,
          label = x_label,
          levels = lapply(levels, function(level) {
            list(value = level, label = level)
          })
        ),
        groups = groups,
        observations = observations,
        y_label = y_label,
        outcome = y_label,
        point_display = list(
          total_n = length(observations),
          included_n = length(observations),
          truncated = FALSE
        )
      )
    ),
    metadata = list(
      source = "r_ggplot",
      available_figures = list(list(
        type = "nonparametric_group_plot",
        label = "Box plot",
        default = TRUE
      ))
    ),
    packages = ms_packages_basic("ggplot2")
  )
}

mellio_ggplot_lm_scatter_payload <- function(plot, title = NULL, .call = NULL,
                                             max_points = 1500L) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("rlang", quietly = TRUE)) {
    return(NULL)
  }

  layers <- plot$layers %||% list()
  if (length(layers) != 2L) return(NULL)
  point_idx <- which(vapply(layers, function(layer) {
    inherits(layer$geom, "GeomPoint") && inherits(layer$stat, "StatIdentity")
  }, logical(1)))
  smooth_idx <- which(vapply(layers, function(layer) {
    inherits(layer$geom, "GeomSmooth") && inherits(layer$stat, "StatSmooth")
  }, logical(1)))
  if (length(point_idx) != 1L || length(smooth_idx) != 1L) return(NULL)

  smooth_layer <- layers[[smooth_idx]]
  method <- smooth_layer$stat_params$method
  method_is_lm <- identical(method, "lm") ||
    (is.function(method) && identical(method, stats::lm))
  if (!isTRUE(method_is_lm)) return(NULL)

  point_layer <- layers[[point_idx]]
  point_data <- mellio_ggplot_layer_data(plot, point_layer)
  if (!is.data.frame(point_data) || nrow(point_data) < 3L) return(NULL)
  point_mapping <- mellio_ggplot_layer_mapping(plot, point_layer)
  smooth_mapping <- mellio_ggplot_layer_mapping(plot, smooth_layer)
  if (is.null(point_mapping$x) || is.null(point_mapping$y)) return(NULL)
  if (!is.null(point_mapping$colour) || !is.null(point_mapping$color) ||
      !is.null(point_mapping$fill) || !is.null(point_mapping$group) ||
      !is.null(smooth_mapping$colour) || !is.null(smooth_mapping$color) ||
      !is.null(smooth_mapping$fill) || !is.null(smooth_mapping$group)) {
    return(NULL)
  }

  x_raw <- mellio_ggplot_eval_mapping(point_mapping$x, point_data)
  y_raw <- mellio_ggplot_eval_mapping(point_mapping$y, point_data)
  if (!is.numeric(x_raw) || !is.numeric(y_raw)) return(NULL)
  x_vals <- as.numeric(x_raw)
  y_vals <- as.numeric(y_raw)
  if (length(x_vals) != nrow(point_data) || length(y_vals) != nrow(point_data)) return(NULL)
  valid <- is.finite(x_vals) & is.finite(y_vals)
  if (sum(valid) < 3L || sum(valid) > max_points) return(NULL)

  built <- tryCatch(
    suppressMessages(ggplot2::ggplot_build(plot)),
    error = function(e) NULL
  )
  if (is.null(built) || length(built$data) != 2L) return(NULL)
  smooth_built <- built$data[[smooth_idx]]
  point_built <- built$data[[point_idx]]
  if (length(unique(as.character(point_built$PANEL))) > 1L ||
      length(unique(as.character(smooth_built$PANEL))) > 1L ||
      length(unique(as.character(smooth_built$group))) > 1L) {
    return(NULL)
  }

  line <- mellio_ggplot_smooth_line(smooth_built)
  if (length(line) < 2L) return(NULL)

  idx <- which(valid)
  observations <- lapply(seq_along(idx), function(i) {
    row <- idx[[i]]
    list(
      id = as.integer(i),
      row = as.integer(row),
      x = x_vals[[row]],
      y = y_vals[[row]]
    )
  })

  x_var <- mellio_ggplot_mapping_label(point_mapping$x)
  y_var <- mellio_ggplot_mapping_label(point_mapping$y)
  x_label <- mellio_ggplot_plot_label(plot, "x", x_var)
  y_label <- mellio_ggplot_plot_label(plot, "y", y_var)
  conf_level <- suppressWarnings(as.numeric(smooth_layer$stat_params$level))
  if (length(conf_level) != 1L || !isTRUE(is.finite(conf_level)) ||
      is.na(conf_level) || conf_level <= 0 || conf_level >= 1) {
    conf_level <- 0.95
  }

  payload_title <- title %||%
    mellio_ggplot_label_value(plot, "title") %||%
    paste0("Scatter plot of ", y_label, " and ", x_label)

  list(
    schema_version = "0.1",
    result_id = paste0("r_ggplot_", format(Sys.time(), "%Y%m%d%H%M%OS6",
                                           tz = "UTC")),
    card_kind = "table",
    type = "r_ggplot_lm_scatter",
    type_label = "ggplot regression scatterplot",
    name = payload_title,
    call = .call %||% "ggplot object",
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    fields = list(table_type = "scatter_plot", x = x_var, y = y_var),
    figure_data = list(
      scatter_plot = list(
        source = "r_ggplot",
        method = "ols",
        x_label = x_label,
        y_label = y_label,
        observations = observations,
        fit = list(
          kind = "ols",
          conf_level = conf_level,
          line = line
        ),
        n = length(observations),
        point_display = list(
          total_n = length(observations),
          included_n = length(observations),
          truncated = FALSE
        )
      )
    ),
    metadata = list(
      source = "r_ggplot",
      available_figures = list(list(
        type = "scatter_plot",
        label = "Regression scatterplot",
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

mellio_ggplot_category_order <- function(values) {
  text <- as.character(values)
  text <- text[!is.na(text) & nzchar(text)]
  if (is.factor(values)) {
    levels_present <- intersect(levels(values), unique(text))
    c(levels_present, setdiff(unique(text), levels_present))
  } else {
    unique(text)
  }
}

mellio_ggplot_count_categories <- function(values) {
  order <- mellio_ggplot_category_order(values)
  text <- as.character(values)
  lapply(order, function(value) {
    list(
      value = value,
      label = value,
      count = as.integer(sum(text == value, na.rm = TRUE))
    )
  })
}

mellio_ggplot_value_categories <- function(values, y) {
  order <- mellio_ggplot_category_order(values)
  text <- as.character(values)
  lapply(order, function(value) {
    list(
      value = value,
      label = value,
      count = sum(y[text == value], na.rm = TRUE)
    )
  })
}

mellio_ggplot_box_stats <- function(values) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  n <- length(values)
  if (n < 1L) return(NULL)
  sorted <- sort(values)
  q1 <- unname(stats::quantile(sorted, 0.25, names = FALSE, type = 7))
  median <- unname(stats::quantile(sorted, 0.5, names = FALSE, type = 7))
  q3 <- unname(stats::quantile(sorted, 0.75, names = FALSE, type = 7))
  iqr <- q3 - q1
  lo <- q1 - 1.5 * iqr
  hi <- q3 + 1.5 * iqr
  in_range <- sorted[sorted >= lo & sorted <= hi]
  outliers <- sorted[sorted < lo | sorted > hi]
  list(
    min = sorted[[1]],
    q1 = q1,
    median = median,
    q3 = q3,
    max = sorted[[n]],
    iqr = iqr,
    whisker_lo = if (length(in_range)) in_range[[1]] else q1,
    whisker_hi = if (length(in_range)) in_range[[length(in_range)]] else q3,
    outliers = outliers
  )
}

mellio_ggplot_group_aes_safe <- function(mapping, aes, data, valid,
                                         target_label = NULL) {
  expr <- mapping[[aes]]
  if (is.null(expr)) return(TRUE)

  label <- mellio_ggplot_mapping_label(expr)
  if (!is.null(target_label) && identical(label, target_label)) return(TRUE)

  value <- mellio_ggplot_eval_mapping(expr, data)
  if (is.null(value) || length(value) != nrow(data)) return(FALSE)
  text <- as.character(value)
  text <- text[valid & !is.na(text) & nzchar(text)]
  length(unique(text)) <= 1L
}

mellio_ggplot_smooth_line <- function(data) {
  if (!is.data.frame(data) || !all(c("x", "y") %in% names(data))) return(list())
  valid <- is.finite(data$x) & is.finite(data$y)
  has_band <- all(c("ymin", "ymax") %in% names(data))
  if (has_band) {
    valid <- valid & is.finite(data$ymin) & is.finite(data$ymax)
  }
  if (sum(valid) < 2L) return(list())
  data <- data[valid, , drop = FALSE]
  data <- data[order(data$x), , drop = FALSE]
  lapply(seq_len(nrow(data)), function(i) {
    point <- list(x = data$x[[i]], y = data$y[[i]])
    if (has_band) {
      point$lower <- data$ymin[[i]]
      point$upper <- data$ymax[[i]]
    }
    point
  })
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
