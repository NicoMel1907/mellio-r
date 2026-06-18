#' Change the citation style
#'
#' @param x A melliotab object
#' @param style Style name: "apa7" or "ieee"
#' @return Modified melliotab object
#' @export
mt_set_style <- function(x, style = c("apa7", "ieee")) {
  style <- match.arg(style)
  sc <- get_style(style)

  x$style <- style
  x$style_config <- sc
  x$options$remove_leading_zeros <- isTRUE(sc$remove_leading_zeros)
  x$options$italic_stat_headers <- isTRUE(sc$italic_stat_headers)
  x$options$bold_section_titles <- isTRUE(sc$bold_section_titles)

  col_types <- detect_column_types(names(x$raw_data))
  lz_cols <- detect_leading_zero_cols(names(x$raw_data))
  x$column_types <- col_types
  x$data <- format_table_data(
    x$raw_data, col_types, lz_cols,
    decimals = x$decimals, p_decimals = x$p_decimals,
    remove_lz = isTRUE(sc$remove_leading_zeros),
    style_config = sc
  )

  if (isTRUE(x$options$is_correlation)) {
    x <- apply_correlation_formatting(x)
  }

  if (isTRUE(x$options$sig_stars)) {
    x <- apply_sig_stars(x)
  }

  x
}

#' Set decimal places
#'
#' @param x A melliotab object
#' @param decimals Decimal places for estimates/statistics
#' @param p_decimals Decimal places for p-values
#' @return Modified melliotab object
#' @export
mt_decimals <- function(x, decimals = 2L, p_decimals = 3L) {
  x$decimals <- as.integer(decimals)
  x$p_decimals <- as.integer(p_decimals)

  lz_cols <- detect_leading_zero_cols(names(x$raw_data))
  x$data <- format_table_data(
    x$raw_data, x$column_types, lz_cols,
    decimals = x$decimals, p_decimals = x$p_decimals,
    remove_lz = x$options$remove_leading_zeros,
    style_config = x$style_config
  )

  if (isTRUE(x$options$is_correlation)) {
    x <- apply_correlation_formatting(x)
  }

  x
}

#' Control leading zero removal
#'
#' @param x A melliotab object
#' @param enabled Whether to remove leading zeros
#' @return Modified melliotab object
#' @export
mt_remove_leading_zeros <- function(x, enabled = TRUE) {
  x$options$remove_leading_zeros <- enabled

  lz_cols <- detect_leading_zero_cols(names(x$raw_data))
  x$data <- format_table_data(
    x$raw_data, x$column_types, lz_cols,
    decimals = x$decimals, p_decimals = x$p_decimals,
    remove_lz = enabled,
    style_config = x$style_config
  )

  if (isTRUE(x$options$is_correlation)) {
    x <- apply_correlation_formatting(x)
  }

  x
}

#' Add significance stars
#'
#' @param x A melliotab object
#' @param target Column to add stars to ("auto" to detect)
#' @param remove_p Whether to remove the p-value column
#' @param levels Named vector of significance thresholds
#' @return Modified melliotab object
#' @export
mt_sig_stars <- function(x, target = "auto", remove_p = TRUE,
                          levels = c("*" = 0.05, "**" = 0.01, "***" = 0.001)) {
  apply_sig_stars(x, target = target, remove_p = remove_p, levels = levels)
}

#' Format confidence intervals
#'
#' @param x A melliotab object
#' @param decimals Decimal places for CI values (NULL uses table default)
#' @param bracket Whether to use bracket notation
#' @return Modified melliotab object
#' @export
mt_format_ci <- function(x, decimals = NULL, bracket = TRUE) {
  dec <- decimals %||% x$decimals
  ci_cols <- detect_ci_cols(names(x$data))

  for (j in which(ci_cols)) {
    x$data[[j]] <- vapply(as.character(x$data[[j]]),
                           function(v) format_ci_value(v, dec),
                           character(1), USE.NAMES = FALSE)
  }

  x$options$format_ci <- TRUE
  x
}

#' Simplify SPSS verbose headers
#'
#' @param x A melliotab object
#' @return Modified melliotab object
#' @export
mt_simplify_headers <- function(x) {
  old_names <- names(x$data)
  new_names <- vapply(old_names, simplify_spss_header, character(1),
                       USE.NAMES = FALSE)
  names(x$data) <- new_names
  names(x$raw_data) <- new_names
  x$column_types <- detect_column_types(new_names)
  x
}

#' Set table title
#'
#' @param x A melliotab object
#' @param title Title text
#' @return Modified melliotab object
#' @export
mt_title <- function(x, title) {
  x$title <- title
  x
}

#' Set table number
#'
#' @param x A melliotab object
#' @param number Table number
#' @return Modified melliotab object
#' @export
mt_number <- function(x, number) {
  x$number <- number
  x
}

#' Set table note
#'
#' @param x A melliotab object
#' @param note Note text
#' @return Modified melliotab object
#' @export
mt_note <- function(x, note) {
  x$note <- note
  x
}

#' Set table source
#'
#' @param x A melliotab object
#' @param source Source text
#' @return Modified melliotab object
#' @export
mt_source <- function(x, source) {
  x$source <- source
  x
}

#' Set correlation matrix display options
#'
#' @param x A melliotab object
#' @param mode Diagonal display: "dash", "one", "blank"
#' @param triangle Triangle to show: "all", "lower", "upper"
#' @return Modified melliotab object
#' @export
mt_diagonal <- function(x, mode = c("dash", "one", "blank"),
                         triangle = c("all", "lower", "upper")) {
  mode <- match.arg(mode)
  triangle <- match.arg(triangle)

  x$options$is_correlation <- TRUE
  x$options$diagonal_mode <- mode
  x$options$triangle <- triangle

  lz_cols <- detect_leading_zero_cols(names(x$raw_data))
  x$data <- format_table_data(
    x$raw_data, x$column_types, lz_cols,
    decimals = x$decimals, p_decimals = x$p_decimals,
    remove_lz = x$options$remove_leading_zeros,
    style_config = x$style_config
  )
  x <- apply_correlation_formatting(x)

  x
}

#' Add a spanning column header
#'
#' @param x A melliotab object
#' @param label Spanner label text
#' @param columns Column indices or names covered by the spanner
#' @param level Spanner level (1 = closest to data)
#' @return Modified melliotab object
#' @export
mt_spanner <- function(x, label, columns, level = 1L) {
  if (is.character(columns)) {
    columns <- match(columns, names(x$data))
    columns <- columns[!is.na(columns)]
  }
  x$spanners <- c(x$spanners, list(list(
    label = label, columns = columns, level = as.integer(level)
  )))
  x
}

#' Add a section title row
#'
#' @param x A melliotab object
#' @param label Section title text
#' @param before Row index to insert before
#' @param after Row index to insert after
#' @return Modified melliotab object
#' @export
mt_section_title <- function(x, label, before = NULL, after = NULL) {
  x$section_titles <- c(x$section_titles, list(list(
    label = label, before = before, after = after
  )))
  x
}

#' Indent row labels
#'
#' @param x A melliotab object
#' @param rows Row indices to indent
#' @param level Indentation level (1-3)
#' @return Modified melliotab object
#' @export
mt_indent <- function(x, rows, level = 1L) {
  for (r in rows) {
    x$indent_levels[as.character(r)] <- as.integer(level)
  }
  x
}
