# Markdown export
# Pipe-delimited table format

#' Convert a melliotab object to Markdown
#'
#' @param x A melliotab object
#' @return Character string of Markdown table
#' @export
mt_as_markdown <- function(x) {
  if (!inherits(x, "melliotab")) {
    cli::cli_abort("{.arg x} must be a melliotab object.")
  }

  data <- x$data
  sc <- x$style_config
  n_cols <- ncol(data)
  n_rows <- nrow(data)

  lines <- character(0)
  label <- NULL

  # Table label and title
  if (!is.null(x$number)) {
    num <- x$number
    if (isTRUE(sc$table_label$numbering == "roman")) {
      num <- to_roman(as.integer(num))
    }
    separator <- sc$table_label$separator %||% ""
    label <- paste0(sc$table_label$prefix, " ", num, separator)
    if (isTRUE(sc$table_label$bold)) {
      label <- paste0("**", label, "**")
    }
  }

  title_text <- NULL
  if (!is.null(x$title)) {
    title_text <- apply_case(x$title, sc$table_title$case)
    if (isTRUE(sc$table_title$italic)) {
      title_text <- paste0("*", title_text, "*")
    }
  }

  if (!is.null(label) && !is.null(title_text) &&
      isFALSE(sc$table_title$separate_line)) {
    joiner <- if (nzchar(sc$table_label$separator %||% "")) " " else ": "
    lines <- c(lines, paste0(label, joiner, title_text))
  } else {
    if (!is.null(label)) lines <- c(lines, label)
    if (!is.null(title_text)) lines <- c(lines, title_text)
  }

  if (length(lines) > 0) {
    lines <- c(lines, "")
  }

  # Calculate column widths
  headers <- names(data)
  if (isTRUE(sc$italic_stat_headers)) {
    headers <- vapply(headers, function(h) {
      if (is_stat_symbol(h)) paste0("*", h, "*") else h
    }, character(1), USE.NAMES = FALSE)
  }

  col_widths <- vapply(seq_len(n_cols), function(j) {
    vals <- c(headers[j], as.character(data[[j]]))
    max(nchar(vals), na.rm = TRUE)
  }, integer(1))
  col_widths <- pmax(col_widths, 3L)

  # Header row
  header_parts <- vapply(seq_len(n_cols), function(j) {
    formatC(headers[j], width = col_widths[j], flag = if (j == 1) "-" else "")
  }, character(1))
  lines <- c(lines, paste0("| ", paste(header_parts, collapse = " | "), " |"))

  # Separator row
  sep_parts <- vapply(seq_len(n_cols), function(j) {
    align_char <- if (j == 1) ":" else "-"
    right_char <- if (j > 1) ":" else "-"
    paste0(align_char, strrep("-", col_widths[j] - 1), right_char)
  }, character(1))
  lines <- c(lines, paste0("| ", paste(sep_parts, collapse = " | "), " |"))

  # Data rows
  for (i in seq_len(n_rows)) {
    row_parts <- vapply(seq_len(n_cols), function(j) {
      val <- as.character(data[i, j])
      if (is.na(val)) val <- ""
      formatC(val, width = col_widths[j], flag = if (j == 1) "-" else "")
    }, character(1))
    lines <- c(lines, paste0("| ", paste(row_parts, collapse = " | "), " |"))
  }

  # Note
  if (!is.null(x$note)) {
    note_label <- sc$notes$general_label
    note_text <- capitalize_note_after_label(x$note, note_label)
    lines <- c(lines, "")
    if (!is.null(note_label) && nzchar(note_label)) {
      if (isTRUE(sc$notes$general_label_italic)) {
        lines <- c(lines, paste0("*", note_label, "* ", note_text))
      } else {
        lines <- c(lines, paste(note_label, note_text))
      }
    } else {
      lines <- c(lines, note_text)
    }
  }

  # Source
  if (!is.null(x$source) && !is.null(sc$notes$source_label)) {
    lines <- c(lines, "")
    lines <- c(lines, paste(sc$notes$source_label, x$source))
  }

  paste(lines, collapse = "\n")
}
