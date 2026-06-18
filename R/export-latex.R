# LaTeX export
# Generates tabular environment with Unicode-to-LaTeX conversion

#' Convert a melliotab object to LaTeX
#'
#' @param x A melliotab object
#' @return Character string of LaTeX code
#' @export
mt_as_latex <- function(x) {
  if (!inherits(x, "melliotab")) {
    cli::cli_abort("{.arg x} must be a melliotab object.")
  }

  data <- x$data
  sc <- x$style_config
  col_types <- x$column_types
  n_cols <- ncol(data)
  n_rows <- nrow(data)

  # Build column alignment spec
  # First column left, rest right
  aligns <- c("l", rep("r", n_cols - 1))
  col_spec <- paste(aligns, collapse = "")

  lines <- character(0)

  # Table label and title (rendered above tabular, outside borders)
  lines <- c(lines, "\\begin{table}[htbp]")

  label_text <- NULL
  if (!is.null(x$number)) {
    num <- x$number
    if (isTRUE(sc$table_label$numbering == "roman")) {
      num <- to_roman(as.integer(num))
    }
    separator <- sc$table_label$separator %||% ""
    label_text <- paste0(sc$table_label$prefix, " ", num, separator)
    if (isTRUE(sc$table_label$bold)) {
      label_text <- paste0("\\textbf{", escape_latex(label_text), "}")
    } else {
      label_text <- escape_latex(label_text)
    }
  }

  title_text <- NULL
  if (!is.null(x$title)) {
    title_text <- apply_case(x$title, sc$table_title$case)
    if (isTRUE(sc$table_title$italic)) {
      title_text <- paste0("\\textit{", escape_latex(title_text), "}")
    } else {
      title_text <- escape_latex(title_text)
    }
  }

  if (!is.null(label_text) && !is.null(title_text) &&
      isFALSE(sc$table_title$separate_line)) {
    joiner <- if (nzchar(sc$table_label$separator %||% "")) " " else ": "
    lines <- c(lines, paste0(label_text, joiner, title_text))
  } else {
    if (!is.null(label_text)) lines <- c(lines, label_text)
    if (!is.null(title_text)) lines <- c(lines, paste0("\\vspace{1mm}", title_text))
  }

  lines <- c(lines, paste0("\\begin{tabular}{", col_spec, "}"))

  # Top rule
  if (isTRUE(sc$borders$top)) {
    lines <- c(lines, "\\toprule")
  }

  # Header row
  header_vals <- vapply(names(data), function(h) {
    h_latex <- unicode_to_latex(escape_latex(h))
    # Italicize stat symbols
    if (isTRUE(sc$italic_stat_headers) && is_stat_symbol(h)) {
      h_latex <- paste0("\\textit{", escape_latex(h), "}")
      h_latex <- unicode_to_latex(h_latex)
    }
    h_latex
  }, character(1), USE.NAMES = FALSE)

  lines <- c(lines, paste0(paste(header_vals, collapse = " & "), " \\\\"))

  # Header bottom rule
  if (isTRUE(sc$borders$header_bottom)) {
    lines <- c(lines, "\\midrule")
  }

  # Data rows
  for (i in seq_len(n_rows)) {
    row_vals <- vapply(seq_len(n_cols), function(j) {
      val <- data[i, j]
      unicode_to_latex(escape_latex(as.character(val)))
    }, character(1))

    row_line <- paste(row_vals, collapse = " & ")

    # Optional internal row borders
    if (isTRUE(sc$borders$internal_rows) && i < n_rows) {
      lines <- c(lines, paste0(row_line, " \\\\ \\hline"))
    } else {
      lines <- c(lines, paste0(row_line, " \\\\"))
    }
  }

  # Bottom rule
  if (isTRUE(sc$borders$bottom)) {
    lines <- c(lines, "\\bottomrule")
  }

  lines <- c(lines, "\\end{tabular}")

  # Note
  if (!is.null(x$note)) {
    note_label <- sc$notes$general_label
    note_text <- capitalize_note_after_label(x$note, note_label)
    note_text <- unicode_to_latex(escape_latex(note_text))
    if (!is.null(note_label) && nzchar(note_label)) {
      if (isTRUE(sc$notes$general_label_italic)) {
        note_label <- paste0("\\textit{", escape_latex(note_label), "}")
      } else {
        note_label <- escape_latex(note_label)
      }
      note_text <- paste(note_label, note_text)
    }
    lines <- c(lines, "")
    lines <- c(lines, paste0("\\vspace{2mm}\\noindent\\small ", note_text))
  }

  # Source
  if (!is.null(x$source) && !is.null(sc$notes$source_label)) {
    source_text <- paste(
      escape_latex(sc$notes$source_label),
      unicode_to_latex(escape_latex(x$source))
    )
    lines <- c(lines, "")
    lines <- c(lines, paste0("\\vspace{1mm}\\noindent\\small ", source_text))
  }

  lines <- c(lines, "\\end{table}")

  paste(lines, collapse = "\n")
}
