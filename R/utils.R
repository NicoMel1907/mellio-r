# Shared utility functions

#' Convert integer to Roman numeral
#' @keywords internal
to_roman <- function(num) {
  as.character(utils::as.roman(num))
}

#' Escape special HTML characters
#' @keywords internal
escape_html <- function(text) {
  text <- gsub("&", "&amp;", text, fixed = TRUE)
  text <- gsub("<", "&lt;", text, fixed = TRUE)
  text <- gsub(">", "&gt;", text, fixed = TRUE)
  text
}

#' Capitalize prose note text after a period-ending note label
#' @keywords internal
capitalize_note_after_label <- function(note, label) {
  note <- trimws(note %||% "")
  label <- trimws(label %||% "")
  if (!nzchar(note) || !nzchar(label) || !grepl("\\.\\s*$", label, perl = TRUE)) {
    return(note)
  }
  if (grepl("^\\*+\\s*p\\s*[<=>\u2264\u2265]", note, perl = TRUE) ||
      grepl("^(p|n|r|t|z|d|b|f)\\s*([<=>\u2264\u2265]|\\()", note, perl = TRUE)) {
    return(note)
  }
  sub("^(\\s*[\"'\\(\\[]*)([a-z])", "\\1\\U\\2", note, perl = TRUE)
}

#' Escape special LaTeX characters
#' @keywords internal
escape_latex <- function(text) {
  text <- gsub("\\", "\\textbackslash{}", text, fixed = TRUE)
  text <- gsub("{", "\\{", text, fixed = TRUE)
  text <- gsub("}", "\\}", text, fixed = TRUE)
  text <- gsub("$", "\\$", text, fixed = TRUE)
  text <- gsub("#", "\\#", text, fixed = TRUE)
  text <- gsub("%", "\\%", text, fixed = TRUE)
  text <- gsub("&", "\\&", text, fixed = TRUE)
  text <- gsub("_", "\\_", text, fixed = TRUE)
  text <- gsub("^", "\\^{}", text, fixed = TRUE)
  text <- gsub("~", "\\~{}", text, fixed = TRUE)
  text
}

#' Convert Unicode statistical symbols to LaTeX equivalents
#' @keywords internal
unicode_to_latex <- function(text) {
  if (is.null(text) || !nzchar(text)) return("")

  # Order matters: longer patterns first to avoid partial matches
  replacements <- list(
    c("\u03B7\u209A\u00B2", "$\\eta_p^2$"),
    c("\u03B7p\u00B2",      "$\\eta_p^2$"),
    c("\u03B7\u00B2",       "$\\eta^2$"),
    c("\u03C9\u00B2",       "$\\omega^2$"),
    c("\u03C7\u00B2",       "$\\chi^2$"),
    c("\u03B5\u00B2",       "$\\varepsilon^2$"),
    c("\u0394R\u00B2",      "$\\Delta R^2$"),
    c("R\u00B2",            "$R^2$"),
    c("f\u00B2",            "$f^2$"),
    c("\u03B2",             "$\\beta$"),
    c("\u03B7",             "$\\eta$"),
    c("\u03C9",             "$\\omega$"),
    c("\u03B1",             "$\\alpha$"),
    c("\u03BA",             "$\\kappa$"),
    c("\u03C1",             "$\\rho$"),
    c("\u03C4",             "$\\tau$"),
    c("\u03C6",             "$\\varphi$"),
    c("\u0394",             "$\\Delta$"),
    c("\u039B",             "$\\Lambda$"),
    c("\u03B8",             "$\\theta$"),
    c("\u00B2",             "$^2$"),
    c("\u2014",             "---"),
    c("\u2013",             "--"),
    c("\u2264",             "$\\leq$"),
    c("\u2265",             "$\\geq$"),
    c("\u2260",             "$\\neq$"),
    c("\u00D7",             "$\\times$")
  )

  result <- text
  for (r in replacements) {
    result <- gsub(r[1], r[2], result, fixed = TRUE)
  }
  result
}

#' Extract trailing annotations (significance stars, superscripts)
#' @keywords internal
extract_annotation <- function(text) {
  pattern <- "([*\u00B2\u00B3\u00B9\u1D43-\u1D6A\u1D9C-\u1DBF\u02B0-\u02FF\u2070-\u209F]+)$"
  m <- regmatches(text, regexpr(pattern, text, perl = TRUE))
  if (length(m) == 0 || !nzchar(m)) {
    return(list(core = text, annotation = ""))
  }
  list(
    core = trimws(sub(pattern, "", text, perl = TRUE)),
    annotation = m
  )
}

#' Convert text to title case
#' @keywords internal
to_title_case <- function(text) {
  tools::toTitleCase(tolower(text))
}

#' Convert text to sentence case
#' @keywords internal
to_sentence_case <- function(text) {
  result <- tolower(text)
  substring(result, 1, 1) <- toupper(substring(result, 1, 1))
  result
}
