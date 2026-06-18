# Citation style definitions
# Ported from citation-styles.js — APA and IEEE citation styles
# Each style includes both table and figure formatting configurations

.styles <- list(
  apa7 = list(
    label = "APA 7th Edition",
    table_label = list(prefix = "Table", numbering = "arabic", bold = TRUE),
    table_title = list(
      position = "above", separate_line = TRUE, italic = TRUE,
      case = "title", align = "left"
    ),
    figure_label = list(
      prefix = "Figure", numbering = "arabic", bold = TRUE, italic = FALSE
    ),
    figure_title = list(
      position = "above", separate_line = TRUE, italic = TRUE,
      case = "title", align = "left"
    ),
    figure_notes = list(
      general_label = "Note.", general_label_italic = TRUE,
      specific_marker = "superscript_lowercase_alpha",
      probability_markers = c("*p < .05", "**p < .01", "***p < .001"),
      source_label = NULL,
      order = c("general", "specific", "probability")
    ),
    label_align = "left",
    borders = list(
      top = TRUE, bottom = TRUE, header_bottom = TRUE,
      vertical = FALSE, internal_rows = FALSE
    ),
    notes = list(
      position = "below", general_label = "Note.",
      general_label_italic = TRUE,
      specific_marker = "superscript_lowercase_alpha",
      probability_markers = c("*p < .05", "**p < .01", "***p < .001"),
      order = c("general", "specific", "probability")
    ),
    spacing = "double",
    font = "Times New Roman",
    font_size = 12,
    cell_padding = c(4, 8),
    italic_stat_headers = TRUE,
    remove_leading_zeros = TRUE,
    bold_section_titles = TRUE,
    p_decimals = 3L,
    stat_decimals = 2L
  ),

  ieee = list(
    label = "IEEE",
    table_label = list(
      prefix = "TABLE", numbering = "roman", bold = TRUE, case = "upper"
    ),
    table_title = list(
      position = "above", separate_line = TRUE, italic = FALSE,
      case = "upper", align = "center"
    ),
    figure_label = list(
      prefix = "Fig.", numbering = "arabic", bold = FALSE, italic = FALSE,
      separator = "."
    ),
    figure_title = list(
      position = "below", separate_line = FALSE, italic = FALSE,
      case = "sentence", align = "center"
    ),
    figure_notes = list(
      general_label = NULL, source_label = NULL
    ),
    label_align = "center",
    borders = list(
      top = TRUE, bottom = TRUE, header_bottom = TRUE,
      vertical = FALSE, internal_rows = FALSE
    ),
    notes = list(
      position = "below", general_label = NULL,
      specific_marker = "superscript_lowercase_alpha"
    ),
    spacing = "single",
    font = "Times New Roman",
    font_size = 10,
    cell_padding = c(3, 6),
    italic_stat_headers = FALSE,
    remove_leading_zeros = FALSE,
    bold_section_titles = FALSE,
    uppercase_p = TRUE,
    p_decimals = 3L,
    stat_decimals = 2L
  )
)

#' Get a citation style configuration
#'
#' @param name Style name: "apa7" or "ieee"
#' @return A named list with all style properties
#' @keywords internal
get_style <- function(name) {
  name <- match.arg(tolower(name), names(.styles))
  .styles[[name]]
}

#' List available citation styles
#'
#' @return Character vector of style names
#' @keywords internal
list_styles <- function() {
  names(.styles)
}
