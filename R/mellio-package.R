#' mellio: Polished, Editable Tables and Statistical Results
#'
#' mellio provides two main workflows:
#'
#' * `mellio_open()` sends supported R objects to the Mellio web app.
#'   Statistical results open as Result Cards, tabular data opens in the
#'   Tables workspace, and supported plots or image files open in Figures.
#' * `melliotab()` creates polished, editable tables directly in R for data
#'   frames, model summaries, correlation matrices, and model comparisons.
#'
#' @section Common `mellio_open()` inputs:
#' `mellio_open()` supports common `htest` objects such as `t.test()`,
#' `cor.test()`, and `chisq.test()`; model objects such as `lm`, `glm`,
#' `aov`, `lme4` mixed models, `lavaan` fits, survival models, ordinal
#' models, and selected optional-package result classes; data frames,
#' matrices, and base `table` objects; `ggplot2`, lattice, htmlwidget,
#' recorded plot, and image-file inputs.
#'
#' @section Common `melliotab()` workflow:
#' Use `melliotab(x)` to create a table, then apply table modifiers such as
#' `mt_title()`, `mt_note()`, `mt_decimals()`, `mt_format_ci()`,
#' `mt_sig_stars()`, `mt_spanner()`, `mt_section_title()`, and
#' `mt_simplify_headers()`. Pass two or more model objects to
#' `melliotab(m1, m2, ...)` to create a side-by-side model comparison table.
#' For file-based handoff, use `mt_copy()` or `mt_save()` to copy or save a
#' finished table as HTML, LaTeX, or Markdown.
#'
#' @section Quick reference:
#' \tabular{lll}{
#' Function or option \tab What it does \tab Common values \cr
#' \code{melliotab()} \tab Creates a formatted table in R \tab \code{style = "apa7"} or \code{"ieee"} \cr
#' \code{mellio_open()} \tab Opens supported objects in Mellio \tab Models, tests, tables, data, plots \cr
#' \code{style} / \code{mt_set_style()} \tab Sets or changes table style \tab \code{"apa7"}, \code{"ieee"} \cr
#' \code{title} / \code{mt_title()} \tab Sets the table title \tab Text \cr
#' \code{number} / \code{mt_number()} \tab Sets the table number \tab Number or text \cr
#' \code{note} / \code{mt_note()} \tab Adds a table note \tab Text \cr
#' \code{source} / \code{mt_source()} \tab Adds source text \tab Text \cr
#' \code{decimals}, \code{p_decimals} / \code{mt_decimals()} \tab Controls rounding \tab \code{decimals = 2}, \code{p_decimals = 3} \cr
#' \code{mt_sig_stars()} \tab Adds significance stars to an existing table \tab \code{remove_p = TRUE} or \code{FALSE} \cr
#' \code{mt_remove_leading_zeros()} \tab Controls leading zeros \tab \code{TRUE}, \code{FALSE} \cr
#' \code{mt_diagonal()} \tab Formats correlation matrices \tab \code{mode = "dash"}; \code{triangle = "lower"} \cr
#' \code{mt_spanner()} \tab Adds a spanning column header \tab Label text and columns \cr
#' \code{mt_section_title()} \tab Adds a section-title row \tab \code{before =} or \code{after =} a row number \cr
#' \code{mt_indent()} \tab Indents selected rows \tab \code{rows =}; \code{level = 1}, \code{2}, or \code{3} \cr
#' \code{mt_copy()} / \code{mt_save()} \tab Copies or saves a table \tab Clipboard, \code{.html}, \code{.tex}, \code{.md}
#' }
#'
#' @section Privacy and provenance:
#' By default, Mellio payloads include R/package-version metadata and data
#' fingerprints where available. Local machine details such as user name,
#' host name, working directory, git state, and script path are not included
#' unless you opt in with `options(mellio.provenance = "full")`. Set
#' `options(mellio.provenance = FALSE)` to omit provenance metadata.
#'
#' @section Citation:
#' Use `citation("mellio")` to get the package citation.
#'
#' @seealso [mellio_open()], [melliotab()], [mt_save()],
#'   [mt_copy()]
"_PACKAGE"

#' @importFrom rlang %||%
#' @importFrom rlang check_installed
#' @importFrom rlang set_names
#' @importFrom cli cli_warn
#' @importFrom cli cli_inform
#' @importFrom stats setNames
#' @importFrom utils head
NULL
