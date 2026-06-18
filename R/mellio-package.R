#' mellio: Publication-Ready Tables and Statistical Results
#'
#' mellio provides two main workflows:
#'
#' * `mellio_open()` sends supported R objects to the Mellio web app.
#'   Statistical results open as Result Cards, tabular data opens in the
#'   Tables workspace, and supported plots or image files open in Figures.
#' * `melliotab()` creates publication-ready tables directly in R for data
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
#' @seealso [mellio_open()], [melliotab()], [mt_compare()], [mt_save()],
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
