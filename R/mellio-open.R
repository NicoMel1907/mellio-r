#' Open a table, statistical result, or plot in Mellio
#'
#' Unified entry point for sending R objects to the Mellio web app.
#' Statistical results and tables open as Stats Result Cards; plots open
#' in the web app's Figures editor.
#'
#' @param x An R object supported by Mellio.
#' @param ... Additional arguments forwarded to the relevant Mellio
#'   constructor or payload method.
#' @param browse Open the browser? `TRUE` by default. Set `FALSE` to
#'   return the URL without launching anything.
#' @param .call Internal captured user call for provenance.
#' @section Supported inputs:
#' `mellio_open()` routes objects by class:
#'
#' * Statistical results such as `htest`, `lm`, `glm`, `aov`, `lavaan`,
#'   mixed-model, survival, ordinal, mediation, and selected optional-package
#'   result classes open as structured Result Cards.
#' * Data frames, matrices, base `table` objects, and non-statistical
#'   `melliotab` objects open in the Tables workspace.
#' * `ggplot2`, lattice, htmlwidget, recorded plot, and supported image-file
#'   inputs open in Figures.
#'
#' Optional model and figure integrations use packages listed in `Suggests`.
#' If an optional package is not installed, the corresponding object class is
#' skipped or reported with an informative message.
#'
#' @section Privacy and destination:
#' `mellio_open()` opens `https://www.mellioapp.com` by default. Advanced
#' users can set `options(mellio.editor_url = "https://...")`, but should only
#' use a trusted Mellio deployment. The R payload is encoded in the URL
#' fragment. URL fragments are not sent as HTTP requests to the server, but the
#' full URL can still be visible to the browser, the opened web app, browser
#' history, extensions, and anyone the URL is shared with.
#'
#' By default, Mellio payloads include R/package-version metadata and data
#' fingerprints where available. Local machine details are opt-in:
#' `options(mellio.provenance = "full")`. To omit provenance metadata, use
#' `options(mellio.provenance = FALSE)`.
#' @section Plot size limits:
#' Plots are transported through the URL hash. Very large images may exceed
#' browser URL limits; reduce `width`, `height`, or `dpi`, or save the image and
#' upload it manually in Mellio.
#' @return Invisibly, the URL string.
#' @export
#' @family R bridge
#' @examples
#' \donttest{
#' if (interactive()) {
#'   mellio_open(t.test(extra ~ group, data = sleep))
#'
#'   m1 <- lm(Ozone ~ Temp, data = airquality)
#'   m2 <- lm(Ozone ~ Temp + Wind, data = airquality)
#'   m3 <- lm(Ozone ~ Temp + Wind + Solar.R, data = airquality)
#'   mellio_open(
#'     m1, m2, m3,
#'     labels = c("Step 1", "Step 2", "Step 3"),
#'     dep.var.labels = "Ozone concentration"
#'   )
#' }
#' }
mellio_open <- function(x, ..., browse = TRUE, .call = NULL) {
  user_env <- parent.frame()
  dots <- list(...)
  dot_names <- names(dots) %||% rep("", length(dots))
  if (length(dots) && !nzchar(dot_names[[1]]) &&
      is.logical(dots[[1]]) && length(dots[[1]]) == 1L &&
      identical(browse, TRUE)) {
    browse <- dots[[1]]
    dots <- dots[-1L]
  }

  if (is.null(.call)) {
    user_call <- match.call()$x
    .call <- if (!is.null(user_call)) {
      paste(deparse(user_call, width.cutoff = 500L), collapse = " ")
    } else NA_character_
  }

  if (mellio_should_build_model_comparison_table(x, dots)) {
    tbl <- mellio_model_comparison_table(x, dots)
    return(send_table_to_mellio(tbl, browse = browse))
  }

  if (inherits(x, "htest") || inherits(x, "pairwise.htest")) {
    payload <- do.call(
      mellio_payload,
      c(list(x, .call = .call, .env = user_env), dots)
    )
    return(send_payload_to_stats(payload, browse = browse))
  }

  if (inherits(x, "summary.aov")) {
    payload <- do.call(
      mellio_payload,
      c(list(x, .call = .call, .env = user_env), dots)
    )
    return(send_payload_to_stats(payload, browse = browse))
  }

  do.call(
    mellio_open_dispatch,
    c(list(x, browse = browse, .call = .call), dots)
  )
}

mellio_open_dispatch <- function(x, browse = TRUE, ..., .call = NULL) {
  UseMethod("mellio_open_dispatch")
}

#' @export
mellio_open_dispatch.default <- function(x, browse = TRUE, ..., .call = NULL) {
  payload <- mellio_payload(x, ..., .call = .call)
  send_payload_to_stats(payload, browse = browse)
}

#' @export
mellio_open_dispatch.dunn_test <- function(x, browse = TRUE, ..., .call = NULL) {
  payload <- mellio_payload(x, ..., .call = .call)
  send_payload_to_stats(payload, browse = browse)
}

#' @export
mellio_open_dispatch.t_test <- function(x, browse = TRUE, ..., .call = NULL) {
  payload <- mellio_payload(x, ..., .call = .call)
  send_payload_to_stats(payload, browse = browse)
}

#' @export
mellio_open_dispatch.wilcox_test <- function(x, browse = TRUE, ..., .call = NULL) {
  payload <- mellio_payload(x, ..., .call = .call)
  send_payload_to_stats(payload, browse = browse)
}

#' @export
mellio_open_dispatch.mellio_payload <- function(x, browse = TRUE, ..., .call = NULL) {
  send_payload_to_stats(x, browse = browse)
}

#' @export
mellio_open_dispatch.anova <- function(x, browse = TRUE, ..., .call = NULL) {
  payload <- tryCatch(
    mellio_payload(x, ..., .call = .call),
    error = function(e) NULL
  )
  if (!is.null(payload)) {
    return(send_payload_to_stats(payload, browse = browse))
  }

  send_table_to_mellio(melliotab(x, ...), browse = browse)
}

#' @export
mellio_open_dispatch.data.frame <- function(x, browse = TRUE, ..., .call = NULL) {
  # Bare data.frames have no statistical-class wrapper — pure tabular
  # data goes straight to the Tables workspace. Statistical results
  # (lm, htest, lavaan, …) route through mellio_open_dispatch.default
  # to Stats, and the melliotab method below covers the
  # melliotab(data.frame) and melliotab(model) cases separately.
  tbl <- melliotab(x, ...)
  send_table_to_mellio(tbl, browse = browse)
}

#' @export
mellio_open_dispatch.melliotab <- function(x, browse = TRUE, ..., .call = NULL) {
  # Discriminator: a melliotab with a `pvalue` column type carries
  # statistical content (built from lm / glm / aov / htest / lavaan /
  # mediate / fa / FitDiff, etc.) and routes to Stats so the result can
  # be narrated. A melliotab without a p-value column is pure tabular
  # (built from a data.frame, matrix, summary, contingency table) and
  # routes to Tables.
  #
  # Why $model isn't the discriminator: melliotab.character,
  # .summaryDefault and .table also populate $model with the input
  # object even though their output is non-statistical. Column-type
  # detection is the reliable signal.
  if ("pvalue" %in% (x$column_types %||% character(0))) {
    payload <- mellio_payload(x, .call = .call)
    send_payload_to_stats(payload, browse = browse)
  } else {
    send_table_to_mellio(x, browse = browse)
  }
}

#' @export
mellio_open_dispatch.matrix <- function(x, browse = TRUE, ..., .call = NULL) {
  # Correlation matrices have statistical semantics and figure affordances,
  # so they route to Stats. Plain numeric matrices remain table inputs.
  if (is_correlation_matrix(x)) {
    payload <- mellio_payload(x, ..., .call = .call)
    if (identical(payload$type, "correlation_matrix")) {
      return(send_payload_to_stats(payload, browse = browse))
    }
  }

  tbl <- melliotab(x, ...)
  send_table_to_mellio(tbl, browse = browse)
}

#' @export
mellio_open_dispatch.table <- function(x, browse = TRUE, ..., .call = NULL) {
  # xtabs() / table() output — raw counts with no test attached. Goes
  # to Tables. Users who want chi-square inference call chisq.test()
  # which has class `htest` and routes to Stats.
  tbl <- melliotab(x, ...)
  send_table_to_mellio(tbl, browse = browse)
}

#' @export
mellio_open_dispatch.gg <- function(x, browse = TRUE, ..., .call = NULL) {
  send_figure_to_mellio(
    mellio_open_figure_from_gg(x, ..., .call = .call),
    browse = browse
  )
}

#' @export
mellio_open_dispatch.ggplot <- mellio_open_dispatch.gg

#' @export
mellio_open_dispatch.recordedplot <- function(x, browse = TRUE, ..., .call = NULL) {
  send_figure_to_mellio(
    mellio_open_figure_from_recordedplot(x, ...),
    browse = browse
  )
}

#' @export
mellio_open_dispatch.trellis <- function(x, browse = TRUE, ..., .call = NULL) {
  send_figure_to_mellio(
    mellio_open_figure_from_printed(x, ...),
    browse = browse
  )
}

#' @export
mellio_open_dispatch.htmlwidget <- function(x, browse = TRUE, ..., .call = NULL) {
  send_figure_to_mellio(
    mellio_open_figure_from_htmlwidget(x, ...),
    browse = browse
  )
}

#' @export
mellio_open_dispatch.ggExtraPlot <- function(x, browse = TRUE, ..., .call = NULL) {
  send_figure_to_mellio(
    mellio_open_figure_from_printed(x, ...),
    browse = browse
  )
}

#' @export
mellio_open_dispatch.character <- function(x, browse = TRUE, ..., .call = NULL) {
  if (mellio_open_is_figure_path(x)) {
    return(send_figure_to_mellio(
      mellio_open_figure_from_file(x, ...),
      browse = browse
    ))
  }

  mellio_open_dispatch.default(x, browse = browse, ..., .call = .call)
}

mellio_open_is_figure_path <- function(x) {
  if (!is.character(x) || length(x) != 1L || !file.exists(x)) return(FALSE)
  ext <- tolower(tools::file_ext(x))
  ext %in% c("png", "jpg", "jpeg", "tiff", "tif", "bmp", "pdf", "svg")
}
