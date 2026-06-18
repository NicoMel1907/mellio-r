# R bridge — payload envelope.
#
# Builds the JSON Result Card sent from R to the Mellio web app.
# See docs/STATS-R-BRIDGE-SCHEMA.md for the wire format spec.
#
# This file holds mellio_payload() (the S3 generic), the envelope assembly,
# and helpers (ID generation, raw-output capture, JSON serialisation).
# Per-class extractors live in bridge-extract-*.R and dispatch via
# mellio_payload().
#
# No network here. mellio_open() handles browser transport in bridge-open helpers.

#' Build a Mellio Result Card payload from an R object
#'
#' Converts a supported R model or test object into a structured JSON
#' payload (a "Result Card") that Mellio renders in the Stats section.
#' Pure and offline — no network calls. Use [mellio_open()] to additionally
#' open the payload in the Mellio web app.
#'
#' Supported in v0: objects of class `htest` (`t.test`, `cor.test`,
#' `wilcox.test`, `chisq.test`, `fisher.test`).
#'
#' Schema: see `docs/STATS-R-BRIDGE-SCHEMA.md` in the repo. The payload
#' is value-only: numeric statistics, p-values, and CIs go through
#' unformatted; the Mellio web app applies citation-style rules
#' (decimals, italics, leading zeros).
#'
#' @param x A supported R object.
#' @param ... Additional arguments passed to the dispatched method.
#' @param .call Internal — captured user call, used to attribute
#'   provenance to the calling script rather than this dispatch.
#' @param controls Optional character vector of model terms to mark as
#'   controls (vs focal predictors) — recognised by lm/aov/anova methods.
#' @param focal Optional character vector of model terms to mark as the
#'   focal predictor(s). Recognised by lm/aov/anova methods.
#' @param name Optional name to attach to the payload (used by the
#'   numeric/vector descriptive method).
#' @param tests Multivariate test statistics to extract (manova method).
#'   Defaults to all four classical tests.
#' @param title Optional title override (psych/character/data.frame
#'   methods).
#' @return A list with class `mellio_payload`, structured according to
#'   the v0.1 schema.
#' @name mellio_payload
#' @export
#' @family R bridge
#' @seealso [mellio_open()] to open the result in Mellio. [mellio_to_json()] to
#'   serialise the payload offline.
#' @examples
#' tt <- t.test(extra ~ group, data = sleep)
#' p <- mellio_payload(tt)
#' p$type
#' p$fields$statistic
mellio_payload <- function(x, ...) {
  UseMethod("mellio_payload")
}

#' @rdname mellio_payload
#' @export
mellio_payload.mellio_payload <- function(x, ..., .call = NULL) {
  x
}

#' @rdname mellio_payload
#' @export
mellio_payload.default <- function(x, ..., .call = NULL) {
  if (is.function(x)) {
    stop(ms_function_hint(x, .call = .call), call. = FALSE)
  }
  if (ms_is_bruce_process(x)) {
    return(ms_bruce_process_payload(x, .call = .call))
  }
  cls <- class(x)
  printed <- ms_capture_output_safe(x)
  call_str <- ms_unsupported_call(.call)

  partial <- ms_partial_payload(x, cls = cls, printed = printed,
                                call_str = call_str)
  if (!is.null(partial)) return(partial)

  suggestions <- ms_unsupported_suggestions(
    cls = cls,
    call_str = call_str,
    printed_text = printed$text
  )
  fields <- list(
    class   = I(as.character(cls)),
    printed = printed$text,
    message = paste(
      "Mellio saved the original R output but could not safely extract",
      "structured fields for this object yet. Copy the part you need into",
      "Stats or Tables, or send us this class name so we can add a",
      "dedicated adapter."
    ),
    truncated = printed$truncated
  )
  if (length(suggestions)) fields$suggestions <- suggestions

  ms_build_envelope(
    type       = "unsupported",
    type_label = "Unrecognized object",
    call       = call_str,
    fields     = fields,
    raw_output = printed$text,
    card_kind  = "unsupported"
  )
}

#' @rdname mellio_payload
#' @export
mellio_payload.NULL <- function(x, ..., .call = NULL) {
  stop(ms_null_hint(.call = .call), call. = FALSE)
}

#' @rdname mellio_payload
#' @export
mellio_payload.table <- function(x, ..., .call = NULL) {
  df <- frequency_table_data(x)
  rows <- ms_rows_from_df(df)
  columns <- ms_table_columns_from_df(df)

  ms_build_envelope(
    type       = "frequency_table",
    type_label = "Frequency table",
    call       = ms_unsupported_call(.call),
    fields     = list(
      table_type   = "frequency_table",
      rows         = rows,
      columns      = columns,
      total        = sum(df$n),
      n_dimensions = length(dim(x))
    ),
    raw_output = ms_capture_output(x),
    card_kind  = "table"
  )
}

#' @rdname mellio_payload
#' @export
mellio_payload.matrix <- function(x, ..., .call = NULL) {
  is_corr <- nrow(x) == ncol(x) &&
    isTRUE(all.equal(unname(diag(x)), rep(1, nrow(x)), tolerance = 1e-8)) &&
    isTRUE(all.equal(unname(x), unname(t(x)), tolerance = 1e-8))

  if (!is_corr) {
    df <- as.data.frame(x, stringsAsFactors = FALSE)
    row_names <- rownames(x)
    if (!is.null(row_names) && length(row_names) == nrow(df)) {
      df <- cbind(Row = row_names, df, stringsAsFactors = FALSE)
    }
    return(ms_build_envelope(
      type       = "matrix_table",
      type_label = "Matrix",
      call       = ms_unsupported_call(.call),
      fields     = list(
        table_type = "matrix",
        rows       = ms_rows_from_df(df),
        columns    = ms_table_columns_from_df(df),
        n_rows     = nrow(x),
        n_cols     = ncol(x)
      ),
      raw_output = ms_capture_output(x),
      card_kind  = "table"
    ))
  }

  var_names <- rownames(x) %||% paste0("V", seq_len(nrow(x)))
  col_names <- colnames(x) %||% var_names
  if (length(var_names) != nrow(x)) var_names <- paste0("V", seq_len(nrow(x)))
  if (length(col_names) != ncol(x)) col_names <- paste0("V", seq_len(ncol(x)))

  rows <- lapply(seq_len(nrow(x)), function(i) {
    values <- as.list(stats::setNames(as.numeric(x[i, ]), col_names))
    c(list(Variable = var_names[[i]]), values)
  })
  columns <- c(
    list(list(key = "Variable", label = "Variable", format = "text")),
    lapply(col_names, function(nm) {
      list(key = nm, label = nm, format = "bounded")
    })
  )

  ms_build_envelope(
    type       = "correlation_matrix",
    type_label = "Correlation matrix",
    call       = ms_unsupported_call(.call),
    fields     = list(
      table_type = "correlation_matrix",
      rows       = rows,
      columns    = columns,
      variables  = var_names,
      n_vars     = nrow(x)
    ),
    raw_output = ms_capture_output(x),
    card_kind  = "table",
    figure_data = list(
      correlation_heatmap = ms_correlation_heatmap_data(
        r = x,
        variables = var_names
      ),
      correlation_forest = ms_correlation_forest_data(
        r = x,
        variables = var_names
      )
    )
  )
}

#' @rdname mellio_payload
#' @export
mellio_payload.melliotab <- function(x, ..., .call = NULL) {
  rows <- ms_rows_from_df(x$data)
  columns <- lapply(names(x$data), function(nm) {
    list(key = nm, label = nm, format = "text")
  })

  ms_build_envelope(
    type       = "styled_table",
    type_label = x$title %||% "Table",
    call       = ms_unsupported_call(.call),
    fields     = list(
      table_type     = "styled_table",
      rows           = rows,
      columns        = columns,
      title          = x$title,
      number         = x$number,
      note           = x$note,
      source_note    = x$source,
      spanners       = x$spanners,
      section_titles = x$section_titles,
      style          = x$style
    ),
    raw_output = table_to_tsv(x$data),
    card_kind  = "table"
  )
}

ms_null_hint <- function(.call = NULL) {
  call_txt <- if (!is.null(.call) && length(.call) > 0L && !is.na(.call)) {
    paste(.call, collapse = " ")
  } else {
    ""
  }
  call_txt <- trimws(call_txt)

  hints <- c(
    "Nothing to send to Mellio.",
    "The R code printed output but returned NULL."
  )

  if (grepl("\\b(plot|hist|barplot|boxplot|pairs|qqplot|curve|image|contour|persp|stripchart|dotchart)\\s*\\(", call_txt)) {
    hints <- c(
      hints,
      "For base R plots, draw the plot first, then call `mellio_capture(...)`, for example: `plot(x, y); mellio_capture(title = \"My figure\")`."
    )
  } else if (grepl("\\baov\\s*\\(|\\banova\\s*\\(", call_txt)) {
    hints <- c(
      hints,
      "For ANOVA, pass an object: `fit <- aov(...); mellio_open(anova(fit))`, or return a data.frame with term, df1, df2, F, and p."
    )
  } else if (grepl("\\bcor\\.test\\s*\\(", call_txt)) {
    hints <- c(
      hints,
      "For several correlations, return a data.frame with variable, r, and p, or use psych::corr.test()."
    )
  } else if (grepl("^for\\s*\\(", call_txt) ||
             grepl("\\bcat\\s*\\(", call_txt) ||
             grepl("\\bprint\\s*\\(", call_txt)) {
    hints <- c(
      hints,
      "Return a supported result object or data.frame; use mellio_open(melliotab(x)) for manuscript tables."
    )
  } else {
    hints <- c(
      hints,
      "Return a supported result object or data.frame."
    )
  }

  paste(hints, collapse = "\n")
}

# ── Internal helpers ──────────────────────────────────────────────────

ms_rows_from_df <- function(df) {
  lapply(seq_len(nrow(df)), function(i) {
    row <- lapply(df[i, , drop = FALSE], function(value) {
      value <- value[[1L]]
      if (is.factor(value)) return(as.character(value))
      if (inherits(value, c("Date", "POSIXct", "POSIXlt"))) return(as.character(value))
      if (is.atomic(value) && length(value) == 1L) return(value)
      as.character(value)
    })
    names(row) <- names(df)
    row
  })
}

ms_table_columns_from_df <- function(df) {
  col_types <- detect_column_types(names(df))
  lapply(seq_along(df), function(i) {
    list(
      key = names(df)[[i]],
      label = names(df)[[i]],
      format = ms_table_column_format(col_types[[i]], df[[i]], names(df)[[i]])
    )
  })
}

ms_table_column_format <- function(col_type, values, name) {
  lname <- tolower(trimws(name %||% ""))
  if (identical(col_type, "pvalue")) return("pvalue")
  if (identical(col_type, "integer") || identical(lname, "n")) return("integer")
  if (identical(col_type, "statistic")) return("statistic")
  if (identical(col_type, "estimate")) return("number")
  if (is.numeric(values) || is.integer(values)) return("number")
  "text"
}

ms_partial_payload <- function(x, cls, printed, call_str) {
  ms_partial_broom_payload(x, cls = cls, printed = printed,
                           call_str = call_str)
}

ms_partial_broom_payload <- function(x, cls, printed, call_str) {
  if (!requireNamespace("broom", quietly = TRUE)) return(NULL)

  tidy_df <- tryCatch(
    suppressWarnings(broom::tidy(x)),
    error = function(e) NULL
  )
  if (!is.data.frame(tidy_df) || nrow(tidy_df) < 1L || ncol(tidy_df) < 1L) {
    return(NULL)
  }

  tidy_df <- ms_partial_clean_df(tidy_df)
  glance_df <- tryCatch(
    suppressWarnings(broom::glance(x)),
    error = function(e) NULL
  )
  glance <- ms_partial_glance_fields(glance_df)

  class_label <- paste(as.character(cls), collapse = "/")
  fields <- list(
    table_type = "partial_broom_summary",
    parse_state = "partial",
    adapter = "broom::tidy",
    class = I(as.character(cls)),
    rows = ms_rows_from_df(tidy_df),
    columns = ms_table_columns_from_df(tidy_df),
    printed = printed$text,
    truncated = printed$truncated,
    message = paste(
      "Mellio did not have a dedicated adapter for this R object, so it",
      "created a partial table with broom::tidy(). Verify the table against",
      "the source output before reporting it."
    ),
    table_note = paste0(
      "Partially parsed from an object of class ", class_label,
      " using broom::tidy(); verify against the attached source output."
    )
  )
  if (length(glance) > 0L) {
    fields$glance <- glance
    fields$glance_adapter <- "broom::glance"
  }

  ms_build_envelope(
    type = "partial_broom_summary",
    type_label = "Partially parsed R result",
    call = call_str,
    fields = fields,
    raw_output = printed$text,
    packages = ms_packages_basic(extras = "broom"),
    card_kind = "table"
  )
}

ms_partial_clean_df <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- ms_partial_clean_names(names(df))
  for (nm in names(df)) {
    value <- df[[nm]]
    if (is.factor(value)) {
      df[[nm]] <- as.character(value)
    } else if (inherits(value, c("Date", "POSIXct", "POSIXlt"))) {
      df[[nm]] <- as.character(value)
    } else if (is.list(value)) {
      df[[nm]] <- vapply(value, ms_partial_cell_text, character(1))
    }
  }
  df
}

ms_partial_clean_names <- function(names) {
  names <- as.character(names %||% character(0))
  missing <- !nzchar(names) | is.na(names)
  if (any(missing)) names[missing] <- paste0("column_", which(missing))
  make.unique(names)
}

ms_partial_cell_text <- function(value) {
  if (is.null(value) || length(value) == 0L) return(NA_character_)
  if (is.atomic(value)) return(paste(as.character(value), collapse = ", "))
  paste(utils::capture.output(utils::str(value, give.attr = FALSE)), collapse = " ")
}

ms_partial_glance_fields <- function(df) {
  if (!is.data.frame(df) || nrow(df) < 1L || ncol(df) < 1L) return(list())
  df <- ms_partial_clean_df(df[1L, , drop = FALSE])
  row <- ms_rows_from_df(df)[[1L]]
  row[!vapply(row, function(x) length(x) == 1L && is.na(x), logical(1))]
}

# Generate a Result Card ID: "rs_" + 8 base36 chars.
# Uses the global RNG; that's fine for v0 — collision risk is ~10^-6
# per million cards. CSPRNG-backed in v1 if/when we need it.
ms_result_id <- function() {
  pool <- c(letters, as.character(0:9))
  paste0("rs_", paste(sample(pool, 8, replace = TRUE), collapse = ""))
}

# ISO 8601 UTC timestamp, no microseconds.
ms_now_iso8601 <- function() {
  format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

# Capture the printed form of an R object for the raw_output field.
ms_capture_output <- function(x) {
  out <- utils::capture.output(print(x))
  paste(out, collapse = "\n")
}

ms_capture_output_safe <- function(x, max_chars = 120000L) {
  cls <- class(x)
  text <- tryCatch(
    paste(utils::capture.output(print(x)), collapse = "\n"),
    error = function(e) paste(
      "(could not print object of class",
      paste(cls, collapse = "/"), ")"
    )
  )
  text <- as.character(text %||% "")
  truncated <- FALSE
  if (nchar(text, type = "chars", allowNA = FALSE, keepNA = FALSE) > max_chars) {
    text <- paste0(
      substr(text, 1L, max_chars),
      "\n\n[Output truncated by Mellio; copy the original R console output if you need the rest.]"
    )
    truncated <- TRUE
  }
  list(text = text, truncated = truncated)
}

ms_unsupported_call <- function(.call = NULL) {
  if (!is.null(.call) && length(.call) > 0L && !is.na(.call)) {
    return(trimws(gsub("\\s+", " ", paste(.call, collapse = " "))))
  }
  NA_character_
}

ms_unsupported_suggestions <- function(cls, call_str, printed_text = "") {
  cls <- tolower(as.character(cls %||% character(0)))
  call_str <- trimws(as.character(call_str %||% ""))
  if (!length(call_str) || is.na(call_str)) call_str <- ""
  printed_text <- paste(as.character(printed_text %||% ""), collapse = "\n")
  if (!length(printed_text) || is.na(printed_text)) printed_text <- ""
  suggestions <- list()

  add <- function(title, detail, code = NULL) {
    item <- list(title = title, detail = detail)
    if (!is.null(code) && nzchar(trimws(as.character(code)))) {
      item$code <- trimws(as.character(code))
    }
    key <- paste(item$title, item$code %||% "", sep = "\r")
    existing <- vapply(suggestions, function(suggestion) {
      paste(suggestion$title, suggestion$code %||% "", sep = "\r")
    }, character(1))
    existing_codes <- vapply(suggestions, function(suggestion) {
      trimws(as.character(suggestion$code %||% ""))
    }, character(1))
    if (!key %in% existing &&
        (!nzchar(item$code %||% "") || !item$code %in% existing_codes)) {
      suggestions[[length(suggestions) + 1L]] <<- item
    }
  }

  summary_target <- ms_summary_call_target(call_str)
  if (nzchar(summary_target)) {
    add(
      "Try the original model object",
      "Mellio can often extract richer structured results from the fitted object than from its printed summary.",
      paste0("mellio_open(", summary_target, ")")
    )
  }

  if (any(cls %in% c(
    "summary.lm", "summary.glm", "summary.polr", "summary.clm",
    "summary.lme", "summary.merMod", "summary.glmerMod", "summary.lmerMod",
    "summary.aovlist"
  ))) {
    add(
      "Pass the fitted model, not summary()",
      "The fitted model keeps formula, data, and model metadata that a summary object may discard.",
      if (nzchar(summary_target)) paste0("mellio_open(", summary_target, ")") else "mellio_open(fit)"
    )
  }

  if (ms_unsupported_looks_like_anova(printed_text)) {
    add(
      "Try the ANOVA/model object",
      "This printed output looks like an ANOVA table. Mellio usually does better with the original aov/lm model or anova(fit).",
      "mellio_open(anova(fit))"
    )
  }

  if (any(cls %in% c("emm_list", "emmgrid")) ||
      ms_unsupported_looks_like_emmeans(printed_text)) {
    add(
      "Try the emmeans contrasts",
      "For estimated marginal means pairwise output, pass the contrasts grid or use pairs() on the emmeans grid.",
      "mellio_open(result$contrasts)"
    )
  }

  suggestions[seq_len(min(length(suggestions), 3L))]
}

ms_summary_call_target <- function(call_str) {
  call_str <- trimws(as.character(call_str %||% ""))
  if (!length(call_str) || is.na(call_str)) return("")
  if (!nzchar(call_str)) return("")
  match <- regexec("^summary\\(([^,\\)]+)(?:[,\\)].*)$", call_str)
  parts <- regmatches(call_str, match)[[1]]
  if (length(parts) < 2L) return("")
  target <- trimws(parts[[2]])
  if (grepl("^[A-Za-z.][A-Za-z0-9._]*$", target)) target else ""
}

ms_unsupported_looks_like_anova <- function(text) {
  text <- paste(as.character(text %||% ""), collapse = "\n")
  if (!length(text) || is.na(text)) return(FALSE)
  has_cols <- isTRUE(grepl("\\bDf\\b", text)) &&
    isTRUE(grepl("\\bSum Sq\\b|\\bSum of Sq\\b", text)) &&
    isTRUE(grepl("\\bPr\\(>F\\)", text))
  has_cols && isTRUE(grepl("\\bResiduals?\\b", text, ignore.case = TRUE))
}

ms_unsupported_looks_like_emmeans <- function(text) {
  text <- paste(as.character(text %||% ""), collapse = "\n")
  if (!length(text) || is.na(text)) return(FALSE)
  isTRUE(grepl("\\$emmeans|\\$contrasts", text)) &&
    isTRUE(grepl("\\bemmean\\b|\\bcontrast\\b|\\bt\\.ratio\\b|\\bp\\.value\\b", text))
}

# Deparse an R call to a single-line string. Falls back to method name
# from htest objects when no call is available (interactive use).
ms_deparse_call <- function(call_obj) {
  if (is.null(call_obj)) return(NA_character_)
  txt <- paste(deparse(call_obj, width.cutoff = 500L), collapse = " ")
  trimws(gsub("\\s+", " ", txt))
}

# Coerce numeric values that JSON can't represent to NA (later -> null).
# Inf, -Inf, NaN all become NA. Used before JSON serialisation.
ms_safe_numeric <- function(x) {
  if (is.null(x)) return(NA_real_)
  if (length(x) == 0) return(NA_real_)
  ifelse(is.finite(x), as.numeric(x), NA_real_)
}

# Build a v0.1 envelope from extractor outputs.
# The extractor supplies type, type_label, fields. The envelope adds
# schema_version, result_id, created_at, card_kind, raw_output, and
# provenance/packages (optional but populated by default in v0.1+).
ms_build_envelope <- function(type, type_label, call, fields,
                              raw_output, packages = NULL,
                              provenance = NULL,
                              card_kind = "inline",
                              figure_data = NULL,
                              available_figures = NULL,
                              available_tables = NULL) {
  payload <- list(
    schema_version = "0.1",
    result_id = ms_result_id(),
    card_kind = card_kind,
    type = type,
    type_label = type_label,
    call = call,
    created_at = ms_now_iso8601(),
    fields = fields,
    raw_output = raw_output,
    provenance = provenance %||% ms_provenance_basic(),
    packages   = packages   %||% ms_packages_basic()
  )
  if (!is.null(figure_data)) payload$figure_data <- figure_data
  if (!is.null(available_figures)) {
    payload$metadata <- payload$metadata %||% list()
    payload$metadata$available_figures <- available_figures
  }
  if (!is.null(available_tables)) {
    payload$metadata <- payload$metadata %||% list()
    payload$metadata$available_tables <- available_tables
  }
  payload <- ms_attach_available_figures(payload)
  class(payload) <- c("mellio_payload", "list")
  payload
}

ms_attach_available_figures <- function(payload) {
  payload <- ms_attach_existing_figure_metadata(payload)
  payload <- ms_attach_interaction_plot_figure(payload)
  payload <- ms_attach_coefficient_plot_figure(payload)
  payload
}

ms_attach_existing_figure_metadata <- function(payload) {
  figure_data <- payload$figure_data %||% list()
  if (!is.null(figure_data$correlation_heatmap)) {
    payload <- ms_add_available_figure(
      payload,
      type = "correlation_heatmap",
      label = "Correlation heatmap"
    )
  }
  if (!is.null(figure_data$correlation_forest)) {
    payload <- ms_add_available_figure(
      payload,
      type = "correlation_forest",
      label = "Correlation forest plot",
      default = FALSE
    )
  }
  if (!is.null(figure_data$adjusted_means)) {
    means_source <- as.character(figure_data$adjusted_means$source %||% "")
    payload <- ms_add_available_figure(
      payload,
      type = "adjusted_means",
      label = if (identical(means_source, "glmer_emmeans")) "Predicted probabilities" else "Means plot"
    )
  }
  # Register AFTER adjusted_means so the canonical means plot stays the default
  # (first-registered figure wins; see ms_add_available_figure). Estimation is opt-in.
  if (!is.null(figure_data$estimation_plot)) {
    payload <- ms_add_available_figure(
      payload,
      type = "estimation_plot",
      label = "Estimation plot"
    )
  }
  if (!is.null(figure_data$scatter_plot)) {
    payload <- ms_add_available_figure(
      payload,
      type = "scatter_plot",
      label = "Scatter plot"
    )
  }
  if (!is.null(figure_data$paired_difference_plot)) {
    payload <- ms_add_available_figure(
      payload,
      type = "paired_difference_plot",
      label = "Paired plot"
    )
  }
  if (!is.null(figure_data$one_sample_mean_plot)) {
    payload <- ms_add_available_figure(
      payload,
      type = "one_sample_mean_plot",
      label = "Mean plot"
    )
  }
  if (!is.null(figure_data$nonparametric_group_plot)) {
    payload <- ms_add_available_figure(
      payload,
      type = "nonparametric_group_plot",
      label = "Distribution plot"
    )
  }
  if (!is.null(figure_data$coefficient_plot)) {
    payload <- ms_add_available_figure(
      payload,
      type = "coefficient_plot",
      label = "Coefficient plot"
    )
  }
  if (!is.null(figure_data$km_curve)) {
    payload <- ms_add_available_figure(
      payload,
      type = "km_curve",
      label = "Kaplan-Meier curve"
    )
  }
  if (!is.null(figure_data$structural_path_diagram)) {
    diagram_type <- as.character(figure_data$structural_path_diagram$model_type %||% "")
    payload <- ms_add_available_figure(
      payload,
      type = "structural_path_diagram",
      label = if (identical(diagram_type, "lavaan_cfa")) "CFA path diagram" else "Path diagram"
    )
  }
  if (!is.null(figure_data$cfa_path_diagram)) {
    payload <- ms_add_available_figure(
      payload,
      type = "cfa_path_diagram",
      label = "CFA path diagram"
    )
  }
  if (!is.null(figure_data$sem_path_diagram)) {
    payload <- ms_add_available_figure(
      payload,
      type = "sem_path_diagram",
      label = "SEM path diagram"
    )
  }
  payload
}

ms_attach_interaction_plot_figure <- function(payload) {
  figure_data <- payload$figure_data %||% list()
  if (is.null(figure_data$interaction_plot)) return(payload)
  kind <- as.character(figure_data$interaction_plot$interaction_kind %||% "")

  payload <- ms_add_available_figure(
    payload,
    type = "interaction_plot",
    label = if (identical(kind, "continuous_main_effect")) "Effect plot" else "Interaction plot",
    default = length(payload$metadata$available_figures %||% list()) == 0L
  )
  payload
}

ms_add_available_figure <- function(payload, type, label, default = NULL) {
  payload$metadata <- payload$metadata %||% list()
  figures <- payload$metadata$available_figures %||% list()
  existing <- vapply(figures, function(figure) {
    is.list(figure) && identical(figure$type %||% NULL, type)
  }, logical(1))
  if (any(existing)) return(payload)

  is_default <- if (is.null(default)) length(figures) == 0L else isTRUE(default)
  figures[[length(figures) + 1L]] <- list(
    type = type,
    label = label,
    default = is_default
  )
  payload$metadata$available_figures <- figures
  payload
}

ms_attach_coefficient_plot_figure <- function(payload) {
  fields <- payload$fields %||% list()
  coefficients <- fields$coefficients
  if (!is.list(coefficients) || length(coefficients) == 0L) {
    return(payload)
  }

  usable <- Filter(function(row) {
    is.list(row) &&
      !is.null(row$term) &&
      !is.null(row$estimate) &&
      !is.na(ms_safe_numeric(row$estimate))
  }, coefficients)
  if (length(usable) == 0L) return(payload)

  coefficient_scale <- ms_coefficient_plot_scale(fields, usable)
  figure_coefficients <- ms_coefficient_plot_rows(fields, usable, coefficient_scale)
  figure_coefficients <- Filter(function(row) {
    is.list(row) &&
      !is.null(row$term) &&
      !is.null(row$estimate) &&
      !is.na(ms_safe_numeric(row$estimate))
  }, figure_coefficients)
  if (length(figure_coefficients) == 0L) return(payload)
  if (ms_coefficient_plot_suppressed_by_diagnostics(
    fields,
    figure_coefficients,
    coefficient_scale
  )) {
    return(payload)
  }

  payload$figure_data <- payload$figure_data %||% list()
  payload$figure_data$coefficient_plot <- list(
    coefficients       = figure_coefficients,
    coefficient_scale  = coefficient_scale,
    estimate_label     = ms_coefficient_plot_estimate_label(fields, figure_coefficients, coefficient_scale),
    statistic_label    = fields$statistic_label %||% NULL,
    outcome            = fields$outcome %||% NULL,
    predictor          = fields$predictor %||% NULL,
    model_family       = fields$model_family %||% NULL,
    model_link         = fields$model_link %||% NULL,
    n                  = fields$n %||% NULL,
    events             = fields$events %||% NULL,
    model_fit          = fields$model_fit %||% NULL,
    groups             = fields$groups %||% NULL,
    default_show_intercept = FALSE
  )
  payload <- ms_add_available_figure(
    payload,
    type = "coefficient_plot",
    label = "Coefficient plot",
    default = length(payload$metadata$available_figures %||% list()) == 0L
  )
  payload
}

ms_coefficient_plot_scale <- function(fields, rows) {
  fields <- fields %||% list()
  scale <- as.character(fields$coefficient_scale %||% "")
  if (nzchar(scale)) return(scale)

  labels <- vapply(rows, function(row) {
    as.character(row$estimate_name %||% "")
  }, character(1))
  if (any(labels == "HR")) return("hazard_ratio")
  if (any(labels == "OR")) return("odds_ratio")
  if (ms_coefficient_plot_is_logistic(fields)) return("odds_ratio")
  "raw"
}

ms_coefficient_plot_is_logistic <- function(fields) {
  family <- tolower(trimws(as.character(fields$model_family %||% "")))
  link <- tolower(trimws(as.character(fields$model_link %||% "")))
  family %in% c("binomial", "quasibinomial") &&
    link %in% c("logit", "logistic")
}

ms_coefficient_plot_rows <- function(fields, rows, scale) {
  if (!ms_coefficient_plot_should_exponentiate(fields, rows, scale)) {
    return(rows)
  }
  lapply(rows, ms_coefficient_plot_exp_row)
}

ms_coefficient_plot_should_exponentiate <- function(fields, rows, scale) {
  if (!identical(scale, "odds_ratio")) return(FALSE)
  existing_scale <- as.character(fields$coefficient_scale %||% "")
  if (existing_scale %in% c("odds_ratio", "proportional_odds", "hazard_ratio")) {
    return(FALSE)
  }
  already_or <- any(vapply(rows, function(row) {
    identical(as.character(row$estimate_name %||% ""), "OR")
  }, logical(1)))
  !already_or && ms_coefficient_plot_is_logistic(fields)
}

ms_coefficient_plot_exp_row <- function(row) {
  estimate <- ms_safe_numeric(row$estimate)
  if (is.na(estimate)) return(row)
  row$log_estimate <- estimate
  row$estimate <- ms_exp_safe_numeric(estimate)
  row$estimate_name <- "OR"
  if (!is.null(row$ci_lower)) row$ci_lower <- ms_exp_safe_numeric(row$ci_lower)
  if (!is.null(row$ci_upper)) row$ci_upper <- ms_exp_safe_numeric(row$ci_upper)
  row$std_estimate <- NULL
  row$std_std_error <- NULL
  row$std_ci_lower <- NULL
  row$std_ci_upper <- NULL
  row
}

ms_coefficient_plot_suppressed_by_diagnostics <- function(fields, rows, scale) {
  if (!as.character(scale %||% "") %in% c("odds_ratio", "proportional_odds", "hazard_ratio")) {
    return(FALSE)
  }
  warnings <- fields$model_warnings %||% list()
  warning_types <- vapply(warnings, function(row) {
    if (!is.list(row)) return("")
    as.character(row$type %||% "")
  }, character(1))
  if (!any(warning_types %in% c("separation_or_boundary", "boundary_fit"))) {
    return(FALSE)
  }
  rows <- Filter(function(row) {
    is.list(row) &&
      !grepl("^\\(?intercept\\)?$", as.character(row$term %||% ""),
             ignore.case = TRUE)
  }, rows)
  if (length(rows) == 0L) return(FALSE)
  any(vapply(rows, ms_coefficient_ratio_row_unstable, logical(1)))
}

ms_coefficient_ratio_row_unstable <- function(row) {
  estimate <- ms_safe_numeric(row$estimate)
  log_estimate <- ms_safe_numeric(row$log_estimate)
  ci_lower <- ms_safe_numeric(row$ci_lower)
  ci_upper <- ms_safe_numeric(row$ci_upper)
  if (is.na(estimate) || estimate <= 0) return(TRUE)
  if (estimate > 1e6 || estimate < 1e-6) return(TRUE)
  if (!is.na(log_estimate) && abs(log_estimate) > 8 &&
      (is.na(ci_lower) || is.na(ci_upper))) {
    return(TRUE)
  }
  FALSE
}

ms_exp_safe_numeric <- function(x) {
  value <- ms_safe_numeric(x)
  if (is.na(value)) return(value)
  ms_safe_numeric(exp(value))
}

ms_coefficient_plot_estimate_label <- function(fields, rows, scale = NULL) {
  for (row in rows) {
    label <- row$estimate_name %||% NULL
    if (!is.null(label) && length(label) == 1L && nzchar(label)) {
      return(as.character(label))
    }
  }
  scale <- scale %||% fields$coefficient_scale %||% ""
  if (identical(scale, "odds_ratio") ||
      identical(scale, "proportional_odds")) {
    return("OR")
  }
  if (identical(scale, "hazard_ratio")) return("HR")
  "B"
}

ms_correlation_heatmap_data <- function(r, p = NULL, n = NULL,
                                        variables = NULL, method = NULL,
                                        missing = NULL, adjust = NULL) {
  r_mat <- ms_square_numeric_matrix(r)
  if (is.null(r_mat)) return(NULL)

  vars <- variables %||% rownames(r_mat) %||% colnames(r_mat)
  vars <- as.character(vars)
  if (length(vars) != nrow(r_mat) || any(!nzchar(vars))) {
    vars <- paste0("V", seq_len(nrow(r_mat)))
  }

  out <- list(
    variables = vars,
    r = ms_numeric_matrix_to_list(r_mat),
    default_triangle = "lower"
  )

  p_mat <- ms_square_numeric_matrix(p, n = nrow(r_mat))
  if (!is.null(p_mat)) out$p <- ms_numeric_matrix_to_list(p_mat)

  n_value <- ms_heatmap_n_value(n, nrow(r_mat))
  if (!is.null(n_value)) out$n <- n_value

  method <- ms_correlation_method_label(method)
  if (!is.null(method)) out$method <- method
  missing <- trimws(as.character(missing %||% ""))
  if (length(missing) == 1L && nzchar(missing)) out$missing <- missing
  adjust <- trimws(as.character(adjust %||% ""))
  if (length(adjust) == 1L && nzchar(adjust) && !identical(adjust, "none")) {
    out$adjust <- adjust
  }
  out
}

ms_correlation_forest_data <- function(r = NULL, p = NULL, n = NULL,
                                       variables = NULL, method = NULL,
                                       adjust = NULL, pairs = NULL) {
  if (is.null(pairs)) {
    r_mat <- ms_square_numeric_matrix(r)
    if (is.null(r_mat) || nrow(r_mat) < 2L) return(NULL)
    vars <- variables %||% rownames(r_mat) %||% colnames(r_mat)
    vars <- as.character(vars)
    if (length(vars) != nrow(r_mat) || any(!nzchar(vars))) {
      vars <- paste0("V", seq_len(nrow(r_mat)))
    }
    p_mat <- ms_square_numeric_matrix(p, n = nrow(r_mat))
    pairs <- list()
    for (i in seq_len(nrow(r_mat) - 1L)) {
      for (j in seq.int(i + 1L, ncol(r_mat))) {
        row <- list(
          x = vars[[i]],
          y = vars[[j]],
          r = ms_safe_numeric(r_mat[i, j])
        )
        if (!is.null(p_mat)) row$p_value <- ms_safe_numeric(p_mat[i, j])
        n_pair <- ms_corr_n_pair_safe(n, i, j)
        if (!is.na(n_pair)) row$n <- n_pair
        pairs[[length(pairs) + 1L]] <- row
      }
    }
  }

  pairs <- Filter(function(row) {
    is.list(row) &&
      !is.null(row$x) &&
      !is.null(row$y) &&
      !is.null(row$r) &&
      !is.na(ms_safe_numeric(row$r))
  }, pairs)
  if (length(pairs) == 0L) return(NULL)

  vars <- unique(unlist(lapply(pairs, function(row) {
    c(as.character(row$x %||% ""), as.character(row$y %||% ""))
  }), use.names = FALSE))
  vars <- vars[nzchar(vars)]

  out <- list(
    variables = vars,
    pairs = pairs
  )
  method <- ms_correlation_method_label(method)
  if (!is.null(method)) out$method <- method
  adjust <- trimws(as.character(adjust %||% ""))
  if (length(adjust) == 1L && nzchar(adjust) && !identical(adjust, "none")) {
    out$adjust <- adjust
  }
  n_value <- ms_heatmap_n_value(n, NA_integer_)
  if (!is.null(n_value) && length(n_value) == 1L) out$n <- n_value
  out
}

ms_square_numeric_matrix <- function(x, n = NULL) {
  if (is.null(x)) return(NULL)
  m <- suppressWarnings(as.matrix(x))
  if (!is.numeric(m) && !is.integer(m)) return(NULL)
  if (nrow(m) == 0L || nrow(m) != ncol(m)) return(NULL)
  if (!is.null(n) && nrow(m) != n) return(NULL)
  storage.mode(m) <- "numeric"
  m[] <- ms_safe_numeric(m)
  if (all(is.na(m))) return(NULL)
  m
}

ms_numeric_matrix_to_list <- function(m) {
  lapply(seq_len(nrow(m)), function(i) {
    as.list(unname(ms_safe_numeric(m[i, ])))
  })
}

ms_heatmap_n_value <- function(n, size) {
  if (is.null(n)) return(NULL)
  if (length(n) == 1L) {
    value <- ms_safe_numeric(n)
    if (is.na(value)) return(NULL)
    return(as.integer(value))
  }
  n_mat <- suppressWarnings(as.matrix(n))
  if (is.na(size)) return(NULL)
  if (nrow(n_mat) != size || ncol(n_mat) != size) return(NULL)
  storage.mode(n_mat) <- "numeric"
  ms_numeric_matrix_to_list(n_mat)
}

ms_corr_n_pair_safe <- function(n, i, j) {
  if (is.null(n)) return(NA_integer_)
  if (length(n) == 1L) {
    value <- ms_safe_numeric(n)
    return(if (is.na(value)) NA_integer_ else as.integer(value))
  }
  m <- suppressWarnings(as.matrix(n))
  if (i > nrow(m) || j > ncol(m)) return(NA_integer_)
  value <- ms_safe_numeric(m[i, j])
  if (is.na(value)) NA_integer_ else as.integer(value)
}

ms_correlation_method_label <- function(method) {
  method <- trimws(tolower(as.character(method %||% "")))
  if (length(method) != 1L || !nzchar(method)) return(NULL)
  if (method %in% c("pearson", "pearson's", "pearsons")) return("Pearson")
  if (method %in% c("spearman", "spearman's", "spearmans")) return("Spearman")
  if (method %in% c("kendall", "kendall's", "kendalls", "kendall_tau")) {
    return("Kendall")
  }
  sentence_case <- paste0(toupper(substr(method, 1L, 1L)), substr(method, 2L, nchar(method)))
  sentence_case
}

ms_provenance_mode <- function() {
  opt <- getOption("mellio.provenance", "standard")
  if (isFALSE(opt) || identical(opt, "none")) return("none")
  if (isTRUE(opt)) return("full")
  opt <- tolower(as.character(opt[[1]] %||% "standard"))
  if (opt %in% c("standard", "minimal")) return("standard")
  if (opt %in% c("full", "local")) return("full")
  "standard"
}

ms_provenance_basic <- function() {
  mode <- ms_provenance_mode()
  if (identical(mode, "none")) return(NULL)

  prov <- list(
    r_version      = R.version.string,
    platform       = R.version$platform,
    mellio_version = ms_package_version("mellio")
  )

  if (identical(mode, "full")) {
    prov$working_dir <- ms_safe_getwd()
    prov$sender <- ms_sender()
    git <- ms_git_state()
    if (!is.null(git)) prov$git <- git
    script <- ms_script_provenance()
    if (!is.null(script)) prov$script <- script
  }

  prov
}

ms_sender <- function() {
  info <- tryCatch(Sys.info(), error = function(e) NULL)
  if (is.null(info)) return(NULL)
  user <- if (!is.null(info[["user"]])) as.character(info[["user"]]) else NA_character_
  host <- if (!is.null(info[["nodename"]])) as.character(info[["nodename"]]) else NA_character_
  if (is.na(user) && is.na(host)) return(NULL)
  list(user = user, host = host)
}

# Get working directory, returning NA_character_ if the call fails.
ms_safe_getwd <- function() {
  tryCatch(getwd(), error = function(e) NA_character_)
}

ms_git_state <- function() {
  if (!ms_git_available()) return(NULL)
  if (!ms_in_git_repo()) return(NULL)
  list(
    commit = ms_git_run(c("rev-parse", "HEAD")),
    branch = ms_git_run(c("rev-parse", "--abbrev-ref", "HEAD")),
    dirty  = ms_git_dirty()
  )
}

ms_git_available <- function() {
  res <- tryCatch(
    suppressWarnings(system2("git", "--version", stdout = FALSE, stderr = FALSE)),
    error = function(e) 127L
  )
  identical(res, 0L)
}

ms_in_git_repo <- function() {
  out <- ms_git_run(c("rev-parse", "--is-inside-work-tree"))
  identical(out, "true")
}

ms_git_run <- function(args) {
  tryCatch({
    res <- suppressWarnings(
      system2("git", args, stdout = TRUE, stderr = FALSE)
    )
    if (length(res) == 0L) return(NA_character_)
    trimws(res[1])
  }, error = function(e) NA_character_)
}

ms_git_dirty <- function() {
  tryCatch({
    res <- suppressWarnings(
      system2("git", c("status", "--porcelain"),
              stdout = TRUE, stderr = FALSE)
    )
    length(res) > 0L
  }, error = function(e) NA)
}

ms_script_provenance <- function() {
  src <- ms_script_source()
  if (is.null(src)) return(NULL)
  out <- list(file = src$file, line = src$line)
  if (file.exists(src$file) &&
      requireNamespace("digest", quietly = TRUE)) {
    content <- tryCatch(
      readLines(src$file, warn = FALSE),
      error = function(e) NULL
    )
    if (!is.null(content)) {
      out$hash <- digest::digest(content, algo = "sha1")
    }
  }
  out
}

ms_mellio_package_dir <- function() {
  norm <- function(p) {
    if (is.null(p) || length(p) != 1L || !nzchar(p) || is.na(p)) {
      return(NA_character_)
    }
    tryCatch(normalizePath(p, winslash = "/", mustWork = FALSE),
             error = function(e) NA_character_)
  }

  ns_path <- tryCatch(
    getNamespaceInfo("mellio", "path"),
    error = function(e) NULL
  )
  d <- norm(ns_path)
  if (!is.na(d)) return(d)

  sys_path <- tryCatch(system.file(package = "mellio"),
                       error = function(e) "")
  d <- norm(sys_path)
  if (!is.na(d)) return(d)

  self_file <- tryCatch({
    sr <- attr(sys.function(), "srcref")
    if (is.null(sr)) NULL else attr(sr, "srcfile")$filename
  }, error = function(e) NULL)
  d <- norm(self_file)
  if (!is.na(d)) return(norm(dirname(dirname(d))))

  NA_character_
}

ms_path_inside <- function(fname, dir) {
  if (is.null(dir) || length(dir) != 1L || is.na(dir) || !nzchar(dir)) {
    return(FALSE)
  }
  norm <- tryCatch(
    normalizePath(fname, winslash = "/", mustWork = FALSE),
    error = function(e) fname
  )
  prefix <- if (endsWith(dir, "/")) dir else paste0(dir, "/")
  startsWith(norm, prefix)
}

ms_script_source <- function() {
  pkg_dir <- ms_mellio_package_dir()

  calls <- sys.calls()
  for (i in seq_along(calls)) {
    sr <- attr(calls[[i]], "srcref")
    if (is.null(sr)) next
    srcfile <- attr(sr, "srcfile")
    if (is.null(srcfile)) next
    fname <- srcfile$filename %||% ""
    if (!nzchar(fname) || fname == "<text>") next
    if (ms_path_inside(fname, pkg_dir)) next
    return(list(file = fname, line = as.integer(sr[1])))
  }

  args <- commandArgs(trailingOnly = FALSE)
  fa <- args[grepl("^--file=", args)]
  if (length(fa) == 1L) {
    fname <- sub("^--file=", "", fa)
    if (nzchar(fname)) {
      norm <- tryCatch(normalizePath(fname, mustWork = FALSE),
                       error = function(e) fname)
      return(list(file = norm, line = NA_integer_))
    }
  }

  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    is_avail <- isTRUE(tryCatch(rstudioapi::isAvailable(),
                                error = function(e) FALSE))
    if (is_avail) {
      path <- tryCatch(
        rstudioapi::getSourceEditorContext()$path,
        error = function(e) ""
      )
      if (nzchar(path %||% "")) {
        ctx <- tryCatch(rstudioapi::getActiveDocumentContext(),
                        error = function(e) NULL)
        line <- NA_integer_
        sel <- tryCatch(ctx$selection %||% list(), error = function(e) list())
        if (length(sel) > 0L) {
          start <- tryCatch(sel[[1]]$range$start, error = function(e) NULL)
          line <- ms_rstudio_selection_line(start)
        }
        norm <- tryCatch(normalizePath(path, mustWork = FALSE),
                         error = function(e) path)
        return(list(file = norm, line = line))
      }
    }
  }

  NULL
}

ms_rstudio_selection_line <- function(start) {
  if (is.null(start)) return(NA_integer_)

  value <- NULL
  if (is.list(start)) {
    value <- start[["row"]]
    if (is.null(value) && length(start) > 0L) value <- start[[1L]]
  } else if (is.atomic(start) && length(start) > 0L) {
    if (!is.null(names(start)) && "row" %in% names(start)) {
      value <- start[["row"]]
    } else {
      value <- start[[1L]]
    }
  }

  out <- suppressWarnings(as.integer(value))
  if (length(out) == 0L || is.na(out[[1L]])) NA_integer_ else out[[1L]]
}

ms_data_provenance <- function(df) {
  if (is.null(df)) return(NULL)
  if (!requireNamespace("digest", quietly = TRUE)) return(NULL)
  # model.frame() attaches `terms`/`na.action` with environment refs
  # that vary between fits, so two identical analyses produce different
  # hashes from the raw object. Hash a sorted column-value list to get
  # a stable fingerprint of the data content alone.
  content <- as.list(df)
  content <- content[order(names(content))]
  list(
    hash = digest::digest(content, algo = "sha1"),
    n    = nrow(df)
  )
}

ms_provenance_add_data <- function(provenance, data_provenance) {
  if (is.null(provenance)) return(NULL)
  if (!is.null(data_provenance)) provenance$data <- data_provenance
  provenance
}

ms_packages_basic <- function(extras = NULL) {
  base_pkgs <- list(
    list(name = "R",      version = paste(R.version$major, R.version$minor, sep = ".")),
    list(name = "mellio", version = ms_package_version("mellio")),
    list(name = "stats",  version = "base")
  )
  if (!is.null(extras) && length(extras) > 0) {
    extras <- lapply(extras, function(p) {
      if (is.character(p)) {
        list(name = p, version = ms_package_version(p))
      } else p
    })
    return(c(base_pkgs, extras))
  }
  base_pkgs
}

# Resolve a package version, returning NA_character_ if not installed.
ms_package_version <- function(pkg) {
  tryCatch(
    as.character(utils::packageVersion(pkg)),
    error = function(e) NA_character_
  )
}

#' Serialise a Mellio payload to JSON
#'
#' @param payload A `mellio_payload` object from [mellio_payload()].
#' @param pretty Pretty-print with line breaks (default `FALSE`).
#' @return A length-1 character vector containing the JSON.
#' @export
#' @family R bridge
#' @seealso [mellio_payload()] to build the payload, [mellio_open()] to send it
#'   directly to the Mellio web app.
#' @examples
#' p <- mellio_payload(t.test(extra ~ group, data = sleep))
#' cat(mellio_to_json(p, pretty = TRUE))
mellio_to_json <- function(payload, pretty = FALSE) {
  rlang::check_installed("jsonlite", reason = "to serialise R bridge payloads")
  jsonlite::toJSON(
    payload,
    auto_unbox = TRUE,
    na = "null",
    null = "null",
    digits = NA,
    pretty = pretty
  )
}

#' Print method for Mellio payloads
#'
#' Compact summary of a `mellio_payload` object. Use [mellio_to_json()] to
#' see the full structure or [str()] to inspect the field tree.
#'
#' @param x A `mellio_payload` object.
#' @param ... Unused.
#' @return Invisibly returns `x`.
#' @export
#' @family R bridge
print.mellio_payload <- function(x, ...) {
  cat("<mellio_payload v", x$schema_version, ">\n", sep = "")
  cat("  id:    ", x$result_id, "\n", sep = "")
  cat("  type:  ", x$type, " (", x$type_label, ")\n", sep = "")
  cat("  call:  ", x$call, "\n", sep = "")
  cat("  saved: ", x$created_at, "\n", sep = "")
  invisible(x)
}
