# R bridge + table extractors for processR summary objects.
#
# processR builds many PROCESS-style models through lavaan, but its printed
# summary objects are already reduced effect tables. Mellio treats those
# summaries as table Result Cards and leaves topology-specific diagrams/prose
# to the underlying lavaan fit, where the regression and defined-parameter
# structure can be checked.

#' @rdname mellio_payload
#' @export
mellio_payload.medSummary <- function(x, ..., .call = NULL) {
  ms_processr_med_summary_payload(x, .call = .call)
}

#' @rdname mellio_payload
#' @export
mellio_payload.medSummary2 <- function(x, ..., .call = NULL) {
  ms_processr_med_summary2_payload(x, .call = .call)
}

#' @rdname mellio_payload
#' @export
mellio_payload.modmedSummary <- function(x, ..., .call = NULL) {
  ms_processr_modmed_summary_payload(x, .call = .call)
}

#' @rdname mellio_payload
#' @export
mellio_payload.modmedSummary2 <- function(x, ..., .call = NULL) {
  ms_processr_modmed_summary2_payload(x, .call = .call)
}

#' @export
mellio_open_dispatch.medSummary <- function(x, browse = TRUE, ..., .call = NULL) {
  send_payload_to_stats(mellio_payload(x, ..., .call = .call), browse = browse)
}

#' @export
mellio_open_dispatch.medSummary2 <- function(x, browse = TRUE, ..., .call = NULL) {
  send_payload_to_stats(mellio_payload(x, ..., .call = .call), browse = browse)
}

#' @export
mellio_open_dispatch.modmedSummary <- function(x, browse = TRUE, ..., .call = NULL) {
  send_payload_to_stats(mellio_payload(x, ..., .call = .call), browse = browse)
}

#' @export
mellio_open_dispatch.modmedSummary2 <- function(x, browse = TRUE, ..., .call = NULL) {
  send_payload_to_stats(mellio_payload(x, ..., .call = .call), browse = browse)
}

#' @rdname melliotab
#' @export
melliotab.medSummary <- function(x, ..., section = NULL) {
  melliotab_from_payload(mellio_payload(x, ...), section = section, ...)
}

#' @rdname melliotab
#' @export
melliotab.medSummary2 <- function(x, ..., section = NULL) {
  melliotab_from_payload(mellio_payload(x, ...), section = section, ...)
}

#' @rdname melliotab
#' @export
melliotab.modmedSummary <- function(x, ..., section = NULL) {
  melliotab_from_payload(mellio_payload(x, ...), section = section, ...)
}

#' @rdname melliotab
#' @export
melliotab.modmedSummary2 <- function(x, ..., section = NULL) {
  melliotab_from_payload(mellio_payload(x, ...), section = section, ...)
}

ms_processr_med_summary_payload <- function(x, .call = NULL) {
  df <- as.data.frame(x, stringsAsFactors = FALSE)
  rows <- lapply(seq_len(nrow(df)), function(i) {
    list(
      effect = ms_processr_value(df, i, c("Effect", "effect", "lhs")),
      equation = ms_processr_value(df, i, c("equation", "Equation", "rhs")),
      estimate = ms_processr_numeric(df, i, c("est", "estimate")),
      ci_lower = ms_processr_numeric(df, i, c("ci.lower", "lower", "conf.low")),
      ci_upper = ms_processr_numeric(df, i, c("ci.upper", "upper", "conf.high")),
      p_value = ms_processr_numeric(df, i, c("pvalue", "p.value", "p"))
    )
  })
  rows <- Filter(function(row) !is.na(row$estimate), rows)

  ms_processr_table_payload(
    type = "processr_mediation_summary",
    type_label = "Mediation summary (processR)",
    call = ms_processr_call(.call, "processR::medSummary(...)"),
    rows = rows,
    columns = ms_processr_effect_columns(include_equation = TRUE),
    note = ms_processr_summary_note(x, moderated = FALSE, rows = rows),
    raw_output = ms_capture_output(x)
  )
}

ms_processr_med_summary2_payload <- function(x, .call = NULL) {
  df <- as.data.frame(x, stringsAsFactors = FALSE)
  effects <- attr(x, "effects", exact = TRUE)
  equations <- attr(x, "equations", exact = TRUE)
  if (is.null(effects) || !length(effects)) {
    effects <- unique(sub("^(?:est|lower|upper|p)\\.", "", names(df)[-1L]))
  }

  rows <- list()
  for (i in seq_len(nrow(df))) {
    interval_type <- ms_processr_value(df, i, c("type", "Type"))
    for (j in seq_along(effects)) {
      effect <- as.character(effects[[j]])
      row <- list(
        interval_type = interval_type,
        effect = effect,
        equation = if (!is.null(equations) && length(equations) >= j) {
          as.character(equations[[j]])
        } else NA_character_,
        estimate = ms_processr_numeric(df, i, paste0("est.", effect)),
        ci_lower = ms_processr_numeric(df, i, paste0("lower.", effect)),
        ci_upper = ms_processr_numeric(df, i, paste0("upper.", effect)),
        p_value = ms_processr_numeric(df, i, paste0("p.", effect))
      )
      if (!is.na(row$estimate)) rows[[length(rows) + 1L]] <- row
    }
  }

  ms_processr_table_payload(
    type = "processr_mediation_summary",
    type_label = "Mediation summary (processR)",
    call = ms_processr_call(.call, "processR::medSummary(..., boot.ci.type = \"all\")"),
    rows = rows,
    columns = c(
      list(list(key = "interval_type", label = "CI type", format = "text")),
      ms_processr_effect_columns(include_equation = TRUE)
    ),
    note = ms_processr_summary_note(x, moderated = FALSE, rows = rows),
    raw_output = ms_capture_output(x)
  )
}

ms_processr_modmed_summary_payload <- function(x, .call = NULL) {
  df <- as.data.frame(x, stringsAsFactors = FALSE)
  moderator <- attr(x, "mod", exact = TRUE) %||% "Moderator"
  has_label <- "label" %in% names(df)
  label_values <- if (has_label) {
    unique(stats::na.omit(as.character(df[["label"]])))
  } else character(0)
  include_label <- has_label && length(label_values) > 1L
  rows <- list()
  for (i in seq_len(nrow(df))) {
    condition <- ms_processr_value(df, i, c("values", "W", moderator))
    label <- if (include_label) ms_processr_value(df, i, "label") else NA_character_
    rows[[length(rows) + 1L]] <- list(
      moderator_value = condition,
      label = label,
      effect = "Conditional indirect effect",
      estimate = ms_processr_numeric(df, i, "indirect"),
      ci_lower = ms_processr_numeric(df, i, "lower"),
      ci_upper = ms_processr_numeric(df, i, "upper"),
      p_value = ms_processr_numeric(df, i, "indirectp")
    )
    rows[[length(rows) + 1L]] <- list(
      moderator_value = condition,
      label = label,
      effect = "Conditional direct effect",
      estimate = ms_processr_numeric(df, i, "direct"),
      ci_lower = ms_processr_numeric(df, i, "lowerd"),
      ci_upper = ms_processr_numeric(df, i, "upperd"),
      p_value = ms_processr_numeric(df, i, "directp")
    )
  }
  rows <- Filter(function(row) !is.na(row$estimate), rows)

  columns <- list(
    list(key = "moderator_value",
         label = ms_processr_moderator_label(moderator),
         format = "text")
  )
  if (include_label) {
    columns <- c(columns, list(list(key = "label", label = "Path", format = "text")))
  }
  columns <- c(columns, ms_processr_effect_columns(include_equation = FALSE))

  ms_processr_table_payload(
    type = "processr_moderated_mediation_summary",
    type_label = "Moderated mediation summary (processR)",
    call = ms_processr_call(.call, "processR::modmedSummary(...)"),
    rows = rows,
    columns = columns,
    note = ms_processr_summary_note(x, moderated = TRUE, rows = rows),
    raw_output = ms_capture_output(x)
  )
}

ms_processr_modmed_summary2_payload <- function(x, .call = NULL) {
  df <- as.data.frame(x, stringsAsFactors = FALSE)
  columns <- ms_table_columns_from_df(df)
  note <- ms_processr_summary_note(x, moderated = TRUE)
  if (nzchar(note)) {
    note <- paste(
      note,
      "This compact processR object does not include confidence intervals; Mellio renders it as a summary table only."
    )
  }
  ms_processr_table_payload(
    type = "processr_moderated_mediation_summary",
    type_label = "Moderated mediation summary (processR)",
    call = ms_processr_call(.call, "processR::modmedSummary(...)"),
    rows = ms_rows_from_df(df),
    columns = columns,
    note = note,
    raw_output = ms_capture_output(x)
  )
}

ms_processr_table_payload <- function(type, type_label, call, rows, columns,
                                      note = "", raw_output = "") {
  fields <- list(
    table_type = type,
    rows = rows,
    columns = columns,
    source = "processR"
  )
  if (nzchar(note)) fields$note <- note

  ms_build_envelope(
    type = type,
    type_label = type_label,
    call = trimws(gsub("\\s+", " ", call)),
    fields = fields,
    raw_output = raw_output,
    packages = ms_packages_basic(extras = c("processR", "lavaan")),
    card_kind = "table"
  )
}

ms_processr_effect_columns <- function(include_equation = TRUE) {
  columns <- list(
    list(key = "effect", label = "Effect", format = "text")
  )
  if (isTRUE(include_equation)) {
    columns <- c(columns, list(
      list(key = "equation", label = "Equation", format = "text")
    ))
  }
  c(columns, list(
    list(key = "estimate", label = "Estimate", format = "number"),
    list(key = "ci", label = "95% CI", format = "ci"),
    list(key = "p_value", label = "p", format = "pvalue")
  ))
}

ms_processr_summary_note <- function(x, moderated = FALSE, rows = NULL) {
  bits <- c(
    if (isTRUE(moderated)) {
      "Conditional effects are from a processR moderated-mediation summary"
    } else {
      "Effects are from a processR mediation summary"
    }
  )
  se <- attr(x, "se", exact = TRUE)
  boot <- attr(x, "boot.ci.type", exact = TRUE)
  if (!is.null(boot) && length(boot) && !is.na(boot[[1]]) && nzchar(boot[[1]])) {
    bits <- c(bits, paste0("boot.ci.type = ", boot[[1]]))
  } else if (!is.null(se) && length(se) && identical(as.character(se[[1]]), "standard")) {
    bits <- c(bits, "standard lavaan CIs")
  }
  bits <- c(bits,
    "Mellio renders this processR summary as a table only; send the underlying lavaan fit for supported mediation diagrams and topology-specific publication output"
  )
  if (ms_processr_has_prop_mediated(rows)) {
    bits <- c(bits,
      "proportion mediated is a ratio and can be unstable when the total effect is small; its CI may fall outside [0, 1]"
    )
  }
  paste0(paste(bits, collapse = "; "), ".")
}

ms_processr_has_prop_mediated <- function(rows) {
  if (is.null(rows) || !length(rows)) return(FALSE)
  any(vapply(rows, function(row) {
    effect <- tolower(as.character(row$effect %||% ""))
    effect <- gsub("[._[:space:]-]+", ".", trimws(effect))
    effect %in% c("prop.mediated", "proportion.mediated")
  }, logical(1)))
}

ms_processr_moderator_label <- function(moderator) {
  moderator <- ms_processr_display_label(moderator)
  if (!nzchar(moderator)) moderator <- "Moderator"
  paste0(moderator, " (W)")
}

ms_processr_display_label <- function(value) {
  value <- trimws(as.character(value %||% ""))
  if (!nzchar(value)) return("")
  value <- sub("^`(.+)`$", "\\1", value)
  value <- sub(
    "(?:[_\\-. ](?:c|centered|centred|center|centre|z|std|standardized|standardised|scaled|scale))$",
    "",
    value,
    perl = TRUE,
    ignore.case = TRUE
  )
  value <- gsub("[_.]+", " ", value)
  value <- gsub("\\s+", " ", value)
  trimws(value)
}

ms_processr_call <- function(.call, fallback) {
  if (!is.null(.call) && !is.na(.call) && nzchar(.call)) return(.call)
  fallback
}

ms_processr_value <- function(df, i, keys) {
  for (key in keys) {
    if (key %in% names(df)) {
      value <- df[[key]][[i]]
      if (is.factor(value)) value <- as.character(value)
      if (is.null(value) || length(value) == 0L || is.na(value)) return(NA_character_)
      return(as.character(value))
    }
  }
  NA_character_
}

ms_processr_numeric <- function(df, i, keys) {
  for (key in keys) {
    if (key %in% names(df)) {
      return(ms_safe_numeric(df[[key]][[i]]))
    }
  }
  NA_real_
}
