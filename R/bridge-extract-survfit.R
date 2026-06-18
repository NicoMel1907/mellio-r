# R bridge -- survival::survfit Kaplan-Meier estimates.
#
# survival::survfit returns step-function estimates of the survival curve,
# optionally split by strata. Mellio renders the curve client-side from
# structured data; the R side extracts per-stratum tracks (time, survival,
# CI, at-risk, events, censor counts) plus a median-survival summary table
# and a log-rank companion test when the original call can be replayed.
#
# Schema: docs/STATS-R-BRIDGE-SCHEMA.md

#' @rdname mellio_payload
#' @export
mellio_payload.survfit <- function(x, ..., .call = NULL, .env = parent.frame()) {
  rlang::check_installed("survival", reason = "to extract Kaplan-Meier estimates")

  # Multi-state survfitms objects dispatch separately; here we only handle
  # plain two-state survival curves.
  if (inherits(x, "survfitms")) {
    stop("Multi-state survfit objects are not supported yet.", call. = FALSE)
  }

  summary_obj <- tryCatch(summary(x), error = function(e) NULL)
  s_table <- summary_obj$table
  if (is.null(s_table)) {
    stop("Could not extract a survfit summary table.", call. = FALSE)
  }

  if (is.null(x$strata)) {
    strata_names <- "Overall"
    strata_lens <- length(x$time)
    has_strata <- FALSE
  } else {
    strata_names <- names(x$strata)
    strata_lens <- as.integer(x$strata)
    has_strata <- TRUE
  }

  build_track <- function(start, len, name) {
    idx <- seq.int(start, length.out = len)
    out <- list(
      name = unname(name),
      label = ms_survfit_strata_label(name),
      time = as.numeric(x$time[idx]),
      surv = as.numeric(x$surv[idx]),
      n_risk = as.integer(x$n.risk[idx]),
      n_event = as.integer(x$n.event[idx]),
      n_censor = as.integer(x$n.censor[idx])
    )
    if (!is.null(x$lower)) out$lower <- as.numeric(x$lower[idx])
    if (!is.null(x$upper)) out$upper <- as.numeric(x$upper[idx])
    out
  }

  tracks <- list()
  cursor <- 1L
  for (i in seq_along(strata_lens)) {
    len <- strata_lens[[i]]
    if (len > 0L) {
      tracks[[length(tracks) + 1L]] <- build_track(cursor, len, strata_names[[i]])
    }
    cursor <- cursor + len
  }

  conf_level <- ms_survfit_conf_level(x)
  tracks <- ms_survfit_attach_medians(tracks, s_table, has_strata)

  log_rank <- ms_survfit_log_rank(x, .env)
  formula_info <- ms_survfit_formula_info(x)
  formula_info$time_unit <- ms_survfit_time_unit(x, formula_info)

  rows <- lapply(tracks, function(tr) {
    list(
      group = tr$label %||% tr$name,
      group_key = tr$name,
      n = tr$n_total,
      events = tr$events_total,
      median = tr$median,
      ci_lower = tr$median_lower,
      ci_upper = tr$median_upper
    )
  })
  columns <- list(
    list(key = "group",  label = "Group",           format = "text"),
    list(key = "n",      label = "n",               format = "integer"),
    list(key = "events", label = "Events",          format = "integer"),
    list(key = "median", label = "Median survival", format = "number"),
    list(key = "ci",     label = ms_survfit_ci_label(conf_level), format = "ci")
  )

  total_n <- ms_survfit_sum_field(tracks, "n_total")
  total_events <- ms_survfit_sum_field(tracks, "events_total")

  fields <- list(
    source            = "survival",
    outcome           = formula_info$outcome %||% "survival",
    predictor         = formula_info$predictor,
    time_label        = formula_info$time_label,
    time_unit         = formula_info$time_unit,
    event_label       = formula_info$event_label,
    n                 = total_n,
    events            = total_events,
    groups_count      = length(tracks),
    has_strata        = has_strata,
    conf_level        = conf_level,
    conf_type         = as.character(x$conf.type %||% "log"),
    rows              = rows,
    columns           = columns,
    table_type        = "kaplan_meier_summary",
    test_name         = if (!is.null(log_rank)) "Log-rank test" else NULL,
    statistic_label   = if (!is.null(log_rank)) "\u03c7\u00b2" else NULL,
    note              = ms_survfit_table_note(log_rank, has_strata, formula_info, conf_level)
  )

  if (!is.null(log_rank)) {
    fields$statistic <- list(
      name  = "\u03c7\u00b2",
      value = ms_safe_numeric(log_rank$chi_sq),
      df    = ms_safe_numeric(log_rank$df)
    )
    fields$p_value <- ms_safe_numeric(log_rank$p_value)
  }

  figure_data <- list(km_curve = list(
    source       = "survfit",
    tracks       = tracks,
    has_strata   = has_strata,
    n_total      = total_n,
    events_total = total_events,
    time_label   = formula_info$time_label,
    time_unit    = formula_info$time_unit,
    event_label  = formula_info$event_label,
    predictor    = formula_info$predictor,
    conf_level   = conf_level,
    conf_type    = as.character(x$conf.type %||% "log"),
    log_rank     = log_rank
  ))

  ms_build_envelope(
    type        = "kaplan_meier_survival",
    type_label  = "Kaplan-Meier survival",
    call        = trimws(gsub("\\s+", " ", ms_model_call_string(x, .call = .call))),
    fields      = fields,
    raw_output  = ms_capture_output(summary_obj),
    packages    = ms_packages_basic(extras = "survival"),
    figure_data = figure_data
  )
}

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------

# Strata are reported by survfit as "var=Level" pairs (or "var1=A, var2=B" for
# multi-strata). The figure should show "Level" / "A, B"; the Group column
# falls back to the full string when the prefix can't be stripped.
ms_survfit_strata_label <- function(name) {
  if (is.null(name) || !nzchar(as.character(name))) return("Overall")
  parts <- strsplit(as.character(name), ",\\s*", fixed = FALSE)[[1L]]
  cleaned <- vapply(parts, function(part) {
    eq <- regexpr("=", part, fixed = TRUE)
    if (eq < 1L) return(trimws(part))
    trimws(substr(part, eq + 1L, nchar(part)))
  }, character(1))
  paste(cleaned, collapse = ", ")
}

# summary(survfit_obj)$table is a numeric vector for single-stratum fits and
# a matrix for multi-stratum fits. Column names also differ across survival
# versions ("records" vs "n", "n.start" vs "n.max"), so be defensive.
ms_survfit_attach_medians <- function(tracks, s_table, has_strata) {
  if (is.null(s_table) || length(tracks) == 0L) return(tracks)

  if (is.matrix(s_table)) {
    col_names <- colnames(s_table)
    pick_row <- function(idx) {
      if (idx > nrow(s_table)) return(NULL)
      stats::setNames(s_table[idx, , drop = TRUE], col_names)
    }
  } else {
    col_names <- names(s_table)
    pick_row <- function(idx) {
      if (idx > 1L) return(NULL)
      stats::setNames(as.numeric(s_table), col_names)
    }
  }

  read_field <- function(row, candidates) {
    if (is.null(row)) return(NA_real_)
    for (cand in candidates) {
      if (cand %in% names(row)) {
        value <- row[[cand]]
        if (!is.null(value) && !is.na(value)) return(ms_safe_numeric(value))
      }
    }
    NA_real_
  }

  for (i in seq_along(tracks)) {
    row <- pick_row(i)
    n_val <- read_field(row, c("records", "n", "n.start", "n.max"))
    evt_val <- read_field(row, c("events"))
    tracks[[i]]$n_total <- if (is.na(n_val)) NA_integer_ else as.integer(round(n_val))
    tracks[[i]]$events_total <- if (is.na(evt_val)) NA_integer_ else as.integer(round(evt_val))
    tracks[[i]]$median <- read_field(row, c("median"))
    tracks[[i]]$median_lower <- read_field(row, ms_survfit_ci_column_candidates(row, "LCL"))
    tracks[[i]]$median_upper <- read_field(row, ms_survfit_ci_column_candidates(row, "UCL"))
  }
  tracks
}

ms_survfit_ci_column_candidates <- function(row, suffix) {
  suffix <- as.character(suffix %||% "")
  fallback <- c(
    paste0("0.95", suffix),
    paste0("0.9", suffix),
    paste0("0.99", suffix)
  )
  if (is.null(row) || !nzchar(suffix)) return(fallback)
  hits <- grep(paste0(suffix, "$"), names(row), value = TRUE)
  unique(c(hits, fallback))
}

ms_survfit_conf_level <- function(x) {
  level <- ms_safe_numeric(x$conf.int %||% NA_real_)
  if (is.na(level) || level <= 0 || level >= 1) 0.95 else level
}

ms_survfit_ci_label <- function(conf_level) {
  paste0(ms_survfit_conf_percent(conf_level), "% CI")
}

# Log-rank companion. We replay the survfit call with survival::survdiff
# substituted in. Bail quietly if the call can't be reconstructed (e.g. the
# data argument has gone out of scope).
ms_survfit_log_rank <- function(x, env) {
  if (is.null(x$strata)) return(NULL)
  if (length(x$strata) < 2L) return(NULL)
  call_obj <- x$call
  if (is.null(call_obj)) return(NULL)

  rho <- if (is.environment(env)) env else parent.frame()
  diff <- tryCatch({
    new_call <- call_obj
    new_call[[1L]] <- quote(survival::survdiff)
    eval(new_call, envir = rho)
  }, error = function(e) NULL)
  if (is.null(diff)) return(NULL)

  chi <- ms_safe_numeric(diff$chisq)
  df_val <- if (!is.null(diff$n)) max(length(diff$n) - 1L, 1L) else NA_real_
  if (is.na(df_val) && !is.null(diff$obs)) {
    df_val <- max(length(diff$obs) - 1L, 1L)
  }
  p_val <- if (!is.na(chi) && !is.na(df_val) && df_val > 0L) {
    ms_safe_numeric(stats::pchisq(chi, df = df_val, lower.tail = FALSE))
  } else {
    NA_real_
  }

  list(
    chi_sq   = chi,
    df       = ms_safe_numeric(df_val),
    p_value  = p_val,
    method   = "log_rank"
  )
}

# Pull readable predictor / time labels from the survfit call's formula. The
# left-hand side of Surv(time, status) gives the time variable; the right-
# hand side gives the stratification predictor(s).
ms_survfit_formula_info <- function(x) {
  out <- list(
    time_label = NULL,
    event_label = NULL,
    predictor = NULL,
    outcome = NULL
  )
  formula_obj <- tryCatch(stats::formula(x), error = function(e) NULL)
  if (is.null(formula_obj)) {
    call_obj <- x$call
    if (!is.null(call_obj)) {
      raw <- call_obj[["formula"]]
      if (!is.null(raw)) {
        formula_obj <- tryCatch(stats::as.formula(raw), error = function(e) NULL)
      }
    }
  }
  if (is.null(formula_obj)) return(out)

  if (length(formula_obj) >= 3L) {
    lhs <- formula_obj[[2L]]
    fn_name <- ms_survfit_call_head(lhs)
    if (is.call(lhs) && identical(fn_name, "Surv") && length(lhs) >= 3L) {
      if (length(lhs) >= 4L) {
        out$time_label <- ms_survfit_expr_label(lhs[[3L]])
        out$event_label <- ms_survfit_expr_label(lhs[[4L]])
      } else {
        out$time_label <- ms_survfit_expr_label(lhs[[2L]])
        out$event_label <- ms_survfit_expr_label(lhs[[3L]])
      }
    } else if (is.symbol(lhs) || is.character(lhs)) {
      out$time_label <- ms_survfit_expr_label(lhs)
    }

    rhs <- formula_obj[[3L]]
    rhs_label <- trimws(paste(deparse(rhs, width.cutoff = 200L), collapse = " "))
    if (nzchar(rhs_label) && rhs_label != "1" && rhs_label != "NULL") {
      out$predictor <- rhs_label
    }
  }

  out
}

ms_survfit_expr_label <- function(expr) {
  label <- paste(deparse(expr, width.cutoff = 200L), collapse = " ")
  trimws(label)
}

ms_survfit_time_unit <- function(x, formula_info) {
  time_label <- tolower(trimws(as.character(formula_info$time_label %||% "")))
  if (!nzchar(time_label)) return(NULL)
  if (grepl("(^|[_.\\s-])days?($|[_.\\s-])", time_label)) return("days")
  if (grepl("(^|[_.\\s-])months?($|[_.\\s-])", time_label)) return("months")
  if (grepl("(^|[_.\\s-])years?($|[_.\\s-])", time_label)) return("years")

  data_expr <- tryCatch(x$call[["data"]], error = function(e) NULL)
  data_label <- if (!is.null(data_expr)) {
    tolower(paste(deparse(data_expr, width.cutoff = 200L), collapse = " "))
  } else {
    ""
  }
  if (identical(time_label, "time") && grepl("lung", data_label, fixed = TRUE)) {
    return("days")
  }
  NULL
}

# Returns the unqualified function name from a call head. Handles bare
# symbols (`Surv`) and namespace-qualified calls (`survival::Surv`).
ms_survfit_call_head <- function(node) {
  if (!is.call(node) && !is.symbol(node)) return(NA_character_)
  head <- if (is.call(node)) node[[1L]] else node
  if (is.symbol(head)) return(as.character(head))
  if (is.call(head)) {
    parts <- as.character(head)
    return(parts[length(parts)])
  }
  NA_character_
}

ms_survfit_sum_field <- function(tracks, field) {
  if (length(tracks) == 0L) return(NA_integer_)
  vals <- vapply(tracks, function(tr) {
    v <- tr[[field]]
    if (is.null(v)) NA_integer_ else as.integer(v)
  }, integer(1))
  if (any(is.na(vals))) return(NA_integer_)
  as.integer(sum(vals))
}

ms_survfit_table_note <- function(log_rank, has_strata, formula_info, conf_level = 0.95) {
  bits <- character(0)
  ci_label <- ms_survfit_ci_label(conf_level)
  if (isTRUE(has_strata)) {
    bits <- c(bits, paste0("Median survival and ", ci_label, " are Kaplan-Meier estimates per ",
                           formula_info$predictor %||% "group", "."))
  } else {
    bits <- c(bits, paste0("Median survival and ", ci_label, " are Kaplan-Meier estimates."))
  }
  if (!is.null(log_rank) && !is.na(log_rank$chi_sq) && !is.na(log_rank$df)) {
    chi <- formatC(log_rank$chi_sq, format = "f", digits = 2)
    df_str <- format(as.integer(log_rank$df))
    p_str <- ms_survfit_format_p(log_rank$p_value)
    bits <- c(bits, paste0("Log-rank test computed with survival::survdiff: \u03c7\u00b2(",
                           df_str, ") = ", chi, ", ", p_str, "."))
  }
  paste(bits, collapse = " ")
}

ms_survfit_conf_percent <- function(conf_level) {
  level <- ms_survfit_conf_level(list(conf.int = conf_level))
  pct <- level * 100
  if (abs(pct - round(pct)) < 1e-8) return(format(as.integer(round(pct))))
  sub("\\.?0+$", "", formatC(pct, format = "f", digits = 1))
}

ms_survfit_format_p <- function(p) {
  if (is.null(p) || is.na(p)) return("p = NA")
  if (p < .001) return("p < .001")
  formatted <- sub("^0", "", formatC(p, format = "f", digits = 3))
  paste0("p = ", formatted)
}
