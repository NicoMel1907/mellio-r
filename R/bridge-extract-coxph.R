# R bridge -- Cox proportional hazards models.
#
# survival::coxph has a solid broom path, but a Cox model needs survival-
# specific editorial semantics: hazard ratios, event counts, and a
# likelihood-ratio chi-square headline. Keep that logic here rather than
# letting the generic broom model renderer describe it as ordinary
# regression.

#' @rdname mellio_payload
#' @export
mellio_payload.coxph <- function(x, ..., .call = NULL) {
  rlang::check_installed("survival", reason = "to extract Cox models")

  s <- tryCatch(summary(x), error = function(e) NULL)
  if (is.null(s) || is.null(s$coefficients)) {
    stop("Could not extract a Cox model summary.", call. = FALSE)
  }

  coef_mat <- as.matrix(s$coefficients)
  conf_mat <- if (!is.null(s$conf.int)) as.matrix(s$conf.int) else NULL
  terms <- rownames(coef_mat)
  if (is.null(terms)) terms <- paste0("term_", seq_len(nrow(coef_mat)))
  assign_map <- tryCatch(ms_broom_coef_assign_map(x), error = function(e) NULL)
  coef_label_map <- tryCatch(ms_coefficient_label_map(x), error = function(e) list())

  col_pick <- function(pattern) {
    hit <- grep(pattern, colnames(coef_mat), ignore.case = TRUE, value = TRUE)
    if (length(hit)) hit[[1L]] else NA_character_
  }
  conf_pick <- function(pattern) {
    if (is.null(conf_mat)) return(NA_character_)
    hit <- grep(pattern, colnames(conf_mat), ignore.case = TRUE, value = TRUE)
    if (length(hit)) hit[[1L]] else NA_character_
  }

  coef_col <- col_pick("^coef$")
  hr_col <- col_pick("^exp\\(coef\\)$")
  se_col <- col_pick("se\\(coef\\)|robust se")
  z_col <- col_pick("^z$")
  p_col <- col_pick("Pr\\(>\\|z\\|\\)")
  lo_col <- conf_pick("^lower")
  hi_col <- conf_pick("^upper")

  coef_value <- function(i, col) {
    if (is.na(col) || !col %in% colnames(coef_mat)) return(NA_real_)
    ms_safe_numeric(coef_mat[i, col])
  }
  conf_value <- function(term, col) {
    if (is.null(conf_mat) || is.na(col) || !term %in% rownames(conf_mat) ||
        !col %in% colnames(conf_mat)) {
      return(NA_real_)
    }
    ms_safe_numeric(conf_mat[term, col])
  }

  coefficients <- lapply(seq_len(nrow(coef_mat)), function(i) {
    term <- terms[[i]]
    log_estimate <- coef_value(i, coef_col)
    hr <- coef_value(i, hr_col)
    if (is.na(hr) && !is.na(log_estimate)) hr <- exp(log_estimate)
    row <- list(
      term = term,
      estimate_name = "HR",
      estimate = hr,
      log_estimate = log_estimate
    )
    if (!is.na(se_col)) row$std_error <- coef_value(i, se_col)
    if (!is.na(z_col)) {
      row$statistic <- coef_value(i, z_col)
      row$statistic_label <- "z"
    }
    if (!is.na(p_col)) row$p_value <- coef_value(i, p_col)
    if (!is.null(conf_mat) && term %in% rownames(conf_mat)) {
      if (!is.na(lo_col)) row$ci_lower <- conf_value(term, lo_col)
      if (!is.na(hi_col)) row$ci_upper <- conf_value(term, hi_col)
      if (!is.null(row$ci_lower) && !is.null(row$ci_upper) &&
          !is.na(row$ci_lower) && !is.na(row$ci_upper)) {
        row$ci_method <- "wald"
      }
    }
    row$term_source <- if (!is.null(assign_map) && term %in% names(assign_map)) {
      assign_map[[term]]
    } else {
      term
    }
    ms_apply_coefficient_label(row, coef_label_map)
  })

  logtest <- s$logtest %||% numeric(0)
  fields <- list(
    coefficients = coefficients,
    coefficient_scale = "hazard_ratio",
    coefficient_ci_method = "wald",
    coefficient_ci_scale = "hazard_ratio",
    coefficient_p_value_method = "wald_z",
    conf_level = 0.95,
    source = "R",
    outcome = "survival",
    statistic_label = "z",
    aic = ms_safe_numeric(stats::AIC(x)),
    bic = ms_safe_numeric(stats::BIC(x)),
    logLik = ms_safe_numeric(as.numeric(stats::logLik(x)))
  )

  if (length(logtest) >= 3L) {
    fields$statistic <- list(
      name = "\u03c7\u00b2",
      value = ms_safe_numeric(logtest[["test"]] %||% logtest[[1L]]),
      df = ms_safe_numeric(logtest[["df"]] %||% logtest[[2L]])
    )
    fields$p_value <- ms_safe_numeric(logtest[["pvalue"]] %||% logtest[[3L]])
    fields$test_name <- "Likelihood ratio test"
  }

  n <- ms_safe_numeric(s$n)
  if (!is.na(n)) fields$n <- as.integer(n)
  events <- ms_safe_numeric(s$nevent)
  if (!is.na(events)) fields$events <- as.integer(events)

  model_warnings <- ms_coxph_model_warnings(coefficients)
  if (length(model_warnings) > 0L) fields$model_warnings <- model_warnings

  model_terms <- tryCatch(attr(stats::terms(x), "term.labels"), error = function(e) character(0))
  if (length(model_terms) > 0L) {
    fields$terms <- lapply(model_terms, function(term) {
      list(name = term, label = ms_model_clean_term(term), role = "term", type = "predictor")
    })
    fields$focal_terms <- model_terms
    fields$predictor <- ms_model_term_phrase(vapply(fields$terms, function(t) {
      t$name %||% t$label
    }, character(1)))
  }

  ms_build_envelope(
    type = "cox_proportional_hazards",
    type_label = "Cox proportional hazards model",
    call = trimws(gsub("\\s+", " ", ms_model_call_string(x, .call = .call))),
    fields = fields,
    raw_output = ms_capture_output(s),
    packages = ms_packages_basic(extras = "survival")
  )
}

ms_coxph_model_warnings <- function(coefficients) {
  rows <- Filter(is.list, coefficients %||% list())
  if (length(rows) == 0L) return(list())

  unstable <- vapply(rows, function(row) {
    estimate <- ms_safe_numeric(row$estimate)
    log_estimate <- ms_safe_numeric(row$log_estimate)
    std_error <- ms_safe_numeric(row$std_error)
    ci_lower <- ms_safe_numeric(row$ci_lower)
    ci_upper <- ms_safe_numeric(row$ci_upper)

    extreme_hr <- !is.na(estimate) && (estimate > 1e6 || estimate < 1e-6)
    extreme_log <- !is.na(log_estimate) && abs(log_estimate) > 8
    huge_se <- !is.na(std_error) && std_error > 100
    missing_ci <- is.na(ci_lower) || is.na(ci_upper)

    (extreme_hr || extreme_log) && (huge_se || missing_ci)
  }, logical(1))

  if (!any(unstable)) return(list())
  list(list(
    type = "separation_or_boundary",
    severity = "warning",
    message = paste(
      "The Cox model produced extreme hazard-ratio estimates or infinite confidence limits;",
      "monotone likelihood, separation, or sparse event patterns may make coefficient estimates unstable."
    )
  ))
}
