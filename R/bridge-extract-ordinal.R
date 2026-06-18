# R bridge -- ordinal regression models.

#' @rdname mellio_payload
#' @export
mellio_payload.polr <- function(x, ..., .call = NULL,
                            conf.int = TRUE,
                            conf.level = 0.95) {
  rlang::check_installed("broom", reason = "to extract ordinal regression models")
  rlang::check_installed("MASS", reason = "to extract MASS::polr models")

  tidy_log <- ms_polr_tidy(x, conf.int = conf.int, conf.level = conf.level,
                           exponentiate = FALSE)
  tidy_or <- ms_polr_tidy(x, conf.int = conf.int, conf.level = conf.level,
                          exponentiate = TRUE)

  coef_rows <- ms_ordinal_polr_coefficients(tidy_log, tidy_or)
  if (!length(coef_rows)) {
    stop("Could not extract ordinal regression coefficients.", call. = FALSE)
  }
  separation <- ms_ordinal_empty_cell_diagnostics(x, coef_rows)
  coef_rows <- separation$coefficients
  coef_rows <- ms_ordinal_add_wald_fallback_ci(
    coef_rows,
    conf.int = conf.int,
    conf.level = conf.level
  )
  inference <- ms_ordinal_missing_inference_diagnostics(coef_rows)
  coef_rows <- inference$coefficients

  fields <- ms_ordinal_base_fields(x, link = ms_ordinal_link_name(x$method %||% "logit"))
  fields$coefficients <- coef_rows
  fields$thresholds <- ms_ordinal_polr_thresholds(tidy_log)
  fields$threshold_count <- length(fields$thresholds)
  model_term_tests <- ms_ordinal_model_term_tests(x)
  if (length(model_term_tests) > 0L) fields$model_term_tests <- model_term_tests
  if (isTRUE(conf.int) && ms_ordinal_any_ci(coef_rows)) {
    fields$conf_level <- conf.level
    fields$coefficient_ci_method <- ms_ordinal_coefficient_ci_method(coef_rows)
  }
  fields <- ms_ordinal_add_probability_curve_warning(fields, x)
  fields <- ms_ordinal_add_warning_rows(fields, separation$warnings)
  fields <- ms_ordinal_add_warning_rows(fields, inference$warnings)
  fields <- ms_ordinal_add_unstable_coefficient_warning(fields, coef_rows)

  glance_df <- tryCatch(as.data.frame(broom::glance(x)), error = function(e) data.frame())
  fields <- c(fields, ms_broom_glance_fields(glance_df))

  probability_curve <- tryCatch(
    ms_ordinal_probability_curve_figure_data(x, conf.level = conf.level),
    error = function(e) NULL
  )
  figure_data <- list()
  available_figures <- NULL
  if (!is.null(probability_curve)) {
    figure_data$interaction_plot <- probability_curve
    available_figures <- list(list(
      type = "interaction_plot",
      label = ms_ordinal_probability_curve_label(probability_curve),
      default = TRUE
    ))
  }
  if (!length(figure_data)) figure_data <- NULL

  ms_ordinal_payload(
    x = x,
    fields = fields,
    call = ms_model_call_string(x, .call = .call),
    raw_output = ms_capture_output(summary(x)),
    packages = ms_packages_basic(extras = c("MASS", "broom")),
    figure_data = figure_data,
    available_figures = available_figures
  )
}

#' @rdname mellio_payload
#' @export
mellio_payload.clm <- function(x, ..., .call = NULL,
                               conf.int = TRUE,
                               conf.level = 0.95) {
  rlang::check_installed("ordinal", reason = "to extract ordinal::clm models")

  s <- tryCatch(summary(x), error = function(e) NULL)
  if (is.null(s) || is.null(s$coefficients)) {
    stop("Could not extract ordinal regression coefficients.", call. = FALSE)
  }

  coef_mat <- as.matrix(s$coefficients)
  beta_names <- names(x$beta %||% numeric(0))
  alpha_names <- names(x$alpha %||% numeric(0))
  if (!length(beta_names)) {
    beta_names <- rownames(coef_mat)[!grepl("\\|", rownames(coef_mat) %||% character(0))]
  }
  profile_ci <- ms_ordinal_clm_profile_ci(
    x,
    beta_names,
    conf.int = conf.int,
    conf.level = conf.level
  )
  coef_rows <- ms_ordinal_clm_coefficients(
    coef_mat,
    beta_names,
    profile_ci = profile_ci,
    conf.int = conf.int,
    conf.level = conf.level
  )
  if (!length(coef_rows)) {
    stop("Could not extract ordinal regression coefficients.", call. = FALSE)
  }
  separation <- ms_ordinal_empty_cell_diagnostics(x, coef_rows)
  coef_rows <- separation$coefficients
  coef_rows <- ms_ordinal_add_wald_fallback_ci(
    coef_rows,
    conf.int = conf.int,
    conf.level = conf.level
  )
  inference <- ms_ordinal_missing_inference_diagnostics(coef_rows)
  coef_rows <- inference$coefficients

  fields <- ms_ordinal_base_fields(x, link = ms_ordinal_link_name(x$link %||% ""))
  fields$coefficients <- coef_rows
  fields$thresholds <- ms_ordinal_clm_thresholds(coef_mat, alpha_names)
  fields$threshold_count <- length(fields$thresholds)
  model_term_tests <- ms_ordinal_model_term_tests(x)
  if (length(model_term_tests) > 0L) fields$model_term_tests <- model_term_tests
  if (isTRUE(conf.int) && ms_ordinal_any_ci(coef_rows)) {
    fields$conf_level <- conf.level
    fields$coefficient_ci_method <- ms_ordinal_coefficient_ci_method(coef_rows)
  }
  fields <- ms_ordinal_add_probability_curve_warning(fields, x)
  fields <- ms_ordinal_add_warning_rows(fields, separation$warnings)
  fields <- ms_ordinal_add_warning_rows(fields, inference$warnings)
  fields <- ms_ordinal_add_unstable_coefficient_warning(fields, coef_rows)

  n <- tryCatch(stats::nobs(x), error = function(e) NA_real_)
  if (!is.na(ms_safe_numeric(n))) fields$n <- as.integer(ms_safe_numeric(n))
  aic <- tryCatch(stats::AIC(x), error = function(e) NA_real_)
  if (!is.na(ms_safe_numeric(aic))) fields$aic <- ms_safe_numeric(aic)
  bic <- tryCatch(stats::BIC(x), error = function(e) NA_real_)
  if (!is.na(ms_safe_numeric(bic))) fields$bic <- ms_safe_numeric(bic)
  loglik <- tryCatch(as.numeric(stats::logLik(x)), error = function(e) NA_real_)
  if (!is.na(ms_safe_numeric(loglik))) fields$logLik <- ms_safe_numeric(loglik)
  residual_df <- tryCatch(stats::df.residual(x), error = function(e) NA_real_)
  if (!is.na(ms_safe_numeric(residual_df))) fields$residual_df <- ms_safe_numeric(residual_df)

  probability_curve <- tryCatch(
    ms_ordinal_probability_curve_figure_data(x, conf.level = conf.level),
    error = function(e) NULL
  )
  figure_data <- list()
  available_figures <- NULL
  if (!is.null(probability_curve)) {
    figure_data$interaction_plot <- probability_curve
    available_figures <- list(list(
      type = "interaction_plot",
      label = ms_ordinal_probability_curve_label(probability_curve),
      default = TRUE
    ))
  }
  if (!length(figure_data)) figure_data <- NULL

  ms_ordinal_payload(
    x = x,
    fields = fields,
    call = ms_model_call_string(x, .call = .call),
    raw_output = ms_capture_output(s),
    packages = ms_packages_basic(extras = "ordinal"),
    figure_data = figure_data,
    available_figures = available_figures
  )
}

ms_ordinal_payload <- function(x, fields, call, raw_output, packages,
                               figure_data = NULL, available_figures = NULL) {
  ms_build_envelope(
    type = "ordinal_regression",
    type_label = "Ordinal regression",
    call = trimws(gsub("\\s+", " ", call %||% NA_character_)),
    fields = fields,
    raw_output = raw_output,
    packages = packages,
    figure_data = figure_data,
    available_figures = available_figures
  )
}

ms_ordinal_base_fields <- function(x, link) {
  fields <- list(
    p_value = NA_real_,
    source = "R",
    model_kind = "ordinal_regression",
    model_family = "ordinal",
    model_link = as.character(link %||% ""),
    link_function = as.character(link %||% ""),
    coefficient_scale = "proportional_odds",
    coefficient_ci_scale = "proportional_odds",
    coefficient_p_value_method = "wald_z",
    statistic_label = "z"
  )
  term_roles <- tryCatch(ms_lm_term_roles(x), error = function(e) NULL)
  if (!is.null(term_roles)) {
    if (!is.null(term_roles$outcome)) fields$outcome <- term_roles$outcome
    if (length(term_roles$terms) > 0L) fields$terms <- term_roles$terms
    if (length(term_roles$focal_terms) > 0L) fields$focal_terms <- term_roles$focal_terms
    if (length(term_roles$control_terms) > 0L) fields$control_terms <- term_roles$control_terms
    if (!is.null(term_roles$predictor)) fields$predictor <- term_roles$predictor
    if (!is.null(term_roles$model_kind)) fields$model_kind <- "ordinal_regression"
  }
  outcome_levels <- ms_ordinal_outcome_levels(x)
  if (length(outcome_levels)) {
    fields$outcome_levels <- as.list(outcome_levels)
    fields$outcome_level_count <- length(outcome_levels)
  }
  fields
}

ms_ordinal_model_term_tests <- function(x) {
  tf <- tryCatch(stats::terms(x), error = function(e) NULL)
  mf <- tryCatch(stats::model.frame(x), error = function(e) NULL)
  if (is.null(tf) || is.null(mf)) return(list())

  term_labels <- attr(tf, "term.labels") %||% character(0)
  orders <- attr(tf, "order") %||% integer(0)
  if (!length(term_labels)) return(list())
  if (length(orders) != length(term_labels)) {
    orders <- rep(1L, length(term_labels))
  }

  eligible <- character(0)
  levels_by_term <- list()
  for (i in seq_along(term_labels)) {
    term <- term_labels[[i]]
    if (!nzchar(term) || orders[[i]] != 1L || grepl(":", term, fixed = TRUE)) next
    values <- ms_ordinal_term_values(term, mf)
    if (is.null(values) ||
        !(is.factor(values) || is.character(values) || is.logical(values))) {
      next
    }
    term_levels <- ms_ordinal_factor_levels(values)
    if (length(term_levels) <= 2L) next
    eligible <- c(eligible, term)
    levels_by_term[[term]] <- term_levels
  }
  if (!length(eligible)) return(list())

  drop_rows <- tryCatch(ms_model_term_tests_drop1_chisq(x, eligible),
                        error = function(e) list())
  out <- vector("list", length(eligible))
  for (i in seq_along(eligible)) {
    term <- eligible[[i]]
    row <- drop_rows[[term]]
    if (is.null(row)) next
    row$term_type <- "main"
    row$predictor_type <- "factor"
    row$test_scope <- "factor_omnibus"
    row$method <- "ordinal_factor_lrt"
    row$levels <- as.list(levels_by_term[[term]] %||% character(0))
    out[[i]] <- row
  }
  Filter(Negate(is.null), out)
}

ms_ordinal_term_values <- function(term, mf) {
  term <- as.character(term %||% "")
  if (!nzchar(term) || is.null(mf)) return(NULL)
  if (term %in% names(mf)) return(mf[[term]])
  vars <- tryCatch(all.vars(stats::as.formula(paste0("~", term))),
                   error = function(e) character(0))
  vars <- unique(vars[nzchar(vars)])
  if (length(vars) == 1L && vars[[1L]] %in% names(mf)) return(mf[[vars[[1L]]]])
  NULL
}

ms_ordinal_factor_levels <- function(values) {
  if (is.factor(values)) {
    return(as.character(levels(values)))
  }
  values <- values[!is.na(values)]
  if (is.logical(values)) {
    return(as.character(sort(unique(values))))
  }
  sort(unique(as.character(values)))
}

ms_ordinal_link_name <- function(link) {
  link <- as.character(link %||% "")
  if (identical(tolower(link), "logistic")) return("logit")
  link
}

ms_polr_tidy <- function(x, conf.int, conf.level, exponentiate) {
  out <- tryCatch(
    suppressMessages(suppressWarnings(broom::tidy(
      x,
      conf.int = conf.int,
      conf.level = conf.level,
      exponentiate = exponentiate
    ))),
    error = function(e) NULL
  )
  if (is.null(out) && isTRUE(conf.int)) {
    out <- suppressMessages(suppressWarnings(broom::tidy(
      x,
      conf.int = FALSE,
      exponentiate = exponentiate
    )))
  }
  as.data.frame(out)
}

ms_ordinal_polr_coefficients <- function(tidy_log, tidy_or) {
  if (!nrow(tidy_log) || !nrow(tidy_or)) return(list())
  keep <- if ("coef.type" %in% names(tidy_log)) {
    tidy_log$coef.type %in% "coefficient"
  } else {
    !grepl("\\|", tidy_log$term)
  }
  log_rows <- tidy_log[keep, , drop = FALSE]
  or_rows <- tidy_or[match(log_rows$term, tidy_or$term), , drop = FALSE]
  rows <- list()
  for (i in seq_len(nrow(log_rows))) {
    term <- as.character(log_rows$term[[i]])
    stat <- ms_ordinal_pick(log_rows, i, "statistic")
    p <- if ("p.value" %in% names(log_rows)) {
      ms_ordinal_pick(log_rows, i, "p.value")
    } else if (!is.na(stat)) {
      ms_safe_numeric(2 * stats::pnorm(abs(stat), lower.tail = FALSE))
    } else {
      NA_real_
    }
    row <- list(
      term = term,
      estimate_name = "OR",
      estimate = ms_ordinal_pick(or_rows, i, "estimate"),
      log_estimate = ms_ordinal_pick(log_rows, i, "estimate"),
      std_error = ms_ordinal_pick(log_rows, i, "std.error"),
      statistic = stat,
      statistic_label = "z",
      p_value = p,
      term_source = term
    )
    if ("conf.low" %in% names(or_rows)) row$ci_lower <- ms_ordinal_pick(or_rows, i, "conf.low")
    if ("conf.high" %in% names(or_rows)) row$ci_upper <- ms_ordinal_pick(or_rows, i, "conf.high")
    if (!is.na(row$ci_lower %||% NA_real_) && !is.na(row$ci_upper %||% NA_real_)) {
      row$ci_method <- "profile_likelihood"
    }
    rows[[length(rows) + 1L]] <- row
  }
  rows
}

ms_ordinal_polr_thresholds <- function(tidy_log) {
  if (!nrow(tidy_log)) return(list())
  keep <- if ("coef.type" %in% names(tidy_log)) {
    tidy_log$coef.type %in% "scale"
  } else {
    grepl("\\|", tidy_log$term)
  }
  rows <- tidy_log[keep, , drop = FALSE]
  lapply(seq_len(nrow(rows)), function(i) {
    list(
      threshold = as.character(rows$term[[i]]),
      estimate = ms_ordinal_pick(rows, i, "estimate"),
      std_error = ms_ordinal_pick(rows, i, "std.error"),
      statistic = ms_ordinal_pick(rows, i, "statistic")
    )
  })
}

ms_ordinal_clm_coefficients <- function(coef_mat, beta_names,
                                        profile_ci = NULL,
                                        conf.int = TRUE,
                                        conf.level = 0.95) {
  rn <- rownames(coef_mat) %||% character(0)
  if (!length(beta_names)) {
    beta_names <- rn[!grepl("\\|", rn)]
  }
  beta_names <- intersect(beta_names, rn)
  est_col <- ms_ordinal_matrix_col(coef_mat, "^Estimate$")
  se_col <- ms_ordinal_matrix_col(coef_mat, "Std\\. Error|Std\\.Error")
  stat_col <- ms_ordinal_matrix_col(coef_mat, "z value|t value|statistic")
  p_col <- ms_ordinal_matrix_col(coef_mat, "Pr\\(>|p.value")

  ci_crit <- stats::qnorm((1 + conf.level) / 2)

  lapply(beta_names, function(term) {
    i <- match(term, rn)
    log_est <- ms_ordinal_matrix_cell(coef_mat, i, est_col)
    se <- ms_ordinal_matrix_cell(coef_mat, i, se_col)
    stat <- ms_ordinal_matrix_cell(coef_mat, i, stat_col)
    p <- ms_ordinal_matrix_cell(coef_mat, i, p_col)
    if (is.na(p) && !is.na(stat)) {
      p <- ms_safe_numeric(2 * stats::pnorm(abs(stat), lower.tail = FALSE))
    }
    row <- list(
      term = term,
      estimate_name = "OR",
      estimate = if (!is.na(log_est)) ms_safe_numeric(exp(log_est)) else NA_real_,
      log_estimate = log_est,
      std_error = se,
      statistic = stat,
      statistic_label = "z",
      p_value = p,
      term_source = term
    )
    profile <- ms_ordinal_profile_ci_row(profile_ci, term)
    if (isTRUE(conf.int) && !is.null(profile)) {
      row$ci_lower <- ms_safe_numeric(exp(profile[[1L]]))
      row$ci_upper <- ms_safe_numeric(exp(profile[[2L]]))
      row$ci_method <- "profile_likelihood"
    } else if (isTRUE(conf.int) && is.finite(ci_crit) && !is.na(log_est) && !is.na(se)) {
      row$ci_lower <- ms_safe_numeric(exp(log_est - ci_crit * se))
      row$ci_upper <- ms_safe_numeric(exp(log_est + ci_crit * se))
      row$ci_method <- "wald"
      row$ci_fallback <- TRUE
    }
    row
  })
}

ms_ordinal_clm_profile_ci <- function(x, beta_names, conf.int = TRUE,
                                      conf.level = 0.95) {
  if (!isTRUE(conf.int)) return(NULL)
  beta_names <- as.character(beta_names %||% character(0))
  if (!length(beta_names)) return(NULL)
  ci <- tryCatch(
    suppressMessages(suppressWarnings(stats::confint(
      x,
      parm = beta_names,
      level = conf.level,
      type = "profile"
    ))),
    error = function(e) NULL
  )
  if (is.null(ci)) {
    ci <- tryCatch(
      suppressMessages(suppressWarnings(stats::confint(
        x,
        parm = beta_names,
        level = conf.level
      ))),
      error = function(e) NULL
    )
  }
  ms_ordinal_profile_ci_matrix(ci, beta_names)
}

ms_ordinal_profile_ci_matrix <- function(ci, beta_names) {
  if (is.null(ci)) return(NULL)
  mat <- tryCatch(as.matrix(ci), error = function(e) NULL)
  if (is.null(mat) || !length(mat) || ncol(mat) < 2L) return(NULL)
  suppressWarnings(storage.mode(mat) <- "numeric")
  if (all(!is.finite(mat[, seq_len(2L), drop = FALSE]))) return(NULL)
  if (is.null(rownames(mat))) {
    if (nrow(mat) == length(beta_names)) rownames(mat) <- beta_names
  }
  if (is.null(rownames(mat))) return(NULL)
  mat[, seq_len(2L), drop = FALSE]
}

ms_ordinal_profile_ci_row <- function(profile_ci, term) {
  if (is.null(profile_ci)) return(NULL)
  rn <- rownames(profile_ci) %||% character(0)
  i <- match(term, rn)
  if (is.na(i)) return(NULL)
  values <- ms_safe_numeric(profile_ci[i, seq_len(2L)])
  if (length(values) < 2L || any(is.na(values))) return(NULL)
  values
}

ms_ordinal_coefficient_ci_method <- function(coefficients) {
  methods <- unique(vapply(coefficients %||% list(), function(row) {
    if (!is.list(row)) return("")
    if (is.na(ms_safe_numeric(row$ci_lower %||% NA_real_)) ||
        is.na(ms_safe_numeric(row$ci_upper %||% NA_real_))) {
      return("")
    }
    as.character(row$ci_method %||% "")
  }, character(1)))
  methods <- methods[nzchar(methods)]
  if (!length(methods)) return(NULL)
  if (length(methods) == 1L) return(methods[[1L]])
  "mixed"
}

ms_ordinal_add_wald_fallback_ci <- function(coefficients, conf.int = TRUE,
                                           conf.level = 0.95) {
  if (!isTRUE(conf.int)) return(coefficients)
  ci_crit <- stats::qnorm((1 + conf.level) / 2)
  if (!is.finite(ci_crit)) return(coefficients)
  lapply(coefficients %||% list(), function(row) {
    if (!is.list(row)) return(row)
    if (isTRUE(row$estimate_unstable) || isTRUE(row$unstable)) return(row)
    ci_lower <- ms_safe_numeric(row$ci_lower %||% NA_real_)
    ci_upper <- ms_safe_numeric(row$ci_upper %||% NA_real_)
    if (!is.na(ci_lower) && !is.na(ci_upper)) return(row)
    log_est <- ms_safe_numeric(row$log_estimate %||% NA_real_)
    se <- ms_safe_numeric(row$std_error %||% NA_real_)
    if (is.na(log_est) || is.na(se)) return(row)
    row$ci_lower <- ms_safe_numeric(exp(log_est - ci_crit * se))
    row$ci_upper <- ms_safe_numeric(exp(log_est + ci_crit * se))
    row$ci_method <- "wald"
    row$ci_fallback <- TRUE
    row
  })
}

ms_ordinal_missing_inference_diagnostics <- function(coefficients) {
  rows <- Filter(is.list, coefficients %||% list())
  if (!length(rows)) return(list(coefficients = coefficients, warnings = list()))

  missing <- vapply(rows, ms_ordinal_row_inference_missing, logical(1))
  if (!any(missing)) return(list(coefficients = coefficients, warnings = list()))

  coefficients <- lapply(coefficients %||% list(), function(row) {
    if (!is.list(row) || !ms_ordinal_row_inference_missing(row)) return(row)
    row$inference_unavailable <- TRUE
    row$inference_note <- "singular information matrix"
    if (is.na(ms_ordinal_finite_numeric(row$ci_lower %||% NA_real_)) ||
        is.na(ms_ordinal_finite_numeric(row$ci_upper %||% NA_real_))) {
      if (is.null(row$ci_note) || !nzchar(as.character(row$ci_note))) {
        row$ci_note <- "not available: singular information"
      }
    }
    row
  })

  scope <- if (all(missing)) "all terms in this fit" else "some terms in this fit"
  warning <- list(
    type = "singular_information",
    severity = "warning",
    message = paste(
      "The ordinal model information matrix is singular, so standard errors,",
      "Wald tests, p values, and coefficient confidence intervals are unavailable for",
      paste0(scope, "."),
      "Blank inferential cells are intentional and should not be read as missing output."
    )
  )
  list(coefficients = coefficients, warnings = list(warning))
}

ms_ordinal_row_inference_missing <- function(row) {
  if (!is.list(row)) return(FALSE)
  is.na(ms_ordinal_finite_numeric(row$std_error %||% NA_real_))
}

ms_ordinal_finite_numeric <- function(value) {
  if (is.null(value) || length(value) == 0L) return(NA_real_)
  out <- suppressWarnings(as.numeric(value[[1L]]))
  if (!is.finite(out)) NA_real_ else out
}

ms_ordinal_clm_thresholds <- function(coef_mat, alpha_names) {
  rn <- rownames(coef_mat) %||% character(0)
  if (!length(alpha_names)) {
    alpha_names <- rn[grepl("\\|", rn)]
  }
  alpha_names <- intersect(alpha_names, rn)
  est_col <- ms_ordinal_matrix_col(coef_mat, "^Estimate$")
  se_col <- ms_ordinal_matrix_col(coef_mat, "Std\\. Error|Std\\.Error")
  stat_col <- ms_ordinal_matrix_col(coef_mat, "z value|t value|statistic")
  lapply(alpha_names, function(term) {
    i <- match(term, rn)
    list(
      threshold = term,
      estimate = ms_ordinal_matrix_cell(coef_mat, i, est_col),
      std_error = ms_ordinal_matrix_cell(coef_mat, i, se_col),
      statistic = ms_ordinal_matrix_cell(coef_mat, i, stat_col)
    )
  })
}

ms_ordinal_pick <- function(df, i, col) {
  if (is.null(df) || !nrow(df) || !col %in% names(df) || i > nrow(df)) return(NA_real_)
  ms_safe_numeric(df[[col]][[i]])
}

ms_ordinal_matrix_col <- function(mat, pattern) {
  hit <- grep(pattern, colnames(mat), ignore.case = TRUE, value = TRUE)
  if (length(hit)) hit[[1L]] else NA_character_
}

ms_ordinal_matrix_cell <- function(mat, i, col) {
  if (is.na(col) || is.na(i) || i < 1L || i > nrow(mat) || !col %in% colnames(mat)) {
    return(NA_real_)
  }
  ms_safe_numeric(mat[i, col])
}

ms_ordinal_any_ci <- function(coefficients) {
  any(vapply(coefficients, function(row) {
    is.list(row) &&
      !is.na(ms_safe_numeric(row$ci_lower %||% NA_real_)) &&
      !is.na(ms_safe_numeric(row$ci_upper %||% NA_real_))
  }, logical(1)))
}

ms_ordinal_outcome_levels <- function(x) {
  mf <- tryCatch(stats::model.frame(x), error = function(e) NULL)
  if (is.null(mf)) return(character(0))
  y <- tryCatch(stats::model.response(mf), error = function(e) NULL)
  if (is.null(y)) return(character(0))
  lev <- levels(y)
  if (is.null(lev)) lev <- sort(unique(as.character(y)))
  as.character(lev %||% character(0))
}

ms_ordinal_probability_curve_figure_data <- function(x, conf.level = 0.95,
                                                      grid_points = 80L,
                                                      max_levels = 8L) {
  tf <- tryCatch(stats::terms(x), error = function(e) NULL)
  mf <- tryCatch(stats::model.frame(x), error = function(e) NULL)
  if (is.null(tf) || is.null(mf) || nrow(mf) < 2L) return(NULL)

  term_labels <- attr(tf, "term.labels") %||% character(0)
  orders <- attr(tf, "order") %||% integer(0)
  if (length(term_labels) == 0L || length(orders) != length(term_labels)) return(NULL)
  if (any(orders != 1L)) return(NULL)

  outcome_levels <- ms_ordinal_outcome_levels(x)
  if (length(outcome_levels) < 2L) return(NULL)

  var_map <- ms_model_frame_variable_map(mf)
  term_info <- ms_glm_response_term_info(term_labels, mf, var_map, max_levels = max_levels)
  if (length(term_info) != length(term_labels)) return(NULL)

  term_roles <- tryCatch(ms_lm_term_roles(x), error = function(e) NULL)
  x_info <- ms_glm_response_selected_term(term_info, term_roles = term_roles)
  if (is.null(x_info) || !x_info$type %in% c("numeric", "categorical")) return(NULL)

  x_set <- NULL
  if (identical(x_info$type, "numeric")) {
    x_values <- ms_interaction_numeric_values(x_info$values)
    if (length(x_values) < 2L) return(NULL)
    x_range <- range(x_values, na.rm = TRUE)
    if (!all(is.finite(x_range)) || x_range[[1L]] == x_range[[2L]]) return(NULL)
    x_grid <- seq(x_range[[1L]], x_range[[2L]], length.out = max(12L, as.integer(grid_points)))
    interaction_kind <- "continuous_by_categorical"
  } else {
    x_set <- ms_interaction_categorical_set(x_info, max_levels = max_levels)
    if (is.null(x_set) || length(x_set$values) < 2L) return(NULL)
    x_range <- c(1, length(x_set$values))
    x_grid <- x_set$values
    interaction_kind <- "categorical_by_categorical"
  }

  base <- ms_glm_response_base_values(
    x = x,
    var_map = var_map,
    focal_variable = x_info$variable,
    max_levels = max_levels
  )
  if (is.null(base)) return(NULL)

  prediction_rows <- lapply(x_grid, function(x_value) {
    row <- base$values
    row[[x_info$variable]] <- x_value
    row
  })
  newdata <- ms_glm_response_newdata(prediction_rows, var_map)
  if (is.null(newdata) || nrow(newdata) != length(x_grid)) return(NULL)

  probabilities <- ms_ordinal_predict_probability_matrix(
    x,
    newdata = newdata,
    outcome_levels = outcome_levels
  )
  if (is.null(probabilities) || nrow(probabilities) != length(x_grid)) return(NULL)

  outcome_set <- list(
    id = "outcome_levels",
    label = "Outcome levels",
    rule = "ordinal_outcome_levels",
    values = outcome_levels,
    labels = outcome_levels,
    value_labels = outcome_levels,
    slices = rep(NA_character_, length(outcome_levels))
  )
  outcome_label <- ms_lm_outcome_label(x) %||% ""
  moderator_label <- if (nzchar(outcome_label)) {
    paste(outcome_label, "level")
  } else {
    "Outcome level"
  }

  grid <- list()
  for (level_index in seq_along(outcome_levels)) {
    level <- outcome_levels[[level_index]]
    for (x_index in seq_along(x_grid)) {
      probability <- ms_safe_numeric(probabilities[x_index, level_index])
      if (is.na(probability)) return(NULL)
      row <- list(
        x = if (!is.null(x_set)) x_index else ms_safe_numeric(x_grid[[x_index]]),
        moderator_value = as.character(level),
        moderator_label = as.character(level),
        estimate = probability
      )
      if (!is.null(x_set)) {
        row$x_value <- as.character(x_set$values[[x_index]])
        row$x_label <- x_set$labels[[x_index]]
      }
      grid[[length(grid) + 1L]] <- row
    }
  }
  if (length(grid) < length(x_grid) * length(outcome_levels)) return(NULL)

  link <- ms_ordinal_model_link(x)
  out <- list(
    interaction_term = x_info$term,
    interaction_kind = interaction_kind,
    source = "ordinal_predictions",
    mean_kind = "predicted_probability",
    variables = c(x_info$variable, "outcome_level"),
    x = list(
      variable = x_info$variable,
      term = x_info$term,
      label = x_info$label,
      type = x_info$type,
      range = I(ms_safe_numeric(x_range))
    ),
    moderator = list(
      variable = "outcome_level",
      term = "outcome_level",
      label = moderator_label,
      type = "categorical",
      levels = ms_interaction_set_levels(outcome_set)
    ),
    grid = grid,
    held_constant = base$held_constant,
    outcome = outcome_label %||% NULL,
    y_label = "Predicted probability",
    scale = "response",
    ci_level = conf.level,
    ci_method = "none",
    model_family = "ordinal",
    model_link = link,
    bounded_response = TRUE
  )
  if (!is.null(x_set)) {
    out$x$levels <- ms_interaction_set_levels(x_set)
    out$connect_levels <- FALSE
  }
  out
}

ms_ordinal_predict_probability_matrix <- function(x, newdata, outcome_levels = NULL) {
  predicted <- tryCatch(
    stats::predict(x, newdata = newdata, type = "probs"),
    error = function(e) NULL
  )
  if (is.null(predicted)) {
    predicted <- tryCatch(
      stats::predict(x, newdata = newdata, type = "prob"),
      error = function(e) NULL
    )
  }
  ms_ordinal_probability_matrix(predicted, n = nrow(newdata), outcome_levels = outcome_levels)
}

ms_ordinal_probability_matrix <- function(predicted, n, outcome_levels = NULL) {
  if (is.null(predicted)) return(NULL)
  candidates <- list(predicted)
  if (is.list(predicted) && !is.data.frame(predicted)) {
    for (key in c("fit", "prob", "probs", "probabilities")) {
      if (!is.null(predicted[[key]])) candidates[[length(candidates) + 1L]] <- predicted[[key]]
    }
  }
  outcome_levels <- as.character(outcome_levels %||% character(0))

  for (candidate in candidates) {
    mat <- ms_ordinal_probability_candidate_matrix(candidate, n = n,
                                                   outcome_levels = outcome_levels)
    if (!is.null(mat)) return(mat)
  }
  NULL
}

ms_ordinal_probability_candidate_matrix <- function(candidate, n, outcome_levels = character(0)) {
  if (is.null(candidate)) return(NULL)

  if (is.data.frame(candidate)) {
    if (length(outcome_levels) > 0L && all(outcome_levels %in% names(candidate))) {
      mat <- as.matrix(candidate[, outcome_levels, drop = FALSE])
    } else {
      numeric_cols <- vapply(candidate, is.numeric, logical(1))
      if (!any(numeric_cols)) return(NULL)
      mat <- as.matrix(candidate[, numeric_cols, drop = FALSE])
    }
  } else {
    mat <- tryCatch(as.matrix(candidate), error = function(e) NULL)
    if (is.null(mat)) return(NULL)
  }

  ok_numeric <- tryCatch({
    suppressWarnings(storage.mode(mat) <- "numeric")
    TRUE
  }, error = function(e) FALSE)
  if (!isTRUE(ok_numeric)) return(NULL)
  if (!is.numeric(mat) || !length(mat)) return(NULL)
  if (nrow(mat) != n && ncol(mat) == n) mat <- t(mat)
  if (nrow(mat) != n || ncol(mat) < 2L) return(NULL)

  if (length(outcome_levels) > 0L) {
    if (!is.null(colnames(mat)) && all(outcome_levels %in% colnames(mat))) {
      mat <- mat[, outcome_levels, drop = FALSE]
    } else if (ncol(mat) == length(outcome_levels)) {
      colnames(mat) <- outcome_levels
    } else {
      return(NULL)
    }
  }

  if (any(!is.finite(mat))) return(NULL)
  if (any(mat < -1e-8 | mat > 1 + 1e-8)) return(NULL)
  mat
}

ms_ordinal_model_link <- function(x) {
  if (inherits(x, "polr")) return(ms_ordinal_link_name(x$method %||% "logit"))
  ms_ordinal_link_name(x$link %||% "")
}

ms_ordinal_probability_curve_label <- function(plot) {
  if (identical(as.character(plot$x$type %||% ""), "categorical")) {
    return("Predicted probabilities")
  }
  "Predicted probability curves"
}

ms_ordinal_empty_cell_diagnostics <- function(x, coefficients) {
  coefficients <- coefficients %||% list()
  mf <- tryCatch(stats::model.frame(x), error = function(e) NULL)
  tf <- tryCatch(stats::terms(x), error = function(e) NULL)
  y <- tryCatch(stats::model.response(mf), error = function(e) NULL)
  if (is.null(mf) || is.null(tf) || is.null(y) || ncol(mf) < 2L) {
    return(list(coefficients = coefficients, warnings = list()))
  }

  term_labels <- attr(tf, "term.labels") %||% character(0)
  orders <- attr(tf, "order") %||% integer(0)
  if (!length(term_labels) || length(orders) != length(term_labels)) {
    return(list(coefficients = coefficients, warnings = list()))
  }

  main_terms <- term_labels[orders == 1L]
  main_variables <- unique(unlist(lapply(main_terms, function(term) {
    tryCatch(all.vars(stats::as.formula(paste0("~", term))),
             error = function(e) character(0))
  }), use.names = FALSE))
  if (!length(main_variables)) {
    return(list(coefficients = coefficients, warnings = list()))
  }

  outcome_label <- ms_lm_outcome_label(x) %||% names(mf)[[1L]] %||% "outcome"
  y_fac <- if (is.factor(y)) droplevels(y) else factor(y)
  if (nlevels(y_fac) < 2L) {
    return(list(coefficients = coefficients, warnings = list()))
  }

  warnings <- list()
  predictor_names <- names(mf)[-1L]
  for (variable in predictor_names) {
    if (!variable %in% main_variables) next
    values <- mf[[variable]]
    if (!is.factor(values) && !is.character(values) && !is.logical(values)) next
    x_fac <- if (is.factor(values)) droplevels(values) else factor(values)
    if (nlevels(x_fac) < 2L) next

    tab <- table(x_fac, y_fac)
    zero_index <- which(tab == 0L, arr.ind = TRUE)
    if (!nrow(zero_index)) next

    zero_cells <- lapply(seq_len(nrow(zero_index)), function(i) {
      list(
        predictor = variable,
        predictor_level = rownames(tab)[zero_index[i, 1L]],
        outcome = outcome_label,
        outcome_level = colnames(tab)[zero_index[i, 2L]],
        n = 0L
      )
    })
    affected_terms <- ms_ordinal_affected_terms(coefficients, variable)
    message <- ms_ordinal_empty_cell_message(variable, outcome_label, zero_cells)
    coefficients <- ms_ordinal_mark_empty_cell_coefficients(
      coefficients = coefficients,
      variable = variable,
      message = message
    )
    warnings[[length(warnings) + 1L]] <- list(
      type = "separation_or_boundary",
      severity = "warning",
      message = message,
      affected_terms = as.list(affected_terms),
      zero_cells = zero_cells
    )
  }

  list(coefficients = coefficients, warnings = warnings)
}

ms_ordinal_empty_cell_message <- function(variable, outcome, zero_cells) {
  cells <- vapply(zero_cells, function(cell) {
    paste0(cell$predictor_level %||% "level", " x ",
           cell$outcome_level %||% "outcome level")
  }, character(1))
  if (length(cells) > 4L) {
    cells <- c(cells[seq_len(4L)], paste0(length(cells) - 4L, " more"))
  }
  paste(
    "Empty outcome-by-predictor cells were detected for",
    paste0(variable, " by ", outcome, " (", paste(cells, collapse = "; "), ")."),
    "Maximum-likelihood proportional-odds estimates for affected terms may be",
    "unstable or unbounded; Wald z tests, p values, odds ratios, and predicted",
    "probabilities should be interpreted cautiously. Consider collapsing sparse",
    "levels, collecting more data, or using a penalized/Bayesian ordinal model."
  )
}

ms_ordinal_affected_terms <- function(coefficients, variable) {
  var_key <- ms_term_key(variable)
  terms <- vapply(coefficients %||% list(), function(row) {
    if (!is.list(row)) return(NA_character_)
    term <- as.character(row$term %||% row$term_source %||% "")
    term_source <- as.character(row$term_source %||% term)
    keys <- unique(c(ms_term_key(term), ms_term_key(term_source)))
    if (nzchar(var_key) && any(grepl(var_key, keys, fixed = TRUE))) {
      return(term)
    }
    NA_character_
  }, character(1))
  unique(stats::na.omit(terms))
}

ms_ordinal_mark_empty_cell_coefficients <- function(coefficients, variable, message) {
  var_key <- ms_term_key(variable)
  lapply(coefficients %||% list(), function(row) {
    if (!is.list(row) || !nzchar(var_key)) return(row)
    keys <- unique(c(ms_term_key(row$term %||% ""),
                     ms_term_key(row$term_source %||% "")))
    if (!any(grepl(var_key, keys, fixed = TRUE))) return(row)
    row$unstable <- TRUE
    row$estimate_unstable <- TRUE
    row$estimate_note <- "unstable/unbounded"
    row$diagnostic <- "separation_or_boundary"
    row$diagnostic_label <- "empty outcome-by-predictor cells"
    row$diagnostic_message <- message
    if (!is.null(row$statistic) && !is.na(ms_safe_numeric(row$statistic))) {
      row$wald_statistic <- row$statistic
    }
    if (!is.null(row$p_value) && !is.na(ms_safe_numeric(row$p_value))) {
      row$wald_p_value <- row$p_value
    }
    row$statistic <- "not interpreted"
    row$p_value <- "not interpreted"
    ci_lower <- ms_safe_numeric(row$ci_lower %||% NA_real_)
    ci_upper <- ms_safe_numeric(row$ci_upper %||% NA_real_)
    if (is.na(ci_lower) || is.na(ci_upper)) {
      row$ci_note <- "not estimable: separation"
    } else {
      row$unstable_ci_lower <- row$ci_lower
      row$unstable_ci_upper <- row$ci_upper
      row$ci_lower <- NA_real_
      row$ci_upper <- NA_real_
      row$ci_note <- "not interpreted: separation"
    }
    row
  })
}

ms_ordinal_add_probability_curve_warning <- function(fields, x) {
  reason <- ms_ordinal_probability_curve_unsupported_reason(x)
  if (is.null(reason)) return(fields)
  fields$model_warnings <- c(fields$model_warnings %||% list(), list(list(
    type = "unsupported_ordinal_probability_curve",
    severity = "info",
    message = reason
  )))
  fields
}

ms_ordinal_add_warning_rows <- function(fields, warnings) {
  warnings <- Filter(is.list, warnings %||% list())
  if (!length(warnings)) return(fields)
  fields$model_warnings <- c(fields$model_warnings %||% list(), warnings)
  fields
}

ms_ordinal_add_unstable_coefficient_warning <- function(fields, coefficients) {
  rows <- Filter(is.list, coefficients %||% list())
  if (!length(rows)) return(fields)
  unstable <- vapply(rows, ms_coefficient_ratio_row_unstable, logical(1))
  if (!any(unstable)) return(fields)
  existing <- fields$model_warnings %||% list()
  existing_types <- vapply(existing, function(row) {
    if (!is.list(row)) return("")
    as.character(row$type %||% "")
  }, character(1))
  if (any(existing_types == "separation_or_boundary")) return(fields)
  fields$model_warnings <- c(fields$model_warnings %||% list(), list(list(
    type = "separation_or_boundary",
    severity = "warning",
    message = paste(
      "The ordinal model produced extreme proportional-odds estimates;",
      "sparse outcome-by-predictor patterns or quasi-separation may make",
      "coefficient estimates unstable. Interpret odds ratios and Wald tests cautiously."
    )
  )))
  fields
}

ms_ordinal_probability_curve_unsupported_reason <- function(x) {
  tf <- tryCatch(stats::terms(x), error = function(e) NULL)
  if (is.null(tf)) return(NULL)
  term_labels <- attr(tf, "term.labels") %||% character(0)
  orders <- attr(tf, "order") %||% integer(0)
  if (!length(term_labels) || length(orders) != length(term_labels)) return(NULL)
  if (any(orders != 1L)) {
    return(paste(
      "Ordinal probability curves are currently generated for main-effect",
      "proportional-odds models only; this model includes interaction or",
      "higher-order terms, so Mellio rendered the coefficient table without",
      "a probability curve."
    ))
  }
  NULL
}
