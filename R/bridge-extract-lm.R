# R bridge — lm extractor.
#
# v0.1: emits the overall model F-test + R² as an inline card. Per-
# coefficient detail lives in raw_output for now; the structural
# "model_summary" archetype (coefficient table as a separate section)
# ships in v1.5.
#
# Card shape produced:
#   type:        "lm_model_summary"
#   type_label:  "Linear model (n predictors)"
#   statistic:   { name: "F", value: <F>, df: [df1, df2] }
#   p_value:     <overall model p>
#   estimate:    { name: "R²", value: <r.squared> }
#   raw_output:  output of print(summary(x))
#
# Schema: docs/STATS-R-BRIDGE-SCHEMA.md

#' @rdname mellio_payload
#' @export
mellio_payload.lm <- function(x, focal = NULL, controls = NULL, ..., .call = NULL) {
  # Capture user's call — prefer explicit .call from mellio_open,
  # then x$call from lm() itself (lm preserves it, unlike t.test),
  # then match.call()$x.
  call_str <- if (!is.null(.call)) {
    .call
  } else if (!is.null(x$call)) {
    paste(deparse(x$call, width.cutoff = 500L), collapse = " ")
  } else {
    user_call <- match.call()$x
    if (!is.null(user_call) && !identical(user_call, as.name("x"))) {
      paste(deparse(user_call, width.cutoff = 500L), collapse = " ")
    } else NA_character_
  }

  s <- summary(x)

  # F-statistic: s$fstatistic is c(value, df1, df2). Absent for
  # intercept-only models — treat that as NULL.
  fstat <- s$fstatistic
  if (is.null(fstat) || length(fstat) < 3) {
    f_value <- NA_real_
    df1     <- NA_real_
    df2     <- NA_real_
    p_value <- NA_real_
  } else {
    f_value <- ms_safe_numeric(unname(fstat["value"]))
    df1     <- ms_safe_numeric(unname(fstat["numdf"]))
    df2     <- ms_safe_numeric(unname(fstat["dendf"]))
    p_value <- ms_safe_numeric(stats::pf(f_value, df1, df2, lower.tail = FALSE))
  }

  # Count predictors (excluding intercept).
  n_pred <- max(0L, length(stats::coef(x)) - 1L)

  fields <- list(
    statistic = list(
      name  = "F",
      value = f_value,
      df    = I(c(df1, df2))
    ),
    p_value = p_value,
    estimate = list(
      name  = "R\u00b2",  # R-squared (\u escape: ASCII-safe source)
      value = ms_safe_numeric(s$r.squared)
    ),
    r_squared = ms_safe_numeric(s$r.squared),
    adj_r_squared = ms_safe_numeric(s$adj.r.squared),
    sigma = ms_safe_numeric(s$sigma),
    residual_df = ms_safe_numeric(stats::df.residual(x)),
    aic = ms_safe_numeric(stats::AIC(x)),
    bic = ms_safe_numeric(stats::BIC(x)),
    logLik = ms_safe_numeric(as.numeric(stats::logLik(x))),
    n = ms_safe_numeric(length(x$residuals)),
    conf_level = 0.95,
    coefficient_ci_method = "wald_t"
  )

  term_roles <- ms_lm_term_roles(x, focal = focal, controls = controls)
  if (!is.null(term_roles$outcome)) fields$outcome <- term_roles$outcome
  if (length(term_roles$terms) > 0L) fields$terms <- term_roles$terms
  if (length(term_roles$focal_terms) > 0L) fields$focal_terms <- term_roles$focal_terms
  if (length(term_roles$control_terms) > 0L) fields$control_terms <- term_roles$control_terms
  if (!is.null(term_roles$predictor)) fields$predictor <- term_roles$predictor
  if (!is.null(term_roles$model_kind)) fields$model_kind <- term_roles$model_kind
  model_term_tests <- ms_model_term_tests(x)
  if (length(model_term_tests) > 0L) fields$model_term_tests <- model_term_tests
  interaction_tests <- ms_lm_interaction_tests(x)
  if (length(interaction_tests) > 0L) fields$interaction_tests <- interaction_tests

  # P3: per-coefficient structured list. summary.lm() gives the matrix
  # (Estimate, Std. Error, t value, Pr(>|t|)); confint() gives 95% CIs.
  # confint may fail (e.g. rank-deficient fits) \u2014 degrade silently when
  # it does. Intercept rows are included; the paragraph builder skips
  # them when rendering "X predicted Y" sentences.
  coefs_mat <- tryCatch(stats::coef(s), error = function(e) NULL)
  if (!is.null(coefs_mat) && is.matrix(coefs_mat) && nrow(coefs_mat) > 0L) {
    ci_mat <- tryCatch(stats::confint(x), error = function(e) NULL)
    mm <- tryCatch(stats::model.matrix(x), error = function(e) NULL)
    mf <- tryCatch(stats::model.frame(x), error = function(e) NULL)
    y <- if (!is.null(mf)) tryCatch(stats::model.response(mf), error = function(e) NULL) else NULL
    y_sd <- if (is.numeric(y) && length(y) > 1L) stats::sd(y, na.rm = TRUE) else NA_real_
    std_scales <- NULL
    if (!is.null(mm) && is.matrix(mm) && is.finite(y_sd) && y_sd > 0) {
      x_sds <- apply(mm, 2L, stats::sd, na.rm = TRUE)
      std_scales <- ms_safe_numeric(x_sds / y_sd)
      names(std_scales) <- colnames(mm)
    }
    coef_assign <- if (!is.null(mm)) attr(mm, "assign") else NULL
    coef_label_map <- ms_coefficient_label_map(x, mm = mm, mf = mf)
    term_labels <- attr(stats::terms(x), "term.labels") %||% character(0)
    cn <- colnames(coefs_mat)
    n_col <- ncol(coefs_mat)
    # Bounds-checked cell read — non-standard lm subclasses (rlm, mlm)
    # can return a coefficient matrix without the expected columns;
    # an out-of-range index must degrade to NA, never crash.
    cell <- function(i, idx) {
      if (length(idx) == 1L && !is.na(idx) && idx >= 1L && idx <= n_col) {
        ms_safe_numeric(coefs_mat[i, idx])
      } else {
        NA_real_
      }
    }
    est_idx <- match("Estimate", cn, nomatch = 1L)
    se_idx  <- match("Std. Error", cn, nomatch = 2L)
    t_idx   <- match("t value", cn, nomatch = 3L)
    p_idx   <- match("Pr(>|t|)", cn, nomatch = 0L)
    terms_named <- rownames(coefs_mat)
    fields$coefficients <- lapply(seq_len(nrow(coefs_mat)), function(i) {
      row <- list(
        term            = terms_named[i],
        estimate_name   = "B",
        estimate        = cell(i, est_idx),
        std_error       = cell(i, se_idx),
        statistic       = cell(i, t_idx),
        statistic_label = "t",
        p_value         = cell(i, p_idx)
      )
      if (!is.null(coef_assign) && i <= length(coef_assign) &&
          coef_assign[[i]] > 0L && coef_assign[[i]] <= length(term_labels)) {
        row$term_source <- term_labels[[coef_assign[[i]]]]
      }
      row <- ms_apply_coefficient_label(row, coef_label_map)
      if (!is.null(ci_mat) && i <= nrow(ci_mat) && ncol(ci_mat) >= 2L) {
        row$ci_lower <- ms_safe_numeric(ci_mat[i, 1])
        row$ci_upper <- ms_safe_numeric(ci_mat[i, 2])
        if (!is.na(row$ci_lower) && !is.na(row$ci_upper)) {
          row$ci_method <- "wald_t"
        }
      }
      if (!is.null(std_scales) && !identical(terms_named[i], "(Intercept)")) {
        scale_idx <- match(terms_named[i], names(std_scales), nomatch = 0L)
        scale <- if (scale_idx > 0L) std_scales[[scale_idx]] else NA_real_
        if (ms_should_skip_standardized_coefficient(row, terms_named[i])) {
          row$std_estimate_skipped <- TRUE
          row$std_estimate_skipped_reason <- "interaction"
        } else if (!is.na(scale) && is.finite(scale) && scale > 0) {
          row$std_estimate <- ms_safe_numeric(row$estimate * scale)
          if (!is.null(row$std_error)) {
            row$std_std_error <- ms_safe_numeric(row$std_error * scale)
          }
          if (!is.null(row$ci_lower) && !is.null(row$ci_upper)) {
            row$std_ci_lower <- ms_safe_numeric(row$ci_lower * scale)
            row$std_ci_upper <- ms_safe_numeric(row$ci_upper * scale)
          }
        }
      }
      row
    })
    # Surface the statistic label at the fields level too — the JS
    # table renderer (coefficientStatisticTableLabel) checks rows[i]
    # first, but the fields-level value is the documented fallback
    # and makes the payload more self-describing.
    fields$statistic_label <- "t"
  }

  type_label <- "Linear Regression"

  data_prov <- tryCatch(
    ms_data_provenance(stats::model.frame(x)),
    error = function(e) NULL
  )
  provenance <- ms_provenance_basic()
  provenance <- ms_provenance_add_data(provenance, data_prov)

  figure_data <- list()
  interaction_plot <- ms_interaction_plot_data(x)
  if (!is.null(interaction_plot)) figure_data$interaction_plot <- interaction_plot
  if (is.null(interaction_plot)) {
    adjusted_means <- ms_ancova_adjusted_means_figure_data(x)
    if (is.null(adjusted_means)) {
      adjusted_means <- ms_factorial_main_effect_means_figure_data(
        x,
        focal = focal,
        controls = controls
      )
    }
    if (!is.null(adjusted_means)) figure_data$adjusted_means <- adjusted_means
  }
  if (length(figure_data) == 0L) figure_data <- NULL

  ms_build_envelope(
    type       = "lm_model_summary",
    type_label = type_label,
    call       = trimws(gsub("\\s+", " ", call_str)),
    fields     = fields,
    raw_output = ms_capture_output(s),
    provenance = provenance,
    figure_data = figure_data,
    packages = if (!is.null(figure_data$adjusted_means)) {
      ms_packages_basic(extras = "emmeans")
    } else {
      NULL
    }
  )
}

ms_lm_term_roles <- function(x, focal = NULL, controls = NULL) {
  tf <- tryCatch(stats::terms(x), error = function(e) NULL)
  if (is.null(tf)) {
    return(list(
      terms = list(),
      focal_terms = character(0),
      control_terms = character(0)
    ))
  }

  term_labels <- attr(tf, "term.labels") %||% character(0)
  mf <- tryCatch(stats::model.frame(x), error = function(e) NULL)
  outcome <- ms_lm_outcome_label(x)

  terms <- lapply(term_labels, function(term) {
    list(
      name = term,
      label = ms_model_clean_term(term),
      role = "term",
      type = ms_model_term_type(term, mf)
    )
  })

  roles <- ms_assign_term_roles(
    terms,
    focal = focal,
    controls = controls,
    infer_covariate_roles = FALSE
  )
  terms <- roles$terms

  focal_terms <- vapply(Filter(function(t) identical(t$role, "focal"), terms),
                        function(t) t$name, character(1))
  control_terms <- vapply(Filter(function(t) identical(t$role, "control"), terms),
                          function(t) t$name, character(1))

  predictor_terms <- if (length(focal_terms) > 0L) {
    Filter(function(t) identical(t$role, "focal"), terms)
  } else {
    Filter(function(t) !identical(t$role, "control"), terms)
  }
  # Use $name (raw identifier with underscores) for the field value, not
  # $label (humanised). The extracted-fields panel lists variable names
  # verbatim; prose italicisation is applied in the browser.
  predictor <- ms_model_term_phrase(vapply(predictor_terms, function(t) {
    t$name %||% t$label
  }, character(1)))

  list(
    outcome = outcome,
    terms = terms,
    focal_terms = focal_terms,
    control_terms = control_terms,
    predictor = if (nzchar(predictor)) predictor else NULL,
    model_kind = roles$model_kind
  )
}

ms_lm_interaction_tests <- function(x) {
  if (!inherits(x, "lm") || inherits(x, "glm")) return(list())
  tf <- tryCatch(stats::terms(x), error = function(e) NULL)
  if (is.null(tf)) return(list())

  term_labels <- attr(tf, "term.labels") %||% character(0)
  orders <- attr(tf, "order") %||% integer(0)
  if (!length(term_labels) || length(term_labels) != length(orders)) return(list())

  interaction_terms <- term_labels[orders >= 2L & grepl(":", term_labels, fixed = TRUE)]
  if (!length(interaction_terms)) return(list())

  tests <- tryCatch(stats::drop1(x, test = "F"), error = function(e) NULL)
  if (is.null(tests) || !inherits(tests, "data.frame") || !nrow(tests)) return(list())

  rn <- rownames(tests) %||% character(0)
  df2 <- ms_safe_numeric(stats::df.residual(x))
  df_col <- match("Df", names(tests), nomatch = 0L)
  ss_col <- match("Sum of Sq", names(tests), nomatch = 0L)
  f_col <- match("F value", names(tests), nomatch = 0L)
  p_col <- grep("^Pr\\(>F\\)$", names(tests), perl = TRUE)
  p_col <- if (length(p_col)) p_col[[1L]] else 0L

  out <- lapply(interaction_terms, function(term) {
    idx <- match(term, rn, nomatch = 0L)
    if (!idx) return(NULL)
    df1 <- if (df_col) ms_safe_numeric(tests[[df_col]][idx]) else NA_real_
    f_value <- if (f_col) ms_safe_numeric(tests[[f_col]][idx]) else NA_real_
    p_value <- if (p_col) ms_safe_numeric(tests[[p_col]][idx]) else NA_real_
    if (is.na(df1) || is.na(f_value)) return(NULL)

    row <- list(
      term = term,
      statistic = list(
        name = "F",
        value = f_value,
        df = I(c(df1, df2))
      ),
      p_value = p_value,
      method = "drop1_f"
    )
    if (ss_col) row$sum_sq <- ms_safe_numeric(tests[[ss_col]][idx])
    if (!is.na(f_value) && !is.na(df1) && !is.na(df2) && df2 > 0) {
      row$effect <- list(
        name = "eta_sq_partial",
        value = ms_safe_numeric((f_value * df1) / ((f_value * df1) + df2))
      )
    }
    row
  })

  Filter(Negate(is.null), out)
}

ms_model_term_tests <- function(x) {
  if (!inherits(x, c("lm", "glm"))) return(list())
  if (inherits(x, c("lmerMod", "glmerMod", "merMod", "gam", "rlm"))) return(list())

  tf <- tryCatch(stats::terms(x), error = function(e) NULL)
  if (is.null(tf)) return(list())

  term_labels <- attr(tf, "term.labels") %||% character(0)
  orders <- attr(tf, "order") %||% integer(0)
  if (!length(term_labels)) return(list())
  if (length(orders) != length(term_labels)) {
    orders <- rep(1L, length(term_labels))
  }

  is_glm <- inherits(x, "glm")
  drop_rows <- if (is_glm) {
    ms_model_term_tests_drop1_chisq(x, term_labels)
  } else {
    ms_model_term_tests_drop1_f(x, term_labels)
  }
  seq_rows <- if (is_glm) {
    ms_model_term_tests_anova_chisq(x, term_labels)
  } else {
    ms_model_term_tests_anova_f(x, term_labels)
  }

  out <- vector("list", length(term_labels))
  for (i in seq_along(term_labels)) {
    term <- term_labels[[i]]
    row <- drop_rows[[term]]
    if (is.null(row)) row <- seq_rows[[term]]
    if (!is.null(row)) {
      row$term_type <- if (orders[[i]] >= 2L || grepl(":", term, fixed = TRUE)) {
        "interaction"
      } else {
        "main"
      }
      out[[i]] <- row
    }
  }

  Filter(Negate(is.null), out)
}

ms_model_term_test_base_row <- function(term, statistic_name, statistic_value,
                                        df = NULL, p_value = NA_real_,
                                        method = NULL) {
  if (is.na(statistic_value)) return(NULL)
  statistic <- list(
    name = statistic_name,
    value = ms_safe_numeric(statistic_value)
  )
  if (!is.null(df) && length(df) > 0L && !all(is.na(df))) {
    statistic$df <- I(ms_safe_numeric(df))
  }
  row <- list(
    term = term,
    label = ms_model_term_test_label(term),
    statistic = statistic,
    p_value = ms_safe_numeric(p_value),
    method = method,
    test_type = "omnibus"
  )
  row
}

ms_model_term_test_label <- function(term) {
  term <- trimws(as.character(term %||% ""))
  if (!nzchar(term)) return("")
  if (grepl(":", term, fixed = TRUE)) {
    parts <- strsplit(term, ":", fixed = TRUE)[[1]]
    parts <- vapply(parts, ms_model_clean_term, character(1))
    return(paste(parts[nzchar(parts)], collapse = " \u00d7 "))
  }
  ms_model_clean_term(term)
}

ms_model_term_tests_drop1_f <- function(x, term_labels) {
  tests <- tryCatch(stats::drop1(x, test = "F"), error = function(e) NULL)
  if (is.null(tests) || !inherits(tests, "data.frame") || !nrow(tests)) return(list())

  rn <- rownames(tests) %||% character(0)
  df2 <- ms_safe_numeric(stats::df.residual(x))
  df_col <- match("Df", names(tests), nomatch = 0L)
  ss_col <- match("Sum of Sq", names(tests), nomatch = 0L)
  f_col <- match("F value", names(tests), nomatch = 0L)
  p_col <- grep("^Pr\\(>F\\)$", names(tests), perl = TRUE)
  p_col <- if (length(p_col)) p_col[[1L]] else 0L

  out <- list()
  for (term in term_labels) {
    idx <- match(term, rn, nomatch = 0L)
    if (!idx || !df_col || !f_col) next
    df1 <- ms_safe_numeric(tests[[df_col]][idx])
    f_value <- ms_safe_numeric(tests[[f_col]][idx])
    p_value <- if (p_col) ms_safe_numeric(tests[[p_col]][idx]) else NA_real_
    row <- ms_model_term_test_base_row(
      term = term,
      statistic_name = "F",
      statistic_value = f_value,
      df = c(df1, df2),
      p_value = p_value,
      method = "drop1_f"
    )
    if (is.null(row)) next
    if (ss_col) row$sum_sq <- ms_safe_numeric(tests[[ss_col]][idx])
    if (!is.na(f_value) && !is.na(df1) && !is.na(df2) && df2 > 0) {
      row$effect <- list(
        name = "eta_sq_partial",
        value = ms_safe_numeric((f_value * df1) / ((f_value * df1) + df2))
      )
    }
    out[[term]] <- row
  }
  out
}

ms_model_term_tests_anova_f <- function(x, term_labels) {
  tests <- tryCatch(stats::anova(x), error = function(e) NULL)
  if (is.null(tests) || !inherits(tests, "data.frame") || !nrow(tests)) return(list())

  rn <- rownames(tests) %||% character(0)
  is_resid <- grepl("^Residuals?$", rn, ignore.case = TRUE)
  resid_idx <- which(is_resid)
  df2 <- if (length(resid_idx) && "Df" %in% names(tests)) {
    ms_safe_numeric(tests[["Df"]][resid_idx[[1L]]])
  } else {
    ms_safe_numeric(stats::df.residual(x))
  }
  resid_sum_sq <- if (length(resid_idx) && "Sum Sq" %in% names(tests)) {
    ms_safe_numeric(tests[["Sum Sq"]][resid_idx[[1L]]])
  } else NA_real_

  df_col <- match("Df", names(tests), nomatch = 0L)
  ss_col <- match("Sum Sq", names(tests), nomatch = 0L)
  f_col <- if ("F" %in% names(tests)) match("F", names(tests), nomatch = 0L) else match("F value", names(tests), nomatch = 0L)
  p_col <- grep("^Pr\\(>F\\)$", names(tests), perl = TRUE)
  p_col <- if (length(p_col)) p_col[[1L]] else 0L

  out <- list()
  for (term in term_labels) {
    idx <- match(term, rn, nomatch = 0L)
    if (!idx || !df_col || !f_col) next
    df1 <- ms_safe_numeric(tests[[df_col]][idx])
    f_value <- ms_safe_numeric(tests[[f_col]][idx])
    p_value <- if (p_col) ms_safe_numeric(tests[[p_col]][idx]) else NA_real_
    row <- ms_model_term_test_base_row(
      term = term,
      statistic_name = "F",
      statistic_value = f_value,
      df = c(df1, df2),
      p_value = p_value,
      method = "anova_f_sequential"
    )
    if (is.null(row)) next
    if (ss_col) {
      sum_sq <- ms_safe_numeric(tests[[ss_col]][idx])
      row$sum_sq <- sum_sq
      if (!is.na(sum_sq) && !is.na(resid_sum_sq) && (sum_sq + resid_sum_sq) > 0) {
        row$effect <- list(
          name = "eta_sq_partial",
          value = ms_safe_numeric(sum_sq / (sum_sq + resid_sum_sq))
        )
      }
    }
    row$ss_type <- "type_i_sequential"
    out[[term]] <- row
  }
  out
}

ms_model_term_tests_drop1_chisq <- function(x, term_labels) {
  tests <- tryCatch(stats::drop1(x, test = "Chisq"), error = function(e) NULL)
  if (is.null(tests) || !inherits(tests, "data.frame") || !nrow(tests)) return(list())

  rn <- rownames(tests) %||% character(0)
  df_col <- match("Df", names(tests), nomatch = 0L)
  if (!df_col) df_col <- match("npar", names(tests), nomatch = 0L)
  stat_col <- match("LRT", names(tests), nomatch = 0L)
  if (!stat_col) stat_col <- match("Deviance", names(tests), nomatch = 0L)
  p_col <- grep("^Pr\\(>?Chi\\)$", names(tests), perl = TRUE)
  p_col <- if (length(p_col)) p_col[[1L]] else 0L

  out <- list()
  for (term in term_labels) {
    idx <- match(term, rn, nomatch = 0L)
    if (!idx || !df_col || !stat_col) next
    df <- ms_safe_numeric(tests[[df_col]][idx])
    chisq <- ms_safe_numeric(tests[[stat_col]][idx])
    p_value <- if (p_col) ms_safe_numeric(tests[[p_col]][idx]) else NA_real_
    row <- ms_model_term_test_base_row(
      term = term,
      statistic_name = "chi2",
      statistic_value = chisq,
      df = df,
      p_value = p_value,
      method = "drop1_chisq"
    )
    if (!is.null(row)) out[[term]] <- row
  }
  out
}

ms_model_term_tests_anova_chisq <- function(x, term_labels) {
  tests <- tryCatch(stats::anova(x, test = "Chisq"), error = function(e) NULL)
  if (is.null(tests) || !inherits(tests, "data.frame") || !nrow(tests)) return(list())

  rn <- rownames(tests) %||% character(0)
  df_col <- match("Df", names(tests), nomatch = 0L)
  stat_col <- match("Deviance", names(tests), nomatch = 0L)
  p_col <- grep("^Pr\\(>Chi\\)$", names(tests), perl = TRUE)
  p_col <- if (length(p_col)) p_col[[1L]] else 0L

  out <- list()
  for (term in term_labels) {
    idx <- match(term, rn, nomatch = 0L)
    if (!idx || !df_col || !stat_col) next
    df <- ms_safe_numeric(tests[[df_col]][idx])
    chisq <- ms_safe_numeric(tests[[stat_col]][idx])
    p_value <- if (p_col) ms_safe_numeric(tests[[p_col]][idx]) else NA_real_
    row <- ms_model_term_test_base_row(
      term = term,
      statistic_name = "chi2",
      statistic_value = chisq,
      df = df,
      p_value = p_value,
      method = "anova_chisq_sequential"
    )
    if (is.null(row)) next
    row$ss_type <- "type_i_sequential"
    out[[term]] <- row
  }
  out
}

ms_lm_outcome_label <- function(x) {
  f <- tryCatch(stats::formula(x), error = function(e) NULL)
  if (is.null(f) || length(f) < 3L) return(NULL)
  label <- ms_model_clean_term(paste(deparse(f[[2]], width.cutoff = 500L), collapse = " "))
  if (nzchar(label)) label else NULL
}

ms_assign_term_roles <- function(terms, focal = NULL, controls = NULL,
                                 infer_covariate_roles = FALSE) {
  if (!length(terms)) return(list(terms = terms, model_kind = NULL))

  focal_keys <- ms_role_keys(focal)
  control_keys <- ms_role_keys(controls)
  has_focal <- length(focal_keys) > 0L
  has_controls <- length(control_keys) > 0L

  terms <- lapply(terms, function(term) {
    keys <- ms_term_match_keys(term)
    if (length(intersect(keys, focal_keys)) > 0L) {
      term$role <- "focal"
    } else if (length(intersect(keys, control_keys)) > 0L) {
      term$role <- "control"
    }
    term
  })

  if (has_focal) {
    terms <- lapply(terms, function(term) {
      if (identical(term$role, "term") && !identical(term$type, "interaction")) {
        term$role <- "control"
      }
      term
    })
  } else if (has_controls) {
    terms <- lapply(terms, function(term) {
      if (identical(term$role, "term") && !identical(term$type, "interaction")) {
        term$role <- "focal"
      }
      term
    })
  } else if (isTRUE(infer_covariate_roles)) {
    is_interaction <- vapply(terms, function(t) identical(t$type, "interaction"), logical(1))
    is_factor <- vapply(terms, function(t) identical(t$type, "factor"), logical(1))
    is_supported_factor <- is_factor & vapply(terms, function(t) {
      ms_ancova_supported_factor_term(
        t$name %||% "",
        ms_ancova_single_variable(t$name %||% "")
      )
    }, logical(1))
    is_numeric <- vapply(terms, function(t) identical(t$type, "numeric"), logical(1))
    is_other <- !(is_supported_factor | is_numeric)
    if (!any(is_interaction) &&
        sum(is_supported_factor) == 1L &&
        sum(is_numeric) >= 1L &&
        !any(is_other)) {
      terms <- Map(function(term, factor_term, numeric_term) {
        if (factor_term) term$role <- "focal"
        if (numeric_term) term$role <- "control"
        term
      }, terms, is_supported_factor, is_numeric)
    } else {
      terms <- lapply(terms, function(term) {
        if (identical(term$role, "term") && !identical(term$type, "interaction")) {
          term$role <- "focal"
        }
        term
      })
    }
  } else {
    terms <- lapply(terms, function(term) {
      if (identical(term$role, "term") && !identical(term$type, "interaction")) {
        term$role <- "focal"
      }
      term
    })
  }

  has_control_role <- any(vapply(terms, function(t) identical(t$role, "control"), logical(1)))
  has_focal_role <- any(vapply(terms, function(t) identical(t$role, "focal"), logical(1)))

  model_kind <- NULL
  if (has_focal_role && has_control_role) {
    model_kind <- "controlled_regression"
  }

  list(terms = terms, model_kind = model_kind)
}

ms_role_keys <- function(values) {
  if (is.null(values) || length(values) == 0L) return(character(0))
  unique(vapply(as.character(values), ms_term_key, character(1)))
}

ms_term_match_keys <- function(term) {
  unique(ms_term_key(c(term$name %||% "", term$label %||% "")))
}

ms_term_key <- function(value) {
  value <- ms_model_clean_term(value)
  tolower(gsub("[^[:alnum:]]+", "", value))
}

ms_model_clean_term <- function(value) {
  value <- trimws(as.character(value %||% ""))
  value <- sub("^`(.+)`$", "\\1", value)
  value <- sub("^.*\\$", "", value)
  value <- gsub("\\[[^]]*\\]", "", value)
  value <- gsub(
    "\\b(as\\.factor|factor|scale|I)\\(([^()]+)\\)",
    "\\2",
    value,
    perl = TRUE
  )
  value <- gsub(":", " by ", value, fixed = TRUE)
  value <- gsub("[_.]+", " ", value)
  value <- gsub("\\s+", " ", value)
  trimws(value)
}

ms_model_term_type <- function(term, mf = NULL) {
  if (grepl(":", term, fixed = TRUE)) return("interaction")
  if (is.null(mf)) return("other")
  # Look up the literal term name in the model frame first — model.frame()
  # stores transformed columns under their full call ("factor(cyl)",
  # "cut(wt, 2)", "I(cyl > 4)") not under the bare variable name. Falling
  # back to all.vars() handles bare-name terms.
  x <- if (term %in% names(mf)) mf[[term]] else NULL
  if (is.null(x)) {
    vars <- tryCatch(all.vars(stats::as.formula(paste0("~", term))),
                     error = function(e) character(0))
    if (length(vars) == 1L && vars[[1]] %in% names(mf)) x <- mf[[vars[[1]]]]
  }
  if (is.null(x)) return("other")
  if (is.factor(x) || is.character(x) || is.logical(x)) return("factor")
  if (is.numeric(x) || is.integer(x)) return("numeric")
  "other"
}

ms_coefficient_label_map <- function(x, mm = NULL, mf = NULL) {
  mm <- mm %||% tryCatch(stats::model.matrix(x), error = function(e) NULL)
  if (is.null(mm) || !is.matrix(mm)) return(list())

  cols <- colnames(mm) %||% character(0)
  assign <- attr(mm, "assign")
  if (!length(cols) || is.null(assign) || length(assign) != length(cols)) return(list())

  tf <- tryCatch(stats::terms(x), error = function(e) NULL)
  term_labels <- if (!is.null(tf)) attr(tf, "term.labels") %||% character(0) else character(0)
  if (!length(term_labels)) return(list())

  mf <- mf %||% tryCatch(stats::model.frame(x), error = function(e) NULL)
  if (is.null(mf)) return(list())

  out <- list()
  for (i in seq_along(cols)) {
    idx <- assign[[i]]
    if (is.na(idx) || idx <= 0L || idx > length(term_labels)) next
    term <- term_labels[[idx]]
    meta <- ms_coefficient_factor_contrast_metadata(cols[[i]], term, mf)
    if (!is.null(meta)) out[[cols[[i]]]] <- meta
  }
  out
}

ms_apply_coefficient_label <- function(row, label_map) {
  if (!is.list(row) || !length(label_map)) return(row)
  term <- as.character(row$term %||% "")
  if (!nzchar(term) || is.null(label_map[[term]])) return(row)
  meta <- label_map[[term]]
  if (!is.null(meta$label) && nzchar(meta$label)) row$label <- meta$label
  if (!is.null(meta$contrast_variable)) row$contrast_variable <- meta$contrast_variable
  if (!is.null(meta$contrast_level)) row$contrast_level <- meta$contrast_level
  if (!is.null(meta$contrast_reference)) row$contrast_reference <- meta$contrast_reference
  if (!is.null(meta$contrast_type)) row$contrast_type <- meta$contrast_type
  if (!is.null(meta$contrast_column)) row$contrast_column <- meta$contrast_column
  row
}

ms_should_skip_standardized_coefficient <- function(row, term) {
  candidates <- c(
    as.character(term %||% ""),
    as.character(row$term %||% ""),
    as.character(row$term_source %||% "")
  )
  any(grepl(":", candidates, fixed = TRUE))
}

ms_coefficient_factor_contrast_metadata <- function(coef_name, term, mf) {
  coef_name <- as.character(coef_name %||% "")
  term <- as.character(term %||% "")
  if (!nzchar(coef_name) || !nzchar(term)) {
    return(NULL)
  }
  if (grepl(":", term, fixed = TRUE)) {
    return(ms_coefficient_interaction_contrast_metadata(coef_name, term, mf))
  }

  vars <- tryCatch(all.vars(stats::as.formula(paste0("~", term))),
                   error = function(e) character(0))
  vars <- unique(vars[nzchar(vars)])
  if (length(vars) != 1L) return(NULL)

  values <- if (term %in% names(mf)) {
    mf[[term]]
  } else if (vars[[1L]] %in% names(mf)) {
    mf[[vars[[1L]]]]
  } else {
    NULL
  }
  if (is.null(values)) return(NULL)
  if (!(is.factor(values) || is.character(values) || is.logical(values))) {
    return(NULL)
  }
  factor_values <- if (is.factor(values)) values else factor(values)
  levels <- levels(factor_values)
  if (length(levels) < 2L) return(NULL)

  contrast <- tryCatch(stats::contrasts(factor_values, contrasts = TRUE),
                       error = function(e) NULL)
  if (is.null(contrast) || !is.matrix(contrast) || ncol(contrast) < 1L) {
    return(NULL)
  }
  contrast_cols <- colnames(contrast) %||% as.character(seq_len(ncol(contrast)))
  contrast_col <- ms_coefficient_contrast_column(coef_name, term, vars[[1L]], contrast_cols)
  if (!nzchar(contrast_col)) return(NULL)
  contrast_idx <- match(contrast_col, contrast_cols, nomatch = 0L)
  if (!contrast_idx) return(NULL)

  values <- suppressWarnings(as.numeric(contrast[, contrast_idx]))
  names(values) <- rownames(contrast) %||% levels
  variable_label <- ms_model_clean_term(term)
  if (!nzchar(variable_label)) variable_label <- vars[[1L]]

  treatment <- ms_coefficient_treatment_contrast_label(
    variable_label,
    values,
    contrast
  )
  if (!is.null(treatment)) {
    treatment$contrast_variable <- vars[[1L]]
    treatment$contrast_column <- contrast_col
    return(treatment)
  }

  sum_label <- ms_coefficient_sum_contrast_label(variable_label, values)
  if (!is.null(sum_label)) {
    sum_label$contrast_variable <- vars[[1L]]
    sum_label$contrast_column <- contrast_col
    return(sum_label)
  }

  list(
    label = paste0(variable_label, ": contrast ", contrast_col),
    contrast_variable = vars[[1L]],
    contrast_column = contrast_col,
    contrast_type = "factor_contrast"
  )
}

ms_coefficient_interaction_contrast_metadata <- function(coef_name, term, mf) {
  coef_parts <- strsplit(coef_name, ":", fixed = TRUE)[[1]]
  term_parts <- strsplit(term, ":", fixed = TRUE)[[1]]
  if (length(coef_parts) != length(term_parts)) return(NULL)

  labels <- character(0)
  contrast_types <- character(0)
  for (i in seq_along(term_parts)) {
    part_meta <- ms_coefficient_factor_contrast_metadata(coef_parts[[i]], term_parts[[i]], mf)
    if (!is.null(part_meta)) {
      variable_label <- ms_model_clean_term(part_meta$contrast_variable %||% term_parts[[i]])
      if (!nzchar(variable_label)) variable_label <- ms_model_clean_term(term_parts[[i]])
      if (identical(part_meta$contrast_type, "sum")) {
        labels <- c(labels, paste0(variable_label, " (sum contrast)"))
      } else if (!is.null(part_meta$contrast_level) && !is.null(part_meta$contrast_reference)) {
        labels <- c(labels, paste0(
          variable_label,
          ": ",
          part_meta$contrast_level,
          " vs. ",
          part_meta$contrast_reference
        ))
      } else {
        labels <- c(labels, part_meta$label %||% variable_label)
      }
      contrast_types <- c(contrast_types, part_meta$contrast_type %||% "factor_contrast")
    } else {
      labels <- c(labels, ms_model_clean_term(term_parts[[i]]))
    }
  }

  labels <- labels[nzchar(labels)]
  if (!length(labels)) return(NULL)
  out <- list(label = paste(labels, collapse = " \u00d7 "))
  if (length(contrast_types)) {
    out$contrast_type <- paste(unique(contrast_types[nzchar(contrast_types)]), collapse = "+")
  }
  out
}

ms_coefficient_contrast_column <- function(coef_name, term, variable, contrast_cols) {
  candidates <- unique(c(term, variable, ms_model_clean_term(term)))
  for (prefix in candidates[nzchar(candidates)]) {
    if (identical(coef_name, prefix)) next
    if (startsWith(coef_name, prefix)) {
      suffix <- substring(coef_name, nchar(prefix) + 1L)
      if (suffix %in% contrast_cols) return(suffix)
    }
  }
  if (coef_name %in% contrast_cols) return(coef_name)
  ""
}

ms_coefficient_treatment_contrast_label <- function(variable_label, values, contrast) {
  eps <- 1e-8
  focal_idx <- which(abs(values - 1) <= eps)
  if (length(focal_idx) != 1L) return(NULL)

  col_zero <- apply(abs(contrast) <= eps, 1L, all)
  reference_idx <- which(col_zero)
  if (length(reference_idx) != 1L) return(NULL)
  if (!all(abs(values[-focal_idx]) <= eps)) return(NULL)

  level <- names(values)[[focal_idx]]
  reference <- rownames(contrast)[[reference_idx]]
  list(
    label = paste0(variable_label, ": ", level, " vs. ", reference),
    contrast_level = level,
    contrast_reference = reference,
    contrast_type = "treatment"
  )
}

ms_coefficient_sum_contrast_label <- function(variable_label, values) {
  eps <- 1e-8
  rounded <- round(values)
  if (!all(abs(values - rounded) <= eps) || !all(rounded %in% c(-1, 0, 1))) {
    return(NULL)
  }
  positive_idx <- which(rounded == 1)
  if (length(positive_idx) != 1L || !any(rounded == -1)) return(NULL)

  level <- names(values)[[positive_idx]]
  list(
    label = paste0(variable_label, ": ", level, " (sum contrast)"),
    contrast_level = level,
    contrast_type = "sum"
  )
}

ms_model_term_phrase <- function(items, max_items = 3L, overflow = "other predictors") {
  # Comma-separated identifier list for the extracted-fields "predictors"
  # value. We deliberately drop the English "and" join here: the field is
  # a reference list of variable names, not running prose. The browser
  # prose generators (predictorRoleText / associationVerbText in
  # stats-bridge.js) already detect commas and switch to plural verbs.
  items <- items[nzchar(items)]
  if (!length(items)) return("")
  shown <- head(items, max_items)
  if (length(items) > length(shown) && nzchar(overflow)) shown <- c(shown, overflow)
  paste(shown, collapse = ", ")
}

ms_interaction_plot_data <- function(x, grid_points = 80L, max_levels = 6L,
                                     conf.level = 0.95) {
  if (!inherits(x, c("lm", "glm"))) return(NULL)
  if (inherits(x, c("lmerMod", "glmerMod", "merMod", "gam", "rlm"))) return(NULL)

  tf <- tryCatch(stats::terms(x), error = function(e) NULL)
  mf <- tryCatch(stats::model.frame(x), error = function(e) NULL)
  if (is.null(tf) || is.null(mf) || nrow(mf) < 2L) return(NULL)

  term_labels <- attr(tf, "term.labels") %||% character(0)
  orders <- attr(tf, "order") %||% integer(0)
  if (length(term_labels) == 0L || length(orders) != length(term_labels)) return(NULL)
  if (any(orders > 3L)) return(NULL)

  interaction_terms <- term_labels[orders >= 2L & grepl(":", term_labels, fixed = TRUE)]
  if (length(interaction_terms) == 0L) return(NULL)

  three_way_terms <- term_labels[orders == 3L & grepl(":", term_labels, fixed = TRUE)]
  two_way_terms <- term_labels[orders == 2L & grepl(":", term_labels, fixed = TRUE)]
  if (length(three_way_terms) == 1L) {
    interaction_term <- three_way_terms[[1]]
  } else if (length(three_way_terms) == 0L && length(two_way_terms) == 1L) {
    interaction_term <- two_way_terms[[1]]
  } else {
    return(NULL)
  }

  parts <- strsplit(interaction_term, ":", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  if (!(length(parts) %in% c(2L, 3L)) || any(!nzchar(parts))) return(NULL)
  interaction_variables <- vapply(parts, function(part) {
    ms_interaction_component_variable(part) %||% ""
  }, character(1))
  if (any(!nzchar(interaction_variables)) ||
      length(unique(interaction_variables)) != length(interaction_variables)) {
    return(NULL)
  }
  if (!ms_interaction_terms_within_components(interaction_terms, interaction_variables)) {
    return(NULL)
  }
  if (!all(vapply(parts, ms_interaction_component_supported, logical(1)))) {
    return(NULL)
  }

  var_map <- ms_model_frame_variable_map(mf)
  components <- lapply(parts, ms_interaction_component_info,
                       var_map = var_map, max_levels = max_levels)
  if (any(vapply(components, is.null, logical(1)))) return(NULL)

  is_numeric <- vapply(components, function(info) identical(info$type, "numeric"), logical(1))
  is_categorical <- vapply(components, function(info) identical(info$type, "categorical"), logical(1))
  interaction_kind <- NULL
  x_set <- NULL
  moderator_sets <- NULL
  facet_info <- NULL
  facet_set <- NULL
  prediction_specs <- NULL

  if (length(parts) == 2L) {
    if (sum(is_numeric) == 1L && sum(is_categorical) == 1L) {
      x_info <- components[[which(is_numeric)[[1]]]]
      moderator_info <- components[[which(is_categorical)[[1]]]]
      moderator_set <- ms_interaction_categorical_set(moderator_info, max_levels = max_levels)
      if (is.null(moderator_set)) return(NULL)
      moderator_sets <- list(moderator_set)
      prediction_specs <- list(list(moderator_set = moderator_set, facet_set = NULL))
      interaction_kind <- "continuous_by_categorical"
    } else if (sum(is_numeric) == 2L && sum(is_categorical) == 0L) {
      # For continuous-by-continuous interactions, keep the formula order:
      # first term on the x-axis, second term as the simple-slope moderator.
      x_info <- components[[1L]]
      moderator_info <- components[[2L]]
      moderator_sets <- ms_interaction_numeric_moderator_presets(moderator_info$values)
      if (!length(moderator_sets)) return(NULL)
      prediction_specs <- lapply(moderator_sets, function(set) {
        list(moderator_set = set, facet_set = NULL)
      })
      interaction_kind <- "continuous_by_continuous"
    } else if (sum(is_numeric) == 0L && sum(is_categorical) == 2L) {
      x_info <- components[[1L]]
      moderator_info <- components[[2L]]
      x_set <- ms_interaction_categorical_set(x_info, max_levels = max_levels)
      moderator_set <- ms_interaction_categorical_set(moderator_info, max_levels = max_levels)
      if (is.null(x_set) || is.null(moderator_set)) return(NULL)
      moderator_sets <- list(moderator_set)
      prediction_specs <- list(list(moderator_set = moderator_set, facet_set = NULL))
      interaction_kind <- "categorical_by_categorical"
    } else {
      return(NULL)
    }
  } else if (length(parts) == 3L) {
    if (sum(is_numeric) == 1L && sum(is_categorical) == 2L) {
      x_info <- components[[which(is_numeric)[[1]]]]
      categorical_idx <- which(is_categorical)
      moderator_info <- components[[categorical_idx[[1]]]]
      facet_info <- components[[categorical_idx[[2]]]]
      moderator_set <- ms_interaction_categorical_set(moderator_info, max_levels = max_levels)
      facet_set <- ms_interaction_categorical_set(facet_info, max_levels = max_levels)
      if (is.null(moderator_set) || is.null(facet_set)) return(NULL)
      moderator_sets <- list(moderator_set)
      prediction_specs <- list(list(moderator_set = moderator_set, facet_set = facet_set))
      interaction_kind <- "continuous_by_categorical_by_categorical"
    } else if (sum(is_numeric) == 2L && sum(is_categorical) == 1L) {
      numeric_idx <- which(is_numeric)
      x_info <- components[[numeric_idx[[1]]]]
      moderator_info <- components[[numeric_idx[[2]]]]
      facet_info <- components[[which(is_categorical)[[1]]]]
      moderator_sets <- ms_interaction_numeric_moderator_presets(moderator_info$values)
      facet_set <- ms_interaction_categorical_set(facet_info, max_levels = max_levels)
      if (!length(moderator_sets) || is.null(facet_set)) return(NULL)
      prediction_specs <- lapply(moderator_sets, function(set) {
        list(moderator_set = set, facet_set = facet_set)
      })
      interaction_kind <- "continuous_by_continuous_by_categorical"
    } else if (sum(is_numeric) == 3L && sum(is_categorical) == 0L) {
      x_info <- components[[1L]]
      moderator_info <- components[[2L]]
      facet_info <- components[[3L]]
      moderator_sets <- ms_interaction_numeric_moderator_presets(moderator_info$values)
      facet_sets <- ms_interaction_numeric_moderator_presets(facet_info$values)
      if (!length(moderator_sets) || !length(facet_sets)) return(NULL)
      prediction_specs <- ms_interaction_numeric_numeric_specs(moderator_sets, facet_sets)
      if (!length(prediction_specs)) return(NULL)
      interaction_kind <- "continuous_by_continuous_by_continuous"
    } else if (sum(is_numeric) == 0L && sum(is_categorical) == 3L) {
      x_info <- components[[1L]]
      moderator_info <- components[[2L]]
      facet_info <- components[[3L]]
      x_set <- ms_interaction_categorical_set(x_info, max_levels = max_levels)
      moderator_set <- ms_interaction_categorical_set(moderator_info, max_levels = max_levels)
      facet_set <- ms_interaction_categorical_set(facet_info, max_levels = max_levels)
      if (is.null(x_set) || is.null(moderator_set) || is.null(facet_set)) return(NULL)
      moderator_sets <- list(moderator_set)
      prediction_specs <- list(list(moderator_set = moderator_set, facet_set = facet_set))
      interaction_kind <- "categorical_by_categorical_by_categorical"
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }

  if (identical(x_info$type, "numeric")) {
    x_values <- ms_interaction_numeric_values(x_info$values)
    if (length(x_values) < 2L) return(NULL)
    x_range <- range(x_values, na.rm = TRUE)
    if (!all(is.finite(x_range)) || x_range[[1]] == x_range[[2]]) return(NULL)
    x_grid <- seq(x_range[[1]], x_range[[2]], length.out = max(12L, as.integer(grid_points)))
  } else {
    if (is.null(x_set)) x_set <- ms_interaction_categorical_set(x_info, max_levels = max_levels)
    if (is.null(x_set) || length(x_set$values) < 2L) return(NULL)
    x_range <- c(1, length(x_set$values))
    x_grid <- x_set$values
  }

  base <- ms_interaction_base_values(
    x, mf, var_map,
    x_info$variable,
    moderator_info$variable,
    facet_variable = if (!is.null(facet_info)) facet_info$variable else NULL,
    max_levels = max_levels
  )
  if (is.null(base)) return(NULL)

  prediction_sets <- lapply(prediction_specs, function(spec) {
    ms_interaction_prediction_set(
      x = x,
      x_grid = x_grid,
      moderator_set = spec$moderator_set,
      base_values = base$values,
      var_map = var_map,
      x_info = x_info,
      x_set = x_set,
      moderator_info = moderator_info,
      facet_info = facet_info,
      facet_set = spec$facet_set,
      conf.level = conf.level
    )
  })
  prediction_sets <- Filter(Negate(is.null), prediction_sets)
  if (!length(prediction_sets)) return(NULL)
  default_set <- prediction_sets[[1L]]

  family <- tryCatch(stats::family(x), error = function(e) NULL)
  fam <- if (!is.null(family)) as.character(family$family %||% "") else ""
  link <- if (!is.null(family)) as.character(family$link %||% "") else ""
  bounded <- tolower(fam) %in% c("binomial", "quasibinomial")
  outcome <- ms_lm_outcome_label(x)
  y_label <- if (bounded) {
    "Predicted probability"
  } else if (!is.null(outcome) && nzchar(outcome)) {
    paste("Predicted", outcome)
  } else {
    "Predicted response"
  }

  out <- list(
    interaction_term = interaction_term,
    interaction_kind = interaction_kind,
    variables = c(x_info$variable, moderator_info$variable,
                  if (!is.null(facet_info)) facet_info$variable else character(0)),
    x = list(
      variable = x_info$variable,
      term = x_info$term,
      label = x_info$label,
      type = x_info$type,
      range = I(ms_safe_numeric(x_range))
    ),
    moderator = list(
      variable = moderator_info$variable,
      term = moderator_info$term,
      label = moderator_info$label,
      type = moderator_info$type,
      levels = default_set$levels
    ),
    grid = default_set$grid,
    held_constant = base$held_constant,
    outcome = outcome %||% NULL,
    y_label = y_label,
    scale = "response",
    ci_level = conf.level,
    bounded_response = bounded
  )
  if (!is.null(x_set)) {
    out$x$levels <- ms_interaction_set_levels(x_set)
  }
  if (!is.null(facet_info)) {
    out$facet <- list(
      variable = facet_info$variable,
      term = facet_info$term,
      label = facet_info$label,
      type = facet_info$type,
      levels = default_set$facet_levels %||% ms_interaction_set_levels(facet_set)
    )
  }
  if (identical(moderator_info$type, "numeric") && length(prediction_sets) > 1L) {
    out$moderator_default_preset <- default_set$id
    out$moderator_value_presets <- lapply(prediction_sets, function(set) {
      list(
        id = set$id,
        label = set$label,
        rule = set$rule,
        levels = set$levels,
        grid = set$grid,
        facet_levels = set$facet_levels %||% NULL
      )
    })
  }
  if (nzchar(fam)) out$model_family <- fam
  if (nzchar(link)) out$model_link <- link
  out
}

ms_glm_response_main_effect_figure_data <- function(x, term_roles = NULL,
                                                    grid_points = 80L,
                                                    max_levels = 8L,
                                                    conf.level = 0.95) {
  if (!inherits(x, "glm")) return(NULL)
  if (inherits(x, c("lmerMod", "glmerMod", "merMod", "gam", "rlm"))) return(NULL)

  family <- tryCatch(stats::family(x), error = function(e) NULL)
  if (is.null(family) || !is.function(family$linkinv)) return(NULL)
  family_name <- tolower(as.character(family$family %||% ""))
  family_link <- tolower(as.character(family$link %||% ""))

  tf <- tryCatch(stats::terms(x), error = function(e) NULL)
  mf <- tryCatch(stats::model.frame(x), error = function(e) NULL)
  if (is.null(tf) || is.null(mf) || nrow(mf) < 2L) return(NULL)

  term_labels <- attr(tf, "term.labels") %||% character(0)
  orders <- attr(tf, "order") %||% integer(0)
  if (length(term_labels) == 0L || length(orders) != length(term_labels)) return(NULL)
  if (any(orders != 1L)) return(NULL)

  offset_info <- ms_glm_offset_info(tf = tf, term_labels = term_labels,
                                    family_link = family_link)
  if (!is.null(offset_info) && !isTRUE(offset_info$supported)) return(NULL)

  family_info <- ms_glm_response_family_info(
    family_name = family_name,
    family_link = family_link,
    offset_info = offset_info
  )
  if (is.null(family_info)) return(NULL)

  var_map <- ms_model_frame_variable_map(mf)
  term_info <- ms_glm_response_term_info(term_labels, mf, var_map, max_levels = max_levels)
  if (!length(term_info)) return(NULL)

  focal_info <- ms_glm_response_selected_term(term_info, term_roles = term_roles)
  if (is.null(focal_info)) return(NULL)

  if (identical(focal_info$type, "categorical")) {
    adjusted_means <- ms_glm_response_factor_means_figure_data(
      x = x,
      mf = mf,
      var_map = var_map,
      term_info = term_info,
      factor_info = focal_info,
      family_info = family_info,
      offset_info = offset_info,
      max_levels = max_levels,
      conf.level = conf.level
    )
    if (is.null(adjusted_means)) return(NULL)
    return(list(adjusted_means = adjusted_means))
  }

  if (identical(focal_info$type, "numeric")) {
    effect_plot <- ms_glm_response_continuous_effect_figure_data(
      x = x,
      var_map = var_map,
      x_info = focal_info,
      family_info = family_info,
      offset_info = offset_info,
      grid_points = grid_points,
      max_levels = max_levels,
      conf.level = conf.level
    )
    if (is.null(effect_plot)) return(NULL)
    return(list(interaction_plot = effect_plot))
  }

  NULL
}

ms_glm_response_family_info <- function(family_name = "", family_link = "",
                                        offset_info = NULL) {
  family_name <- tolower(as.character(family_name %||% ""))
  family_link <- tolower(as.character(family_link %||% ""))
  if (family_name %in% c("binomial", "quasibinomial") &&
      family_link %in% c("logit", "logistic")) {
    return(list(
      mean_kind = "predicted_probability",
      y_label = "Predicted probability",
      bounded_response = TRUE,
      model_family = family_name,
      model_link = family_link
    ))
  }

  if (family_name %in% c("poisson", "quasipoisson") &&
      family_link %in% c("log", "identity", "sqrt")) {
    is_rate <- !is.null(offset_info) &&
      isTRUE(offset_info$rate) &&
      identical(family_link, "log")
    return(list(
      mean_kind = if (is_rate) "predicted_rate" else "predicted_count",
      y_label = if (is_rate) "Predicted rate" else "Predicted count",
      bounded_response = FALSE,
      model_family = family_name,
      model_link = family_link
    ))
  }

  if (identical(family_name, "gamma") && identical(family_link, "log")) {
    return(list(
      mean_kind = "predicted_mean",
      y_label = "Predicted mean",
      bounded_response = FALSE,
      model_family = family_name,
      model_link = family_link
    ))
  }

  NULL
}

ms_glm_offset_info <- function(tf, term_labels = character(0),
                               family_link = "") {
  offset_idx <- attr(tf, "offset") %||% integer(0)
  if (!length(offset_idx)) return(NULL)

  variables_call <- attr(tf, "variables")
  if (is.null(variables_call) || length(variables_call) == 0L) {
    return(list(present = TRUE, supported = FALSE))
  }

  term_variables <- unique(unlist(lapply(term_labels, function(term) {
    tryCatch(
      all.vars(stats::as.formula(paste0("~", term))),
      error = function(e) character(0)
    )
  }), use.names = FALSE))
  term_variables <- term_variables[nzchar(term_variables)]

  items <- lapply(offset_idx, function(idx) {
    expr <- tryCatch(variables_call[[idx + 1L]], error = function(e) NULL)
    ms_glm_offset_item(expr, term_variables = term_variables)
  })
  items <- Filter(Negate(is.null), items)
  if (!length(items)) return(list(present = TRUE, supported = FALSE))

  supported <- all(vapply(items, function(item) isTRUE(item$supported), logical(1)))
  if (!supported) return(list(present = TRUE, supported = FALSE))

  prediction_values <- list()
  held <- list()
  for (item in items) {
    if (!isTRUE(item$offset_only)) next
    if (is.null(item$variable) || !nzchar(item$variable)) next
    prediction_values[[item$variable]] <- item$value
    held[[length(held) + 1L]] <- list(
      variable = item$variable,
      term = item$term,
      label = ms_model_clean_term(item$variable),
      value = as.character(item$value),
      value_label = item$value_label,
      rule = item$rule
    )
  }

  list(
    present = TRUE,
    supported = TRUE,
    rate = identical(tolower(as.character(family_link %||% "")), "log") &&
      length(held) > 0L,
    variables = unique(vapply(items, function(item) item$variable %||% "", character(1))),
    offset_only_variables = unique(vapply(Filter(function(item) isTRUE(item$offset_only), items),
                                          function(item) item$variable %||% "", character(1))),
    prediction_values = prediction_values,
    held_constant = held,
    terms = lapply(items, function(item) {
      list(
        term = item$term,
        variable = item$variable,
        label = ms_model_clean_term(item$variable),
        value = as.character(item$value),
        value_label = item$value_label,
        rule = item$rule,
        offset_only = isTRUE(item$offset_only)
      )
    })
  )
}

ms_glm_offset_item <- function(expr, term_variables = character(0)) {
  if (is.null(expr)) return(NULL)
  term <- paste(deparse(expr, width.cutoff = 500L), collapse = " ")
  arg <- if (is.call(expr) &&
             identical(as.character(expr[[1L]] %||% ""), "offset") &&
             length(expr) >= 2L) {
    expr[[2L]]
  } else {
    expr
  }
  vars <- tryCatch(all.vars(arg), error = function(e) character(0))
  vars <- unique(vars[nzchar(vars)])
  if (length(vars) != 1L) {
    variable <- if (length(vars) > 0L) vars[[1L]] else ""
    return(list(term = term, variable = variable, supported = FALSE))
  }

  variable <- vars[[1L]]
  offset_only <- !(variable %in% term_variables)
  if (is.call(arg) &&
      identical(as.character(arg[[1L]] %||% ""), "log") &&
      length(arg) == 2L &&
      is.symbol(arg[[2L]]) &&
      identical(as.character(arg[[2L]]), variable)) {
    return(list(
      term = term,
      variable = variable,
      value = 1,
      value_label = "1",
      rule = "unit_exposure",
      offset_only = offset_only,
      supported = TRUE
    ))
  }

  if (is.symbol(arg)) {
    return(list(
      term = term,
      variable = variable,
      value = 0,
      value_label = "0",
      rule = "zero_offset",
      offset_only = offset_only,
      supported = TRUE
    ))
  }

  list(term = term, variable = variable, supported = FALSE)
}

ms_glm_response_term_info <- function(term_labels, mf, var_map, max_levels = 8L) {
  out <- lapply(term_labels, function(term) {
    if (!ms_interaction_component_supported(term)) return(NULL)
    variable <- ms_interaction_component_variable(term)
    if (is.null(variable) || is.null(var_map[[variable]])) return(NULL)
    values <- var_map[[variable]]$values
    type <- ms_interaction_variable_type(values, max_levels = max_levels)
    if (is.null(type)) return(NULL)
    info <- list(
      term = term,
      variable = variable,
      label = ms_model_clean_term(term),
      type = type,
      values = values
    )
    if (identical(type, "categorical")) {
      info$levels <- ms_interaction_category_levels(values)
    }
    info
  })
  Filter(Negate(is.null), out)
}

ms_glm_response_selected_term <- function(term_info, term_roles = NULL) {
  if (!length(term_info)) return(NULL)
  focal_terms <- if (!is.null(term_roles) && length(term_roles$focal_terms %||% character(0))) {
    as.character(term_roles$focal_terms)
  } else {
    character(0)
  }
  focal_keys <- ms_role_keys(focal_terms)

  candidates <- term_info
  if (length(focal_keys)) {
    is_focal <- vapply(term_info, function(info) {
      keys <- ms_role_keys(c(info$term %||% "", info$variable %||% "", info$label %||% ""))
      length(intersect(keys, focal_keys)) > 0L
    }, logical(1))
    if (any(is_focal)) candidates <- term_info[is_focal]
  }

  categorical_idx <- which(vapply(candidates, function(info) {
    identical(info$type, "categorical")
  }, logical(1)))
  if (length(categorical_idx)) return(candidates[[categorical_idx[[1L]]]])

  numeric_idx <- which(vapply(candidates, function(info) {
    identical(info$type, "numeric")
  }, logical(1)))
  if (length(numeric_idx)) return(candidates[[numeric_idx[[1L]]]])

  NULL
}

ms_glm_response_continuous_effect_figure_data <- function(x, var_map, x_info,
                                                          family_info = NULL,
                                                          offset_info = NULL,
                                                          grid_points = 80L,
                                                          max_levels = 8L,
                                                          conf.level = 0.95) {
  x_values <- ms_interaction_numeric_values(x_info$values)
  if (length(x_values) < 2L) return(NULL)
  x_range <- range(x_values, na.rm = TRUE)
  if (!all(is.finite(x_range)) || x_range[[1L]] == x_range[[2L]]) return(NULL)
  x_grid <- seq(x_range[[1L]], x_range[[2L]], length.out = max(12L, as.integer(grid_points)))

  base <- ms_glm_response_base_values(
    x = x,
    var_map = var_map,
    focal_variable = x_info$variable,
    offset_info = offset_info,
    max_levels = max_levels
  )
  if (is.null(base)) return(NULL)

  rows <- lapply(x_grid, function(x_value) {
    row <- base$values
    row[[x_info$variable]] <- x_value
    row
  })
  predicted <- ms_glm_response_predict_rows(x, rows, var_map, conf.level = conf.level)
  if (is.null(predicted)) return(NULL)

  grid <- lapply(seq_along(rows), function(i) {
    row <- list(
      x = ms_safe_numeric(x_grid[[i]]),
      moderator_value = "estimate",
      moderator_label = "Estimate",
      estimate = predicted$estimate[[i]],
      ci_lower = predicted$ci_lower[[i]],
      ci_upper = predicted$ci_upper[[i]],
      se = predicted$se[[i]]
    )
    if (!is.null(predicted$link_estimate)) row$link_estimate <- predicted$link_estimate[[i]]
    row
  })
  grid <- Filter(function(row) {
    is.list(row) &&
      !is.na(ms_safe_numeric(row$x)) &&
      !is.na(ms_safe_numeric(row$estimate))
  }, grid)
  if (length(grid) < 2L) return(NULL)

  outcome <- ms_lm_outcome_label(x)
  out <- list(
    interaction_term = x_info$term,
    interaction_kind = "continuous_main_effect",
    source = "glm_predictions",
    mean_kind = family_info$mean_kind %||% "predicted",
    variables = c(x_info$variable),
    x = list(
      variable = x_info$variable,
      term = x_info$term,
      label = x_info$label,
      type = "numeric",
      range = I(ms_safe_numeric(x_range))
    ),
    moderator = list(
      variable = "estimate",
      term = "estimate",
      label = "Estimate",
      type = "categorical",
      levels = list(list(value = "estimate", label = "Estimate"))
    ),
    grid = grid,
    held_constant = base$held_constant,
    outcome = outcome %||% NULL,
    y_label = family_info$y_label %||% "Predicted value",
    scale = "response",
    ci_level = conf.level,
    ci_method = "wald_link",
    bounded_response = isTRUE(family_info$bounded_response),
    model_family = family_info$model_family %||% "",
    model_link = family_info$model_link %||% ""
  )
  if (!is.null(offset_info) && length(offset_info$terms %||% list())) out$offset <- offset_info
  out
}

ms_glm_response_factor_means_figure_data <- function(x, mf, var_map, term_info,
                                                     factor_info,
                                                     family_info = NULL,
                                                     offset_info = NULL,
                                                     max_levels = 8L,
                                                     conf.level = 0.95) {
  factor_set <- ms_interaction_categorical_set(factor_info, max_levels = max_levels)
  if (is.null(factor_set)) return(NULL)

  base <- ms_glm_response_base_values(
    x = x,
    var_map = var_map,
    focal_variable = factor_info$variable,
    offset_info = offset_info,
    max_levels = max_levels
  )
  if (is.null(base)) return(NULL)

  rows <- lapply(factor_set$values, function(level) {
    row <- base$values
    row[[factor_info$variable]] <- level
    row
  })
  predicted <- ms_glm_response_predict_rows(x, rows, var_map, conf.level = conf.level)
  if (is.null(predicted)) return(NULL)

  counts <- table(as.character(factor_info$values[!is.na(factor_info$values)]))
  groups <- lapply(seq_along(factor_set$values), function(i) {
    level <- factor_set$values[[i]]
    n_value <- if (as.character(level) %in% names(counts)) {
      as.integer(counts[[as.character(level)]])
    } else {
      NA_integer_
    }
    row <- list(
      level = as.character(level),
      label = factor_set$labels[[i]],
      n = n_value,
      mean = predicted$estimate[[i]],
      ci_lower = predicted$ci_lower[[i]],
      ci_upper = predicted$ci_upper[[i]],
      se = predicted$se[[i]]
    )
    if (!is.null(predicted$link_estimate)) row$link_estimate <- predicted$link_estimate[[i]]
    row
  })
  groups <- Filter(function(row) {
    is.list(row) &&
      !is.null(row$level) &&
      !is.na(ms_safe_numeric(row$mean))
  }, groups)
  if (length(groups) < 2L) return(NULL)

  held_numeric <- Filter(function(row) identical(row$type %||% "", "mean"), base$held_constant)
  held_reference <- Filter(function(row) !identical(row$type %||% "", "mean"), base$held_constant)
  covariates <- lapply(held_numeric, function(row) {
    value <- suppressWarnings(as.numeric(row$value))
    if (is.na(value) || !is.finite(value)) return(NULL)
    list(
      variable = row$variable,
      term = row$variable,
      label = row$label,
      value = value,
      value_label = row$value_label,
      rule = "sample_mean"
    )
  })
  covariates <- Filter(Negate(is.null), covariates)

  outcome <- ms_lm_outcome_label(x)
  out <- list(
    mean_kind = family_info$mean_kind %||% "predicted",
    source = "glm_predictions",
    factor = list(
      variable = factor_info$variable,
      term = factor_info$term,
      label = factor_info$label,
      levels = ms_interaction_set_levels(factor_set)
    ),
    groups = groups,
    marginalized_terms = list(),
    outcome = outcome %||% NULL,
    y_label = family_info$y_label %||% "Predicted value",
    ci_level = conf.level,
    ci_method = "wald_link",
    bounded_response = isTRUE(family_info$bounded_response),
    model_family = family_info$model_family %||% "",
    model_link = family_info$model_link %||% ""
  )
  if (length(covariates) > 0L) {
    out$covariates <- covariates
    out$adjustment <- list(rule = "sample_mean", label = "sample means")
  }
  if (length(held_reference) > 0L) out$held_constant <- held_reference
  if (!is.null(offset_info) && length(offset_info$terms %||% list())) out$offset <- offset_info
  out
}

ms_glm_response_base_values <- function(x, var_map, focal_variable,
                                        offset_info = NULL,
                                        max_levels = 8L) {
  tf <- tryCatch(stats::delete.response(stats::terms(x)), error = function(e) NULL)
  predictor_vars <- if (!is.null(tf)) all.vars(tf) else character(0)
  predictor_vars <- unique(predictor_vars[nzchar(predictor_vars)])
  offset_only <- as.character(offset_info$offset_only_variables %||% character(0))
  offset_only <- offset_only[nzchar(offset_only)]
  if (length(offset_only) > 0L) {
    predictor_vars <- predictor_vars[!(predictor_vars %in% offset_only)]
  }
  if (length(predictor_vars) == 0L) return(NULL)

  values <- list()
  held_constant <- list()
  for (variable in predictor_vars) {
    if (identical(variable, focal_variable)) next
    entry <- var_map[[variable]]
    if (is.null(entry)) return(NULL)
    hold <- ms_interaction_hold_value(entry$values, max_levels = max_levels)
    if (is.null(hold)) return(NULL)
    values[[variable]] <- hold$value
    held_constant[[length(held_constant) + 1L]] <- list(
      variable = variable,
      label = ms_model_clean_term(variable),
      value = as.character(hold$value),
      value_label = hold$label,
      type = hold$type
    )
  }
  prediction_values <- offset_info$prediction_values %||% list()
  if (length(prediction_values) > 0L) {
    for (variable in names(prediction_values)) {
      if (!nzchar(variable)) next
      values[[variable]] <- prediction_values[[variable]]
    }
  }
  list(values = values, held_constant = held_constant)
}

ms_glm_response_predict_rows <- function(x, rows, var_map, conf.level = 0.95) {
  newdata <- ms_glm_response_newdata(rows, var_map)
  if (is.null(newdata) || nrow(newdata) == 0L) return(NULL)
  predicted <- ms_interaction_predict_grid(x, newdata, conf.level = conf.level)
  if (is.null(predicted) || length(predicted$estimate) != nrow(newdata)) return(NULL)
  predicted
}

ms_glm_response_newdata <- function(rows, var_map) {
  out <- tryCatch(
    do.call(rbind, lapply(rows, function(row) {
      as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
    })),
    error = function(e) NULL
  )
  if (is.null(out)) return(NULL)
  rownames(out) <- NULL
  for (variable in names(out)) {
    values <- var_map[[variable]]$values
    if (is.factor(values)) {
      out[[variable]] <- factor(as.character(out[[variable]]), levels = levels(values))
    } else if (is.logical(values)) {
      out[[variable]] <- as.logical(out[[variable]])
    } else if (is.numeric(values) || is.integer(values)) {
      out[[variable]] <- suppressWarnings(as.numeric(out[[variable]]))
    } else {
      out[[variable]] <- as.character(out[[variable]])
    }
  }
  out
}

ms_interaction_component_supported <- function(part) {
  part <- trimws(as.character(part %||% ""))
  if (!nzchar(part)) return(FALSE)
  if (!grepl("\\(", part)) return(TRUE)
  grepl("^\\s*(factor|as\\.factor)\\s*\\(", part)
}

ms_interaction_component_variable <- function(part) {
  vars <- tryCatch(
    all.vars(stats::as.formula(paste0("~", part))),
    error = function(e) character(0)
  )
  if (length(vars) != 1L || !nzchar(vars[[1]])) return(NULL)
  vars[[1]]
}

ms_interaction_terms_within_components <- function(terms, component_variables) {
  allowed <- unique(as.character(component_variables %||% character(0)))
  allowed <- allowed[nzchar(allowed)]
  if (!length(allowed)) return(FALSE)
  all(vapply(terms, function(term) {
    parts <- strsplit(term, ":", fixed = TRUE)[[1]]
    vars <- vapply(parts, function(part) {
      ms_interaction_component_variable(part) %||% ""
    }, character(1))
    length(vars) >= 2L && all(nzchar(vars)) && all(vars %in% allowed)
  }, logical(1)))
}

ms_interaction_component_info <- function(part, var_map, max_levels = 6L) {
  variable <- ms_interaction_component_variable(part)
  if (is.null(variable) || is.null(var_map[[variable]])) return(NULL)
  values <- var_map[[variable]]$values
  type <- ms_interaction_variable_type(values, max_levels = max_levels)
  if (is.null(type)) return(NULL)
  levels <- if (identical(type, "categorical")) {
    ms_interaction_category_levels(values)
  } else {
    NULL
  }
  list(
    term = part,
    variable = variable,
    label = ms_model_clean_term(part),
    values = values,
    type = type,
    levels = levels
  )
}

ms_model_frame_variable_map <- function(mf) {
  out <- list()
  if (is.null(mf) || ncol(mf) < 2L) return(out)
  cols <- names(mf)[-1L]
  for (col in cols) {
    vars <- tryCatch(
      all.vars(stats::as.formula(paste0("~", col))),
      error = function(e) character(0)
    )
    if (length(vars) == 1L && nzchar(vars[[1]]) && is.null(out[[vars[[1]]]])) {
      out[[vars[[1]]]] <- list(variable = vars[[1]], column = col, values = mf[[col]])
    }
    if (nzchar(col) && is.null(out[[col]])) {
      out[[col]] <- list(variable = col, column = col, values = mf[[col]])
    }
  }
  out
}

ms_interaction_variable_type <- function(values, max_levels = 6L) {
  if (is.factor(values) || is.character(values) || is.logical(values)) {
    return("categorical")
  }
  if (is.numeric(values) || is.integer(values)) {
    unique_values <- unique(ms_interaction_numeric_values(values))
    if (length(unique_values) >= 2L && length(unique_values) <= max_levels) {
      return("categorical")
    }
    return("numeric")
  }
  NULL
}

ms_interaction_numeric_values <- function(values) {
  values <- suppressWarnings(as.numeric(values))
  values[is.finite(values)]
}

ms_interaction_category_levels <- function(values) {
  if (is.factor(values)) {
    observed <- unique(as.character(values[!is.na(values)]))
    lev <- levels(values)
    lev[lev %in% observed]
  } else if (is.logical(values)) {
    lev <- unique(values[!is.na(values)])
    lev <- sort(lev)
    as.logical(lev)
  } else if (is.numeric(values) || is.integer(values)) {
    sort(unique(ms_interaction_numeric_values(values)))
  } else {
    sort(unique(as.character(values[!is.na(values)])))
  }
}

ms_interaction_categorical_set <- function(info, max_levels = 6L) {
  if (is.null(info) || !identical(info$type, "categorical")) return(NULL)
  levels <- info$levels
  if (length(levels) < 2L || length(levels) > max_levels) return(NULL)
  list(
    id = "observed",
    label = "Observed levels",
    rule = "observed_levels",
    values = levels,
    labels = ms_interaction_level_labels(levels),
    value_labels = as.character(levels),
    slices = rep(NA_character_, length(levels))
  )
}

ms_interaction_numeric_moderator_presets <- function(values) {
  numeric_values <- ms_interaction_numeric_values(values)
  if (length(numeric_values) < 2L) return(NULL)
  center <- mean(numeric_values, na.rm = TRUE)
  spread <- stats::sd(numeric_values, na.rm = TRUE)
  presets <- list()
  if (is.finite(center) && is.finite(spread) && spread > 0) {
    presets[[length(presets) + 1L]] <- ms_interaction_numeric_moderator_set(
      id = "sd",
      label = "M +/- 1 SD",
      rule = "mean_plus_minus_1_sd",
      values = c(center - spread, center, center + spread),
      slices = c("low", "mean", "high"),
      label_terms = c("Low", "Mean", "High")
    )
  }

  quartiles <- ms_interaction_numeric_quantiles(numeric_values, c(0.25, 0.5, 0.75))
  presets[[length(presets) + 1L]] <- ms_interaction_numeric_moderator_set(
    id = "quartiles",
    label = "Quartiles",
    rule = "quartiles",
    values = quartiles,
    slices = c("q1", "median", "q3"),
    label_terms = c("Q1", "Median", "Q3")
  )

  percentiles <- ms_interaction_numeric_quantiles(numeric_values, c(0.1, 0.5, 0.9))
  presets[[length(presets) + 1L]] <- ms_interaction_numeric_moderator_set(
    id = "percentiles",
    label = "10/50/90 percentiles",
    rule = "percentiles_10_50_90",
    values = percentiles,
    slices = c("p10", "p50", "p90"),
    label_terms = c("P10", "P50", "P90")
  )

  presets <- Filter(Negate(is.null), presets)
  if (!length(presets)) return(NULL)
  presets
}

ms_interaction_numeric_quantiles <- function(values, probs) {
  out <- tryCatch(
    stats::quantile(values, probs = probs, na.rm = TRUE, names = FALSE, type = 7),
    error = function(e) rep(NA_real_, length(probs))
  )
  ms_safe_numeric(out)
}

ms_interaction_numeric_moderator_set <- function(id, label, rule, values, slices,
                                                 label_terms) {
  values <- ms_safe_numeric(values)
  if (length(values) != 3L || any(is.na(values)) || any(!is.finite(values))) {
    return(NULL)
  }
  if (length(unique(round(values, 10))) < 3L) return(NULL)
  value_labels <- vapply(values, ms_interaction_format_number, character(1))
  list(
    id = id,
    label = label,
    rule = rule,
    values = values,
    labels = paste0(label_terms, " (", value_labels, ")"),
    value_labels = value_labels,
    slices = slices
  )
}

ms_interaction_numeric_numeric_specs <- function(moderator_sets, facet_sets) {
  if (!length(moderator_sets) || !length(facet_sets)) return(list())
  out <- list()
  for (moderator_set in moderator_sets) {
    id <- moderator_set$id %||% ""
    if (!nzchar(id)) next
    facet_idx <- match(id, vapply(facet_sets, function(set) set$id %||% "", character(1)),
                       nomatch = 0L)
    if (facet_idx <= 0L) next
    out[[length(out) + 1L]] <- list(
      moderator_set = moderator_set,
      facet_set = facet_sets[[facet_idx]]
    )
  }
  out
}

ms_interaction_format_number <- function(value) {
  value <- ms_safe_numeric(value)
  if (is.na(value) || !is.finite(value)) return("")
  format(round(value, 3), trim = TRUE, scientific = FALSE)
}

ms_interaction_prediction_set <- function(x, x_grid, moderator_set, base_values,
                                          var_map, x_info, moderator_info, x_set = NULL,
                                          facet_info = NULL, facet_set = NULL,
                                          conf.level = 0.95) {
  if (is.null(moderator_set) || !length(moderator_set$values)) return(NULL)
  prediction_grid <- ms_interaction_prediction_grid(
    x_grid = x_grid,
    moderator_levels = moderator_set$values,
    base_values = base_values,
    var_map = var_map,
    x_variable = x_info$variable,
    moderator_variable = moderator_info$variable,
    facet_variable = if (!is.null(facet_info)) facet_info$variable else NULL,
    facet_levels = if (!is.null(facet_set)) facet_set$values else NULL
  )
  if (is.null(prediction_grid) || nrow(prediction_grid) == 0L) return(NULL)

  predicted <- ms_interaction_predict_grid(x, prediction_grid, conf.level = conf.level)
  if (is.null(predicted)) return(NULL)

  n_x <- length(x_grid)
  n_moderator <- length(moderator_set$values)
  has_facet <- !is.null(facet_info) && !is.null(facet_set) && length(facet_set$values)
  has_categorical_x <- identical(x_info$type, "categorical") && !is.null(x_set)
  grid <- lapply(seq_len(nrow(prediction_grid)), function(i) {
    x_idx <- ((i - 1L) %% n_x) + 1L
    level_idx <- (((i - 1L) %/% n_x) %% n_moderator) + 1L
    x_value <- prediction_grid[[x_info$variable]][[i]]
    row <- list(
      x = if (has_categorical_x) x_idx else ms_safe_numeric(x_value),
      moderator_value = as.character(moderator_set$values[[level_idx]]),
      moderator_label = moderator_set$labels[[level_idx]],
      estimate = predicted$estimate[[i]],
      ci_lower = predicted$ci_lower[[i]],
      ci_upper = predicted$ci_upper[[i]],
      se = predicted$se[[i]]
    )
    if (has_categorical_x) {
      row$x_value <- as.character(x_set$values[[x_idx]])
      row$x_label <- x_set$labels[[x_idx]]
    }
    if (has_facet) {
      facet_idx <- ((i - 1L) %/% (n_x * n_moderator)) + 1L
      row$facet_value <- as.character(facet_set$values[[facet_idx]])
      row$facet_label <- facet_set$labels[[facet_idx]]
    }
    if (!is.null(predicted$link_estimate)) row$link_estimate <- predicted$link_estimate[[i]]
    row
  })

  grid <- Filter(function(row) {
    is.list(row) &&
      !is.na(ms_safe_numeric(row$x)) &&
      !is.na(ms_safe_numeric(row$estimate))
  }, grid)
  if (length(grid) == 0L) return(NULL)

  list(
    id = moderator_set$id,
    label = moderator_set$label,
    rule = moderator_set$rule,
    levels = ms_interaction_set_levels(moderator_set),
    facet_levels = if (has_facet) ms_interaction_set_levels(facet_set) else NULL,
    grid = grid
  )
}

ms_interaction_set_levels <- function(set) {
  if (is.null(set) || !length(set$values)) return(list())
  lapply(seq_along(set$values), function(i) {
    row <- list(
      value = as.character(set$values[[i]]),
      label = as.character(set$labels[[i]] %||% set$values[[i]])
    )
    value_label <- as.character(set$value_labels[[i]] %||% "")
    slice <- as.character(set$slices[[i]] %||% "")
    rule <- as.character(set$rule %||% "")
    if (!is.na(value_label) && nzchar(value_label)) row$value_label <- value_label
    if (!is.na(slice) && nzchar(slice)) row$slice <- slice
    if (!is.na(rule) && nzchar(rule)) row$rule <- rule
    row
  })
}

ms_interaction_base_values <- function(x, mf, var_map, x_variable, moderator_variable,
                                       facet_variable = NULL, max_levels = 6L) {
  tf <- tryCatch(stats::delete.response(stats::terms(x)), error = function(e) NULL)
  predictor_vars <- if (!is.null(tf)) all.vars(tf) else character(0)
  predictor_vars <- predictor_vars[nzchar(predictor_vars)]
  predictor_vars <- unique(predictor_vars)
  if (length(predictor_vars) == 0L) return(NULL)

  values <- list()
  held_constant <- list()
  for (variable in predictor_vars) {
    entry <- var_map[[variable]]
    if (is.null(entry)) return(NULL)
    hold <- ms_interaction_hold_value(entry$values, max_levels = max_levels)
    if (is.null(hold)) return(NULL)
    values[[variable]] <- hold$value
    if (!(variable %in% c(x_variable, moderator_variable, facet_variable))) {
      held_constant[[length(held_constant) + 1L]] <- list(
        variable = variable,
        label = ms_model_clean_term(variable),
        value = as.character(hold$value),
        value_label = hold$label,
        type = hold$type
      )
    }
  }
  list(values = values, held_constant = held_constant)
}

ms_interaction_hold_value <- function(values, max_levels = 6L) {
  type <- ms_interaction_variable_type(values, max_levels = max_levels)
  if (is.null(type)) return(NULL)
  if (identical(type, "numeric")) {
    numeric_values <- ms_interaction_numeric_values(values)
    if (!length(numeric_values)) return(NULL)
    value <- mean(numeric_values, na.rm = TRUE)
    return(list(
      value = value,
      label = ms_interaction_format_number(value),
      type = "mean"
    ))
  }
  levels <- ms_interaction_category_levels(values)
  if (!length(levels)) return(NULL)
  value <- levels[[1]]
  list(value = value, label = as.character(value), type = "reference")
}

ms_interaction_prediction_grid <- function(x_grid, moderator_levels, base_values, var_map,
                                           x_variable, moderator_variable,
                                           facet_variable = NULL, facet_levels = NULL) {
  has_facet <- !is.null(facet_variable) && length(facet_levels)
  facet_values <- if (has_facet) facet_levels else list(NULL)
  rows <- vector("list", length(x_grid) * length(moderator_levels) * length(facet_values))
  idx <- 0L
  for (facet_level in facet_values) {
    for (level in moderator_levels) {
      for (x_value in x_grid) {
        idx <- idx + 1L
        row <- base_values
        row[[x_variable]] <- x_value
        row[[moderator_variable]] <- level
        if (has_facet) row[[facet_variable]] <- facet_level
        rows[[idx]] <- row
      }
    }
  }
  out <- tryCatch(
    do.call(rbind, lapply(rows, function(row) {
      as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
    })),
    error = function(e) NULL
  )
  if (is.null(out)) return(NULL)
  rownames(out) <- NULL
  for (variable in names(out)) {
    values <- var_map[[variable]]$values
    if (is.factor(values)) {
      out[[variable]] <- factor(as.character(out[[variable]]), levels = levels(values))
    } else if (is.logical(values)) {
      out[[variable]] <- as.logical(out[[variable]])
    } else if (is.numeric(values) || is.integer(values)) {
      out[[variable]] <- suppressWarnings(as.numeric(out[[variable]]))
    } else {
      out[[variable]] <- as.character(out[[variable]])
    }
  }
  out
}

ms_interaction_predict_grid <- function(x, newdata, conf.level = 0.95) {
  alpha <- (1 - conf.level) / 2
  if (inherits(x, "glm")) {
    pred <- tryCatch(
      stats::predict(x, newdata = newdata, se.fit = TRUE, type = "link"),
      error = function(e) NULL
    )
    if (is.null(pred) || is.null(pred$fit) || is.null(pred$se.fit)) return(NULL)
    family <- tryCatch(stats::family(x), error = function(e) NULL)
    if (is.null(family) || !is.function(family$linkinv)) return(NULL)
    critical <- stats::qnorm(1 - alpha)
    fit <- ms_safe_numeric(pred$fit)
    se <- ms_safe_numeric(pred$se.fit)
    return(list(
      estimate = ms_safe_numeric(family$linkinv(fit)),
      ci_lower = ms_safe_numeric(family$linkinv(fit - critical * se)),
      ci_upper = ms_safe_numeric(family$linkinv(fit + critical * se)),
      se = se,
      link_estimate = fit
    ))
  }

  pred <- tryCatch(
    stats::predict(x, newdata = newdata, se.fit = TRUE),
    error = function(e) NULL
  )
  if (is.null(pred) || is.null(pred$fit) || is.null(pred$se.fit)) return(NULL)
  df <- tryCatch(stats::df.residual(x), error = function(e) NA_real_)
  critical <- if (!is.na(df) && is.finite(df) && df > 0) {
    stats::qt(1 - alpha, df = df)
  } else {
    stats::qnorm(1 - alpha)
  }
  fit <- ms_safe_numeric(pred$fit)
  se <- ms_safe_numeric(pred$se.fit)
  list(
    estimate = fit,
    ci_lower = ms_safe_numeric(fit - critical * se),
    ci_upper = ms_safe_numeric(fit + critical * se),
    se = se
  )
}

ms_interaction_level_labels <- function(levels) {
  vapply(levels, function(level) {
    text <- as.character(level)
    if (identical(text, "TRUE")) return("TRUE")
    if (identical(text, "FALSE")) return("FALSE")
    text
  }, character(1))
}
