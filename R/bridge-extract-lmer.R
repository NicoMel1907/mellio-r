# R bridge — lme4 lmerMod extractor.
#
# Linear mixed-effects models don't yield a single F/p the way lm does
# (the canonical APA report is per-fixed-effect, not per-model). Emit the
# model-fit headline plus fixed-effect estimates and random-effect variance
# components so the web report can be useful without inventing p-values.
#
# Card shape:
#   type:        "lmer_model_summary"
#   type_label:  "Linear mixed model" + REML/ML hint
#   statistic:   { name: "AIC", value: <num>, df: NA }
#   p_value:     NA_real_   (no global model p)
#   estimate:    { name: "logLik", value: <num> }
#   n:           nobs(x)
#   coefficients: list(term, estimate, std_error, statistic)
#   random_effects: list(group, name, variance, std_dev, corr)
#   raw_output:  print(summary(x))
#
# Requires lme4 (Suggests) — fails with a clear message via
# rlang::check_installed if absent.

#' @rdname mellio_payload
#' @export
mellio_payload.lmerMod <- function(x, ..., .call = NULL) {
  rlang::check_installed("lme4", reason = "to extract lmer model summaries")

  call_str <- if (!is.null(.call)) {
    .call
  } else if (!is.null(x@call)) {
    paste(deparse(x@call, width.cutoff = 500L), collapse = " ")
  } else {
    user_call <- match.call()$x
    if (!is.null(user_call) && !identical(user_call, as.name("x"))) {
      paste(deparse(user_call, width.cutoff = 500L), collapse = " ")
    } else NA_character_
  }

  is_reml <- isTRUE(tryCatch(lme4::isREML(x), error = function(e) FALSE))
  type_label <- paste0(
    "Linear mixed model (",
    if (is_reml) "REML" else "ML",
    ")"
  )

  satterthwaite <- ms_lmer_satterthwaite_info(x)
  s <- if (!is.null(satterthwaite$summary)) satterthwaite$summary else summary(x)
  aic_val    <- ms_safe_numeric(stats::AIC(x))
  bic_val    <- ms_safe_numeric(stats::BIC(x))
  loglik_val <- ms_safe_numeric(as.numeric(stats::logLik(x)))
  n_obs      <- ms_safe_numeric(stats::nobs(x))

  fields <- list(
    statistic = list(name = "AIC", value = aic_val),
    p_value   = NA_real_,
    estimate  = list(name = "logLik", value = loglik_val),
    bic       = bic_val,
    logLik    = loglik_val,
    n         = n_obs,
    model_fit = if (is_reml) "REML" else "ML",
    conf_level = 0.95,
    coefficient_ci_method = "wald"
  )
  model_warnings <- ms_model_diagnostics(x)
  if (length(model_warnings) > 0L) fields$model_warnings <- model_warnings
  if (!is.null(satterthwaite$summary)) {
    fields$fixed_effect_df_method <- "satterthwaite"
    fields$coefficient_ci_method <- "satterthwaite_t"
    fields$coefficient_p_value_method <- "satterthwaite_t"
  }
  r2_fields <- ms_lmer_r2_fields(x)
  if (!is.null(r2_fields)) {
    fields$r2_marginal <- r2_fields$r2_marginal
    fields$r2_conditional <- r2_fields$r2_conditional
    fields$r2_method <- r2_fields$r2_method
  }

  term_roles <- ms_lm_term_roles(x)
  if (is.null(term_roles$model_kind)) {
    tf <- tryCatch(stats::terms(x), error = function(e) NULL)
    mf <- tryCatch(stats::model.frame(x), error = function(e) NULL)
    term_labels <- attr(tf, "term.labels") %||% character(0)
    terms <- lapply(term_labels, function(term) {
      list(
        name = term,
        label = ms_model_clean_term(term),
        role = "term",
        type = ms_model_term_type(term, mf)
      )
    })
    inferred <- ms_assign_term_roles(
      terms,
      infer_covariate_roles = TRUE
    )
    if (!is.null(inferred$model_kind)) {
      term_roles$terms <- inferred$terms
      term_roles$focal_terms <- vapply(
        Filter(function(t) identical(t$role, "focal"), inferred$terms),
        function(t) t$name,
        character(1)
      )
      term_roles$control_terms <- vapply(
        Filter(function(t) identical(t$role, "control"), inferred$terms),
        function(t) t$name,
        character(1)
      )
      term_roles$predictor <- ms_model_term_phrase(term_roles$focal_terms)
      term_roles$model_kind <- inferred$model_kind
    }
  }
  if (!is.null(term_roles$outcome)) fields$outcome <- term_roles$outcome
  if (length(term_roles$terms) > 0L) fields$terms <- term_roles$terms
  if (length(term_roles$focal_terms) > 0L) fields$focal_terms <- term_roles$focal_terms
  if (length(term_roles$control_terms) > 0L) fields$control_terms <- term_roles$control_terms
  if (length(term_roles$predictor) > 0L) fields$predictor <- term_roles$predictor
  if (identical(term_roles$model_kind, "controlled_regression")) {
    fields$model_kind <- "controlled_mixed_model"
  }

  model_term_tests <- if (length(satterthwaite$model_term_tests %||% list())) {
    satterthwaite$model_term_tests
  } else {
    ms_lmer_model_term_tests(x)
  }
  if (length(model_term_tests) > 0L) fields$model_term_tests <- model_term_tests
  if (length(satterthwaite$model_term_tests %||% list())) {
    fields$model_term_test_method <- "satterthwaite_f"
    fields$model_term_ss_type <- "type_iii"
  }

  coefs_mat <- tryCatch(stats::coef(s), error = function(e) NULL)
  if (!is.null(coefs_mat) && is.matrix(coefs_mat) && nrow(coefs_mat) > 0L) {
    fields$coefficients <- ms_lmer_coefficients(x, coefs_mat)
    if (!length(fields$focal_terms %||% character(0))) {
      fields$focal_terms <- ms_lmer_focal_terms(fields$coefficients)
    }
  }

  groups <- ms_lmer_groups(x)
  if (length(groups) > 0L) fields$groups <- groups

  random_terms <- ms_lmer_random_terms(x)
  if (length(random_terms) > 0L) fields$random_terms <- random_terms

  random_effects <- ms_lmer_random_effects(x)
  if (length(random_effects) > 0L) fields$random_effects <- random_effects

  role_prefers_categorical_means <- ms_mixed_role_prefers_categorical_means(term_roles)
  focal_factor_terms <- ms_mixed_focal_factor_terms(term_roles)

  interaction_figure <- ms_lmer_emmeans_interaction_plot_data(x)
  means_figure <- if (is.null(interaction_figure)) {
    if (role_prefers_categorical_means) {
      ms_lmer_emmeans_main_effect_means_figure_data(
        x,
        preferred_terms = focal_factor_terms,
        require_preferred = TRUE
      )
    } else {
      ms_lmer_emmeans_main_effect_means_figure_data(x, require_within = TRUE)
    }
  } else NULL
  effect_figure <- if (is.null(interaction_figure) &&
                       (is.null(means_figure) || role_prefers_categorical_means)) {
    ms_lmer_emmeans_continuous_main_effect_plot_data(x)
  } else NULL
  if (is.null(interaction_figure) && is.null(means_figure) && is.null(effect_figure)) {
    means_figure <- ms_lmer_emmeans_main_effect_means_figure_data(x)
  }
  figure_data <- list()
  if (!is.null(interaction_figure)) figure_data$interaction_plot <- interaction_figure
  if (!is.null(effect_figure)) figure_data$interaction_plot <- effect_figure
  if (!is.null(means_figure)) figure_data$adjusted_means <- means_figure
  if (!length(figure_data)) figure_data <- NULL
  package_extras <- c(
    "lme4",
    if (!is.null(satterthwaite$summary) ||
        length(satterthwaite$model_term_tests %||% list())) "lmerTest" else character(0),
    if (!is.null(r2_fields)) "performance" else character(0),
    if (!is.null(interaction_figure) || !is.null(means_figure) || !is.null(effect_figure)) "emmeans" else character(0)
  )

  prov <- ms_provenance_basic()
  data_prov <- tryCatch(ms_data_provenance(x@frame), error = function(e) NULL)
  prov <- ms_provenance_add_data(prov, data_prov)

  ms_build_envelope(
    type       = "lmer_model_summary",
    type_label = type_label,
    call       = trimws(gsub("\\s+", " ", call_str)),
    fields     = fields,
    raw_output = ms_capture_output(s),
    packages   = ms_packages_basic(extras = package_extras),
    figure_data = figure_data,
    provenance = prov
  )
}

ms_model_diagnostics <- function(x) {
  out <- list()

  add <- function(type, message, severity = "warning", detail = NULL) {
    type <- as.character(type %||% "")
    message <- ms_model_diagnostic_text(message)
    severity <- as.character(severity %||% "warning")
    if (!nzchar(type) || !nzchar(message)) return(invisible(NULL))
    duplicate <- any(vapply(out, function(row) {
      identical(row$type %||% "", type) && identical(row$message %||% "", message)
    }, logical(1)))
    if (duplicate) return(invisible(NULL))
    row <- list(type = type, severity = severity, message = message)
    detail <- ms_model_diagnostic_text(detail)
    if (nzchar(detail)) row$detail <- detail
    out[[length(out) + 1L]] <<- row
    invisible(NULL)
  }

  if (inherits(x, c("merMod", "lmerMod", "glmerMod"))) {
    singular <- isTRUE(tryCatch(lme4::isSingular(x, tol = 1e-4),
                                error = function(e) FALSE))
    opt <- tryCatch(x@optinfo, error = function(e) NULL)
    lme4_messages <- ms_model_diagnostic_vector(
      tryCatch(opt$conv$lme4$messages, error = function(e) NULL)
    )
    opt_warnings <- ms_model_diagnostic_vector(
      tryCatch(opt$warnings, error = function(e) NULL)
    )

    if (singular || any(grepl("singular|boundary", lme4_messages, ignore.case = TRUE))) {
      add(
        "singular_fit",
        "The model returned a singular-fit warning; random-effect variance estimates should be interpreted cautiously."
      )
    }

    convergence_messages <- unique(c(lme4_messages, opt_warnings))
    convergence_messages <- convergence_messages[
      nzchar(convergence_messages) &
        !grepl("singular|boundary", convergence_messages, ignore.case = TRUE)
    ]
    opt_code <- tryCatch(opt$conv$opt, error = function(e) NULL)
    has_nonzero_opt_code <- length(opt_code) > 0L &&
      any(!is.na(suppressWarnings(as.numeric(opt_code))) &
            suppressWarnings(as.numeric(opt_code)) != 0)
    if (length(convergence_messages) > 0L) {
      add(
        "convergence_warning",
        paste0("The model reported convergence warnings: ",
               paste(convergence_messages, collapse = "; "), ".")
      )
    } else if (has_nonzero_opt_code) {
      add(
        "optimizer_convergence",
        paste0("The optimizer reported a non-zero convergence code (",
               paste(opt_code, collapse = ", "),
               "); estimates should be checked before publication.")
      )
    }
  }

  if (inherits(x, "glm")) {
    converged <- tryCatch(x$converged, error = function(e) NULL)
    if (identical(converged, FALSE)) {
      add(
        "convergence_warning",
        "The model did not report successful convergence; estimates should be checked before publication."
      )
    }
    if (isTRUE(tryCatch(x$boundary, error = function(e) FALSE))) {
      add(
        "boundary_fit",
        "The model fit reached a boundary solution; coefficient estimates may be unstable."
      )
    }
    family <- tryCatch(stats::family(x), error = function(e) NULL)
    family_name <- tolower(as.character(family$family %||% ""))
    if (family_name %in% c("binomial", "quasibinomial")) {
      fitted <- tryCatch(stats::fitted(x), error = function(e) numeric(0))
      fitted <- suppressWarnings(as.numeric(fitted))
      fitted <- fitted[is.finite(fitted)]
      if (length(fitted) > 0L && any(fitted <= 1e-8 | fitted >= 1 - 1e-8)) {
        add(
          "separation_or_boundary",
          "The model produced fitted probabilities very close to 0 or 1; separation or sparse outcome patterns may affect estimates."
        )
      }
    }
  }

  out
}

ms_model_diagnostic_vector <- function(value) {
  if (is.null(value)) return(character(0))
  value <- unlist(value, recursive = TRUE, use.names = FALSE)
  value <- as.character(value %||% character(0))
  value <- trimws(value)
  value[nzchar(value)]
}

ms_model_diagnostic_text <- function(value, max_chars = 260L) {
  if (is.null(value) || length(value) == 0L) return("")
  value <- paste(as.character(value), collapse = " ")
  value <- trimws(gsub("\\s+", " ", value))
  if (!nzchar(value)) return("")
  if (nchar(value) > max_chars) {
    value <- paste0(substr(value, 1L, max_chars - 1L), "\u2026")
  }
  value
}

ms_lmer_r2_fields <- function(x) {
  if (!requireNamespace("performance", quietly = TRUE)) return(NULL)
  singular <- isTRUE(tryCatch(lme4::isSingular(x, tol = 1e-4),
                              error = function(e) FALSE))
  if (singular) return(NULL)
  r2_fun <- tryCatch(getExportedValue("performance", "r2"), error = function(e) NULL)
  if (!is.function(r2_fun)) return(NULL)
  r2 <- tryCatch(suppressMessages(suppressWarnings(r2_fun(x))), error = function(e) NULL)
  if (is.null(r2)) return(NULL)
  r2_df <- tryCatch(as.data.frame(r2), error = function(e) NULL)
  if (is.null(r2_df) || !nrow(r2_df)) return(NULL)

  marginal_col <- grep("^R2_?marginal$", names(r2_df), ignore.case = TRUE, value = TRUE)
  conditional_col <- grep("^R2_?conditional$", names(r2_df), ignore.case = TRUE, value = TRUE)
  marginal <- if (length(marginal_col)) ms_safe_numeric(r2_df[[marginal_col[[1L]]]][[1L]]) else NA_real_
  conditional <- if (length(conditional_col)) ms_safe_numeric(r2_df[[conditional_col[[1L]]]][[1L]]) else NA_real_
  if (is.na(marginal) && is.na(conditional)) return(NULL)

  out <- list(r2_method = "nakagawa")
  if (!is.na(marginal)) out$r2_marginal <- marginal
  if (!is.na(conditional)) out$r2_conditional <- conditional
  out
}

ms_lmer_satterthwaite_info <- function(x) {
  if (!inherits(x, "lmerMod") || inherits(x, "glmerMod")) return(NULL)
  if (!requireNamespace("lmerTest", quietly = TRUE)) return(NULL)

  converter <- tryCatch(
    getExportedValue("lmerTest", "as_lmerModLmerTest"),
    error = function(e) NULL
  )
  if (!is.function(converter)) return(NULL)

  test_model <- if (inherits(x, "lmerModLmerTest")) {
    x
  } else {
    tryCatch(
      suppressMessages(suppressWarnings(converter(x))),
      error = function(e) NULL
    )
  }
  if (is.null(test_model)) return(NULL)

  s <- tryCatch(
    suppressMessages(suppressWarnings(summary(test_model, ddf = "Satterthwaite"))),
    error = function(e) NULL
  )
  coefs <- tryCatch(stats::coef(s), error = function(e) NULL)
  if (is.null(coefs) || !is.matrix(coefs) ||
      !all(c("df", "Pr(>|t|)") %in% colnames(coefs))) {
    s <- NULL
  }

  term_tests <- ms_lmer_model_term_tests_satterthwaite(test_model)
  if (is.null(s) && !length(term_tests)) return(NULL)

  list(
    model = test_model,
    summary = s,
    model_term_tests = term_tests
  )
}

ms_lmer_model_term_tests_satterthwaite <- function(x) {
  tf <- tryCatch(stats::terms(x), error = function(e) NULL)
  if (is.null(tf)) return(list())

  term_labels <- attr(tf, "term.labels") %||% character(0)
  orders <- attr(tf, "order") %||% integer(0)
  if (!length(term_labels)) return(list())
  if (length(orders) != length(term_labels)) {
    orders <- rep(1L, length(term_labels))
  }

  tests <- tryCatch(
    suppressMessages(suppressWarnings(stats::anova(
      x,
      type = 3,
      ddf = "Satterthwaite"
    ))),
    error = function(e) NULL
  )
  if (is.null(tests) || !inherits(tests, "data.frame") || !nrow(tests)) {
    return(list())
  }

  rn <- rownames(tests) %||% character(0)
  df1_col <- match("NumDF", names(tests), nomatch = 0L)
  df2_col <- match("DenDF", names(tests), nomatch = 0L)
  f_col <- match("F value", names(tests), nomatch = 0L)
  p_col <- grep("^Pr\\(>F\\)$", names(tests), perl = TRUE)
  p_col <- if (length(p_col)) p_col[[1L]] else 0L
  ss_col <- match("Sum Sq", names(tests), nomatch = 0L)
  ms_col <- match("Mean Sq", names(tests), nomatch = 0L)
  if (!df1_col || !df2_col || !f_col) return(list())

  is_reml <- isTRUE(tryCatch(lme4::isREML(x), error = function(e) FALSE))
  out <- vector("list", length(term_labels))
  for (i in seq_along(term_labels)) {
    term <- term_labels[[i]]
    idx <- match(term, rn, nomatch = 0L)
    if (!idx) next

    df1 <- ms_safe_numeric(tests[[df1_col]][idx])
    df2 <- ms_safe_numeric(tests[[df2_col]][idx])
    f_value <- ms_safe_numeric(tests[[f_col]][idx])
    p_value <- if (p_col) ms_safe_numeric(tests[[p_col]][idx]) else NA_real_
    row <- ms_model_term_test_base_row(
      term = term,
      statistic_name = "F",
      statistic_value = f_value,
      df = c(df1, df2),
      p_value = p_value,
      method = "anova_f_satterthwaite"
    )
    if (is.null(row)) next
    row$term_type <- if (orders[[i]] >= 2L || grepl(":", term, fixed = TRUE)) {
      "interaction"
    } else {
      "main"
    }
    row$test_scope <- "fixed_effect"
    row$model_fit <- if (is_reml) "REML" else "ML"
    row$ddf_method <- "satterthwaite"
    row$ss_type <- "type_iii"
    if (ss_col) row$sum_sq <- ms_safe_numeric(tests[[ss_col]][idx])
    if (ms_col) row$mean_sq <- ms_safe_numeric(tests[[ms_col]][idx])
    if (!is.na(f_value) && !is.na(df1) && !is.na(df2) && df2 > 0) {
      row$effect <- list(
        name = "eta_sq_partial",
        value = ms_safe_numeric((f_value * df1) / ((f_value * df1) + df2))
      )
    }
    out[[i]] <- row
  }

  Filter(Negate(is.null), out)
}

ms_lmer_model_term_tests <- function(x) {
  if (!inherits(x, "lmerMod") || inherits(x, "glmerMod")) return(list())

  tf <- tryCatch(stats::terms(x), error = function(e) NULL)
  if (is.null(tf)) return(list())

  term_labels <- attr(tf, "term.labels") %||% character(0)
  orders <- attr(tf, "order") %||% integer(0)
  if (!length(term_labels)) return(list())
  if (length(orders) != length(term_labels)) {
    orders <- rep(1L, length(term_labels))
  }

  is_reml <- isTRUE(tryCatch(lme4::isREML(x), error = function(e) FALSE))
  test_model <- if (is_reml) {
    tryCatch(lme4::refitML(x), error = function(e) NULL)
  } else {
    x
  }
  if (is.null(test_model)) return(list())

  drop_rows <- ms_model_term_tests_drop1_chisq(test_model, term_labels)
  out <- vector("list", length(term_labels))
  for (i in seq_along(term_labels)) {
    term <- term_labels[[i]]
    row <- drop_rows[[term]]
    if (is.null(row)) next
    row$term_type <- if (orders[[i]] >= 2L || grepl(":", term, fixed = TRUE)) {
      "interaction"
    } else {
      "main"
    }
    row$test_scope <- "fixed_effect"
    row$model_fit <- "ML"
    if (is_reml) {
      row$method <- "drop1_chisq_ml_refit"
      row$refit <- "ML"
      row$refit_from <- "REML"
    }
    out[[i]] <- row
  }

  Filter(Negate(is.null), out)
}

ms_lmer_coefficients <- function(x, coefs_mat) {
  cn <- colnames(coefs_mat)
  est_idx <- match("Estimate", cn, nomatch = 1L)
  se_idx  <- match("Std. Error", cn, nomatch = 2L)
  t_idx   <- match("t value", cn, nomatch = 3L)
  # lmerTest::lmer summaries add `df` and `Pr(>|t|)` (Satterthwaite);
  # plain lme4 lmer has neither. Read them only when present.
  df_idx  <- match("df", cn, nomatch = 0L)
  p_idx   <- match("Pr(>|t|)", cn, nomatch = 0L)

  mm <- tryCatch(stats::model.matrix(x), error = function(e) NULL)
  coef_assign <- if (!is.null(mm)) attr(mm, "assign") else NULL
  mf <- tryCatch(stats::model.frame(x), error = function(e) NULL)
  coef_label_map <- ms_coefficient_label_map(x, mm = mm, mf = mf)
  term_labels <- attr(stats::terms(x), "term.labels") %||% character(0)
  terms_named <- rownames(coefs_mat)

  lapply(seq_len(nrow(coefs_mat)), function(i) {
    estimate <- ms_safe_numeric(coefs_mat[i, est_idx])
    std_error <- ms_safe_numeric(coefs_mat[i, se_idx])
    row <- list(
      term = terms_named[i],
      estimate = estimate,
      estimate_name = "B",
      std_error = std_error,
      statistic = ms_safe_numeric(coefs_mat[i, t_idx]),
      statistic_label = "t"
    )
    if (df_idx > 0L) row$df <- ms_safe_numeric(coefs_mat[i, df_idx])
    if (p_idx > 0L) {
      row$p_value <- ms_safe_numeric(coefs_mat[i, p_idx])
      row$p_value_method <- "satterthwaite_t"
    }
    if (!is.na(estimate) && !is.na(std_error)) {
      crit <- if (!is.na(row$df %||% NA_real_) && row$df > 0) {
        stats::qt(0.975, df = row$df)
      } else {
        stats::qnorm(0.975)
      }
      row$ci_lower <- estimate - crit * std_error
      row$ci_upper <- estimate + crit * std_error
      row$ci_method <- if (!is.na(row$df %||% NA_real_) && row$df > 0) {
        "satterthwaite_t"
      } else {
        "wald"
      }
    }
    if (!is.null(coef_assign) && i <= length(coef_assign) &&
        coef_assign[[i]] > 0L && coef_assign[[i]] <= length(term_labels)) {
      row$term_source <- term_labels[[coef_assign[[i]]]]
    }
    row <- ms_apply_coefficient_label(row, coef_label_map)
    row
  })
}

ms_lmer_focal_terms <- function(coefficients) {
  terms <- vapply(coefficients, function(coef) {
    if (!is.null(coef$term_source)) {
      coef$term_source
    } else {
      coef$term %||% NA_character_
    }
  }, character(1))
  terms <- terms[!is.na(terms)]
  terms <- terms[!tolower(terms) %in% c("(intercept)", "intercept")]
  unique(terms)
}

ms_mixed_role_prefers_categorical_means <- function(term_roles) {
  terms <- term_roles$terms %||% list()
  if (!length(terms)) return(FALSE)
  has_interaction <- any(vapply(terms, function(term) {
    identical(term$type, "interaction")
  }, logical(1)))
  if (has_interaction) return(FALSE)
  focal_terms <- Filter(function(term) identical(term$role, "focal"), terms)
  control_terms <- Filter(function(term) identical(term$role, "control"), terms)
  has_focal_factor <- any(vapply(focal_terms, function(term) {
    identical(term$type, "factor")
  }, logical(1)))
  has_control_numeric <- any(vapply(control_terms, function(term) {
    identical(term$type, "numeric")
  }, logical(1)))
  has_focal_factor && has_control_numeric
}

ms_mixed_focal_factor_terms <- function(term_roles) {
  terms <- term_roles$terms %||% list()
  focal_factors <- Filter(function(term) {
    identical(term$role, "focal") && identical(term$type, "factor")
  }, terms)
  unique(vapply(focal_factors, function(term) term$name %||% "", character(1)))
}

ms_mixed_preferred_term_index <- function(term_info, candidate_idx, preferred_terms) {
  preferred_terms <- unique(as.character(preferred_terms %||% character(0)))
  preferred_terms <- preferred_terms[nzchar(preferred_terms)]
  if (!length(candidate_idx) || !length(preferred_terms)) return(NULL)
  for (candidate in preferred_terms) {
    for (idx in candidate_idx) {
      info <- term_info[[idx]]
      keys <- unique(c(info$term %||% "", info$variable %||% "", info$label %||% ""))
      if (candidate %in% keys) return(idx)
    }
  }
  NULL
}

ms_lmer_groups <- function(x) {
  flist <- tryCatch(lme4::getME(x, "flist"), error = function(e) NULL)
  if (is.null(flist) || !length(flist)) return(list())
  lapply(names(flist), function(name) {
    list(
      label = name,
      n = length(unique(flist[[name]]))
    )
  })
}

ms_lmer_random_terms <- function(x) {
  bars <- tryCatch(
    suppressWarnings(lme4::findbars(stats::formula(x))),
    error = function(e) NULL
  )
  if (is.null(bars) || !length(bars)) return(character(0))
  vapply(bars, function(term) {
    paste0("(", paste(deparse(term, width.cutoff = 500L), collapse = " "), ")")
  }, character(1))
}

ms_lmer_random_effects <- function(x) {
  vc <- tryCatch(lme4::VarCorr(x), error = function(e) NULL)
  if (is.null(vc)) return(list())

  rows <- list()
  for (group in names(vc)) {
    mat <- vc[[group]]
    vars <- diag(mat)
    stddev <- attr(mat, "stddev")
    corr <- attr(mat, "correlation")
    term_names <- names(vars)
    for (i in seq_along(vars)) {
      row <- list(
        group = group,
        name = term_names[[i]],
        variance = ms_safe_numeric(vars[[i]]),
        std_dev = ms_safe_numeric(stddev[[i]])
      )
      if (!is.null(corr) && i > 1L && ncol(corr) >= 1L) {
        row$corr <- ms_safe_numeric(corr[i, 1])
      }
      rows[[length(rows) + 1L]] <- row
    }
  }

  residual_sd <- attr(vc, "sc")
  if (!is.null(residual_sd) && is.finite(residual_sd)) {
    rows[[length(rows) + 1L]] <- list(
      group = "Residual",
      name = "Residual",
      variance = ms_safe_numeric(residual_sd^2),
      std_dev = ms_safe_numeric(residual_sd)
    )
  }

  rows
}

ms_lmer_emmeans_interaction_plot_data <- function(x, conf.level = 0.95,
                                                  max_levels = 8L,
                                                  grid_points = 80L) {
  if (!requireNamespace("emmeans", quietly = TRUE)) return(NULL)

  tf <- tryCatch(stats::terms(x), error = function(e) NULL)
  mf <- tryCatch(stats::model.frame(x), error = function(e) NULL)
  if (is.null(tf) || is.null(mf) || nrow(mf) < 2L) return(NULL)

  term_labels <- attr(tf, "term.labels") %||% character(0)
  orders <- attr(tf, "order") %||% integer(0)
  if (length(term_labels) == 0L || length(orders) != length(term_labels)) return(NULL)
  if (any(orders > 2L)) return(NULL)

  interaction_terms <- term_labels[orders == 2L & grepl(":", term_labels, fixed = TRUE)]
  if (length(interaction_terms) != 1L) return(NULL)
  interaction_term <- interaction_terms[[1L]]
  parts <- strsplit(interaction_term, ":", fixed = TRUE)[[1L]]
  parts <- trimws(parts)
  if (length(parts) != 2L || any(!nzchar(parts))) return(NULL)
  if (!all(vapply(parts, ms_interaction_component_supported, logical(1)))) {
    return(NULL)
  }

  interaction_variables <- vapply(parts, function(part) {
    ms_interaction_component_variable(part) %||% ""
  }, character(1))
  if (any(!nzchar(interaction_variables)) || length(unique(interaction_variables)) != 2L) {
    return(NULL)
  }
  if (!ms_interaction_terms_within_components(interaction_terms, interaction_variables)) {
    return(NULL)
  }

  var_map <- ms_model_frame_variable_map(mf)
  components <- lapply(parts, ms_interaction_component_info,
                       var_map = var_map, max_levels = max_levels)
  if (any(vapply(components, is.null, logical(1)))) return(NULL)
  is_numeric <- vapply(components, function(info) identical(info$type, "numeric"), logical(1))
  is_categorical <- vapply(components, function(info) identical(info$type, "categorical"), logical(1))
  if (sum(is_numeric) == 1L && sum(is_categorical) == 1L) {
    return(ms_lmer_emmeans_continuous_interaction_plot_data(
      x = x,
      mf = mf,
      interaction_term = interaction_term,
      components = components,
      conf.level = conf.level,
      max_levels = max_levels,
      grid_points = grid_points
    ))
  }
  if (!all(vapply(components, function(info) identical(info$type, "categorical"), logical(1)))) {
    return(NULL)
  }

  ordered <- ms_lmer_interaction_role_order(x, mf, components)
  x_info <- ordered$x
  moderator_info <- ordered$moderator
  x_set <- ms_interaction_categorical_set(x_info, max_levels = max_levels)
  moderator_set <- ms_interaction_categorical_set(moderator_info, max_levels = max_levels)
  if (is.null(x_set) || is.null(moderator_set)) return(NULL)

  emm_formula <- stats::as.formula(paste("~", paste(c(x_info$term, moderator_info$term), collapse = " * ")))
  emm <- tryCatch(
    emmeans::emmeans(x, specs = emm_formula),
    error = function(e) NULL
  )
  if (is.null(emm)) return(NULL)
  emm_df <- tryCatch(
    as.data.frame(summary(emm, infer = c(TRUE, FALSE), level = conf.level)),
    error = function(e) NULL
  )
  if (is.null(emm_df) || !nrow(emm_df)) return(NULL)

  x_col <- ms_lmer_emmeans_col(emm_df, c(x_info$variable, x_info$term))
  moderator_col <- ms_lmer_emmeans_col(emm_df, c(moderator_info$variable, moderator_info$term))
  estimate_col <- ms_lmer_emmeans_col(emm_df, c("emmean", "response", "prob", "rate", "estimate"))
  se_col <- ms_lmer_emmeans_col(emm_df, c("SE", "std.error", "std_error"))
  df_col <- ms_lmer_emmeans_col(emm_df, c("df"))
  lower_col <- ms_lmer_emmeans_col(emm_df, c("lower.CL", "asymp.LCL", "lower.HPD", "LCL"))
  upper_col <- ms_lmer_emmeans_col(emm_df, c("upper.CL", "asymp.UCL", "upper.HPD", "UCL"))
  if (is.na(x_col) || is.na(moderator_col) || is.na(estimate_col)) return(NULL)

  grid <- list()
  for (mod_index in seq_along(moderator_set$values)) {
    mod_value <- moderator_set$values[[mod_index]]
    for (x_index in seq_along(x_set$values)) {
      x_value <- x_set$values[[x_index]]
      row_index <- which(
        as.character(emm_df[[x_col]]) == as.character(x_value) &
          as.character(emm_df[[moderator_col]]) == as.character(mod_value)
      )
      if (!length(row_index)) next
      i <- row_index[[1L]]
      row <- list(
        x = x_index,
        x_value = as.character(x_value),
        x_label = x_set$labels[[x_index]],
        moderator_value = as.character(mod_value),
        moderator_label = moderator_set$labels[[mod_index]],
        estimate = ms_safe_numeric(emm_df[[estimate_col]][[i]])
      )
      if (!is.na(se_col)) row$se <- ms_safe_numeric(emm_df[[se_col]][[i]])
      if (!is.na(df_col)) row$df <- ms_safe_numeric(emm_df[[df_col]][[i]])
      if (!is.na(lower_col)) row$ci_lower <- ms_safe_numeric(emm_df[[lower_col]][[i]])
      if (!is.na(upper_col)) row$ci_upper <- ms_safe_numeric(emm_df[[upper_col]][[i]])
      grid[[length(grid) + 1L]] <- row
    }
  }
  grid <- Filter(function(row) {
    is.list(row) &&
      !is.null(row$x) &&
      !is.na(ms_safe_numeric(row$x)) &&
      !is.null(row$estimate) &&
      !is.na(ms_safe_numeric(row$estimate))
  }, grid)
  if (length(grid) < length(x_set$values) * length(moderator_set$values)) return(NULL)

  outcome <- ms_lm_outcome_label(x)
  subject <- ms_lmer_primary_group_summary(x, mf)
  out <- list(
    interaction_term = interaction_term,
    interaction_kind = "categorical_by_categorical",
    source = "lmer_emmeans",
    mean_kind = "estimated_marginal",
    variables = c(x_info$variable, moderator_info$variable),
    x = list(
      variable = x_info$variable,
      term = x_info$term,
      label = x_info$label,
      type = "categorical",
      range = I(c(1, length(x_set$values))),
      levels = ms_interaction_set_levels(x_set)
    ),
    moderator = list(
      variable = moderator_info$variable,
      term = moderator_info$term,
      label = moderator_info$label,
      type = "categorical",
      levels = ms_interaction_set_levels(moderator_set)
    ),
    grid = grid,
    outcome = outcome %||% NULL,
    y_label = if (!is.null(outcome) && nzchar(outcome)) {
      paste("Estimated marginal mean", outcome)
    } else {
      "Estimated marginal mean"
    },
    scale = "response",
    ci_level = conf.level,
    ci_method = "emmeans",
    model_family = "gaussian",
    model_link = "identity",
    model_fit = if (isTRUE(tryCatch(lme4::isREML(x), error = function(e) FALSE))) "REML" else "ML",
    bounded_response = FALSE
  )
  if (!is.null(subject)) out$subject <- subject
  out
}

ms_lmer_emmeans_continuous_interaction_plot_data <- function(x, mf, interaction_term,
                                                            components,
                                                            conf.level = 0.95,
                                                            max_levels = 8L,
                                                            grid_points = 80L) {
  x_info <- components[[which(vapply(components, function(info) identical(info$type, "numeric"), logical(1)))[[1L]]]]
  moderator_info <- components[[which(vapply(components, function(info) identical(info$type, "categorical"), logical(1)))[[1L]]]]
  moderator_set <- ms_interaction_categorical_set(moderator_info, max_levels = max_levels)
  if (is.null(moderator_set)) return(NULL)

  x_values <- ms_interaction_numeric_values(x_info$values)
  if (length(x_values) < 2L) return(NULL)
  x_range <- range(x_values, na.rm = TRUE)
  if (!all(is.finite(x_range)) || x_range[[1L]] == x_range[[2L]]) return(NULL)
  x_grid <- seq(x_range[[1L]], x_range[[2L]], length.out = max(12L, as.integer(grid_points)))

  emm_formula <- stats::as.formula(paste("~", paste(c(x_info$term, moderator_info$term), collapse = " * ")))
  at_values <- setNames(list(x_grid), x_info$variable)
  emm <- tryCatch(
    emmeans::emmeans(x, specs = emm_formula, at = at_values),
    error = function(e) NULL
  )
  if (is.null(emm)) return(NULL)
  emm_df <- tryCatch(
    as.data.frame(summary(emm, infer = c(TRUE, FALSE), level = conf.level)),
    error = function(e) NULL
  )
  if (is.null(emm_df) || !nrow(emm_df)) return(NULL)

  x_col <- ms_lmer_emmeans_col(emm_df, c(x_info$variable, x_info$term))
  moderator_col <- ms_lmer_emmeans_col(emm_df, c(moderator_info$variable, moderator_info$term))
  estimate_col <- ms_lmer_emmeans_col(emm_df, c("emmean", "response", "prob", "rate", "estimate"))
  se_col <- ms_lmer_emmeans_col(emm_df, c("SE", "std.error", "std_error"))
  df_col <- ms_lmer_emmeans_col(emm_df, c("df"))
  lower_col <- ms_lmer_emmeans_col(emm_df, c("lower.CL", "asymp.LCL", "lower.HPD", "LCL"))
  upper_col <- ms_lmer_emmeans_col(emm_df, c("upper.CL", "asymp.UCL", "upper.HPD", "UCL"))
  if (is.na(x_col) || is.na(moderator_col) || is.na(estimate_col)) return(NULL)

  x_observed <- suppressWarnings(as.numeric(emm_df[[x_col]]))
  tolerance <- max(1e-8, diff(x_range) * 1e-8)
  grid <- list()
  for (mod_index in seq_along(moderator_set$values)) {
    mod_value <- moderator_set$values[[mod_index]]
    for (x_value in x_grid) {
      row_index <- which(
        abs(x_observed - x_value) <= tolerance &
          as.character(emm_df[[moderator_col]]) == as.character(mod_value)
      )
      if (!length(row_index)) next
      i <- row_index[[1L]]
      row <- list(
        x = ms_safe_numeric(x_value),
        moderator_value = as.character(mod_value),
        moderator_label = moderator_set$labels[[mod_index]],
        estimate = ms_safe_numeric(emm_df[[estimate_col]][[i]])
      )
      if (!is.na(se_col)) row$se <- ms_safe_numeric(emm_df[[se_col]][[i]])
      if (!is.na(df_col)) row$df <- ms_safe_numeric(emm_df[[df_col]][[i]])
      if (!is.na(lower_col)) row$ci_lower <- ms_safe_numeric(emm_df[[lower_col]][[i]])
      if (!is.na(upper_col)) row$ci_upper <- ms_safe_numeric(emm_df[[upper_col]][[i]])
      grid[[length(grid) + 1L]] <- row
    }
  }
  grid <- Filter(function(row) {
    is.list(row) &&
      !is.null(row$x) &&
      !is.na(ms_safe_numeric(row$x)) &&
      !is.null(row$estimate) &&
      !is.na(ms_safe_numeric(row$estimate))
  }, grid)
  if (length(grid) < length(x_grid) * length(moderator_set$values)) return(NULL)

  outcome <- ms_lm_outcome_label(x)
  subject <- ms_lmer_primary_group_summary(x, mf)
  out <- list(
    interaction_term = interaction_term,
    interaction_kind = "continuous_by_categorical",
    source = "lmer_emmeans",
    mean_kind = "estimated_marginal",
    variables = c(x_info$variable, moderator_info$variable),
    x = list(
      variable = x_info$variable,
      term = x_info$term,
      label = x_info$label,
      type = "numeric",
      range = I(ms_safe_numeric(x_range))
    ),
    moderator = list(
      variable = moderator_info$variable,
      term = moderator_info$term,
      label = moderator_info$label,
      type = "categorical",
      levels = ms_interaction_set_levels(moderator_set)
    ),
    grid = grid,
    outcome = outcome %||% NULL,
    y_label = if (!is.null(outcome) && nzchar(outcome)) {
      paste("Marginal", outcome)
    } else {
      "Marginal response"
    },
    scale = "response",
    ci_level = conf.level,
    ci_method = "emmeans",
    model_family = "gaussian",
    model_link = "identity",
    model_fit = if (isTRUE(tryCatch(lme4::isREML(x), error = function(e) FALSE))) "REML" else "ML",
    bounded_response = FALSE
  )
  if (!is.null(subject)) out$subject <- subject
  out
}

ms_lmer_emmeans_continuous_main_effect_plot_data <- function(x, conf.level = 0.95,
                                                             grid_points = 80L,
                                                             max_levels = 8L) {
  if (!requireNamespace("emmeans", quietly = TRUE)) return(NULL)

  tf <- tryCatch(stats::terms(x), error = function(e) NULL)
  mf <- tryCatch(stats::model.frame(x), error = function(e) NULL)
  if (is.null(tf) || is.null(mf) || nrow(mf) < 2L) return(NULL)

  term_labels <- attr(tf, "term.labels") %||% character(0)
  orders <- attr(tf, "order") %||% integer(0)
  if (length(term_labels) == 0L || length(orders) != length(term_labels)) return(NULL)
  if (any(orders != 1L)) return(NULL)

  var_map <- ms_model_frame_variable_map(mf)
  term_info <- lapply(term_labels, function(term) {
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
  if (any(vapply(term_info, is.null, logical(1)))) return(NULL)

  numeric_idx <- which(vapply(term_info, function(info) identical(info$type, "numeric"), logical(1)))
  if (!length(numeric_idx)) return(NULL)
  x_info <- term_info[[numeric_idx[[1L]]]]
  x_values <- ms_interaction_numeric_values(x_info$values)
  if (length(x_values) < 2L) return(NULL)
  x_range <- range(x_values, na.rm = TRUE)
  if (!all(is.finite(x_range)) || x_range[[1L]] == x_range[[2L]]) return(NULL)
  x_grid <- seq(x_range[[1L]], x_range[[2L]], length.out = max(12L, as.integer(grid_points)))

  other_numeric <- term_info[numeric_idx[-1L]]
  covariates <- lapply(other_numeric, function(info) {
    values <- ms_interaction_numeric_values(info$values)
    if (length(values) < 2L) return(NULL)
    value <- ms_safe_numeric(mean(values, na.rm = TRUE))
    if (is.na(value) || !is.finite(value)) return(NULL)
    list(
      variable = info$variable,
      term = info$term,
      label = info$label,
      value = value,
      value_label = ms_interaction_format_number(value),
      type = "mean"
    )
  })
  covariates <- Filter(Negate(is.null), covariates)

  at_values <- c(
    setNames(list(x_grid), x_info$variable),
    setNames(
      lapply(covariates, function(covariate) covariate$value),
      vapply(covariates, function(covariate) covariate$variable, character(1))
    )
  )
  emm_formula <- tryCatch(stats::as.formula(paste("~", x_info$term)),
                          error = function(e) NULL)
  if (is.null(emm_formula)) return(NULL)
  emm <- tryCatch(
    emmeans::emmeans(x, specs = emm_formula, at = at_values),
    error = function(e) NULL
  )
  if (is.null(emm)) return(NULL)
  emm_df <- tryCatch(
    as.data.frame(summary(emm, infer = c(TRUE, FALSE), level = conf.level)),
    error = function(e) NULL
  )
  if (is.null(emm_df) || nrow(emm_df) < 2L) return(NULL)

  x_col <- ms_lmer_emmeans_col(emm_df, c(x_info$variable, x_info$term))
  estimate_col <- ms_lmer_emmeans_col(emm_df, c("emmean", "response", "prob", "rate", "estimate"))
  se_col <- ms_lmer_emmeans_col(emm_df, c("SE", "std.error", "std_error"))
  df_col <- ms_lmer_emmeans_col(emm_df, c("df"))
  lower_col <- ms_lmer_emmeans_col(emm_df, c("lower.CL", "asymp.LCL", "lower.HPD", "LCL"))
  upper_col <- ms_lmer_emmeans_col(emm_df, c("upper.CL", "asymp.UCL", "upper.HPD", "UCL"))
  if (is.na(x_col) || is.na(estimate_col)) return(NULL)

  x_observed <- suppressWarnings(as.numeric(emm_df[[x_col]]))
  tolerance <- max(1e-8, diff(x_range) * 1e-8)
  grid <- lapply(x_grid, function(x_value) {
    row_index <- which(abs(x_observed - x_value) <= tolerance)
    if (!length(row_index)) return(NULL)
    i <- row_index[[1L]]
    row <- list(
      x = ms_safe_numeric(x_value),
      moderator_value = "estimate",
      moderator_label = "Estimate",
      estimate = ms_safe_numeric(emm_df[[estimate_col]][[i]])
    )
    if (!is.na(se_col)) row$se <- ms_safe_numeric(emm_df[[se_col]][[i]])
    if (!is.na(df_col)) row$df <- ms_safe_numeric(emm_df[[df_col]][[i]])
    if (!is.na(lower_col)) row$ci_lower <- ms_safe_numeric(emm_df[[lower_col]][[i]])
    if (!is.na(upper_col)) row$ci_upper <- ms_safe_numeric(emm_df[[upper_col]][[i]])
    row
  })
  grid <- Filter(function(row) {
    is.list(row) &&
      !is.null(row$x) &&
      !is.na(ms_safe_numeric(row$x)) &&
      !is.null(row$estimate) &&
      !is.na(ms_safe_numeric(row$estimate))
  }, grid)
  if (length(grid) < length(x_grid)) return(NULL)

  categorical_idx <- which(vapply(term_info, function(info) identical(info$type, "categorical"), logical(1)))
  other_factors <- term_info[categorical_idx]
  outcome <- ms_lm_outcome_label(x)
  subject <- ms_lmer_primary_group_summary(x, mf)
  out <- list(
    interaction_term = x_info$term,
    interaction_kind = "continuous_main_effect",
    source = "lmer_emmeans",
    mean_kind = "estimated_marginal",
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
    marginalized_terms = lapply(other_factors, function(info) {
      list(
        variable = info$variable,
        term = info$term,
        label = info$label,
        rule = "estimated_marginal"
      )
    }),
    held_constant = lapply(covariates, function(covariate) {
      list(
        variable = covariate$variable,
        label = covariate$label,
        value = as.character(covariate$value),
        value_label = covariate$value_label,
        type = covariate$type
      )
    }),
    outcome = outcome %||% NULL,
    y_label = if (!is.null(outcome) && nzchar(outcome)) {
      paste("Marginal", outcome)
    } else {
      "Marginal response"
    },
    scale = "response",
    ci_level = conf.level,
    ci_method = "emmeans",
    model_family = "gaussian",
    model_link = "identity",
    model_fit = if (isTRUE(tryCatch(lme4::isREML(x), error = function(e) FALSE))) "REML" else "ML",
    bounded_response = FALSE
  )
  if (!is.null(subject)) out$subject <- subject
  out
}

ms_lmer_emmeans_main_effect_means_figure_data <- function(x, conf.level = 0.95,
                                                          max_levels = 12L,
                                                          require_within = FALSE,
                                                          preferred_terms = NULL,
                                                          require_preferred = FALSE) {
  if (!requireNamespace("emmeans", quietly = TRUE)) return(NULL)

  tf <- tryCatch(stats::terms(x), error = function(e) NULL)
  mf <- tryCatch(stats::model.frame(x), error = function(e) NULL)
  if (is.null(tf) || is.null(mf) || nrow(mf) < 2L) return(NULL)

  term_labels <- attr(tf, "term.labels") %||% character(0)
  orders <- attr(tf, "order") %||% integer(0)
  if (length(term_labels) == 0L || length(orders) != length(term_labels)) return(NULL)
  if (any(orders != 1L)) return(NULL)

  var_map <- ms_model_frame_variable_map(mf)
  term_info <- lapply(term_labels, function(term) {
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
  if (any(vapply(term_info, is.null, logical(1)))) return(NULL)

  categorical_idx <- which(vapply(term_info, function(info) {
    identical(info$type, "categorical")
  }, logical(1)))
  if (!length(categorical_idx)) return(NULL)

  selected_idx <- ms_lmer_main_effect_factor_index(
    x,
    mf,
    term_info,
    categorical_idx,
    preferred_terms = preferred_terms,
    require_preferred = require_preferred
  )
  if (is.null(selected_idx)) return(NULL)
  factor_info <- term_info[[selected_idx]]
  connect_levels <- ms_lmer_factor_varies_within_any_group(x, mf, factor_info)
  if (isTRUE(require_within) && !connect_levels) return(NULL)
  factor_set <- ms_interaction_categorical_set(factor_info, max_levels = max_levels)
  if (is.null(factor_set)) return(NULL)

  numeric_info <- term_info[vapply(term_info, function(info) identical(info$type, "numeric"), logical(1))]
  covariates <- lapply(numeric_info, function(info) {
    values <- ms_interaction_numeric_values(info$values)
    if (length(values) < 2L) return(NULL)
    value <- ms_safe_numeric(mean(values, na.rm = TRUE))
    if (is.na(value) || !is.finite(value)) return(NULL)
    list(
      variable = info$variable,
      term = info$term,
      label = info$label,
      value = value
    )
  })
  covariates <- Filter(Negate(is.null), covariates)

  at_values <- setNames(
    lapply(covariates, function(covariate) covariate$value),
    vapply(covariates, function(covariate) covariate$variable, character(1))
  )
  emm_formula <- tryCatch(stats::as.formula(paste("~", factor_info$term)),
                          error = function(e) NULL)
  if (is.null(emm_formula)) return(NULL)
  emm_args <- list(object = x, specs = emm_formula)
  if (length(at_values) > 0L) emm_args$at <- at_values
  emm <- tryCatch(
    do.call(emmeans::emmeans, emm_args),
    error = function(e) NULL
  )
  if (is.null(emm)) return(NULL)
  emm_df <- tryCatch(
    as.data.frame(summary(emm, infer = c(TRUE, FALSE), level = conf.level)),
    error = function(e) NULL
  )
  if (is.null(emm_df) || nrow(emm_df) < 2L) return(NULL)

  factor_col <- ms_lmer_emmeans_col(emm_df, c(factor_info$variable, factor_info$term))
  estimate_col <- ms_lmer_emmeans_col(emm_df, c("emmean", "response", "prob", "rate", "estimate"))
  se_col <- ms_lmer_emmeans_col(emm_df, c("SE", "std.error", "std_error"))
  df_col <- ms_lmer_emmeans_col(emm_df, c("df"))
  lower_col <- ms_lmer_emmeans_col(emm_df, c("lower.CL", "asymp.LCL", "lower.HPD", "LCL"))
  upper_col <- ms_lmer_emmeans_col(emm_df, c("upper.CL", "asymp.UCL", "upper.HPD", "UCL"))
  if (is.na(factor_col) || is.na(estimate_col)) return(NULL)

  counts <- table(as.character(factor_info$values[!is.na(factor_info$values)]))
  groups <- lapply(seq_along(factor_set$values), function(level_index) {
    level <- factor_set$values[[level_index]]
    row_index <- which(as.character(emm_df[[factor_col]]) == as.character(level))
    if (!length(row_index)) return(NULL)
    i <- row_index[[1L]]
    n_value <- if (as.character(level) %in% names(counts)) {
      as.integer(counts[[as.character(level)]])
    } else {
      NA_integer_
    }
    row <- list(
      level = as.character(level),
      label = factor_set$labels[[level_index]],
      n = n_value,
      mean = ms_safe_numeric(emm_df[[estimate_col]][[i]])
    )
    if (!is.na(se_col)) row$se <- ms_safe_numeric(emm_df[[se_col]][[i]])
    if (!is.na(df_col)) row$df <- ms_safe_numeric(emm_df[[df_col]][[i]])
    if (!is.na(lower_col)) row$ci_lower <- ms_safe_numeric(emm_df[[lower_col]][[i]])
    if (!is.na(upper_col)) row$ci_upper <- ms_safe_numeric(emm_df[[upper_col]][[i]])
    row
  })
  groups <- Filter(function(row) {
    is.list(row) &&
      !is.null(row$level) &&
      !is.null(row$mean) &&
      !is.na(ms_safe_numeric(row$mean))
  }, groups)
  if (length(groups) < 2L) return(NULL)

  other_factors <- term_info[categorical_idx]
  other_factors <- other_factors[vapply(other_factors, function(info) {
    !identical(info$term, factor_info$term)
  }, logical(1))]
  outcome <- ms_lm_outcome_label(x)
  out <- list(
    mean_kind = "estimated_marginal",
    source = "lmer_emmeans",
    factor = list(
      variable = factor_info$variable,
      term = factor_info$term,
      label = factor_info$label,
      levels = ms_interaction_set_levels(factor_set)
    ),
    groups = groups,
    marginalized_terms = lapply(other_factors, function(info) {
      list(
        variable = info$variable,
        term = info$term,
        label = info$label,
        rule = "estimated_marginal"
      )
    }),
    outcome = outcome %||% NULL,
    y_label = if (!is.null(outcome) && nzchar(outcome)) {
      paste("Estimated marginal mean", outcome)
    } else {
      "Estimated marginal mean"
    },
    ci_level = conf.level,
    ci_method = "emmeans",
    model_fit = if (isTRUE(tryCatch(lme4::isREML(x), error = function(e) FALSE))) "REML" else "ML"
  )
  if (connect_levels) out$connect_levels <- TRUE
  if (length(covariates) > 0L) {
    out$covariates <- lapply(covariates, function(covariate) {
      list(
        variable = covariate$variable,
        term = covariate$term,
        label = covariate$label,
        value = covariate$value,
        value_label = ms_interaction_format_number(covariate$value),
        rule = "sample_mean"
      )
    })
    out$adjustment <- list(rule = "sample_mean", label = "sample means")
  }
  subject <- ms_lmer_primary_group_summary(x, mf)
  if (!is.null(subject)) out$subject <- subject
  out
}

ms_lmer_main_effect_factor_index <- function(x, mf, term_info, categorical_idx,
                                             preferred_terms = NULL,
                                             require_preferred = FALSE) {
  preferred_idx <- ms_mixed_preferred_term_index(term_info, categorical_idx, preferred_terms)
  if (!is.null(preferred_idx)) return(preferred_idx)
  if (isTRUE(require_preferred)) return(NULL)

  groups <- tryCatch(lme4::getME(x, "flist"), error = function(e) NULL)
  if (!is.null(groups) && length(groups)) {
    for (group_name in names(groups)) {
      group_values <- if (group_name %in% names(mf)) mf[[group_name]] else groups[[group_name]]
      if (length(group_values) != nrow(mf)) next
      varies <- vapply(term_info[categorical_idx], function(info) {
        ms_lmer_varies_within_group(info$values, group_values)
      }, logical(1))
      if (any(varies)) return(categorical_idx[which(varies)[[1L]]])
    }
  }
  categorical_idx[[1L]]
}

ms_lmer_factor_varies_within_any_group <- function(x, mf, info) {
  groups <- tryCatch(lme4::getME(x, "flist"), error = function(e) NULL)
  if (is.null(groups) || !length(groups) || is.null(info)) return(FALSE)
  for (group_name in names(groups)) {
    group_values <- if (group_name %in% names(mf)) mf[[group_name]] else groups[[group_name]]
    if (length(group_values) != nrow(mf)) next
    if (ms_lmer_varies_within_group(info$values, group_values)) return(TRUE)
  }
  FALSE
}

ms_lmer_interaction_role_order <- function(x, mf, components) {
  groups <- tryCatch(lme4::getME(x, "flist"), error = function(e) NULL)
  if (!is.null(groups) && length(groups)) {
    for (group_name in names(groups)) {
      group_values <- if (group_name %in% names(mf)) mf[[group_name]] else groups[[group_name]]
      if (length(group_values) != nrow(mf)) next
      varies <- vapply(components, function(info) {
        ms_lmer_varies_within_group(info$values, group_values)
      }, logical(1))
      if (sum(varies) == 1L) {
        x_idx <- which(varies)[[1L]]
        return(list(
          x = components[[x_idx]],
          moderator = components[[setdiff(seq_along(components), x_idx)[[1L]]]]
        ))
      }
    }
  }
  list(x = components[[1L]], moderator = components[[2L]])
}

ms_lmer_varies_within_group <- function(values, group_values) {
  if (length(values) != length(group_values)) return(FALSE)
  ok <- !is.na(values) & !is.na(group_values)
  if (!any(ok)) return(FALSE)
  counts <- tapply(as.character(values[ok]), as.character(group_values[ok]), function(value) {
    length(unique(value))
  })
  counts <- ms_safe_numeric(counts)
  counts <- counts[!is.na(counts)]
  length(counts) > 0L && mean(counts > 1L) >= 0.5
}

ms_lmer_primary_group_summary <- function(x, mf) {
  groups <- tryCatch(lme4::getME(x, "flist"), error = function(e) NULL)
  if (is.null(groups) || !length(groups)) return(NULL)
  group_name <- names(groups)[[1L]]
  values <- if (group_name %in% names(mf)) mf[[group_name]] else groups[[group_name]]
  list(
    variable = group_name,
    label = ms_model_clean_term(group_name),
    n = length(unique(values[!is.na(values)]))
  )
}

ms_lmer_emmeans_col <- function(df, candidates) {
  if (is.null(df) || !length(names(df))) return(NA_integer_)
  norm <- function(value) tolower(gsub("[^A-Za-z0-9]+", "", value))
  names_norm <- norm(names(df))
  cand_norm <- norm(candidates)
  hit <- match(cand_norm, names_norm, nomatch = 0L)
  hit <- hit[hit > 0L]
  if (length(hit)) hit[[1L]] else NA_integer_
}
