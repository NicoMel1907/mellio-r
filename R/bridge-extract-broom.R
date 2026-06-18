# R bridge -- broom-backed generic model extractor.
#
# This is the shared extraction path for model classes that have reliable
# broom / broom.mixed methods. Specialist extractors can keep owning model
# families where they add important editorial semantics, but lightweight S3
# methods (glm, glmerMod, etc.) should delegate here rather than drift.

#' @rdname mellio_payload
#' @param focal Optional character vector of coefficient or term names to
#'   surface as focal rows in model Result Cards.
#' @param controls Optional character vector of coefficient or term names to
#'   treat as controls/covariates in model Result Cards.
#' @param exponentiate Logical. For GLM-style models, exponentiate
#'   coefficient estimates (for example, odds ratios in logistic models).
#' @param conf.int Logical. Include coefficient confidence intervals when
#'   supported by the broom method.
#' @param conf.level Confidence level for coefficient intervals.
#' @export
mellio_payload.glm <- function(x, ..., .call = NULL,
                           focal = NULL,
                           controls = NULL,
                           exponentiate = FALSE,
                           conf.int = TRUE,
                           conf.level = 0.95) {
  mellio_payload_broom_model(
    x,
    ...,
    .call = .call,
    focal = focal,
    controls = controls,
    exponentiate = exponentiate,
    conf.int = conf.int,
    conf.level = conf.level,
    mixed = FALSE
  )
}

#' @rdname mellio_payload
#' @export
mellio_payload.rlm <- function(x, ..., .call = NULL) {
  mellio_payload_broom_model(x, ..., .call = .call, mixed = FALSE)
}

#' @rdname mellio_payload
#' @export
mellio_payload.gam <- function(x, ..., .call = NULL) {
  mellio_payload_broom_model(x, ..., .call = .call, mixed = FALSE)
}

#' @rdname mellio_payload
#' @export
mellio_payload.lme <- function(x, ..., .call = NULL) {
  mellio_payload_broom_model(x, ..., .call = .call, mixed = TRUE)
}

#' @rdname mellio_payload
#' @export
mellio_payload.glmmTMB <- function(x, ..., .call = NULL) {
  mellio_payload_broom_model(x, ..., .call = .call, mixed = TRUE)
}

mellio_payload_broom_model <- function(x, ..., .call = NULL,
                                   focal = NULL,
                                   controls = NULL,
                                   exponentiate = FALSE,
                                   conf.int = TRUE,
                                   conf.level = 0.95,
                                   mixed = NULL) {
  rlang::check_installed("broom", reason = "to extract model summaries")

  is_mixed <- isTRUE(mixed) || inherits(x, c("glmerMod", "merMod", "lmerMod", "lme", "glmmTMB"))
  if (is_mixed) {
    rlang::check_installed("broom.mixed", reason = "to extract mixed model summaries")
  }

  call_str <- ms_model_call_string(x, .call = .call)
  tidy_df <- ms_broom_tidy_model(
    x,
    mixed = is_mixed,
    exponentiate = exponentiate,
    conf.int = conf.int,
    conf.level = conf.level
  )
  glance_df <- ms_broom_glance_model(x, mixed = is_mixed)

  fixed_df <- ms_broom_fixed_rows(tidy_df)
  if (is.null(fixed_df) || nrow(fixed_df) < 1L) {
    stop("broom did not return fixed-effect coefficients for this model.",
         call. = FALSE)
  }

  family <- tryCatch(stats::family(x), error = function(e) NULL)
  type_info <- ms_broom_type_info(x, family = family, mixed = is_mixed)
  stat_label <- ms_broom_statistic_label(x, mixed = is_mixed)

  fields <- ms_broom_glance_fields(glance_df)
  fields$p_value <- NA_real_
  fields$coefficients <- ms_broom_coefficients(
    x,
    fixed_df,
    statistic_label = stat_label,
    estimate_name = if (isTRUE(exponentiate)) "OR" else NULL
  )
  fields$source <- "R"
  if (!is.null(stat_label)) fields$statistic_label <- stat_label
  if (isTRUE(exponentiate)) fields$coefficient_scale <- "odds_ratio"
  if (isTRUE(conf.int)) fields$conf_level <- conf.level
  has_coefficient_cis <- ms_coefficients_have_ci(fields$coefficients)
  if (inherits(x, "glmerMod") && isTRUE(conf.int) && has_coefficient_cis) {
    fields$coefficient_ci_method <- "wald"
    fields$coefficient_ci_scale <- if (isTRUE(exponentiate)) "odds_ratio" else "link"
    fields$coefficients <- lapply(fields$coefficients, function(row) {
      if (!is.null(row$ci_lower) && !is.null(row$ci_upper)) row$ci_method <- "wald"
      row
    })
  } else if (!is_mixed && inherits(x, "glm") && !inherits(x, "gam") &&
             isTRUE(conf.int) && has_coefficient_cis) {
    fields$coefficient_ci_method <- "profile_likelihood"
    fields$coefficient_ci_scale <- if (isTRUE(exponentiate)) "odds_ratio" else "link"
    fields$coefficients <- lapply(fields$coefficients, function(row) {
      if (!is.null(row$ci_lower) && !is.null(row$ci_upper)) {
        row$ci_method <- "profile_likelihood"
      }
      row
    })
  }

  if (!is.null(family)) {
    fields$model_family <- as.character(family$family %||% "")
    fields$model_link <- as.character(family$link %||% "")
  }
  model_warnings <- tryCatch(ms_model_diagnostics(x), error = function(e) list())
  if (length(model_warnings) > 0L) fields$model_warnings <- model_warnings

  model_term_tests <- if (inherits(x, "glmerMod")) {
    tryCatch(ms_glmer_model_term_tests(x), error = function(e) list())
  } else if (!is_mixed && inherits(x, "glm") && !inherits(x, "gam")) {
    tryCatch(ms_model_term_tests(x), error = function(e) list())
  } else {
    list()
  }
  if (length(model_term_tests) > 0L) fields$model_term_tests <- model_term_tests

  term_roles <- tryCatch({
    roles <- ms_lm_term_roles(x, focal = focal, controls = controls)
    if (is.null(roles$model_kind)) {
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
        focal = focal,
        controls = controls,
        infer_covariate_roles = TRUE
      )
      if (!is.null(inferred$model_kind)) {
        roles$terms <- inferred$terms
        roles$focal_terms <- vapply(
          Filter(function(t) identical(t$role, "focal"), inferred$terms),
          function(t) t$name,
          character(1)
        )
        roles$control_terms <- vapply(
          Filter(function(t) identical(t$role, "control"), inferred$terms),
          function(t) t$name,
          character(1)
        )
        roles$predictor <- ms_model_term_phrase(roles$focal_terms)
        roles$model_kind <- inferred$model_kind
      }
    }
    roles
  }, error = function(e) NULL)
  if (!is.null(term_roles)) {
    if (!is.null(term_roles$outcome)) fields$outcome <- term_roles$outcome
    if (length(term_roles$terms) > 0L) fields$terms <- term_roles$terms
    if (length(term_roles$focal_terms) > 0L) fields$focal_terms <- term_roles$focal_terms
    if (length(term_roles$control_terms) > 0L) fields$control_terms <- term_roles$control_terms
    if (!is.null(term_roles$predictor)) fields$predictor <- term_roles$predictor
    if (!is.null(term_roles$model_kind)) {
      fields$model_kind <- ms_broom_model_kind(term_roles$model_kind, family = family, mixed = is_mixed)
    }
  }

  glm_response_figures <- if (!is_mixed && inherits(x, "glm") && !inherits(x, "gam")) {
    tryCatch(ms_glm_response_main_effect_figure_data(
      x,
      term_roles = term_roles,
      conf.level = conf.level
    ), error = function(e) NULL)
  } else {
    NULL
  }

  if (is_mixed) {
    groups <- tryCatch(ms_lmer_groups(x), error = function(e) list())
    if (length(groups) > 0L) fields$groups <- groups
    random_terms <- tryCatch(ms_lmer_random_terms(x), error = function(e) character(0))
    if (length(random_terms) > 0L) fields$random_terms <- random_terms
    random_effects <- ms_broom_random_effects(tidy_df)
    if (!length(random_effects)) {
      random_effects <- tryCatch(ms_lmer_random_effects(x), error = function(e) list())
    }
    if (length(random_effects) > 0L) fields$random_effects <- random_effects
  }

  prov <- ms_provenance_basic()
  data_prov <- tryCatch(ms_data_provenance(stats::model.frame(x)), error = function(e) NULL)
  prov <- ms_provenance_add_data(prov, data_prov)

  extras <- c("broom", if (is_mixed) "broom.mixed" else NULL)
  if (inherits(x, c("glmerMod", "lmerMod", "merMod"))) extras <- c(extras, "lme4")
  if (inherits(x, "lme")) extras <- c(extras, "nlme")
  if (inherits(x, "glmmTMB")) extras <- c(extras, "glmmTMB")
  if (inherits(x, "gam")) extras <- c(extras, "mgcv")
  if (inherits(x, "rlm")) extras <- c(extras, "MASS")

  glmer_role_prefers_categorical_means <- inherits(x, "glmerMod") &&
    !is.null(term_roles) &&
    ms_mixed_role_prefers_categorical_means(term_roles)
  adjusted_means <- if (glmer_role_prefers_categorical_means) {
    tryCatch(ms_glmer_probability_means_figure_data(
      x,
      conf.level = conf.level,
      preferred_terms = ms_mixed_focal_factor_terms(term_roles),
      require_preferred = TRUE
    ), error = function(e) NULL)
  } else {
    glm_response_figures$adjusted_means %||% NULL
  }
  interaction_plot <- if (inherits(x, "glmerMod")) {
    plot <- tryCatch(ms_glmer_probability_interaction_plot_data(x, conf.level = conf.level),
                     error = function(e) NULL)
    if (is.null(plot)) {
      plot <- tryCatch(ms_glmer_probability_effect_plot_data(x, conf.level = conf.level),
                       error = function(e) NULL)
    }
    plot
  } else if (!is_mixed && inherits(x, "glm") && !inherits(x, "gam")) {
    plot <- tryCatch(ms_interaction_plot_data(x), error = function(e) NULL)
    if (is.null(plot)) {
      plot <- glm_response_figures$interaction_plot %||% NULL
    }
    plot
  } else {
    NULL
  }
  if ((!is.null(adjusted_means) || !is.null(interaction_plot)) && inherits(x, "glmerMod")) {
    extras <- c(extras, "emmeans")
  }
  figure_data <- list()
  if (!is.null(adjusted_means)) figure_data$adjusted_means <- adjusted_means
  if (!is.null(interaction_plot)) figure_data$interaction_plot <- interaction_plot
  if (length(figure_data) == 0L) figure_data <- NULL

  available_figures <- NULL
  adjusted_mean_kind <- as.character(adjusted_means$mean_kind %||% "")
  if (!is.null(adjusted_means) &&
      adjusted_mean_kind %in% c("predicted_probability", "predicted_count", "predicted_rate",
                                "predicted_mean")) {
    adjusted_label <- switch(
      adjusted_mean_kind,
      predicted_count = "Predicted counts",
      predicted_rate = "Predicted rates",
      predicted_mean = "Predicted means",
      "Predicted probabilities"
    )
    available_figures <- c(available_figures, list(list(
      type = "adjusted_means",
      label = adjusted_label,
      default = TRUE
    )))
  }
  interaction_mean_kind <- as.character(interaction_plot$mean_kind %||% "")
  if (!is.null(interaction_plot) &&
      interaction_mean_kind %in% c("predicted_probability", "predicted_count", "predicted_rate",
                                   "predicted_mean") &&
      identical(as.character(interaction_plot$interaction_kind %||% ""), "continuous_main_effect")) {
    interaction_label <- switch(
      interaction_mean_kind,
      predicted_count = "Predicted count curve",
      predicted_rate = "Predicted rate curve",
      predicted_mean = "Predicted mean curve",
      "Predicted probability curve"
    )
    available_figures <- c(available_figures, list(list(
      type = "interaction_plot",
      label = interaction_label,
      default = TRUE
    )))
  }

  ms_build_envelope(
    type       = type_info$type,
    type_label = type_info$type_label,
    call       = trimws(gsub("\\s+", " ", call_str)),
    fields     = fields,
    raw_output = ms_model_raw_output(x),
    packages   = ms_packages_basic(extras = unique(extras)),
    provenance = prov,
    figure_data = figure_data,
    available_figures = available_figures
  )
}

ms_coefficients_have_ci <- function(coefficients) {
  if (!length(coefficients)) return(FALSE)
  any(vapply(coefficients, function(row) {
    if (is.null(row$ci_lower) || is.null(row$ci_upper)) return(FALSE)
    lower <- ms_safe_numeric(row$ci_lower)
    upper <- ms_safe_numeric(row$ci_upper)
    !is.na(lower) && !is.na(upper)
  }, logical(1)))
}

ms_broom_model_kind <- function(model_kind, family = NULL, mixed = FALSE) {
  model_kind <- as.character(model_kind %||% "")
  if (!model_kind %in% c("ancova", "controlled_regression")) return(model_kind)

  family_name <- tolower(as.character(family$family %||% ""))
  family_link <- tolower(as.character(family$link %||% ""))
  is_gaussian_identity <- family_name %in% c("", "gaussian") &&
    family_link %in% c("", "identity")
  if (isTRUE(mixed)) {
    return(if (is_gaussian_identity) "controlled_mixed_model" else "controlled_glmm")
  }
  if (!is_gaussian_identity) return("controlled_glm")
  "controlled_regression"
}

ms_model_call_string <- function(x, .call = NULL) {
  if (!is.null(.call)) return(.call)
  call_obj <- tryCatch(stats::getCall(x), error = function(e) NULL)
  if (!is.null(call_obj)) return(ms_deparse_call(call_obj))
  user_call <- match.call()$x
  if (!is.null(user_call) && !identical(user_call, as.name("x"))) {
    return(paste(deparse(user_call, width.cutoff = 500L), collapse = " "))
  }
  NA_character_
}

ms_model_raw_output <- function(x) {
  s <- tryCatch(summary(x), error = function(e) NULL)
  if (!is.null(s)) return(ms_capture_output(s))
  ms_capture_output(x)
}

ms_broom_tidy_model <- function(x, mixed, exponentiate, conf.int, conf.level) {
  out <- tryCatch({
    if (isTRUE(mixed)) {
      suppressMessages(suppressWarnings(broom.mixed::tidy(
        x,
        effects = c("fixed", "ran_pars"),
        conf.int = conf.int,
        conf.level = conf.level,
        exponentiate = exponentiate
      )))
    } else {
      suppressMessages(suppressWarnings(broom::tidy(
        x,
        conf.int = conf.int,
        conf.level = conf.level,
        exponentiate = exponentiate
      )))
    }
  }, error = function(e) NULL)

  if (is.null(out) && isTRUE(conf.int)) {
    out <- tryCatch({
      if (isTRUE(mixed)) {
        suppressMessages(suppressWarnings(broom.mixed::tidy(
          x,
          effects = c("fixed", "ran_pars"),
          conf.int = FALSE,
          exponentiate = exponentiate
        )))
      } else {
        suppressMessages(suppressWarnings(broom::tidy(x, conf.int = FALSE, exponentiate = exponentiate)))
      }
    }, error = function(e) NULL)
  }
  # gam etc.: the default tidy() can return a non-coefficient table
  # (smooth terms). Retry for the parametric coefficient table.
  if (!is.null(out) && !isTRUE(mixed) && !"estimate" %in% names(out)) {
    parametric <- tryCatch(
      suppressMessages(suppressWarnings(broom::tidy(x, parametric = TRUE))),
      error = function(e) NULL
    )
    if (!is.null(parametric) && "estimate" %in% names(parametric)) {
      out <- parametric
    }
  }
  if (is.null(out)) {
    stop("broom could not tidy this model object.", call. = FALSE)
  }
  as.data.frame(out)
}

ms_broom_glance_model <- function(x, mixed) {
  out <- tryCatch({
    if (isTRUE(mixed)) broom.mixed::glance(x) else broom::glance(x)
  }, error = function(e) NULL)
  if (is.null(out)) return(data.frame())
  as.data.frame(out)
}

ms_broom_fixed_rows <- function(tidy_df) {
  if (is.null(tidy_df) || !nrow(tidy_df)) return(NULL)
  if ("effect" %in% names(tidy_df)) {
    tidy_df <- tidy_df[is.na(tidy_df$effect) | tidy_df$effect %in% c("", "fixed"), , drop = FALSE]
  }
  if (!nrow(tidy_df) || !"term" %in% names(tidy_df) || !"estimate" %in% names(tidy_df)) {
    return(NULL)
  }
  tidy_df
}

ms_broom_type_info <- function(x, family = NULL, mixed = FALSE) {
  fam <- if (!is.null(family)) as.character(family$family %||% "") else ""
  link <- if (!is.null(family)) as.character(family$link %||% "") else ""
  fam_txt <- if (nzchar(fam) || nzchar(link)) {
    paste0(if (nzchar(fam)) fam else "unknown", "; ", if (nzchar(link)) link else "unknown")
  } else {
    ""
  }

  if (inherits(x, "glmerMod")) {
    return(list(
      type = "glmer_model_summary",
      type_label = paste0("Generalized mixed model (", fam_txt, ")")
    ))
  }

  if (isTRUE(mixed)) {
    return(list(
      type = "mixed_model_summary",
      type_label = if (nzchar(fam_txt)) {
        paste0("Generalized mixed model (", fam_txt, ")")
      } else {
        "Linear mixed model"
      }
    ))
  }

  if (inherits(x, "gam")) {
    return(list(
      type = "generalized_additive_model",
      type_label = if (nzchar(fam_txt)) {
        paste0("Generalized additive model (", fam_txt, ")")
      } else {
        "Generalized additive model"
      }
    ))
  }

  if (inherits(x, "glm")) {
    is_logistic <- identical(tolower(fam), "binomial") && identical(tolower(link), "logit")
    return(list(
      type = if (is_logistic) "logistic_regression" else "generalized_linear_model",
      type_label = if (is_logistic) {
        paste0("Logistic regression (", fam_txt, ")")
      } else if (nzchar(fam_txt)) {
        paste0("Generalized linear model (", fam_txt, ")")
      } else {
        "Generalized linear model"
      }
    ))
  }

  list(
    type = "broom_model_summary",
    type_label = paste0("Model summary (", class(x)[[1]], ")")
  )
}

ms_broom_statistic_label <- function(x, mixed = FALSE) {
  cn <- tryCatch(colnames(summary(x)$coefficients), error = function(e) character(0))
  if (any(grepl("^z value$", cn, ignore.case = TRUE))) return("z")
  if (any(grepl("^t value$", cn, ignore.case = TRUE))) return("t")
  if (any(grepl("^Wald", cn, ignore.case = TRUE))) return("Wald")
  if (isTRUE(mixed) && inherits(x, "glmerMod")) return("z")
  NULL
}

ms_broom_coefficients <- function(x, fixed_df, statistic_label = NULL,
                                  estimate_name = NULL) {
  assign_map <- ms_broom_coef_assign_map(x)
  coef_label_map <- ms_coefficient_label_map(x)
  rows <- lapply(seq_len(nrow(fixed_df)), function(i) {
    term <- as.character(fixed_df$term[i])
    row <- list(
      term = term,
      estimate = ms_safe_numeric(fixed_df$estimate[i])
    )
    if (!is.null(estimate_name)) row$estimate_name <- estimate_name
    if ("std.error" %in% names(fixed_df)) {
      row$std_error <- ms_safe_numeric(fixed_df$std.error[i])
    }
    if ("statistic" %in% names(fixed_df)) {
      row$statistic <- ms_safe_numeric(fixed_df$statistic[i])
      if (!is.null(statistic_label)) row$statistic_label <- statistic_label
    }
    if ("p.value" %in% names(fixed_df)) {
      row$p_value <- ms_safe_numeric(fixed_df$p.value[i])
    }
    if ("conf.low" %in% names(fixed_df)) {
      row$ci_lower <- ms_safe_numeric(fixed_df$conf.low[i])
    }
    if ("conf.high" %in% names(fixed_df)) {
      row$ci_upper <- ms_safe_numeric(fixed_df$conf.high[i])
    }
    if (!is.null(assign_map) && term %in% names(assign_map)) {
      row$term_source <- assign_map[[term]]
    }
    row <- ms_apply_coefficient_label(row, coef_label_map)
    row
  })
  rows
}

ms_broom_coef_assign_map <- function(x) {
  mm <- tryCatch(stats::model.matrix(x), error = function(e) NULL)
  if (is.null(mm)) return(NULL)
  assign <- attr(mm, "assign")
  term_labels <- attr(stats::terms(x), "term.labels") %||% character(0)
  cols <- colnames(mm)
  if (is.null(assign) || is.null(cols)) return(NULL)
  out <- list()
  for (i in seq_along(cols)) {
    idx <- assign[[i]]
    if (idx > 0L && idx <= length(term_labels)) out[[cols[[i]]]] <- term_labels[[idx]]
  }
  out
}

ms_broom_glance_fields <- function(glance_df) {
  fields <- list()
  if (is.null(glance_df) || !nrow(glance_df)) return(fields)
  row <- glance_df[1, , drop = FALSE]

  pick <- function(name) {
    if (name %in% names(row)) ms_safe_numeric(row[[name]][[1]]) else NA_real_
  }

  aic <- pick("AIC")
  if (!is.na(aic)) fields$aic <- aic
  bic <- pick("BIC")
  if (!is.na(bic)) fields$bic <- bic
  loglik <- pick("logLik")
  if (!is.na(loglik)) fields$logLik <- loglik
  nobs <- pick("nobs")
  if (!is.na(nobs)) fields$n <- as.integer(nobs)
  null_dev <- pick("null.deviance")
  if (!is.na(null_dev)) fields$null_deviance <- null_dev
  null_df <- pick("df.null")
  if (!is.na(null_df)) fields$null_df <- null_df
  dev <- pick("deviance")
  if (!is.na(dev)) fields$residual_deviance <- dev
  resid_df <- pick("df.residual")
  if (!is.na(resid_df)) fields$residual_df <- resid_df
  sigma <- pick("sigma")
  if (!is.na(sigma)) fields$sigma <- sigma
  r_sq <- pick("r.squared")
  if (!is.na(r_sq)) fields$r_squared <- r_sq
  adj_r_sq <- pick("adj.r.squared")
  if (!is.na(adj_r_sq)) fields$adj_r_squared <- adj_r_sq

  fields
}

ms_broom_random_effects <- function(tidy_df) {
  if (is.null(tidy_df) || !"effect" %in% names(tidy_df)) return(list())
  ran <- tidy_df[tidy_df$effect %in% "ran_pars", , drop = FALSE]
  if (!nrow(ran)) return(list())

  rows <- list()
  for (i in seq_len(nrow(ran))) {
    term <- as.character(ran$term[i])
    group <- if ("group" %in% names(ran)) as.character(ran$group[i]) else NA_character_
    estimate <- ms_safe_numeric(ran$estimate[i])
    if (is.na(estimate)) next
    if (grepl("^sd_+", term)) {
      name <- sub("^sd_+", "", term)
      rows[[length(rows) + 1L]] <- list(
        group = group,
        name = name,
        variance = ms_safe_numeric(estimate^2),
        std_dev = estimate
      )
    } else if (grepl("^cor_+", term)) {
      name <- sub("^cor_+", "", term)
      rows[[length(rows) + 1L]] <- list(
        group = group,
        name = name,
        corr = estimate
      )
    }
  }
  rows
}
