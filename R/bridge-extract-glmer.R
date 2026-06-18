# R bridge — lme4 glmerMod extractor.
#
# Thin S3 routing shim. The actual extraction is broom.mixed-backed so
# glmer stays aligned with glm and other model families instead of
# drifting as a hand-built extractor.

#' @rdname mellio_payload
#' @export
mellio_payload.glmerMod <- function(x, ..., .call = NULL,
                                exponentiate = FALSE,
                                conf.int = TRUE,
                                conf.level = 0.95) {
  rlang::check_installed("lme4", reason = "to extract glmer model summaries")
  mellio_payload_broom_model(
    x,
    ...,
    .call = .call,
    exponentiate = exponentiate,
    conf.int = conf.int,
    conf.level = conf.level,
    mixed = TRUE
  )
}

ms_glmer_model_term_tests <- function(x) {
  if (!inherits(x, "glmerMod")) return(list())

  tf <- tryCatch(stats::terms(x), error = function(e) NULL)
  if (is.null(tf)) return(list())

  term_labels <- attr(tf, "term.labels") %||% character(0)
  orders <- attr(tf, "order") %||% integer(0)
  if (!length(term_labels)) return(list())
  if (length(orders) != length(term_labels)) {
    orders <- rep(1L, length(term_labels))
  }

  drop_rows <- ms_model_term_tests_drop1_chisq(x, term_labels)
  out <- vector("list", length(term_labels))
  for (i in seq_along(term_labels)) {
    term <- term_labels[[i]]
    row <- drop_rows[[term]]
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

ms_glmer_probability_means_figure_data <- function(x, conf.level = 0.95,
                                                   max_levels = 12L,
                                                   preferred_terms = NULL,
                                                   require_preferred = FALSE) {
  if (!requireNamespace("emmeans", quietly = TRUE)) return(NULL)

  family <- tryCatch(stats::family(x), error = function(e) NULL)
  family_name <- tolower(as.character(family$family %||% ""))
  family_link <- as.character(family$link %||% "")
  if (!family_name %in% c("binomial", "quasibinomial")) return(NULL)

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

  categorical_idx <- which(vapply(term_info, function(info) identical(info$type, "categorical"), logical(1)))
  if (!length(categorical_idx)) return(NULL)
  selected_idx <- ms_mixed_preferred_term_index(term_info, categorical_idx, preferred_terms)
  if (is.null(selected_idx)) {
    if (isTRUE(require_preferred)) return(NULL)
    selected_idx <- categorical_idx[[1L]]
  }
  factor_info <- term_info[[selected_idx]]
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
  emm_args <- list(object = x, specs = emm_formula, type = "response")
  if (length(at_values) > 0L) emm_args$at <- at_values
  emm <- tryCatch(
    do.call(emmeans::emmeans, emm_args),
    error = function(e) NULL
  )
  if (is.null(emm)) return(NULL)
  emm_df <- tryCatch(
    as.data.frame(summary(emm, infer = c(TRUE, FALSE), level = conf.level, type = "response")),
    error = function(e) NULL
  )
  if (is.null(emm_df) || nrow(emm_df) < 2L) return(NULL)

  factor_col <- ms_lmer_emmeans_col(emm_df, c(factor_info$variable, factor_info$term))
  estimate_col <- ms_lmer_emmeans_col(emm_df, c("prob", "response", "emmean", "rate", "estimate"))
  se_col <- ms_lmer_emmeans_col(emm_df, c("SE", "std.error", "std_error"))
  df_col <- ms_lmer_emmeans_col(emm_df, c("df"))
  lower_col <- ms_lmer_emmeans_col(emm_df, c("asymp.LCL", "lower.CL", "lower.HPD", "LCL"))
  upper_col <- ms_lmer_emmeans_col(emm_df, c("asymp.UCL", "upper.CL", "upper.HPD", "UCL"))
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
    mean_kind = "predicted_probability",
    source = "glmer_emmeans",
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
    covariates = lapply(covariates, function(covariate) {
      list(
        variable = covariate$variable,
        term = covariate$term,
        label = covariate$label,
        value = covariate$value,
        value_label = ms_interaction_format_number(covariate$value),
        rule = "sample_mean"
      )
    }),
    adjustment = list(rule = "sample_mean", label = "sample means"),
    outcome = outcome %||% NULL,
    y_label = "Predicted probability",
    ci_level = conf.level,
    ci_method = "emmeans",
    model_family = family_name,
    model_link = family_link,
    bounded_response = TRUE
  )
  subject <- ms_lmer_primary_group_summary(x, mf)
  if (!is.null(subject)) out$subject <- subject
  out
}

ms_glmer_probability_interaction_plot_data <- function(x, conf.level = 0.95,
                                                       max_levels = 8L,
                                                       grid_points = 80L) {
  if (!requireNamespace("emmeans", quietly = TRUE)) return(NULL)

  family <- tryCatch(stats::family(x), error = function(e) NULL)
  family_name <- tolower(as.character(family$family %||% ""))
  family_link <- as.character(family$link %||% "")
  if (!family_name %in% c("binomial", "quasibinomial")) return(NULL)

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
  if (sum(is_numeric) == 0L && sum(is_categorical) == 2L) {
    return(ms_glmer_probability_categorical_interaction_plot_data(
      x = x,
      mf = mf,
      term_labels = term_labels,
      interaction_term = interaction_term,
      interaction_variables = interaction_variables,
      components = components,
      var_map = var_map,
      conf.level = conf.level,
      max_levels = max_levels,
      family_name = family_name,
      family_link = family_link
    ))
  }
  if (sum(is_numeric) != 1L || sum(is_categorical) != 1L) return(NULL)

  x_info <- components[[which(is_numeric)[[1L]]]]
  moderator_info <- components[[which(is_categorical)[[1L]]]]
  moderator_set <- ms_interaction_categorical_set(moderator_info, max_levels = max_levels)
  if (is.null(moderator_set)) return(NULL)

  x_values <- ms_interaction_numeric_values(x_info$values)
  if (length(x_values) < 2L) return(NULL)
  x_range <- range(x_values, na.rm = TRUE)
  if (!all(is.finite(x_range)) || x_range[[1L]] == x_range[[2L]]) return(NULL)
  x_grid <- seq(x_range[[1L]], x_range[[2L]], length.out = max(12L, as.integer(grid_points)))
  adjustments <- ms_glmer_probability_interaction_adjustments(
    term_labels = term_labels,
    var_map = var_map,
    exclude_variables = interaction_variables,
    max_levels = max_levels
  )

  emm_formula <- stats::as.formula(paste("~", paste(c(x_info$term, moderator_info$term), collapse = " * ")))
  at_values <- c(setNames(list(x_grid), x_info$variable), adjustments$at_values)
  emm <- tryCatch(
    emmeans::emmeans(x, specs = emm_formula, at = at_values, type = "response"),
    error = function(e) NULL
  )
  if (is.null(emm)) return(NULL)
  emm_df <- tryCatch(
    as.data.frame(summary(emm, infer = c(TRUE, FALSE), level = conf.level, type = "response")),
    error = function(e) NULL
  )
  if (is.null(emm_df) || !nrow(emm_df)) return(NULL)

  x_col <- ms_lmer_emmeans_col(emm_df, c(x_info$variable, x_info$term))
  moderator_col <- ms_lmer_emmeans_col(emm_df, c(moderator_info$variable, moderator_info$term))
  estimate_col <- ms_lmer_emmeans_col(emm_df, c("prob", "response", "emmean", "rate", "estimate"))
  se_col <- ms_lmer_emmeans_col(emm_df, c("SE", "std.error", "std_error"))
  df_col <- ms_lmer_emmeans_col(emm_df, c("df"))
  lower_col <- ms_lmer_emmeans_col(emm_df, c("asymp.LCL", "lower.CL", "lower.HPD", "LCL"))
  upper_col <- ms_lmer_emmeans_col(emm_df, c("asymp.UCL", "upper.CL", "upper.HPD", "UCL"))
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
    source = "glmer_emmeans",
    mean_kind = "predicted_probability",
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
    held_constant = adjustments$held_constant,
    marginalized_terms = adjustments$marginalized_terms,
    outcome = outcome %||% NULL,
    y_label = "Predicted probability",
    scale = "response",
    ci_level = conf.level,
    ci_method = "emmeans",
    model_family = family_name,
    model_link = family_link,
    bounded_response = TRUE
  )
  if (!is.null(subject)) out$subject <- subject
  out
}

ms_glmer_probability_categorical_interaction_plot_data <- function(x, mf, term_labels,
                                                                   interaction_term,
                                                                   interaction_variables,
                                                                   components,
                                                                   var_map,
                                                                   conf.level = 0.95,
                                                                   max_levels = 8L,
                                                                   family_name = "",
                                                                   family_link = "") {
  ordered <- ms_lmer_interaction_role_order(x, mf, components)
  x_info <- ordered$x
  moderator_info <- ordered$moderator
  x_set <- ms_interaction_categorical_set(x_info, max_levels = max_levels)
  moderator_set <- ms_interaction_categorical_set(moderator_info, max_levels = max_levels)
  if (is.null(x_set) || is.null(moderator_set)) return(NULL)

  adjustments <- ms_glmer_probability_interaction_adjustments(
    term_labels = term_labels,
    var_map = var_map,
    exclude_variables = interaction_variables,
    max_levels = max_levels
  )
  emm_formula <- stats::as.formula(paste("~", paste(c(x_info$term, moderator_info$term), collapse = " * ")))
  emm_args <- list(object = x, specs = emm_formula, type = "response")
  if (length(adjustments$at_values) > 0L) emm_args$at <- adjustments$at_values
  emm <- tryCatch(
    do.call(emmeans::emmeans, emm_args),
    error = function(e) NULL
  )
  if (is.null(emm)) return(NULL)
  emm_df <- tryCatch(
    as.data.frame(summary(emm, infer = c(TRUE, FALSE), level = conf.level, type = "response")),
    error = function(e) NULL
  )
  if (is.null(emm_df) || !nrow(emm_df)) return(NULL)

  x_col <- ms_lmer_emmeans_col(emm_df, c(x_info$variable, x_info$term))
  moderator_col <- ms_lmer_emmeans_col(emm_df, c(moderator_info$variable, moderator_info$term))
  estimate_col <- ms_lmer_emmeans_col(emm_df, c("prob", "response", "emmean", "rate", "estimate"))
  se_col <- ms_lmer_emmeans_col(emm_df, c("SE", "std.error", "std_error"))
  df_col <- ms_lmer_emmeans_col(emm_df, c("df"))
  lower_col <- ms_lmer_emmeans_col(emm_df, c("asymp.LCL", "lower.CL", "lower.HPD", "LCL"))
  upper_col <- ms_lmer_emmeans_col(emm_df, c("asymp.UCL", "upper.CL", "upper.HPD", "UCL"))
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
    source = "glmer_emmeans",
    mean_kind = "predicted_probability",
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
    held_constant = adjustments$held_constant,
    marginalized_terms = adjustments$marginalized_terms,
    outcome = outcome %||% NULL,
    y_label = "Predicted probability",
    scale = "response",
    ci_level = conf.level,
    ci_method = "emmeans",
    model_family = family_name,
    model_link = family_link,
    bounded_response = TRUE
  )
  if (!is.null(subject)) out$subject <- subject
  out
}

ms_glmer_probability_interaction_adjustments <- function(term_labels, var_map,
                                                         exclude_variables,
                                                         max_levels = 8L) {
  term_labels <- as.character(term_labels %||% character(0))
  exclude_variables <- unique(as.character(exclude_variables %||% character(0)))
  other_info <- lapply(term_labels, function(term) {
    if (grepl(":", term, fixed = TRUE)) return(NULL)
    if (!ms_interaction_component_supported(term)) return(NULL)
    variable <- ms_interaction_component_variable(term)
    if (is.null(variable) || variable %in% exclude_variables || is.null(var_map[[variable]])) {
      return(NULL)
    }
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
  other_info <- Filter(Negate(is.null), other_info)
  numeric_info <- other_info[vapply(other_info, function(info) identical(info$type, "numeric"), logical(1))]
  covariates <- lapply(numeric_info, function(info) {
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

  categorical_info <- other_info[vapply(other_info, function(info) identical(info$type, "categorical"), logical(1))]
  list(
    at_values = setNames(
      lapply(covariates, function(covariate) covariate$value),
      vapply(covariates, function(covariate) covariate$variable, character(1))
    ),
    held_constant = lapply(covariates, function(covariate) {
      list(
        variable = covariate$variable,
        label = covariate$label,
        value = as.character(covariate$value),
        value_label = covariate$value_label,
        type = covariate$type
      )
    }),
    marginalized_terms = lapply(categorical_info, function(info) {
      list(
        variable = info$variable,
        term = info$term,
        label = info$label,
        rule = "estimated_marginal"
      )
    })
  )
}

ms_glmer_probability_effect_plot_data <- function(x, conf.level = 0.95,
                                                  grid_points = 80L,
                                                  max_levels = 8L) {
  if (!requireNamespace("emmeans", quietly = TRUE)) return(NULL)

  family <- tryCatch(stats::family(x), error = function(e) NULL)
  family_name <- tolower(as.character(family$family %||% ""))
  family_link <- as.character(family$link %||% "")
  if (!family_name %in% c("binomial", "quasibinomial")) return(NULL)

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
    emmeans::emmeans(x, specs = emm_formula, at = at_values, type = "response"),
    error = function(e) NULL
  )
  if (is.null(emm)) return(NULL)
  emm_df <- tryCatch(
    as.data.frame(summary(emm, infer = c(TRUE, FALSE), level = conf.level, type = "response")),
    error = function(e) NULL
  )
  if (is.null(emm_df) || nrow(emm_df) < 2L) return(NULL)

  x_col <- ms_lmer_emmeans_col(emm_df, c(x_info$variable, x_info$term))
  estimate_col <- ms_lmer_emmeans_col(emm_df, c("prob", "response", "emmean", "rate", "estimate"))
  se_col <- ms_lmer_emmeans_col(emm_df, c("SE", "std.error", "std_error"))
  df_col <- ms_lmer_emmeans_col(emm_df, c("df"))
  lower_col <- ms_lmer_emmeans_col(emm_df, c("asymp.LCL", "lower.CL", "lower.HPD", "LCL"))
  upper_col <- ms_lmer_emmeans_col(emm_df, c("asymp.UCL", "upper.CL", "upper.HPD", "UCL"))
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
    source = "glmer_emmeans",
    mean_kind = "predicted_probability",
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
    y_label = "Predicted probability",
    scale = "response",
    ci_level = conf.level,
    ci_method = "emmeans",
    model_family = family_name,
    model_link = family_link,
    bounded_response = TRUE
  )
  if (!is.null(subject)) out$subject <- subject
  out
}
