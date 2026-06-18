# R bridge — hierarchical / nested model comparison cards.
#
# anova(m1, m2) only carries the comparison table, so it cannot recover
# model-level R2 or the terms added at each step. mellio_compare() keeps the
# original model objects and emits a richer Result Card for hierarchical
# regression reporting.

#' Compare nested regression models as a Mellio Stats card
#'
#' Builds a hierarchical-regression comparison payload from two or more
#' nested model objects. Unlike `mellio_payload(anova(m1, m2))`, this keeps the
#' original models available so Mellio can report model `R^2`, adjusted `R^2`,
#' `Delta R^2`, and the nested-model F-change.
#'
#' @param ... Nested model objects, usually `lm` objects ordered from the
#'   baseline model to the final model.
#' @param labels Optional model labels. Defaults to `"Model 1"`, `"Model 2"`,
#'   and so on.
#' @param .call Optional call string for provenance / display.
#' @return A `mellio_payload` object.
#' @export
#' @family R bridge
#' @examples
#' m1 <- lm(mpg ~ wt + hp, data = mtcars)
#' m2 <- lm(mpg ~ wt + hp + cyl, data = mtcars)
#' p <- mellio_compare(m1, m2)
#' p$type
mellio_compare <- function(..., labels = NULL, .call = NULL) {
  rlang::check_installed("broom", reason = "to compare model summaries")

  dots <- list(...)
  if (length(dots) == 1L && is.list(dots[[1]]) &&
      !inherits(dots[[1]], c("lm", "glm"))) {
    dots <- dots[[1]]
  }
  models <- dots
  n_models <- length(models)
  if (n_models < 2L) {
    stop("mellio_compare() requires at least two nested model objects.", call. = FALSE)
  }

  if (is.null(.call)) {
    .call <- paste(deparse(match.call(), width.cutoff = 500L), collapse = " ")
  }

  if (!is.null(labels) && length(labels) != n_models) {
    stop("labels must have one value per model.", call. = FALSE)
  }
  if (is.null(labels)) {
    model_names <- names(models)
    labels <- if (!is.null(model_names) && all(nzchar(model_names))) {
      model_names
    } else {
      paste("Model", seq_len(n_models))
    }
  }

  ms_payload_model_comparison(models, labels = labels, .call = .call)
}

ms_payload_model_comparison <- function(models, labels, .call = NULL) {
  n_models <- length(models)
  glances <- lapply(seq_len(n_models), function(i) {
    g <- tryCatch(broom::glance(models[[i]]), error = function(e) NULL)
    if (is.null(g)) {
      stop("broom::glance() failed for model ", i, ".", call. = FALSE)
    }
    as.data.frame(g)
  })

  model_rows <- lapply(seq_len(n_models), function(i) {
    ms_compare_model_row(models[[i]], glances[[i]], labels[[i]], i)
  })
  comparisons <- ms_compare_rows(models, glances, labels)
  final_comparison <- if (length(comparisons)) comparisons[[length(comparisons)]] else NULL
  final_model <- model_rows[[length(model_rows)]]
  first_model <- model_rows[[1]]

  outcome <- ms_compare_common_outcome(models)
  all_terms <- unique(unlist(lapply(models, ms_compare_term_labels), use.names = FALSE))

  fields <- list(
    source = "R",
    model_kind = "hierarchical_regression",
    model_count = n_models,
    outcome = outcome,
    terms = I(all_terms),
    models = model_rows,
    comparisons = comparisons,
    final_model = final_model,
    table_note = "Predictors list the cumulative terms in each model. Terms entered list the newly added terms at each step. Change statistics compare each model with the previous nested model."
  )

  if (!is.null(final_comparison)) {
    fields$final_comparison <- final_comparison
    if (!is.null(final_comparison$r_squared_change)) {
      fields$r_squared_change <- final_comparison$r_squared_change
      fields$estimate <- list(
        name = "\u0394R\u00B2",
        value = final_comparison$r_squared_change
      )
    }
    if (!is.null(final_comparison$f_change)) {
      fields$statistic <- list(
        name = "F",
        value = final_comparison$f_change,
        df = I(c(final_comparison$df_change, final_comparison$residual_df))
      )
      fields$f_change <- final_comparison$f_change
    }
    if (!is.null(final_comparison$p_change)) {
      fields$p_value <- final_comparison$p_change
      fields$p_change <- final_comparison$p_change
    }
    if (!is.null(final_comparison$df_change)) fields$df_change <- final_comparison$df_change
    if (!is.null(final_comparison$residual_df)) fields$residual_df <- final_comparison$residual_df
  }

  if (!is.null(final_model$r_squared)) fields$r_squared <- final_model$r_squared
  if (!is.null(final_model$adj_r_squared)) fields$adj_r_squared <- final_model$adj_r_squared
  if (!is.null(final_model$n)) fields$n <- final_model$n
  if (!is.null(first_model$r_squared)) fields$baseline_r_squared <- first_model$r_squared

  prov <- ms_provenance_basic()
  data_prov <- tryCatch(ms_data_provenance(stats::model.frame(models[[n_models]])),
                        error = function(e) NULL)
  prov <- ms_provenance_add_data(prov, data_prov)

  ms_build_envelope(
    type = "hierarchical_regression_comparison",
    type_label = paste0("Hierarchical regression comparison (", n_models, " models)"),
    call = trimws(gsub("\\s+", " ", .call %||% NA_character_)),
    fields = fields,
    raw_output = ms_compare_raw_output(models),
    packages = ms_packages_basic(extras = "broom"),
    provenance = prov
  )
}

ms_compare_model_row <- function(model, glance_df, label, index) {
  pick <- function(name) {
    if (name %in% names(glance_df)) ms_safe_numeric(glance_df[[name]][[1]]) else NA_real_
  }
  formula_txt <- tryCatch(
    paste(deparse(stats::formula(model), width.cutoff = 500L), collapse = " "),
    error = function(e) NA_character_
  )
  terms <- ms_compare_term_labels(model)

  row <- list(
    model = index,
    label = as.character(label),
    formula = formula_txt,
    terms = I(terms)
  )

  nobs <- pick("nobs")
  if (!is.na(nobs)) row$n <- as.integer(nobs)
  r2 <- pick("r.squared")
  if (!is.na(r2)) row$r_squared <- r2
  adj_r2 <- pick("adj.r.squared")
  if (!is.na(adj_r2)) row$adj_r_squared <- adj_r2
  sigma <- pick("sigma")
  if (!is.na(sigma)) row$sigma <- sigma
  aic <- pick("AIC")
  if (!is.na(aic)) row$aic <- aic
  bic <- pick("BIC")
  if (!is.na(bic)) row$bic <- bic
  loglik <- pick("logLik")
  if (!is.na(loglik)) row$logLik <- loglik
  rdf <- pick("df.residual")
  if (!is.na(rdf)) row$residual_df <- rdf

  row
}

ms_compare_rows <- function(models, glances, labels) {
  anova_tbl <- tryCatch(do.call(stats::anova, models), error = function(e) NULL)
  n_models <- length(models)
  rows <- list()

  for (i in seq.int(2L, n_models)) {
    prev_terms <- ms_compare_term_labels(models[[i - 1L]])
    current_terms <- ms_compare_term_labels(models[[i]])
    added_terms <- setdiff(current_terms, prev_terms)

    prev_r2 <- ms_compare_glance_value(glances[[i - 1L]], "r.squared")
    current_r2 <- ms_compare_glance_value(glances[[i]], "r.squared")
    r2_change <- if (!is.na(prev_r2) && !is.na(current_r2)) current_r2 - prev_r2 else NA_real_

    f_change <- df_change <- residual_df <- p_change <- sum_of_sq <- residual_sum_sq <- NA_real_
    if (!is.null(anova_tbl) && nrow(anova_tbl) >= i) {
      f_col <- if ("F" %in% colnames(anova_tbl)) "F" else if ("F value" %in% colnames(anova_tbl)) "F value" else NA_character_
      p_col <- if ("Pr(>F)" %in% colnames(anova_tbl)) "Pr(>F)" else NA_character_
      if (!is.na(f_col)) f_change <- ms_safe_numeric(anova_tbl[[f_col]][i])
      if ("Df" %in% colnames(anova_tbl)) df_change <- ms_safe_numeric(anova_tbl[["Df"]][i])
      if ("Res.Df" %in% colnames(anova_tbl)) residual_df <- ms_safe_numeric(anova_tbl[["Res.Df"]][i])
      if (!is.na(p_col)) p_change <- ms_safe_numeric(anova_tbl[[p_col]][i])
      if ("Sum of Sq" %in% colnames(anova_tbl)) sum_of_sq <- ms_safe_numeric(anova_tbl[["Sum of Sq"]][i])
      if ("RSS" %in% colnames(anova_tbl)) residual_sum_sq <- ms_safe_numeric(anova_tbl[["RSS"]][i])
    }

    row <- list(
      from = i - 1L,
      to = i,
      label = as.character(labels[[i]]),
      added_terms = I(added_terms)
    )
    if (!is.na(r2_change)) row$r_squared_change <- r2_change
    if (!is.na(f_change)) row$f_change <- f_change
    if (!is.na(df_change)) row$df_change <- df_change
    if (!is.na(residual_df)) row$residual_df <- residual_df
    if (!is.na(p_change)) row$p_change <- p_change
    if (!is.na(sum_of_sq)) row$sum_of_sq <- sum_of_sq
    if (!is.na(residual_sum_sq)) row$residual_sum_sq <- residual_sum_sq

    rows[[length(rows) + 1L]] <- row
  }

  rows
}

ms_compare_glance_value <- function(glance_df, name) {
  if (is.null(glance_df) || !(name %in% names(glance_df))) return(NA_real_)
  ms_safe_numeric(glance_df[[name]][[1]])
}

ms_compare_term_labels <- function(model) {
  terms <- tryCatch(attr(stats::terms(model), "term.labels"), error = function(e) character(0))
  as.character(terms %||% character(0))
}

ms_compare_common_outcome <- function(models) {
  outcomes <- vapply(models, function(model) {
    f <- tryCatch(stats::formula(model), error = function(e) NULL)
    if (is.null(f) || length(f) < 3L) return(NA_character_)
    paste(deparse(f[[2]], width.cutoff = 500L), collapse = " ")
  }, character(1))
  outcomes <- outcomes[!is.na(outcomes) & nzchar(outcomes)]
  if (length(outcomes) && length(unique(outcomes)) == 1L) outcomes[[1]] else NA_character_
}

ms_compare_raw_output <- function(models) {
  anova_tbl <- tryCatch(do.call(stats::anova, models), error = function(e) NULL)
  lines <- character(0)
  if (!is.null(anova_tbl)) {
    lines <- c(lines, "Nested model comparison:", ms_capture_output(anova_tbl))
  }
  lines <- c(lines, "", "Models:")
  for (i in seq_along(models)) {
    formula_txt <- tryCatch(
      paste(deparse(stats::formula(models[[i]]), width.cutoff = 500L), collapse = " "),
      error = function(e) NA_character_
    )
    lines <- c(lines, paste0("Model ", i, ": ", formula_txt))
  }
  paste(lines, collapse = "\n")
}
