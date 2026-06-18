# R bridge -- nnet::multinom models.
#
# Multinomial logistic regression has one coefficient equation per
# non-reference outcome category. Represent it as a table so the web app
# can preserve the comparison/term structure instead of flattening the
# model into ordinary binary-logistic prose.

#' @rdname mellio_payload
#' @export
mellio_payload.multinom <- function(x, ..., .call = NULL) {
  s <- tryCatch(summary(x), error = function(e) NULL)
  if (is.null(s) || is.null(s$coefficients) || is.null(s$standard.errors)) {
    stop("Could not extract a multinomial logistic regression summary.",
         call. = FALSE)
  }

  coef_mat <- ms_multinom_matrix(s$coefficients, x)
  se_mat <- ms_multinom_matrix(s$standard.errors, x)
  conf.level <- 0.95
  rows <- ms_multinom_rows(coef_mat, se_mat, conf.level = conf.level)
  if (!length(rows)) {
    stop("Could not extract multinomial coefficient rows.", call. = FALSE)
  }

  columns <- list(
    list(key = "comparison", label = "Comparison", format = "text"),
    list(key = "term", label = "Term", format = "text"),
    list(key = "estimate", label = "B", format = "number"),
    list(key = "std_error", label = "SE", format = "number"),
    list(key = "statistic", label = "z", format = "statistic"),
    list(key = "p_value", label = "p", format = "pvalue"),
    list(key = "ci", label = "95% CI", format = "ci")
  )

  fields <- list(
    table_type = "multinomial_coefficients",
    columns = columns,
    rows = rows,
    source = "R",
    statistic_label = "z",
    model_family = "multinomial logistic",
    model_link = "logit",
    coefficient_scale = "multinomial_logit",
    coefficient_ci_method = "wald",
    coefficient_ci_scale = "link",
    coefficient_p_value_method = "wald_z",
    conf_level = conf.level,
    p_value_note = "Approximate two-sided p-values computed from z = B / SE."
  )

  aic <- ms_safe_numeric(stats::AIC(x))
  if (!is.na(aic)) fields$aic <- aic
  bic <- tryCatch(ms_safe_numeric(stats::BIC(x)), error = function(e) NA_real_)
  if (!is.na(bic)) fields$bic <- bic
  loglik <- tryCatch(ms_safe_numeric(as.numeric(stats::logLik(x))), error = function(e) NA_real_)
  if (!is.na(loglik)) fields$logLik <- loglik
  dev <- ms_safe_numeric(x$deviance)
  if (!is.na(dev)) fields$residual_deviance <- dev
  n <- ms_multinom_nobs(x)
  if (!is.na(n)) fields$n <- as.integer(n)
  outcome_levels <- tryCatch(as.character(x$lev %||% character(0)), error = function(e) character(0))
  if (length(outcome_levels)) {
    fields$outcome_levels <- as.list(outcome_levels)
    fields$outcome_level_count <- length(outcome_levels)
    fields$reference_level <- outcome_levels[[1L]]
  }

  f <- tryCatch(stats::formula(x), error = function(e) NULL)
  if (!is.null(f) && length(f) >= 3L) {
    fields$outcome <- ms_model_clean_term(paste(deparse(f[[2]], width.cutoff = 500L), collapse = " "))
  }
  model_terms <- tryCatch(attr(stats::terms(x), "term.labels"), error = function(e) character(0))
  if (length(model_terms) > 0L) {
    fields$terms <- lapply(model_terms, function(term) {
      list(name = term, label = ms_model_clean_term(term), role = "focal", type = "predictor")
    })
    fields$focal_terms <- model_terms
    fields$predictor <- ms_model_term_phrase(vapply(fields$terms, function(t) {
      t$name %||% t$label
    }, character(1)))
  }
  comparisons <- unique(vapply(rows, function(row) row$comparison %||% "", character(1)))
  comparisons <- comparisons[nzchar(comparisons)]
  if (length(comparisons)) fields$comparisons <- as.list(comparisons)

  ms_build_envelope(
    type = "multinomial_logistic_regression",
    type_label = "Multinomial logistic regression",
    call = trimws(gsub("\\s+", " ", ms_model_call_string(x, .call = .call))),
    fields = fields,
    raw_output = ms_capture_output(s),
    packages = ms_packages_basic(extras = "nnet"),
    card_kind = "table"
  )
}

ms_multinom_matrix <- function(value, x) {
  if (!is.null(dim(value))) return(as.matrix(value))
  terms <- names(value) %||% paste0("term_", seq_along(value))
  comparisons <- ms_multinom_comparisons(x, 1L)
  out <- matrix(as.numeric(value), nrow = 1L,
                dimnames = list(comparisons[[1L]], terms))
  out
}

ms_multinom_comparisons <- function(x, n_rows) {
  lev <- tryCatch(x$lev, error = function(e) character(0))
  if (length(lev) >= n_rows + 1L) return(lev[seq_len(n_rows) + 1L])
  paste0("comparison_", seq_len(n_rows))
}

ms_multinom_rows <- function(coef_mat, se_mat, conf.level = 0.95) {
  comparisons <- rownames(coef_mat)
  terms <- colnames(coef_mat)
  if (is.null(comparisons)) comparisons <- paste0("comparison_", seq_len(nrow(coef_mat)))
  if (is.null(terms)) terms <- paste0("term_", seq_len(ncol(coef_mat)))

  rows <- list()
  ci_crit <- stats::qnorm((1 + conf.level) / 2)
  for (i in seq_len(nrow(coef_mat))) {
    for (j in seq_len(ncol(coef_mat))) {
      estimate <- ms_safe_numeric(coef_mat[i, j])
      se <- if (i <= nrow(se_mat) && j <= ncol(se_mat)) ms_safe_numeric(se_mat[i, j]) else NA_real_
      if (is.na(estimate) || is.na(se)) next
      z <- if (!is.na(se) && se != 0) ms_safe_numeric(estimate / se) else NA_real_
      row <- list(
        comparison = comparisons[[i]],
        term = terms[[j]],
        estimate_name = "B",
        estimate = estimate,
        std_error = se
      )
      if (!is.na(z)) {
        row$statistic <- z
        row$statistic_label <- "z"
        row$p_value <- ms_safe_numeric(2 * stats::pnorm(abs(z), lower.tail = FALSE))
      }
      if (is.finite(ci_crit)) {
        row$ci_lower <- ms_safe_numeric(estimate - ci_crit * se)
        row$ci_upper <- ms_safe_numeric(estimate + ci_crit * se)
        row$ci_method <- "wald"
      }
      rows[[length(rows) + 1L]] <- row
    }
  }
  rows
}

ms_multinom_nobs <- function(x) {
  fv <- tryCatch(x$fitted.values, error = function(e) NULL)
  if (is.matrix(fv) || is.data.frame(fv)) return(ms_safe_numeric(nrow(fv)))
  if (!is.null(fv) && length(fv)) return(ms_safe_numeric(length(fv)))
  NA_real_
}
