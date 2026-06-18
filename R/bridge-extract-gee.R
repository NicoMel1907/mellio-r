# R bridge -- Generalized estimating equations models.
#
# geepack::geeglm coefficient summaries look GLM-like, but the editorial
# contract is different: clustered observations, a working correlation
# structure, and Wald tests. Keep this as a specialist extractor instead
# of letting the generic broom path flatten it into ordinary regression.

#' @rdname mellio_payload
#' @export
mellio_payload.geeglm <- function(x, ..., .call = NULL) {
  s <- tryCatch(summary(x), error = function(e) NULL)
  if (is.null(s) || is.null(s$coefficients)) {
    stop("Could not extract a GEE model summary.", call. = FALSE)
  }

  raw <- ms_capture_output(s)
  coef_mat <- as.matrix(s$coefficients)
  terms <- rownames(coef_mat)
  if (is.null(terms)) terms <- paste0("term_", seq_len(nrow(coef_mat)))

  col_pick <- function(pattern) {
    hit <- grep(pattern, colnames(coef_mat), ignore.case = TRUE, value = TRUE)
    if (length(hit)) hit[[1L]] else NA_character_
  }
  coef_value <- function(i, col) {
    if (is.na(col) || !col %in% colnames(coef_mat)) return(NA_real_)
    ms_safe_numeric(coef_mat[i, col])
  }

  est_col <- col_pick("^Estimate$")
  se_col <- col_pick("^Std\\.err$|std\\.?\\s*err|robust")
  wald_col <- col_pick("^Wald$")
  p_col <- col_pick("Pr\\(>\\|W\\|\\)")
  conf_level <- 0.95
  ci_crit <- stats::qnorm(1 - (1 - conf_level) / 2)

  coefficients <- lapply(seq_len(nrow(coef_mat)), function(i) {
    term <- terms[[i]]
    estimate <- coef_value(i, est_col)
    std_error <- coef_value(i, se_col)
    row <- list(
      term = term,
      estimate_name = "B",
      estimate = estimate,
      term_source = term
    )
    if (!is.na(se_col)) row$std_error <- std_error
    if (!is.na(estimate) && !is.na(std_error)) {
      row$ci_lower <- ms_safe_numeric(estimate - ci_crit * std_error)
      row$ci_upper <- ms_safe_numeric(estimate + ci_crit * std_error)
      row$ci_method <- "wald_robust"
    }
    if (!is.na(wald_col)) {
      row$statistic <- coef_value(i, wald_col)
      row$statistic_label <- "Wald \u03c7\u00b2"
      row$df <- 1
    }
    if (!is.na(p_col)) row$p_value <- coef_value(i, p_col)
    row
  })

  fields <- list(
    coefficients = coefficients,
    coefficient_ci_method = "wald_robust",
    coefficient_p_value_method = "robust_wald",
    std_error_method = "robust_sandwich",
    conf_level = conf_level,
    source = "R",
    statistic_label = "Wald \u03c7\u00b2",
    model_kind = "gee"
  )

  family <- tryCatch(stats::family(x), error = function(e) NULL)
  if (!is.null(family)) {
    fields$model_family <- as.character(family$family %||% "")
    fields$model_link <- as.character(family$link %||% "")
    fields$coefficient_ci_scale <- if (identical(tolower(fields$model_link), "identity")) {
      "response"
    } else {
      "link"
    }
  }

  call_obj <- tryCatch(stats::getCall(x), error = function(e) NULL)
  id_var <- ms_gee_id_variable(call_obj)
  if (nzchar(id_var)) fields$id_variable <- id_var
  n_clusters <- ms_gee_cluster_count(x)
  if (!is.na(n_clusters)) fields$n_clusters <- as.integer(n_clusters)

  corr <- ms_gee_correlation_structure(x, raw)
  if (nzchar(corr)) fields$correlation_structure <- corr

  scale <- ms_gee_scale_parameter(raw)
  if (!is.na(scale$estimate)) fields$scale_parameter <- scale$estimate
  if (!is.na(scale$std_error)) fields$scale_std_error <- scale$std_error

  working_corr <- ms_gee_working_correlation_parameter(x, raw)
  if (!is.na(working_corr$estimate)) {
    fields$working_correlation_parameter <- working_corr$estimate
    fields$working_correlation_parameter_name <- working_corr$name
  }
  if (!is.na(working_corr$std_error)) {
    fields$working_correlation_std_error <- working_corr$std_error
  }

  n <- tryCatch(stats::nobs(x), error = function(e) NA_real_)
  n <- ms_safe_numeric(n)
  if (!is.na(n)) fields$n <- as.integer(n)

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

  ms_build_envelope(
    type = "gee_model_summary",
    type_label = "Generalized estimating equations model",
    call = trimws(gsub("\\s+", " ", ms_model_call_string(x, .call = .call))),
    fields = fields,
    raw_output = raw,
    packages = ms_packages_basic(extras = "geepack")
  )
}

ms_gee_cluster_count <- function(x) {
  id <- tryCatch(x$id, error = function(e) NULL)
  if (is.null(id)) id <- tryCatch(x$geese$id, error = function(e) NULL)
  if (is.null(id)) return(NA_real_)
  ms_safe_numeric(length(unique(id[!is.na(id)])))
}

ms_gee_id_variable <- function(call_obj) {
  if (is.null(call_obj) || is.null(call_obj$id)) return("")
  ms_deparse_call(call_obj$id)
}

ms_gee_correlation_structure <- function(x, raw) {
  corstr <- tryCatch(x$geese$corstr, error = function(e) NULL)
  if (!is.null(corstr) && length(corstr) && !is.na(corstr[[1]])) {
    return(ms_gee_correlation_label(as.character(corstr[[1]])))
  }
  m <- regexec("Correlation structure\\s*=\\s*([^\\r\\n]+)", raw, ignore.case = TRUE)
  hit <- regmatches(raw, m)[[1]]
  if (length(hit) >= 2L) ms_gee_correlation_label(trimws(hit[[2]])) else ""
}

ms_gee_correlation_label <- function(value) {
  value <- trimws(as.character(value %||% ""))
  lower <- tolower(value)
  if (lower %in% c("excha", "exchangeable")) return("exchangeable")
  if (lower %in% c("independence", "independence working", "ind")) return("independence")
  if (lower %in% c("ar1", "ar-1", "ar(1)")) return("AR(1)")
  if (lower %in% c("unstructured", "unstr")) return("unstructured")
  value
}

ms_gee_scale_parameter <- function(raw) {
  rx <- paste0(
    "Estimated Scale Parameters:[\\s\\S]*?Estimate\\s+Std\\.err\\s*",
    "\\n\\s*\\(Intercept\\)\\s+",
    "(-?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)(?:e[+-]?[0-9]+)?)\\s+",
    "(-?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)(?:e[+-]?[0-9]+)?)"
  )
  m <- regexec(rx, raw, ignore.case = TRUE, perl = TRUE)
  hit <- regmatches(raw, m)[[1]]
  if (length(hit) < 3L) {
    return(list(estimate = NA_real_, std_error = NA_real_))
  }
  list(
    estimate = ms_safe_numeric(as.numeric(hit[[2]])),
    std_error = ms_safe_numeric(as.numeric(hit[[3]]))
  )
}

ms_gee_working_correlation_parameter <- function(x, raw) {
  alpha <- tryCatch(x$geese$alpha, error = function(e) NULL)
  estimate <- NA_real_
  name <- "alpha"
  if (!is.null(alpha) && length(alpha) > 0L) {
    estimate <- ms_safe_numeric(as.numeric(alpha[[1L]]))
    alpha_name <- names(alpha)[[1L]] %||% ""
    if (nzchar(alpha_name)) name <- alpha_name
  }

  std_error <- NA_real_
  valpha <- tryCatch(x$geese$valpha, error = function(e) NULL)
  if (!is.null(valpha) && length(valpha) > 0L) {
    v <- if (is.matrix(valpha)) valpha[1L, 1L] else valpha[[1L]]
    v <- ms_safe_numeric(as.numeric(v))
    if (!is.na(v) && v >= 0) std_error <- sqrt(v)
  }

  if (is.na(estimate)) {
    rx <- paste0(
      "Estimated Correlation Parameters:[\\s\\S]*?Estimate\\s+Std\\.err\\s*",
      "\\n\\s*([^\\s]+)\\s+",
      "(-?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)(?:e[+-]?[0-9]+)?)\\s+",
      "(-?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)(?:e[+-]?[0-9]+)?)"
    )
    m <- regexec(rx, raw, ignore.case = TRUE, perl = TRUE)
    hit <- regmatches(raw, m)[[1]]
    if (length(hit) >= 4L) {
      name <- hit[[2]]
      estimate <- ms_safe_numeric(as.numeric(hit[[3]]))
      std_error <- ms_safe_numeric(as.numeric(hit[[4]]))
    }
  }

  list(name = name, estimate = estimate, std_error = std_error)
}
