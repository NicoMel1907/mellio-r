# R bridge -- Bayesian model summary extractors.

#' @rdname mellio_payload
#' @export
mellio_payload.brmsfit <- function(x, ..., .call = NULL) {
  s <- if (is.list(x) && (!is.null(x$fixed) || !is.null(x$fixed_effects))) {
    x
  } else {
    tryCatch(summary(x), error = function(e) NULL)
  }
  if (is.null(s) || !is.list(s)) {
    stop("Could not summarize brmsfit object.", call. = FALSE)
  }

  fixed <- s$fixed %||% s$fixed_effects %||% NULL
  if (is.null(fixed)) {
    stop("brmsfit summary does not contain fixed-effect parameters.", call. = FALSE)
  }

  rows <- ms_bayesian_parameter_rows(fixed)
  if (!length(rows)) {
    stop("No Bayesian parameter rows could be extracted.", call. = FALSE)
  }

  ms_bayesian_payload(
    x = x,
    rows = rows,
    method = "brms",
    source = "brms::brm",
    packages = "brms",
    .call = .call
  )
}

#' @rdname mellio_payload
#' @export
mellio_payload.stanreg <- function(x, ..., .call = NULL) {
  summary_obj <- ms_stanreg_summary_matrix(x)
  rows <- ms_bayesian_parameter_rows(summary_obj)
  if (!length(rows)) {
    stop("No rstanarm parameter rows could be extracted.", call. = FALSE)
  }

  ms_bayesian_payload(
    x = x,
    rows = rows,
    method = "rstanarm",
    source = "rstanarm::stan_glm",
    packages = "rstanarm",
    .call = .call
  )
}

ms_bayesian_payload <- function(x, rows, method, source, packages, .call = NULL) {
  ms_build_envelope(
    type = "bayesian_model_summary",
    type_label = "Bayesian model summary",
    call = ms_bayesian_call(x, .call),
    fields = list(
      table_type = "bayesian_parameters",
      method = method,
      interval_type = "credible_interval",
      columns = ms_bayesian_columns(rows),
      rows = rows,
      n_parameters = length(rows),
      source = source
    ),
    raw_output = ms_capture_output_safe(x)$text,
    packages = ms_packages_basic(extras = packages),
    card_kind = "table"
  )
}

ms_stanreg_summary_matrix <- function(x) {
  if (is.list(x) && !is.null(x$posterior_summary)) return(x$posterior_summary)
  if (is.list(x) && !is.null(x$summary)) return(x$summary)

  ps <- tryCatch({
    if (requireNamespace("rstanarm", quietly = TRUE)) {
      rstanarm::posterior_summary(x, prob = 0.95)
    } else {
      NULL
    }
  }, error = function(e) NULL)
  if (!is.null(ps)) return(ps)

  s <- tryCatch(summary(x), error = function(e) NULL)
  if (is.null(s)) {
    stop("Could not summarize stanreg object.", call. = FALSE)
  }
  s
}

ms_bayesian_parameter_rows <- function(x) {
  df <- as.data.frame(x)
  if (!nrow(df)) return(list())
  lower <- ms_df_names(df)
  rn <- rownames(df)
  if (is.null(rn) || !length(rn)) rn <- paste0("parameter_", seq_len(nrow(df)))

  est_col <- ms_df_col(lower, c("estimate", "mean", "median"))
  se_col <- ms_df_col(lower, c("est_error", "std_error", "standard_error", "se", "sd"))
  low_col <- ms_df_col(lower, c("l_95_ci", "q2_5", "x2_5", "2_5", "ci_lower", "lower", "q025"))
  high_col <- ms_df_col(lower, c("u_95_ci", "q97_5", "x97_5", "97_5", "ci_upper", "upper", "q975"))
  rhat_col <- ms_df_col(lower, c("rhat", "r_hat"))
  ess_bulk_col <- ms_df_col(lower, c("bulk_ess", "ess_bulk", "n_eff", "neff", "ess"))
  ess_tail_col <- ms_df_col(lower, c("tail_ess", "ess_tail"))

  if (is.na(est_col)) {
    numeric_cols <- which(vapply(df, is.numeric, logical(1)))
    if (length(numeric_cols)) est_col <- numeric_cols[[1]]
  }
  if (is.na(est_col)) return(list())

  lapply(seq_len(nrow(df)), function(i) {
    row <- list(
      term = rn[[i]],
      estimate = ms_safe_numeric(df[[est_col]][[i]])
    )
    if (!is.na(se_col)) row$std_error <- ms_safe_numeric(df[[se_col]][[i]])
    if (!is.na(low_col)) row$ci_lower <- ms_safe_numeric(df[[low_col]][[i]])
    if (!is.na(high_col)) row$ci_upper <- ms_safe_numeric(df[[high_col]][[i]])
    if (!is.na(rhat_col)) row$rhat <- ms_safe_numeric(df[[rhat_col]][[i]])
    if (!is.na(ess_bulk_col)) row$ess_bulk <- ms_safe_numeric(df[[ess_bulk_col]][[i]])
    if (!is.na(ess_tail_col)) row$ess_tail <- ms_safe_numeric(df[[ess_tail_col]][[i]])
    row
  })
}

ms_bayesian_columns <- function(rows) {
  has <- function(key) {
    any(vapply(rows, function(row) !is.null(row[[key]]) && !is.na(row[[key]]), logical(1)))
  }
  cols <- list(
    list(key = "term", label = "Term", format = "text"),
    list(key = "estimate", label = "Estimate", format = "number")
  )
  if (has("std_error")) cols <- c(cols, list(list(key = "std_error", label = "SE", format = "number")))
  if (has("ci_lower") && has("ci_upper")) {
    cols <- c(cols, list(
      list(key = "ci_lower", label = "95% CrI lower", format = "number"),
      list(key = "ci_upper", label = "95% CrI upper", format = "number")
    ))
  }
  if (has("rhat")) cols <- c(cols, list(list(key = "rhat", label = "Rhat", format = "number")))
  if (has("ess_bulk")) cols <- c(cols, list(list(key = "ess_bulk", label = "Bulk ESS", format = "integer")))
  if (has("ess_tail")) cols <- c(cols, list(list(key = "ess_tail", label = "Tail ESS", format = "integer")))
  cols
}

ms_bayesian_call <- function(x, .call = NULL) {
  if (!is.null(.call)) return(trimws(gsub("\\s+", " ", .call)))
  call_obj <- tryCatch(stats::getCall(x), error = function(e) NULL)
  if (!is.null(call_obj)) return(ms_deparse_call(call_obj))
  call_obj <- x$call %||% NULL
  if (!is.null(call_obj)) return(ms_deparse_call(call_obj))
  NA_character_
}
