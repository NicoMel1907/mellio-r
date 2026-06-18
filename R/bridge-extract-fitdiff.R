# R bridge -- semTools FitDiff model comparisons.

#' @rdname mellio_payload
#' @param what Optional FitDiff table to extract: `"comparison"`, `"fit"`, or
#'   `"diff"`. Defaults to the first available table.
#' @export
mellio_payload.FitDiff <- function(x, what = NULL, ..., .call = NULL) {
  rlang::check_installed("semTools", reason = "to extract FitDiff model comparisons")

  info <- ms_fitdiff_table(x, what = what)
  ms_build_envelope(
    type       = "model_comparison",
    type_label = info$type_label,
    call       = ms_unsupported_call(.call),
    fields     = list(
      table_type = info$table_type,
      rows       = ms_rows_from_df(info$data),
      columns    = ms_table_columns_from_df(info$data),
      source     = "semTools::compareFit"
    ),
    raw_output = ms_capture_output(x),
    packages   = ms_packages_basic(extras = c("semTools", "lavaan")),
    card_kind  = "table"
  )
}

ms_fitdiff_table <- function(x, what = NULL) {
  nested <- tryCatch(methods::slot(x, "nested"), error = function(e) NULL)
  fit <- tryCatch(methods::slot(x, "fit"), error = function(e) NULL)
  fit_diff <- tryCatch(methods::slot(x, "fit.diff"), error = function(e) NULL)

  available <- character(0)
  if (!is.null(nested) && nrow(nested) > 0L) available <- c(available, "comparison")
  if (!is.null(fit) && nrow(fit) > 0L) available <- c(available, "fit")
  if (!is.null(fit_diff) && nrow(fit_diff) > 0L) available <- c(available, "diff")
  if (length(available) == 0L) {
    cli::cli_abort("FitDiff object does not contain any reportable tables.")
  }

  if (is.null(what)) {
    what <- available[[1L]]
  } else {
    what <- match.arg(what, c("comparison", "fit", "diff"))
    if (!what %in% available) {
      cli::cli_abort("FitDiff object does not contain a {.val {what}} table.")
    }
  }

  switch(what,
    comparison = ms_fitdiff_comparison_df(nested),
    fit = ms_fitdiff_fit_df(fit),
    diff = ms_fitdiff_diff_df(fit_diff)
  )
}

ms_fitdiff_comparison_df <- function(nested) {
  key_cols <- c("Df", "AIC", "BIC", "Chisq", "Chisq diff", "Df diff", "Pr(>Chisq)")
  df <- ms_fitdiff_named_df(nested, key_cols, first_col = "Model")
  names(df) <- ms_fitdiff_rename(names(df), c(
    "Df" = "df",
    "Chisq" = "chi_sq",
    "Chisq diff" = "delta_chi_sq",
    "Df diff" = "delta_df",
    "Pr(>Chisq)" = "p"
  ))
  list(
    data = df,
    table_type = "fitdiff_comparison",
    type_label = "Model fit comparison"
  )
}

ms_fitdiff_fit_df <- function(fit) {
  key_cols <- c("chisq", "df", "pvalue", "rmsea", "cfi", "tli", "srmr", "aic", "bic")
  df <- ms_fitdiff_named_df(fit, key_cols, first_col = "Model")
  names(df) <- ms_fitdiff_rename(names(df), c(
    chisq = "chi_sq",
    pvalue = "p",
    rmsea = "RMSEA",
    cfi = "CFI",
    tli = "TLI",
    srmr = "SRMR",
    aic = "AIC",
    bic = "BIC"
  ))
  list(
    data = df,
    table_type = "fitdiff_fit",
    type_label = "Model fit indices"
  )
}

ms_fitdiff_diff_df <- function(fit_diff) {
  key_cols <- c("df", "rmsea", "cfi", "tli", "srmr", "aic", "bic")
  df <- ms_fitdiff_named_df(fit_diff, key_cols, first_col = "Comparison")
  names(df) <- ms_fitdiff_rename(names(df), c(
    df = "delta_df",
    rmsea = "delta_RMSEA",
    cfi = "delta_CFI",
    tli = "delta_TLI",
    srmr = "delta_SRMR",
    aic = "delta_AIC",
    bic = "delta_BIC"
  ))
  list(
    data = df,
    table_type = "fitdiff_diff",
    type_label = "Model fit differences"
  )
}

ms_fitdiff_named_df <- function(x, key_cols, first_col) {
  avail <- key_cols[key_cols %in% names(x)]
  x <- x[, avail, drop = FALSE]
  labels <- rownames(x) %||% paste0("Model ", seq_len(nrow(x)))
  rownames(x) <- NULL
  cbind(
    data.frame(stats::setNames(list(labels), first_col), stringsAsFactors = FALSE),
    as.data.frame(x, check.names = FALSE)
  )
}

ms_fitdiff_rename <- function(names, map) {
  out <- names
  hit <- match(out, names(map), nomatch = 0L)
  replace <- hit > 0L
  out[replace] <- unname(map[hit[replace]])
  out
}
