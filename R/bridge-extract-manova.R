# R bridge -- MANOVA multivariate tests.
#
# Base R's manova object has class c("manova", "maov", "aov", "mlm", "lm").
# Without explicit methods it falls through to the univariate aov bridge,
# which loses the multivariate criteria. These methods emit a table card:
# one row per effect/test criterion (Pillai, Wilks, Hotelling-Lawley, Roy).

#' @rdname mellio_payload
#' @export
mellio_payload.manova <- function(x, tests = c("Pillai", "Wilks", "Hotelling-Lawley", "Roy"),
                              ..., .call = NULL) {
  call_str <- if (!is.null(.call)) {
    .call
  } else {
    formula_txt <- tryCatch(
      paste(deparse(stats::formula(x), width.cutoff = 500L), collapse = " "),
      error = function(e) NA_character_
    )
    if (!is.na(formula_txt) && nzchar(formula_txt)) formula_txt else ms_model_call_string(x)
  }

  summaries <- lapply(tests, function(test) {
    tryCatch(stats::summary.manova(x, test = test), error = function(e) NULL)
  })
  summaries <- Filter(Negate(is.null), summaries)
  if (length(summaries) == 0L) {
    stop("Could not extract MANOVA multivariate tests.", call. = FALSE)
  }

  info <- ms_manova_table_from_summaries(summaries)
  fields <- info$fields
  vars <- ms_anova_vars_from_call(call_str)
  if (!is.null(vars$outcome)) fields$outcome <- vars$outcome
  if (!is.null(vars$predictor)) fields$predictor <- vars$predictor

  raw_output <- paste(vapply(summaries, function(s) {
    paste(utils::capture.output(print(s)), collapse = "\n")
  }, character(1)), collapse = "\n\n")

  ms_build_envelope(
    type = "manova_multivariate_tests",
    type_label = "MANOVA multivariate tests",
    call = trimws(gsub("\\s+", " ", call_str)),
    fields = fields,
    raw_output = raw_output,
    packages = ms_packages_basic(),
    card_kind = "table"
  )
}

#' @rdname mellio_payload
#' @export
mellio_payload.summary.manova <- function(x, ..., .call = NULL) {
  call_str <- if (!is.null(.call)) .call else "summary(manova(...))"
  info <- ms_manova_table_from_summaries(list(x))
  ms_build_envelope(
    type = "manova_multivariate_tests",
    type_label = "MANOVA multivariate tests",
    call = trimws(gsub("\\s+", " ", call_str)),
    fields = info$fields,
    raw_output = ms_capture_output(x),
    packages = ms_packages_basic(),
    card_kind = "table"
  )
}

ms_manova_table_from_summaries <- function(summaries) {
  rows <- unlist(lapply(summaries, ms_manova_rows_from_summary), recursive = FALSE)
  rows <- Filter(function(row) {
    !is.null(row) &&
      !is.na(row$criterion_value %||% NA_real_) &&
      !is.na(row$f %||% NA_real_)
  }, rows)
  if (length(rows) == 0L) {
    stop("MANOVA summary did not contain reportable multivariate tests.", call. = FALSE)
  }

  list(fields = list(
    table_type = "manova_multivariate_tests",
    columns = ms_manova_columns(),
    rows = rows,
    source = "R"
  ))
}

ms_manova_rows_from_summary <- function(x) {
  stats <- x$stats
  if (is.null(stats) || !is.matrix(stats) || nrow(stats) == 0L) return(list())

  cn <- colnames(stats)
  criterion_col <- setdiff(cn, c("Df", "approx F", "num Df", "den Df", "Pr(>F)"))
  criterion_col <- criterion_col[!is.na(criterion_col) & nzchar(criterion_col)]
  if (length(criterion_col) == 0L) return(list())
  criterion <- criterion_col[[1L]]

  rn <- rownames(stats)
  if (is.null(rn)) rn <- paste0("Effect ", seq_len(nrow(stats)))

  lapply(seq_len(nrow(stats)), function(i) {
    effect <- rn[[i]]
    if (grepl("^Residuals?$", effect, ignore.case = TRUE)) return(NULL)
    f_val <- ms_safe_numeric(stats[i, "approx F"])
    h_df <- ms_safe_numeric(stats[i, "num Df"])
    e_df <- ms_safe_numeric(stats[i, "den Df"])
    p_val <- ms_safe_numeric(stats[i, "Pr(>F)"])
    eta <- ms_manova_partial_eta_sq(f_val, h_df, e_df)
    row <- list(
      effect = effect,
      test = ms_manova_test_label(criterion),
      criterion_value = ms_safe_numeric(stats[i, criterion]),
      f = f_val,
      hypothesis_df = h_df,
      error_df = e_df,
      p_value = p_val
    )
    if (!is.na(eta)) row$eta_sq_partial <- eta
    row
  })
}

ms_manova_columns <- function() {
  list(
    list(key = "effect", label = "Effect", format = "text"),
    list(key = "test", label = "Test", format = "text"),
    list(key = "criterion_value", label = "Value", format = "bounded"),
    list(key = "f", label = "F", format = "number"),
    list(key = "hypothesis_df", label = "Hypothesis df", format = "number"),
    list(key = "error_df", label = "Error df", format = "number"),
    list(key = "p_value", label = "p", format = "pvalue"),
    list(key = "eta_sq_partial", label = "partial \u03b7\u00b2", format = "bounded")
  )
}

ms_manova_test_label <- function(x) {
  x <- as.character(x %||% "")
  switch(
    x,
    "Pillai" = "Pillai's Trace",
    "Wilks" = "Wilks' Lambda",
    "Hotelling-Lawley" = "Hotelling-Lawley Trace",
    "Roy" = "Roy's Largest Root",
    x
  )
}

ms_manova_partial_eta_sq <- function(f_value, hypothesis_df, error_df) {
  if (is.na(f_value) || is.na(hypothesis_df) || is.na(error_df)) return(NA_real_)
  denom <- f_value * hypothesis_df + error_df
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  ms_safe_numeric((f_value * hypothesis_df) / denom)
}
