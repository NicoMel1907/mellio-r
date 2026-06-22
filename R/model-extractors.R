# Model object extractors
# S3 methods for common R model objects

#' @rdname melliotab
#' @param conf.int Include confidence intervals (default TRUE)
#' @param conf.level Confidence level (default 0.95)
#' @export
melliotab.lm <- function(x, style = "apa7", title = NULL,
                          number = NULL, note = NULL,
                          source = NULL, conf.int = TRUE,
                          conf.level = 0.95, decimals = 2L,
                          p_decimals = 3L, ...) {
  rlang::check_installed("broom", reason = "to extract model coefficients")

  tidy_df <- broom::tidy(x, conf.int = conf.int, conf.level = conf.level)

  # Rename columns to APA standard
  col_map <- c(
    term = " ", estimate = "B", std.error = "SE",
    statistic = "t", p.value = "p",
    conf.low = "Lower CI", conf.high = "Upper CI"
  )
  current_names <- names(tidy_df)
  new_names <- vapply(current_names, function(n) {
    if (n %in% names(col_map)) col_map[n] else n
  }, character(1), USE.NAMES = FALSE)
  names(tidy_df) <- new_names

  result <- melliotab.data.frame(
    tidy_df, style = style, title = title,
    number = number, note = note,
    source = source, decimals = decimals,
    p_decimals = p_decimals, ...
  )

  # Store model for auto-note generation
  result$model <- x
  glance_df <- broom::glance(x)
  result$model_summary <- glance_df

  # Auto-generate note if none provided
  if (is.null(note)) {
    n <- stats::nobs(x)
    r2 <- glance_df$r.squared
    adj_r2 <- glance_df$adj.r.squared
    result$note <- sprintf("N = %d. R\u00B2 = %.3f, Adjusted R\u00B2 = %.3f.",
                           n, r2, adj_r2)
  }

  result
}

#' @rdname melliotab
#' @param exponentiate Show exponentiated coefficients (odds ratios for logistic)
#' @export
melliotab.glm <- function(x, style = "apa7", title = NULL,
                           number = NULL, note = NULL,
                           source = NULL, conf.int = TRUE,
                           exponentiate = FALSE, decimals = 2L,
                           p_decimals = 3L, ...) {
  rlang::check_installed("broom", reason = "to extract model coefficients")

  tidy_df <- broom::tidy(x, conf.int = conf.int, exponentiate = exponentiate)

  # Determine statistic name based on family
  stat_name <- if (x$family$family %in% c("binomial", "poisson")) "z" else "t"
  est_name <- if (exponentiate) "OR" else "B"

  col_map <- c(
    term = " ", estimate = est_name, std.error = "SE",
    statistic = stat_name, p.value = "p",
    conf.low = "Lower CI", conf.high = "Upper CI"
  )
  current_names <- names(tidy_df)
  new_names <- vapply(current_names, function(n) {
    if (n %in% names(col_map)) col_map[n] else n
  }, character(1), USE.NAMES = FALSE)
  names(tidy_df) <- new_names

  result <- melliotab.data.frame(
    tidy_df, style = style, title = title,
    number = number, note = note,
    source = source, decimals = decimals,
    p_decimals = p_decimals, ...
  )

  result$model <- x
  result$model_summary <- broom::glance(x)

  if (is.null(note)) {
    n <- stats::nobs(x)
    result$note <- sprintf("N = %d.", n)
  }

  result
}

#' @rdname melliotab
#' @param effect_size Include effect size (eta-squared) for ANOVA
#' @export
melliotab.aov <- function(x, style = "apa7", title = NULL,
                           number = NULL, note = NULL,
                           source = NULL, decimals = 2L,
                           p_decimals = 3L, effect_size = TRUE, ...) {
  rlang::check_installed("broom", reason = "to extract ANOVA results")

  tidy_df <- broom::tidy(x)

  col_map <- c(
    term = "Source", df = "df", sumsq = "SS",
    meansq = "MS", statistic = "F", p.value = "p"
  )
  current_names <- names(tidy_df)
  new_names <- vapply(current_names, function(n) {
    if (n %in% names(col_map)) col_map[n] else n
  }, character(1), USE.NAMES = FALSE)
  names(tidy_df) <- new_names

  # Add eta-squared if requested
  if (effect_size) {
    ss_vals <- tidy_df$SS
    total_ss <- sum(as.numeric(ss_vals), na.rm = TRUE)
    if (total_ss > 0) {
      eta2 <- as.numeric(ss_vals) / total_ss
      tidy_df[["\u03B7\u00B2"]] <- eta2
    }
  }

  result <- melliotab.data.frame(
    tidy_df, style = style, title = title,
    number = number, note = note,
    source = source, decimals = decimals,
    p_decimals = p_decimals, ...
  )

  result$model <- x

  result
}

#' @rdname melliotab
#' @export
melliotab.htest <- function(x, style = "apa7", title = NULL,
                             number = NULL, note = NULL,
                             source = NULL, decimals = 2L,
                             p_decimals = 3L, ...) {
  # Build a summary data frame from the htest object
  rows <- list()

  # Test name as first column
  test_name <- x$method %||% "Test"

  # Statistic
  if (!is.null(x$statistic)) {
    stat_name <- names(x$statistic)
    rows$statistic <- x$statistic
    stat_label <- stat_name
  }

  # Degrees of freedom
  if (!is.null(x$parameter)) {
    rows$df <- x$parameter
  }

  # P-value
  if (!is.null(x$p.value)) {
    rows$p <- x$p.value
  }

  # Estimate
  if (!is.null(x$estimate)) {
    for (i in seq_along(x$estimate)) {
      nm <- names(x$estimate)[i]
      rows[[nm]] <- x$estimate[i]
    }
  }

  # Confidence interval
  if (!is.null(x$conf.int)) {
    rows[["Lower CI"]] <- x$conf.int[1]
    rows[["Upper CI"]] <- x$conf.int[2]
  }

  # Build data frame
  df <- as.data.frame(lapply(rows, function(v) unname(v)),
                       stringsAsFactors = FALSE, check.names = FALSE)

  if (is.null(title)) title <- test_name

  result <- melliotab.data.frame(
    df, style = style, title = title,
    number = number, note = note,
    source = source, decimals = decimals,
    p_decimals = p_decimals, ...
  )

  result$model <- x
  result
}

#' @rdname melliotab
#' @param diagonal Diagonal display: "dash" (em-dash), "one" (keep 1.00), "blank"
#' @param triangle Which triangle to show: "all", "lower", "upper"
#' @export
melliotab.matrix <- function(x, style = "apa7", title = NULL,
                              number = NULL, note = NULL,
                              source = NULL, decimals = 2L,
                              p_decimals = 3L, diagonal = "dash",
                              triangle = "all", ...) {
  # Check if it looks like a correlation matrix
  is_corr <- is_correlation_matrix(x)

  # Convert to data frame with row names as first column
  df <- as.data.frame(x)
  if (!is.null(rownames(x))) {
    df <- cbind(data.frame(Variable = rownames(x), stringsAsFactors = FALSE), df)
  }

  result <- melliotab.data.frame(
    df, style = style, title = title,
    number = number, note = note,
    source = source, decimals = decimals,
    p_decimals = p_decimals, ...
  )

  # For correlation matrices, force all numeric columns to "estimate" type
  # so they get properly formatted with decimals and leading zero removal
  if (is_corr) {
    for (j in seq_along(result$column_types)) {
      if (result$column_types[j] == "default") {
        result$column_types[j] <- "estimate"
      }
    }
    # Re-format with corrected types
    sc <- result$style_config
    lz_cols <- detect_leading_zero_cols(names(result$raw_data))
    # For correlation values, all numeric columns should have leading zeros removed
    if (isTRUE(sc$remove_leading_zeros)) {
      for (j in seq_along(result$column_types)) {
        if (result$column_types[j] == "estimate") lz_cols[j] <- TRUE
      }
    }
    result$data <- format_table_data(
      result$raw_data, result$column_types, lz_cols,
      decimals = decimals, p_decimals = p_decimals,
      remove_lz = isTRUE(sc$remove_leading_zeros),
      style_config = sc
    )
  }

  if (is_corr) {
    result$options$is_correlation <- TRUE
    result$options$diagonal_mode <- diagonal
    result$options$triangle <- triangle

    # Apply correlation formatting
    result <- apply_correlation_formatting(result)
  }

  result
}

#' @rdname melliotab
#' @param standardized Include standardized estimates (default TRUE)
#' @export
melliotab.lavaan <- function(x, style = "apa7", title = NULL,
                              number = NULL, note = NULL,
                              source = NULL, what = NULL, section = NULL,
                              standardized = TRUE, decimals = 2L,
                              p_decimals = 3L, ...) {
  rlang::check_installed("lavaan", reason = "to extract lavaan model results")

  section <- mellio_resolve_section(section = section, what = what)
  payload <- mellio_payload(x, standardized = standardized, ...)
  result <- melliotab_from_payload(
    payload,
    section = section,
    style = style,
    title = title,
    number = number,
    note = note,
    source = source,
    decimals = decimals,
    p_decimals = p_decimals,
    ...
  )
  result$model <- payload
  result
}

#' @rdname melliotab
#' @export
melliotab.FitDiff <- function(x, style = "apa7", title = NULL,
                               number = NULL, note = NULL,
                               source = NULL, section = NULL, what = NULL,
                               decimals = 2L, p_decimals = 3L, ...) {
  rlang::check_installed("semTools",
    reason = "to extract FitDiff model comparisons"
  )

  section <- mellio_resolve_section(section = section, what = what)

  if (is.null(section)) {
    # Build a preview of what's available
    sections <- character(0)
    nested <- methods::slot(x, "nested")
    fit <- methods::slot(x, "fit")
    fit_diff <- methods::slot(x, "fit.diff")

    if (!is.null(nested) && nrow(nested) > 0) {
      sections <- c(sections, paste0(
        "{.val comparison} (chi-squared difference test, ",
        nrow(nested), " rows)"
      ))
    }
    if (!is.null(fit) && nrow(fit) > 0) {
      sections <- c(sections, paste0(
        "{.val fit} (model fit indices, ",
        nrow(fit), " models)"
      ))
    }
    if (!is.null(fit_diff) && nrow(fit_diff) > 0) {
      sections <- c(sections, paste0(
        "{.val diff} (differences in fit indices, ",
        nrow(fit_diff), " rows)"
      ))
    }

    cli::cli_abort(c(
      "FitDiff objects contain multiple tables.",
      "i" = "Specify {.arg section} to choose which table:",
      set_names(sections, rep("*", length(sections))),
      "i" = 'Example: {.code melliotab(x, section = "fit")}'
    ))
  }

  section <- match.arg(section, c("comparison", "fit", "diff"))

  switch(section,
    comparison = .fitdiff_comparison(x, style, title, number, note,
      source, decimals, p_decimals, ...),
    fit = .fitdiff_fit(x, style, title, number, note,
      source, decimals, p_decimals, ...),
    diff = .fitdiff_diff(x, style, title, number, note,
      source, decimals, p_decimals, ...)
  )
}

#' Build chi-squared difference test table from FitDiff
#' @keywords internal
.fitdiff_comparison <- function(x, style, title, number, note,
                                 source, decimals, p_decimals, ...) {
  nested <- methods::slot(x, "nested")
  if (is.null(nested) || nrow(nested) == 0) {
    cli::cli_abort(
      "No nested model comparison data in this FitDiff object."
    )
  }

  # Keep only the key comparison columns
  key_cols <- c("Df", "AIC", "BIC", "Chisq",
                "Chisq diff", "Df diff", "Pr(>Chisq)")
  avail <- key_cols[key_cols %in% names(nested)]
  nested <- nested[, avail, drop = FALSE]

  # Add model names as first column
  df <- cbind(
    data.frame(Model = rownames(nested), stringsAsFactors = FALSE),
    as.data.frame(nested, check.names = FALSE)
  )
  rownames(df) <- NULL

  # Rename columns to cleaner names
  col_renames <- c(
    "Df" = "df",
    "Chisq" = "\u03C7\u00B2",
    "Chisq diff" = "\u0394\u03C7\u00B2",
    "Df diff" = "\u0394df",
    "Pr(>Chisq)" = "p"
  )
  for (i in seq_along(col_renames)) {
    idx <- which(names(df) == names(col_renames)[i])
    if (length(idx) > 0) names(df)[idx] <- col_renames[i]
  }

  if (is.null(title)) title <- "Chi-Squared Difference Test"

  result <- melliotab.data.frame(
    df, style = style, title = title,
    number = number, note = note,
    source = source, decimals = decimals,
    p_decimals = p_decimals, ...
  )
  result
}

#' Build model fit indices table from FitDiff
#' @keywords internal
.fitdiff_fit <- function(x, style, title, number, note,
                          source, decimals, p_decimals, ...) {
  fit <- methods::slot(x, "fit")
  if (is.null(fit) || nrow(fit) == 0) {
    cli::cli_abort("No fit indices data in this FitDiff object.")
  }

  # Keep only the key fit indices (matching summary output)
  key_cols <- c("chisq", "df", "pvalue", "rmsea", "cfi", "tli",
                "srmr", "aic", "bic")
  avail <- key_cols[key_cols %in% names(fit)]
  fit <- fit[, avail, drop = FALSE]

  # Add model names as first column
  df <- cbind(
    data.frame(Model = rownames(fit), stringsAsFactors = FALSE),
    as.data.frame(fit, check.names = FALSE)
  )
  rownames(df) <- NULL

  # Rename to display names
  col_renames <- c(
    "chisq" = "\u03C7\u00B2",
    "df" = "df",
    "pvalue" = "p",
    "rmsea" = "RMSEA",
    "cfi" = "CFI",
    "tli" = "TLI",
    "srmr" = "SRMR",
    "aic" = "AIC",
    "bic" = "BIC"
  )
  for (i in seq_along(col_renames)) {
    idx <- which(names(df) == names(col_renames)[i])
    if (length(idx) > 0) names(df)[idx] <- col_renames[i]
  }

  if (is.null(title)) title <- "Model Fit Indices"

  result <- melliotab.data.frame(
    df, style = style, title = title,
    number = number, note = note,
    source = source, decimals = decimals,
    p_decimals = p_decimals, ...
  )
  result
}

#' Build differences in fit indices table from FitDiff
#' @keywords internal
.fitdiff_diff <- function(x, style, title, number, note,
                           source, decimals, p_decimals, ...) {
  fit_diff <- methods::slot(x, "fit.diff")
  if (is.null(fit_diff) || nrow(fit_diff) == 0) {
    cli::cli_abort("No fit differences data in this FitDiff object.")
  }

  # Keep only the key difference columns
  key_cols <- c("df", "rmsea", "cfi", "tli", "srmr", "aic", "bic")
  avail <- key_cols[key_cols %in% names(fit_diff)]
  fit_diff <- fit_diff[, avail, drop = FALSE]

  # Add comparison names as first column
  df <- cbind(
    data.frame(Comparison = rownames(fit_diff),
               stringsAsFactors = FALSE),
    as.data.frame(fit_diff, check.names = FALSE)
  )
  rownames(df) <- NULL

  # Rename to delta symbols
  col_renames <- c(
    "df" = "\u0394df",
    "rmsea" = "\u0394RMSEA",
    "cfi" = "\u0394CFI",
    "tli" = "\u0394TLI",
    "srmr" = "\u0394SRMR",
    "aic" = "\u0394AIC",
    "bic" = "\u0394BIC"
  )
  for (i in seq_along(col_renames)) {
    idx <- which(names(df) == names(col_renames)[i])
    if (length(idx) > 0) names(df)[idx] <- col_renames[i]
  }

  if (is.null(title)) title <- "Differences in Fit Indices"

  result <- melliotab.data.frame(
    df, style = style, title = title,
    number = number, note = note,
    source = source, decimals = decimals,
    p_decimals = p_decimals, ...
  )
  result
}

#' @rdname melliotab
#' @param cut Minimum absolute loading to display (default 0, show all).
#'   Loadings below this threshold are shown as blank cells.
#' @param sort Sort items by their primary loading (default FALSE)
#' @export
melliotab.fa <- function(x, style = "apa7", title = NULL,
                          number = NULL, note = NULL,
                          source = NULL, section = NULL, what = NULL,
                          cut = 0, sort = FALSE,
                          decimals = 2L, p_decimals = 3L, ...) {
  section <- mellio_resolve_section(
    section = section,
    what = what,
    default = "loadings",
    choices = c("loadings", "variance", "fit")
  )

  switch(section,
    loadings = .fa_loadings_table(x, style, title, number, note,
      source, cut, sort, decimals, p_decimals, ...),
    variance = .fa_variance_table(x, style, title, number, note,
      source, decimals, p_decimals, ...),
    fit = .fa_fit_table(x, style, title, number, note,
      source, decimals, p_decimals, ...)
  )
}

#' Build factor loadings table from psych::fa
#' @keywords internal
.fa_loadings_table <- function(x, style, title, number, note, source,
                                cut, sort, decimals, p_decimals, ...) {
  # Extract loadings matrix
  loads <- as.matrix(unclass(x$loadings))
  n_factors <- ncol(loads)

  # Rename factor columns to "Factor 1", "Factor 2", etc.
  colnames(loads) <- paste("Factor", seq_len(n_factors))

  # Sort by primary loading if requested
  if (sort) {
    primary <- apply(abs(loads), 1, which.max)
    secondary <- apply(abs(loads), 1, max)
    ord <- order(primary, -secondary)
    loads <- loads[ord, , drop = FALSE]
  }

  # Build data frame: Variable + loadings + h2 + u2
  df <- data.frame(
    Variable = rownames(loads),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  for (j in seq_len(n_factors)) {
    df[[colnames(loads)[j]]] <- loads[, j]
  }
  df[["h\u00B2"]] <- x$communality[rownames(loads)]
  df[["u\u00B2"]] <- x$uniquenesses[rownames(loads)]

  # Apply cutoff: blank out loadings below threshold
  if (cut > 0) {
    for (j in seq_len(n_factors)) {
      col_name <- colnames(loads)[j]
      vals <- df[[col_name]]
      df[[col_name]] <- ifelse(abs(vals) < cut, NA, vals)
    }
  }

  if (is.null(title)) title <- "Factor Loadings"

  result <- melliotab.data.frame(
    df, style = style, title = title,
    number = number, note = note,
    source = source, decimals = decimals,
    p_decimals = p_decimals, ...
  )

  # Replace NA strings with empty cells
  for (j in seq_len(ncol(result$data))) {
    result$data[[j]][is.na(result$data[[j]]) |
                     trimws(result$data[[j]]) == "NA"] <- ""
  }
  for (j in seq_len(ncol(result$raw_data))) {
    result$raw_data[[j]][is.na(result$raw_data[[j]])] <- ""
  }

  result$model <- x

  # Auto-generate note
  if (is.null(note)) {
    parts <- character(0)
    n <- x$n.obs
    if (!is.null(n) && !is.na(n)) {
      parts <- c(parts, sprintf("N = %d", as.integer(n)))
    }
    parts <- c(parts, sprintf("Factors = %d", n_factors))
    rot <- x$rotation %||% "none"
    parts <- c(parts, sprintf("Rotation = %s", rot))
    fm <- x$fm %||% "unknown"
    fm_label <- switch(fm,
      pa = "Principal Axis", minres = "Minimum Residual",
      ml = "Maximum Likelihood", wls = "Weighted Least Squares",
      gls = "Generalized Least Squares", fm)
    parts <- c(parts, sprintf("Method = %s", fm_label))
    if (cut > 0) {
      parts <- c(parts, sprintf("Loadings < %.2f suppressed", cut))
    }
    result$note <- paste(parts, collapse = ". ") |> paste0(".")
  }

  result
}

#' Build variance accounted for table from psych::fa
#' @keywords internal
.fa_variance_table <- function(x, style, title, number, note, source,
                                decimals, p_decimals, ...) {
  vacc <- x$Vaccounted
  if (is.null(vacc)) {
    cli::cli_abort("No variance-accounted-for data in this fa object.")
  }

  # Vaccounted is a matrix: rows = metrics, cols = factors
  # Rename columns
  n_factors <- ncol(vacc)
  col_names <- paste("Factor", seq_len(n_factors))

  df <- data.frame(
    Measure = rownames(vacc),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  for (j in seq_len(n_factors)) {
    df[[col_names[j]]] <- vacc[, j]
  }

  if (is.null(title)) title <- "Variance Accounted For"

  result <- melliotab.data.frame(
    df, style = style, title = title,
    number = number, note = note,
    source = source, decimals = decimals,
    p_decimals = p_decimals, ...
  )
  result$model <- x
  result
}

#' Build fit indices table from psych::fa
#' @keywords internal
.fa_fit_table <- function(x, style, title, number, note, source,
                           decimals, p_decimals, ...) {
  indices <- list()

  if (!is.null(x$n.obs) && !is.na(x$n.obs)) {
    indices[["Observations"]] <- as.integer(x$n.obs)
  }
  indices[["Factors"]] <- x$factors
  if (!is.null(x$rotation)) indices[["Rotation"]] <- x$rotation
  if (!is.null(x$fm)) indices[["Method"]] <- x$fm
  if (!is.null(x$fit) && !is.na(x$fit)) {
    indices[["Fit Index"]] <- x$fit
  }
  if (!is.null(x$rms) && !is.na(x$rms)) {
    indices[["RMSR"]] <- x$rms
  }
  if (!is.null(x$RMSEA) && length(x$RMSEA) > 0 && !is.na(x$RMSEA[1])) {
    indices[["RMSEA"]] <- x$RMSEA[1]
  }
  if (!is.null(x$TLI) && !is.na(x$TLI)) {
    indices[["TLI"]] <- x$TLI
  }
  if (!is.null(x$BIC) && !is.na(x$BIC)) {
    indices[["BIC"]] <- x$BIC
  }

  df <- data.frame(
    Measure = names(indices),
    Value = unlist(indices, use.names = FALSE),
    stringsAsFactors = FALSE, check.names = FALSE
  )

  if (is.null(title)) title <- "Factor Analysis Fit Indices"

  result <- melliotab.data.frame(
    df, style = style, title = title,
    number = number, note = note,
    source = source, decimals = decimals,
    p_decimals = p_decimals, ...
  )
  result$model <- x
  result
}

#' Build a fit indices table from a lavaan model
#' @keywords internal
.lavaan_fit_table <- function(x, style, title, number, note, source,
                               decimals, p_decimals, ...) {
  fm <- lavaan::fitMeasures(x)

  # Key indices to display
  keys <- c("chisq", "df", "pvalue", "cfi", "tli",
            "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr")
  display_names <- c(
    chisq = "\u03C7\u00B2",
    df = "df",
    pvalue = "p",
    cfi = "CFI",
    tli = "TLI",
    rmsea = "RMSEA",
    rmsea.ci.lower = "RMSEA CI Lower",
    rmsea.ci.upper = "RMSEA CI Upper",
    srmr = "SRMR"
  )

  # Keep only available indices
  avail <- keys[keys %in% names(fm)]
  df <- data.frame(
    Measure = unname(display_names[avail]),
    Value = unname(fm[avail]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (is.null(title)) title <- "Model Fit Indices"

  result <- melliotab.data.frame(
    df, style = style, title = title,
    number = number, note = note,
    source = source, decimals = decimals,
    p_decimals = p_decimals, ...
  )

  result$model <- x
  result
}
