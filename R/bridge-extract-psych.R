# R bridge -- psych extractors and custom result-table intake.
#
# These methods keep Mellio out of the statistics business: users run
# well-known R functions, then Mellio structures the returned object for
# editing, saving, and export.

#' @rdname mellio_payload
#' @export
mellio_payload.alpha <- function(x, ..., .call = NULL) {
  total <- x$total
  if (is.null(total) || !is.data.frame(total) || nrow(total) < 1L) {
    stop("psych::alpha() object is missing its total summary.", call. = FALSE)
  }

  raw_alpha <- ms_pick_numeric(total, "raw_alpha")
  if (is.na(raw_alpha)) {
    stop("psych::alpha() object does not contain raw_alpha.", call. = FALSE)
  }

  call_str <- ms_psych_call(x, .call, "psych::alpha(...)")
  fields <- list(
    statistic = list(name = "\u03b1", value = raw_alpha),
    p_value   = NA_real_
  )

  # P3c: 95% CI for \u03b1 via the analytic (Feldt) normal approximation
  # raw_alpha \u00b1 qnorm(.975) * ase. The ase column is in $total
  # alongside raw_alpha and is what psych itself uses for inference on
  # \u03b1. Skip when ase is missing or non-positive (degenerate fits).
  ase <- ms_pick_numeric(total, "ase")
  if (!is.na(ase) && ase > 0) {
    z <- stats::qnorm(0.975)
    ci_lo <- ms_safe_numeric(raw_alpha - z * ase)
    ci_hi <- ms_safe_numeric(raw_alpha + z * ase)
    if (!is.na(ci_lo) && !is.na(ci_hi)) {
      fields$statistic$ci <- I(c(ci_lo, ci_hi))
      fields$conf_level <- 0.95
    }
  }

  std_alpha <- ms_pick_numeric(total, "std.alpha")
  if (!is.na(std_alpha)) {
    fields$estimate <- list(name = "std_alpha", value = std_alpha)
  }

  avg_r <- ms_pick_numeric(total, "average_r")
  if (!is.na(avg_r)) {
    fields$effect <- list(name = "average_r", value = avg_r)
  }

  n <- ms_alpha_n(x)
  if (!is.na(n)) fields$n <- n

  k <- x$nvar %||% nrow(x$item.stats %||% data.frame())
  if (!is.null(k) && length(k) > 0L && !is.na(k)) {
    fields$item_count <- as.integer(k)
  }

  if (!is.null(x$item.stats) && nrow(x$item.stats) > 0L) {
    items <- rownames(x$item.stats)
    if (!is.null(items) && any(nzchar(items))) fields$items <- items
  }

  ms_build_envelope(
    type       = "cronbach_alpha",
    type_label = "Cronbach's alpha",
    call       = trimws(gsub("\\s+", " ", call_str)),
    fields     = fields,
    raw_output = ms_capture_output(x),
    packages   = ms_packages_basic("psych")
  )
}

#' @rdname mellio_payload
#' @export
mellio_payload.corr.test <- function(x, ..., .call = NULL) {
  rows <- ms_corr_test_rows(x)
  if (length(rows) == 0L) {
    stop("psych::corr.test() object does not contain reportable correlations.",
         call. = FALSE)
  }

  adjusted <- !identical(x$adjust %||% "none", "none")

  # Three rendering states based on what was passed to corr.test():
  #   A. corr.test(data)        -> sym=TRUE, full matrix
  #                                "Variable 1" / "Variable 2"
  #   B. corr.test(x, y) where y is one variable
  #      sym=FALSE, ncol(r)=1
  #      Drop the y column entirely; the outcome name moves to section
  #      meta ("vs. <name>") since it's redundant on every row.
  #   C. corr.test(x, y) where both multivariate
  #      sym=FALSE, ncol(r)>1
  #      "x-variable" / "y-variable" — mirrors psych's API; avoids
  #      claiming a directional predictor/outcome relationship the
  #      math doesn't actually impose.
  r_mat <- if (is.null(x$r)) matrix(numeric(0)) else as.matrix(x$r)
  sym <- isTRUE(x$sym) && nrow(r_mat) == ncol(r_mat)
  single_y_var <- NULL
  if (!sym && ncol(r_mat) == 1L) {
    cn <- colnames(r_mat)
    if (length(cn) >= 1L && nzchar(cn[[1]])) single_y_var <- cn[[1]]
  }

  if (sym) {
    x_label <- "Variable 1"
    y_label <- "Variable 2"
  } else {
    x_label <- "x-variable"
    y_label <- "y-variable"
  }

  columns <- list(
    list(key = "x", label = x_label, format = "text")
  )
  if (is.null(single_y_var)) {
    columns[[length(columns) + 1L]] <-
      list(key = "y", label = y_label, format = "text")
  }
  columns[[length(columns) + 1L]] <- list(key = "r",       label = "r",      format = "bounded")
  columns[[length(columns) + 1L]] <- list(key = "ci",      label = "95% CI", format = "ci")
  columns[[length(columns) + 1L]] <- list(
    key    = "p_value",
    label  = if (adjusted) "p (adjusted)" else "p",
    format = "pvalue"
  )
  columns[[length(columns) + 1L]] <- list(key = "n", label = "n", format = "integer")
  if (adjusted) {
    # Insert "p (raw)" before "p (adjusted)" — index depends on whether
    # the y column is present.
    insert_after <- if (is.null(single_y_var)) 3L else 2L  # after "ci"
    columns <- append(
      columns,
      list(list(key = "p_raw", label = "p (raw)", format = "pvalue")),
      after = insert_after
    )
  }

  fields <- list(
    table_type   = "correlations",
    columns      = columns,
    rows         = rows,
    adjust       = x$adjust %||% "none",
    sym          = sym
  )
  if (!is.null(single_y_var)) fields$single_y <- single_y_var

  n <- ms_corr_n_summary(x$n)
  if (!is.na(n)) fields$n <- n
  heatmap <- if (sym) ms_corr_test_heatmap_data(x, r_mat) else NULL
  forest <- if (sym) ms_corr_test_forest_data(x, rows) else NULL
  figure_data <- list()
  if (!is.null(heatmap)) figure_data$correlation_heatmap <- heatmap
  if (!is.null(forest)) figure_data$correlation_forest <- forest
  if (length(figure_data) == 0L) figure_data <- NULL

  ms_build_envelope(
    type       = "psych_corr_test",
    type_label = "Correlation table (psych::corr.test)",
    call       = trimws(gsub("\\s+", " ", ms_psych_call(x, .call, "psych::corr.test(...)"))),
    fields     = fields,
    raw_output = ms_capture_output(x),
    packages   = ms_packages_basic("psych"),
    card_kind  = "table",
    figure_data = figure_data
  )
}

#' @rdname mellio_payload
#' @export
mellio_payload.character <- function(x, title = NULL, ..., .call = NULL) {
  df <- console_key_value_table(x)
  if (is.null(df)) {
    stop(ms_character_hint(x), call. = FALSE)
  }

  call_str <- if (!is.null(.call)) .call else "captured console output"
  payload <- mellio_payload.data.frame(df, title = title, .call = call_str)
  payload$raw_output <- paste(x, collapse = "\n")
  payload$fields$source <- "captured_console"
  payload
}

#' @rdname mellio_payload
#' @export
mellio_payload.data.frame <- function(x, title = NULL, ..., .call = NULL) {
  detected <- ms_detect_corr_data_frame(x)
  if (is.null(detected)) detected <- ms_detect_anova_data_frame(x)
  if (is.null(detected)) detected <- ms_detect_effectsize_data_frame(x)
  if (is.null(detected)) detected <- ms_detect_regression_data_frame(x)
  if (is.null(detected)) detected <- ms_detect_descriptive_data_frame(x)

  if (!is.null(detected)) {
    call_str <- if (!is.null(.call)) .call else ms_deparse_call(match.call()$x)
    type_label <- if (!is.null(title) && nzchar(title)) {
      title
    } else {
      detected$type_label
    }
    return(ms_build_envelope(
      type       = detected$type,
      type_label = type_label,
      call       = trimws(gsub("\\s+", " ", call_str %||% "data.frame")),
      fields     = detected$fields,
      raw_output = ms_capture_output(x),
      card_kind  = "table"
    ))
  }

  stop(ms_data_frame_hint(x), call. = FALSE)
}

# ── psych helpers ────────────────────────────────────────────────────

ms_psych_call <- function(x, .call, fallback) {
  if (!is.null(.call)) return(.call)
  call_obj <- x$Call %||% x$call
  if (!is.null(call_obj)) return(ms_deparse_call(call_obj))
  fallback
}

ms_pick_numeric <- function(df, col) {
  if (is.null(df) || !is.data.frame(df) || !(col %in% names(df))) {
    return(NA_real_)
  }
  ms_safe_numeric(unname(df[[col]][1]))
}

ms_alpha_n <- function(x) {
  stats <- x$item.stats
  if (is.null(stats) || !is.data.frame(stats) || !("n" %in% names(stats))) {
    return(NA_real_)
  }
  n <- suppressWarnings(max(stats$n, na.rm = TRUE))
  if (!is.finite(n)) NA_real_ else as.integer(n)
}

ms_corr_test_rows <- function(x) {
  r <- x$r
  if (is.null(r)) return(list())
  r <- as.matrix(r)
  nr <- nrow(r)
  nc <- ncol(r)
  rn <- rownames(r) %||% paste0("x", seq_len(nr))
  cn <- colnames(r) %||% paste0("y", seq_len(nc))
  sym <- isTRUE(x$sym) && nr == nc
  adjusted <- !identical(x$adjust %||% "none", "none")

  pairs <- if (sym) {
    do.call(rbind, lapply(seq_len(max(0L, nr - 1L)), function(i) {
      if (i >= nc) return(NULL)
      cbind(i = i, j = seq.int(i + 1L, nc))
    }))
  } else {
    expand.grid(i = seq_len(nr), j = seq_len(nc))
  }
  if (is.null(pairs) || nrow(pairs) == 0L) return(list())
  pairs <- as.data.frame(pairs)

  rows <- vector("list", nrow(pairs))
  for (idx in seq_len(nrow(pairs))) {
    i <- pairs$i[idx]
    j <- pairs$j[idx]
    x_name <- rn[[i]]
    y_name <- cn[[j]]

    p_raw <- if (sym) {
      ms_matrix_pick(x$p, j, i)
    } else {
      ms_matrix_pick(x$p, i, j)
    }
    p_value <- if (adjusted) {
      if (sym) ms_matrix_pick(x$p, i, j) else ms_matrix_pick(x$p.adj, i, j)
    } else {
      p_raw
    }
    if (is.na(p_value)) p_value <- p_raw

    ci <- ms_corr_ci(x$ci, x_name, y_name)
    row <- list(
      x        = x_name,
      y        = y_name,
      r        = ms_safe_numeric(r[i, j]),
      p_value  = ms_safe_numeric(p_value),
      n        = ms_corr_n_pair(x$n, i, j)
    )
    if (!is.na(p_raw) && adjusted) row$p_raw <- ms_safe_numeric(p_raw)
    if (!is.null(ci)) {
      row$ci_lower <- ci[[1]]
      row$ci_upper <- ci[[2]]
    }
    rows[[idx]] <- row
  }
  rows
}

ms_matrix_pick <- function(x, i, j) {
  if (is.null(x)) return(NA_real_)
  m <- as.matrix(x)
  if (i > nrow(m) || j > ncol(m)) return(NA_real_)
  ms_safe_numeric(m[i, j])
}

ms_corr_n_pair <- function(n, i, j) {
  if (is.null(n)) return(NA_real_)
  if (length(n) == 1L) return(as.integer(n))
  m <- as.matrix(n)
  if (i > nrow(m) || j > ncol(m)) return(NA_real_)
  as.integer(ms_safe_numeric(m[i, j]))
}

ms_corr_n_summary <- function(n) {
  if (is.null(n)) return(NA_real_)
  vals <- suppressWarnings(as.numeric(n))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L) return(NA_real_)
  if (length(unique(vals)) == 1L) as.integer(vals[[1]]) else NA_real_
}

ms_corr_test_heatmap_data <- function(x, r_mat) {
  r_mat <- ms_square_numeric_matrix(r_mat)
  if (is.null(r_mat)) return(NULL)
  p_mat <- ms_corr_test_heatmap_p_matrix(x, nrow(r_mat))
  ms_correlation_heatmap_data(
    r = r_mat,
    p = p_mat,
    n = x$n,
    variables = rownames(r_mat) %||% colnames(r_mat),
    method = ms_corr_test_method(x),
    adjust = x$adjust %||% "none"
  )
}

ms_corr_test_forest_data <- function(x, rows) {
  ms_correlation_forest_data(
    pairs = rows,
    method = ms_corr_test_method(x),
    adjust = x$adjust %||% "none",
    n = x$n
  )
}

ms_corr_test_heatmap_p_matrix <- function(x, size) {
  p <- ms_square_numeric_matrix(x$p, n = size)
  if (is.null(p)) return(NULL)
  adjusted <- !identical(x$adjust %||% "none", "none")
  if (isTRUE(x$sym) && adjusted) {
    for (i in seq_len(size)) {
      for (j in seq_len(size)) {
        if (i > j && !is.na(p[j, i])) p[i, j] <- p[j, i]
      }
    }
  }
  p
}

ms_corr_test_method <- function(x) {
  call_obj <- x$Call %||% x$call
  if (is.call(call_obj) && !is.null(call_obj$method)) {
    method <- gsub("^['\"]|['\"]$", "", paste(deparse(call_obj$method), collapse = " "))
    method <- trimws(method)
    if (tolower(method) %in% c("pearson", "spearman", "kendall")) return(method)
    return(NULL)
  }
  "pearson"
}

ms_corr_ci <- function(ci, x_name, y_name) {
  if (is.null(ci) || !is.data.frame(ci)) return(NULL)
  if (!all(c("lower", "upper") %in% names(ci))) return(NULL)
  rn <- rownames(ci)
  keys <- c(paste0(x_name, "-", y_name), paste0(y_name, "-", x_name))
  hit <- match(keys, rn, nomatch = 0L)
  hit <- hit[hit > 0L]
  if (length(hit) == 0L) return(NULL)
  list(
    ms_safe_numeric(ci$lower[[hit[[1]]]]),
    ms_safe_numeric(ci$upper[[hit[[1]]]])
  )
}

# ── Custom data-frame intake ─────────────────────────────────────────

ms_detect_corr_data_frame <- function(x) {
  if (!is.data.frame(x) || nrow(x) == 0L) return(NULL)
  lower <- ms_df_names(x)

  r_idx <- ms_df_col(lower, c("r", "cor", "corr", "correlation",
                              "correlation_coefficient", "rho",
                              "spearman_rho", "tau", "kendall_tau",
                              "kendalls_tau"))
  p_idx <- ms_df_col(lower, ms_df_p_cols())
  if (is.na(r_idx) || is.na(p_idx)) return(NULL)
  if (!is.numeric(x[[r_idx]]) || !is.numeric(x[[p_idx]])) return(NULL)

  low_idx <- ms_df_col(lower, ms_df_low_cols())
  high_idx <- ms_df_col(lower, ms_df_high_cols())
  n_idx <- ms_df_col(lower, ms_df_n_cols())
  x_idx <- ms_df_col(lower, c("x", "x_var", "var1", "variable1",
                              "variable_1", "parameter1", "parameter_1",
                              "term1", "term_1", "predictor"))
  y_idx <- ms_df_col(lower, c("y", "y_var", "var2", "variable2",
                              "variable_2", "parameter2", "parameter_2",
                              "term2", "term_2", "outcome", "response"))

  label_candidates <- which(!seq_along(lower) %in% c(r_idx, p_idx, low_idx,
                                                     high_idx, n_idx, x_idx,
                                                     y_idx))
  text_candidates <- label_candidates[
    vapply(x[label_candidates], function(col) {
      is.character(col) || is.factor(col)
    }, logical(1))
  ]

  if (is.na(x_idx) && length(text_candidates) > 0L) {
    x_idx <- text_candidates[[1]]
  }
  if (is.na(y_idx) && length(text_candidates) > 1L) {
    y_idx <- text_candidates[[2]]
  }

  columns <- list(
    list(key = "x", label = if (!is.na(y_idx)) "Variable 1" else "Variable",
         format = "text")
  )
  if (!is.na(y_idx)) {
    columns <- c(columns, list(list(key = "y", label = "Variable 2",
                                    format = "text")))
  }
  columns <- c(columns, list(
    list(key = "r", label = "r", format = "bounded")
  ))
  if (!is.na(low_idx) && !is.na(high_idx)) {
    columns <- c(columns, list(list(key = "ci", label = "95% CI",
                                    format = "ci")))
  }
  p_label <- if (lower[[p_idx]] %in% ms_df_adjusted_p_cols()) "p (adjusted)" else "p"
  columns <- c(columns, list(list(key = "p_value", label = p_label,
                                  format = "pvalue")))
  if (!is.na(n_idx)) {
    columns <- c(columns, list(list(key = "n", label = "n",
                                    format = "integer")))
  }

  rows <- lapply(seq_len(nrow(x)), function(i) {
    row <- list(
      x       = if (!is.na(x_idx)) ms_df_text_value(x[[x_idx]][[i]]) else rownames(x)[[i]],
      r       = ms_safe_numeric(x[[r_idx]][[i]]),
      p_value = ms_safe_numeric(x[[p_idx]][[i]])
    )
    if (!is.na(y_idx)) row$y <- ms_df_text_value(x[[y_idx]][[i]])
    if (!is.na(low_idx) && !is.na(high_idx)) {
      row$ci_lower <- ms_safe_numeric(x[[low_idx]][[i]])
      row$ci_upper <- ms_safe_numeric(x[[high_idx]][[i]])
    }
    if (!is.na(n_idx)) row$n <- as.integer(ms_safe_numeric(x[[n_idx]][[i]]))
    row
  })

  fields <- list(
    table_type = "correlations",
    columns = columns,
    rows = rows,
    source = "custom_data_frame"
  )
  list(
    type = "custom_correlation_table",
    type_label = "Correlation table",
    fields = fields
  )
}

ms_detect_anova_data_frame <- function(x) {
  if (!is.data.frame(x) || nrow(x) == 0L) return(NULL)
  lower <- ms_df_names(x)

  f_idx <- ms_df_col(lower, c("f", "f_value", "f_stat", "f_statistic",
                              "f_ratio", "fvalue", "statistic"))
  p_idx <- ms_df_col(lower, c(ms_df_p_cols(), "pr_f", "pr_gt_f",
                              "pr_f_", "prob_f"))
  df1_idx <- ms_df_col(lower, c("df1", "df", "num_df", "numerator_df",
                                "df_num", "effect_df", "term_df"))

  if (is.na(f_idx) || is.na(p_idx) || is.na(df1_idx)) return(NULL)
  if (!is.numeric(x[[f_idx]]) || !is.numeric(x[[p_idx]]) ||
      !is.numeric(x[[df1_idx]])) return(NULL)

  df2_idx <- ms_df_col(lower, c("df2", "res_df", "residual_df",
                                "residuals_df", "den_df",
                                "denominator_df", "df_den", "error_df",
                                "df_error"))
  if (!is.na(df2_idx) && !is.numeric(x[[df2_idx]])) df2_idx <- NA_integer_

  term_idx <- ms_df_col(lower, c("term", "effect", "predictor", "factor",
                                 "source", "source_of_variation",
                                 "parameter"))
  outcome_idx <- ms_df_col(lower, c("outcome", "dependent",
                                    "dependent_variable", "dv", "response",
                                    "measure", "scale"))
  sumsq_idx <- ms_df_col(lower, c("sum_sq", "sum_squares", "sumsq",
                                  "sum_of_squares", "ss"))
  meansq_idx <- ms_df_col(lower, c("mean_sq", "mean_square", "meansq", "ms"))
  eta_idx <- ms_df_col(lower, c("eta_sq", "eta2", "etasq", "eta_squared",
                                "partial_eta_sq", "partial_eta2",
                                "partial_eta_squared"))

  if (!is.na(sumsq_idx) && !is.numeric(x[[sumsq_idx]])) sumsq_idx <- NA_integer_
  if (!is.na(meansq_idx) && !is.numeric(x[[meansq_idx]])) meansq_idx <- NA_integer_
  if (!is.na(eta_idx) && !is.numeric(x[[eta_idx]])) eta_idx <- NA_integer_

  f_vals <- suppressWarnings(as.numeric(x[[f_idx]]))
  p_vals <- suppressWarnings(as.numeric(x[[p_idx]]))
  keep <- which(is.finite(f_vals) & is.finite(p_vals))
  if (length(keep) == 0L) return(NULL)

  residual_df <- NA_real_
  if (is.na(df2_idx)) {
    residual_df <- ms_anova_residual_df(x, term_idx, df1_idx)
  }
  f_is_explicit <- lower[[f_idx]] %in% c("f", "f_value", "f_stat",
                                         "f_statistic", "f_ratio",
                                         "fvalue")
  has_anova_context <- f_is_explicit || !is.na(df2_idx) ||
    !is.na(sumsq_idx) || !is.na(meansq_idx) || is.finite(residual_df)
  if (!has_anova_context) return(NULL)

  if (is.na(term_idx) && !ms_has_informative_rownames(x)) return(NULL)

  columns <- list()
  if (!is.na(outcome_idx)) {
    columns <- c(columns, list(list(key = "outcome", label = "Outcome",
                                    format = "text")))
  }
  columns <- c(columns, list(list(key = "term", label = "Term",
                                  format = "text")))
  if (!is.na(sumsq_idx)) {
    columns <- c(columns, list(list(key = "sum_sq", label = "SS",
                                    format = "number")))
  }
  if (!is.na(meansq_idx)) {
    columns <- c(columns, list(list(key = "mean_sq", label = "MS",
                                    format = "number")))
  }
  columns <- c(columns, list(list(key = "df1", label = "df1",
                                  format = "integer")))
  if (!is.na(df2_idx) || is.finite(residual_df)) {
    columns <- c(columns, list(list(key = "df2", label = "df2",
                                    format = "integer")))
  }
  columns <- c(columns, list(
    list(key = "f", label = "F", format = "statistic"),
    list(key = "p_value", label = "p", format = "pvalue")
  ))
  if (!is.na(eta_idx)) {
    columns <- c(columns, list(list(key = "eta_sq", label = "eta^2",
                                    format = "bounded")))
  }

  rows <- lapply(keep, function(i) {
    row <- list(
      term = if (!is.na(term_idx)) ms_df_text_value(x[[term_idx]][[i]]) else rownames(x)[[i]],
      df1 = as.integer(ms_safe_numeric(x[[df1_idx]][[i]])),
      f = ms_safe_numeric(x[[f_idx]][[i]]),
      p_value = ms_safe_numeric(x[[p_idx]][[i]])
    )
    if (!is.na(outcome_idx)) {
      row$outcome <- ms_df_text_value(x[[outcome_idx]][[i]])
    }
    if (!is.na(sumsq_idx)) row$sum_sq <- ms_safe_numeric(x[[sumsq_idx]][[i]])
    if (!is.na(meansq_idx)) row$mean_sq <- ms_safe_numeric(x[[meansq_idx]][[i]])
    if (!is.na(df2_idx)) {
      row$df2 <- as.integer(ms_safe_numeric(x[[df2_idx]][[i]]))
    } else if (is.finite(residual_df)) {
      row$df2 <- as.integer(residual_df)
    }
    if (!is.na(eta_idx)) row$eta_sq <- ms_safe_numeric(x[[eta_idx]][[i]])
    row
  })

  fields <- list(
    table_type = "anova",
    columns = columns,
    rows = rows,
    source = "custom_data_frame"
  )
  list(
    type = "custom_anova_table",
    type_label = "ANOVA table",
    fields = fields
  )
}

ms_detect_regression_data_frame <- function(x) {
  if (!is.data.frame(x) || nrow(x) == 0L) return(NULL)
  lower <- ms_df_names(x)

  term_idx <- ms_df_col(lower, c("term", "parameter", "predictor",
                                 "variable", "effect"))
  estimate_idx <- ms_df_col(lower, c("estimate", "est", "b", "beta",
                                     "coefficient", "coef", "odds_ratio",
                                     "or"))
  p_idx <- ms_df_col(lower, ms_df_p_cols())
  if (is.na(estimate_idx) || is.na(p_idx)) return(NULL)
  if (!is.numeric(x[[estimate_idx]]) || !is.numeric(x[[p_idx]])) return(NULL)
  if (is.na(term_idx) && !ms_has_informative_rownames(x)) return(NULL)

  se_idx <- ms_df_col(lower, c("std_error", "standard_error", "se",
                               "std_err", "stderr"))
  stat_idx <- ms_df_col(lower, c("statistic", "t", "t_value", "t_stat",
                                 "t_statistic", "z", "z_value", "z_stat",
                                 "z_statistic", "f", "f_value"))
  low_idx <- ms_df_col(lower, ms_df_low_cols())
  high_idx <- ms_df_col(lower, ms_df_high_cols())

  columns <- list(
    list(key = "term", label = "Term", format = "text"),
    list(key = "estimate", label = "Estimate", format = "number")
  )
  if (!is.na(se_idx) && is.numeric(x[[se_idx]])) {
    columns <- c(columns, list(list(key = "std_error", label = "SE",
                                    format = "number")))
  } else {
    se_idx <- NA_integer_
  }
  if (!is.na(stat_idx) && is.numeric(x[[stat_idx]])) {
    stat_label <- ms_df_stat_label(lower[[stat_idx]])
    columns <- c(columns, list(list(key = "statistic",
                                    label = stat_label,
                                    format = "statistic")))
  } else {
    stat_idx <- NA_integer_
    stat_label <- NULL
  }
  if (!is.na(low_idx) && !is.na(high_idx) &&
      is.numeric(x[[low_idx]]) && is.numeric(x[[high_idx]])) {
    columns <- c(columns, list(list(key = "ci", label = "95% CI",
                                    format = "ci")))
  } else {
    low_idx <- high_idx <- NA_integer_
  }

  p_label <- if (lower[[p_idx]] %in% ms_df_adjusted_p_cols()) "p (adjusted)" else "p"
  columns <- c(columns, list(list(key = "p_value", label = p_label,
                                  format = "pvalue")))

  rows <- lapply(seq_len(nrow(x)), function(i) {
    row <- list(
      term = if (!is.na(term_idx)) ms_df_text_value(x[[term_idx]][[i]]) else rownames(x)[[i]],
      estimate = ms_safe_numeric(x[[estimate_idx]][[i]]),
      p_value = ms_safe_numeric(x[[p_idx]][[i]])
    )
    if (!is.na(se_idx)) row$std_error <- ms_safe_numeric(x[[se_idx]][[i]])
    if (!is.na(stat_idx)) {
      row$statistic <- ms_safe_numeric(x[[stat_idx]][[i]])
      row$statistic_label <- stat_label
    }
    if (!is.na(low_idx) && !is.na(high_idx)) {
      row$ci_lower <- ms_safe_numeric(x[[low_idx]][[i]])
      row$ci_upper <- ms_safe_numeric(x[[high_idx]][[i]])
    }
    row
  })

  fields <- list(
    table_type = "coefficients",
    columns = columns,
    rows = rows,
    source = "custom_data_frame"
  )
  if (!is.null(stat_label)) fields$statistic_label <- stat_label
  list(
    type = "custom_regression_table",
    type_label = "Coefficient table",
    fields = fields
  )
}

ms_detect_descriptive_data_frame <- function(x) {
  if (!is.data.frame(x) || nrow(x) == 0L) return(NULL)
  lower <- ms_df_names(x)

  variable_idx <- ms_df_col(lower, c("variable", "var", "name", "term",
                                     "parameter", "measure", "scale"))
  mean_idx <- ms_df_col(lower, c("mean", "m", "average", "avg"))
  sd_idx <- ms_df_col(lower, c("sd", "std_dev", "standard_deviation"))
  if (is.na(mean_idx) || is.na(sd_idx)) return(NULL)
  if (!is.numeric(x[[mean_idx]]) || !is.numeric(x[[sd_idx]])) return(NULL)
  if (is.na(variable_idx) && !ms_has_informative_rownames(x)) return(NULL)

  n_idx <- ms_df_col(lower, ms_df_n_cols())
  median_idx <- ms_df_col(lower, c("median", "mdn"))
  min_idx <- ms_df_col(lower, c("min", "minimum"))
  max_idx <- ms_df_col(lower, c("max", "maximum"))
  skew_idx <- ms_df_col(lower, c("skew", "skewness"))
  kurtosis_idx <- ms_df_col(lower, c("kurtosis", "kurt", "kurtosi"))

  columns <- list(
    list(key = "variable", label = "Variable", format = "text"),
    list(key = "mean", label = "M", format = "number"),
    list(key = "sd", label = "SD", format = "number")
  )
  if (!is.na(n_idx) && is.numeric(x[[n_idx]])) {
    columns <- c(columns, list(list(key = "n", label = "n",
                                    format = "integer")))
  } else {
    n_idx <- NA_integer_
  }
  if (!is.na(median_idx) && is.numeric(x[[median_idx]])) {
    columns <- c(columns, list(list(key = "median", label = "Mdn",
                                    format = "number")))
  } else {
    median_idx <- NA_integer_
  }
  if (!is.na(min_idx) && is.numeric(x[[min_idx]])) {
    columns <- c(columns, list(list(key = "min", label = "Min",
                                    format = "number")))
  } else {
    min_idx <- NA_integer_
  }
  if (!is.na(max_idx) && is.numeric(x[[max_idx]])) {
    columns <- c(columns, list(list(key = "max", label = "Max",
                                    format = "number")))
  } else {
    max_idx <- NA_integer_
  }
  if (!is.na(skew_idx) && is.numeric(x[[skew_idx]])) {
    columns <- c(columns, list(list(key = "skew", label = "Skew",
                                    format = "number")))
  } else {
    skew_idx <- NA_integer_
  }
  if (!is.na(kurtosis_idx) && is.numeric(x[[kurtosis_idx]])) {
    columns <- c(columns, list(list(key = "kurtosis", label = "Kurtosis",
                                    format = "number")))
  } else {
    kurtosis_idx <- NA_integer_
  }

  rows <- lapply(seq_len(nrow(x)), function(i) {
    row <- list(
      variable = if (!is.na(variable_idx)) ms_df_text_value(x[[variable_idx]][[i]]) else rownames(x)[[i]],
      mean = ms_safe_numeric(x[[mean_idx]][[i]]),
      sd = ms_safe_numeric(x[[sd_idx]][[i]])
    )
    if (!is.na(n_idx)) row$n <- as.integer(ms_safe_numeric(x[[n_idx]][[i]]))
    if (!is.na(median_idx)) row$median <- ms_safe_numeric(x[[median_idx]][[i]])
    if (!is.na(min_idx)) row$min <- ms_safe_numeric(x[[min_idx]][[i]])
    if (!is.na(max_idx)) row$max <- ms_safe_numeric(x[[max_idx]][[i]])
    if (!is.na(skew_idx)) row$skew <- ms_safe_numeric(x[[skew_idx]][[i]])
    if (!is.na(kurtosis_idx)) row$kurtosis <- ms_safe_numeric(x[[kurtosis_idx]][[i]])
    row
  })

  fields <- list(
    table_type = "descriptives",
    columns = columns,
    rows = rows,
    source = "custom_data_frame"
  )
  list(
    type = "custom_descriptive_table",
    type_label = "Descriptive statistics table",
    fields = fields
  )
}

ms_data_frame_hint <- function(x) {
  paste(
    "`x` is a data.frame, but Mellio could not recognize it as a",
    "statistical result table.",
    "Supported table shapes include correlations (r, p), ANOVA",
    "(term, df1, df2, F, p), coefficients (term, estimate, p.value),",
    "effect sizes (Cohen's d, eta-squared, CI), and descriptives (M, SD, n).",
    "For manuscript-ready tables, use mellio_open(melliotab(x)).",
    sep = "\n"
  )
}

ms_character_hint <- function(x) {
  txt <- paste(x, collapse = "\n")
  if (grepl("\\bANOVA\\b|\\bDf\\b.*\\bF value\\b|Pr\\s*\\(>F\\)",
            txt, ignore.case = TRUE)) {
    return(paste(
      "Captured text was not recognized as a Stats result.",
      "For ANOVA, pass `anova(fit)` or a data.frame with term, df1, df2, F, and p.",
      sep = "\n"
    ))
  }

  paste(
    "Captured text was not recognized as a Stats result.",
    "Pass a supported R object or a result-shaped data.frame.",
    sep = "\n"
  )
}

ms_anova_residual_df <- function(x, term_idx, df_idx) {
  if (is.na(term_idx)) {
    rn <- rownames(x)
    if (is.null(rn)) return(NA_real_)
    hit <- grep("residual|error", rn, ignore.case = TRUE)
  } else {
    terms <- as.character(x[[term_idx]])
    hit <- grep("residual|error", terms, ignore.case = TRUE)
  }
  if (length(hit) == 0L) return(NA_real_)
  val <- ms_safe_numeric(x[[df_idx]][[hit[[length(hit)]]]])
  if (is.na(val)) NA_real_ else val
}

ms_df_names <- function(x) {
  nms <- names(x)
  nms <- tolower(gsub("[^A-Za-z0-9]+", "_", nms))
  gsub("^_+|_+$", "", nms)
}

ms_df_col <- function(lower_names, candidates) {
  hit <- which(lower_names %in% candidates)
  if (length(hit) == 0L) NA_integer_ else hit[[1]]
}

ms_df_p_cols <- function() {
  c("p", "p_value", "p_val", "pvalue", "p_adjust", "p_adjusted",
    "p_adj", "padj", "p_fdr", "p_holm")
}

ms_df_adjusted_p_cols <- function() {
  c("p_adjust", "p_adjusted", "p_adj", "padj", "p_fdr", "p_holm")
}

ms_df_low_cols <- function() {
  c("ci_low", "ci_lower", "lower", "lower_ci", "conf_low",
    "conf_lower", "conf_int_low", "confint_low", "ci_low_95",
    "lower_95")
}

ms_df_high_cols <- function() {
  c("ci_high", "ci_upper", "upper", "upper_ci", "conf_high",
    "conf_upper", "conf_int_high", "confint_high", "ci_high_95",
    "upper_95")
}

ms_df_n_cols <- function() {
  c("n", "n_obs", "nobs", "sample_size", "observations", "obs",
    "count", "valid_n")
}

ms_df_text_value <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x)) return(NA_character_)
  as.character(x)
}

ms_has_informative_rownames <- function(x) {
  rn <- rownames(x)
  !is.null(rn) && !identical(rn, as.character(seq_len(nrow(x))))
}

ms_df_stat_label <- function(name) {
  if (name %in% c("t", "t_value", "t_stat", "t_statistic")) return("t")
  if (name %in% c("z", "z_value", "z_stat", "z_statistic")) return("z")
  if (name %in% c("f", "f_value")) return("F")
  "Statistic"
}

ms_function_hint <- function(fn, .call = NULL) {
  fml <- formals(fn)
  required <- names(fml)[vapply(fml, function(arg) {
    identical(arg, quote(expr = ))
  }, logical(1))]
  required <- setdiff(required, "...")
  body_txt <- paste(deparse(body(fn), width.cutoff = 500L), collapse = "\n")
  looks_corr <- grepl("\\bcor\\.test\\s*\\(", body_txt)
  looks_alpha <- grepl("psych::alpha\\s*\\(|\\balpha\\s*\\(", body_txt)
  fn_label <- if (!is.null(.call) && nzchar(.call)) {
    .call
  } else {
    "x"
  }

  msg <- c(
    paste0("`", fn_label, "` is a function, not a result.")
  )
  if (length(required) > 0L) {
    msg <- c(msg, paste0(
      "It still needs ",
      if (length(required) == 1L) "argument " else "arguments ",
      paste(required, collapse = ", "), "."
    ))
  }
  if (looks_corr) {
    msg <- c(msg,
      "For correlation tables, pass psych::corr.test() or a result",
      "table with r and p columns to mellio_open()."
    )
  } else if (looks_alpha) {
    msg <- c(msg,
      "For reliability, run psych::alpha() and pass that object to",
      "mellio_open()."
    )
  } else {
    msg <- c(msg,
      "Run the function first, then pass the returned result to mellio_open()."
    )
  }
  paste(msg, collapse = "\n")
}
