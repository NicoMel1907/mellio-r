# Convert a mellio_payload object to a melliotab object where a tabular
# projection is unambiguous.
melliotab_from_payload <- function(payload, section = NULL,
                                   style = "apa7", title = NULL,
                                   number = NULL, note = NULL,
                                   source = NULL,
                                   decimals = 2L, p_decimals = 3L,
                                   what = NULL, ...) {
  section <- mellio_resolve_section(section = section, what = what)
  card_kind <- payload$card_kind %||% "inline"

  if (identical(payload$type, "hierarchical_regression_comparison")) {
    df <- payload_hierarchical_comparison_to_df(payload)
    table_title <- payload$type_label
  } else if (is.list(payload$fields$rows) && length(payload$fields$rows) &&
             is.list(payload$fields$columns) && length(payload$fields$columns)) {
    df <- payload_table_to_df(payload)
    table_title <- payload$type_label
  } else {
    projection <- switch(card_kind,
      table = list(data = payload_table_to_df(payload), title = payload$type_label),
      inline = list(data = payload_inline_to_df(payload), title = payload$type_label),
      structural = payload_structural_to_df(payload, section = section),
      cli::cli_abort(c(
        "Cannot convert {.val {card_kind}} payload to a melliotab.",
        "i" = "Card kinds {.val raw_text} and {.val unsupported} have no tabular projection."
      ))
    )
    df <- projection$data
    table_title <- projection$title %||% payload$type_label
  }

  if (is.null(title)) title <- table_title
  if (is.null(note)) note <- payload$fields$note %||% payload$fields$table_note
  if (is.null(source)) source <- payload$call

  result <- melliotab.data.frame(
    df,
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

mellio_resolve_section <- function(section = NULL, what = NULL,
                                   default = NULL, choices = NULL) {
  if (!is.null(section) && !is.null(what) &&
      !identical(as.character(section), as.character(what))) {
    cli::cli_abort(c(
      "{.arg section} and {.arg what} specify different table sections.",
      "i" = "Use {.arg section}; {.arg what} is kept only as a compatibility alias."
    ))
  }

  value <- section %||% what %||% default
  if (is.null(value)) return(NULL)
  value <- as.character(value[[1L]])
  if (!is.null(choices)) value <- match.arg(value, choices)
  value
}

payload_hierarchical_comparison_to_df <- function(payload) {
  fields <- payload$fields %||% list()
  models <- fields$models
  if (!is.list(models) || length(models) == 0L) {
    cli::cli_abort("Hierarchical comparison payload has no model rows.")
  }

  comparisons <- fields$comparisons
  if (!is.list(comparisons)) comparisons <- list()
  by_model <- list()
  for (change in comparisons) {
    key <- as.character(payload_scalar(change$to))
    if (length(key) == 1L && !is.na(key) && nzchar(key)) {
      by_model[[key]] <- change
    }
  }

  rows <- lapply(seq_along(models), function(i) {
    model <- models[[i]] %||% list()
    model_id <- payload_scalar(model$model)
    if (length(model_id) != 1L || is.na(model_id) || !nzchar(as.character(model_id))) {
      model_id <- i
    }
    change <- by_model[[as.character(model_id)]] %||% list()

    predictors <- payload_scalar(model$terms)
    if (length(predictors) == 1L && is.na(predictors)) predictors <- ""
    terms_entered <- payload_scalar(change$added_terms)
    if (length(terms_entered) == 1L && is.na(terms_entered)) {
      terms_entered <- predictors
    }

    df_pair <- ""
    df1 <- payload_scalar(change$df_change)
    df2 <- payload_scalar(change$residual_df)
    if (!(length(df1) == 1L && is.na(df1)) && !(length(df2) == 1L && is.na(df2))) {
      df_pair <- paste0("(", paste(payload_format_df(df1), payload_format_df(df2), sep = ", "), ")")
    }

    list(
      Model = payload_scalar(model$label %||% paste("Model", model_id)),
      Predictors = predictors,
      terms_entered = terms_entered,
      r_squared = payload_scalar(model$r_squared),
      adj_r_squared = payload_scalar(model$adj_r_squared),
      r_squared_change = payload_scalar(change$r_squared_change),
      f_change = payload_scalar(change$f_change),
      df = df_pair,
      p = payload_scalar(change$p_change),
      n = payload_scalar(model$n)
    )
  })

  row_frames <- lapply(rows, function(row) {
    as.data.frame(row, check.names = FALSE, stringsAsFactors = FALSE)
  })
  df <- as.data.frame(do.call(rbind, row_frames),
                      check.names = FALSE, stringsAsFactors = FALSE)
  # Apply Unicode column names AFTER the data.frame is built. Doing it
  # before (inside the list passed to data.frame) triggers an attempted
  # native-encoding conversion in C/POSIX locale sessions, which replaces
  # characters like \u00B2 and \u0394 with literal escape text (`R<U+00B2>`). Post-
  # construction `names(df) <-` doesn't trigger that path, so the bytes
  # stay correct in any locale. The web app, HTML/LaTeX exports, and
  # everything else further downstream see the same UTF-8 bytes.
  names(df) <- c(
    "Model", "Predictors", "Terms entered", "R\u00B2", "Adjusted R\u00B2",
    "\u0394R\u00B2", "F change", "df", "p", "n"
  )
  df
}

payload_format_df <- function(value) {
  num <- suppressWarnings(as.numeric(value))
  if (!is.na(num)) {
    if (abs(num - round(num)) < 1e-6) return(as.character(as.integer(round(num))))
    return(sub("\\.?0+$", "", formatC(num, digits = 2L, format = "f")))
  }
  as.character(value)
}

payload_table_to_df <- function(payload) {
  fields <- payload$fields %||% list()
  rows <- fields$rows
  if (!is.list(rows) || length(rows) == 0L) {
    cli::cli_abort("Payload has card_kind = {.val table} but no rows.")
  }

  columns <- fields$columns
  if (!is.list(columns) || length(columns) == 0L) {
    columns <- payload_columns_from_rows(rows)
  }

  keys <- vapply(columns, function(col) as.character(col$key %||% ""), character(1))
  keep <- nzchar(keys)
  keys <- keys[keep]
  columns <- columns[keep]
  if (length(keys) == 0L) {
    cli::cli_abort("Payload has card_kind = {.val table} but no usable columns.")
  }

  out <- lapply(seq_along(keys), function(i) {
    values <- lapply(rows, function(row) {
      payload_cell_value(row, keys[[i]], columns[[i]])
    })
    payload_simplify_column(values)
  })

  labels <- vapply(columns, function(col) {
    as.character(col$label %||% col$key %||% "")
  }, character(1))
  labels[!nzchar(labels)] <- keys[!nzchar(labels)]
  names(out) <- make.unique(labels)

  as.data.frame(out, check.names = FALSE, stringsAsFactors = FALSE)
}

payload_inline_to_df <- function(payload) {
  if (!payload$type %in% payload_inline_projection_types()) {
    cli::cli_abort(c(
      "{.fn melliotab} cannot yet project {.val {payload$type}} inline payloads.",
      "i" = "Use {.fn mellio_open} to send this to the Mellio web app instead.",
      "i" = "A class-specific melliotab projection can be added once the table shape is clear."
    ))
  }

  fields <- payload$fields %||% list()
  row <- list()

  row <- payload_add_scalar(row, "Term", fields$term)
  row <- payload_add_scalar(row, "Comparison", fields$comparison)
  row <- payload_add_scalar(row, "Outcome", fields$outcome)
  row <- payload_add_scalar(row, "Predictor", fields$predictor)
  row <- payload_add_scalar(row, "Model fit", fields$model_fit)

  stat <- fields$statistic
  if (is.list(stat)) {
    stat_label <- payload_inline_label(stat$name %||% "Statistic")
    row <- payload_add_scalar(row, stat_label, stat$value)
    row <- payload_add_df(row, stat$df)
    row <- payload_add_ci(row, stat$ci, fields$conf_level)
  }

  estimate <- fields$estimate
  if (is.list(estimate)) {
    est_label <- payload_inline_label(estimate$name %||% "Estimate")
    row <- payload_add_scalar(row, est_label, estimate$value)
    row <- payload_add_ci(row, estimate$ci, fields$conf_level, prefix = "Estimate ")
  }

  effect <- fields$effect %||% fields$effect_size
  if (is.list(effect)) {
    effect_label <- payload_inline_label(effect$name %||% "Effect")
    row <- payload_add_scalar(row, effect_label, effect$value)
  }

  row <- payload_add_scalar(row, "p", fields$p_value)
  row <- payload_add_scalar(row, "N", fields$n %||% fields$sample_size)
  row <- payload_add_scalar(row, "df change", fields$df_change)
  row <- payload_add_scalar(row, "Residual df", fields$residual_df)
  row <- payload_add_scalar(row, "df effect", fields$df_effect)
  row <- payload_add_scalar(row, "df error", fields$df_error)
  row <- payload_add_scalar(row, "Sum Sq", fields$sum_sq %||% fields$sum_of_sq)
  row <- payload_add_scalar(row, "Mean Sq", fields$mean_sq)
  row <- payload_add_scalar(row, "Residual Sum Sq", fields$residual_sum_sq)
  row <- payload_add_scalar(row, "Residual Mean Sq", fields$residual_mean_sq)

  if (length(row) == 0L) {
    cli::cli_abort("Payload has card_kind = {.val inline} but no reportable fields.")
  }

  as.data.frame(row, check.names = FALSE, stringsAsFactors = FALSE)
}

payload_inline_projection_types <- function() {
  c(
    "welch_t_test",
    "students_t_test",
    "paired_t_test",
    "one_sample_t_test",
    "pearson_correlation",
    "spearman_correlation",
    "kendall_correlation",
    "wilcoxon_rank_sum",
    "wilcoxon_signed_rank",
    "chi_squared_test",
    "fisher_exact_test",
    "htest_other",
    "anova_single_model",
    "anova_model_comparison"
  )
}

payload_structural_to_df <- function(payload, section = NULL) {
  tables <- payload_structural_table_options(payload)
  if (length(tables) == 0L) {
    cli::cli_abort(c(
      "{.fn melliotab} cannot find table sections in this {.val {payload$type}} payload.",
      "i" = "Use {.fn mellio_open} to inspect the result in the Mellio web app."
    ))
  }

  if (is.null(section)) {
    choices <- vapply(tables, function(tab) {
      paste0("{.val ", tab$key, "} (", tab$label, ", ",
             payload_row_count_label(tab$n), ")")
    }, character(1))
    cli::cli_abort(c(
      "{.fn melliotab} can make several tables from this result.",
      "i" = "Choose one with {.arg section}:",
      stats::setNames(choices, rep("*", length(choices))),
      "i" = 'Example: {.code melliotab(x, section = "loadings")}'
    ))
  }

  key <- payload_normalize_section(section)
  aliases <- c(
    all = "parameters",
    parameter = "parameters",
    parameter_estimates = "parameters",
    fit_indices = "fit",
    model_fit = "fit",
    regression = "paths",
    regressions = "paths",
    indirect = "defined",
    effect = "defined",
    effects = "defined",
    defined_effects = "defined",
    covariance = "covariances",
    residual_covariances = "covariances",
    variance = "variances",
    residual_variances = "variances",
    r2 = "r_squared",
    rsquared = "r_squared",
    r_square = "r_squared",
    r_squares = "r_squared",
    reliability_estimates = "reliability",
    modification = "modification_indices",
    mi = "modification_indices",
    modification_index = "modification_indices",
    observed = "observed_variables"
  )
  if (key %in% names(aliases)) key <- unname(aliases[[key]])
  table_keys <- vapply(tables, function(tab) tab$key, character(1))
  if (!key %in% table_keys) {
    cli::cli_abort(c(
      "Unknown table section {.val {section}}.",
      "i" = "Available sections: {.val {table_keys}}."
    ))
  }

  selected <- tables[[match(key, table_keys)]]
  list(
    data = payload_structural_section_df(payload, key),
    title = selected$label
  )
}

payload_row_count_label <- function(n) {
  n <- as.integer(n %||% 0L)
  paste0(n, " row", if (identical(n, 1L)) "" else "s")
}

payload_structural_table_options <- function(payload) {
  fields <- payload$fields %||% list()
  rz <- fields$report_zone %||% list()
  iz <- fields$inspection_zone %||% list()
  params <- iz$parameters %||% list()
  out <- list()

  add <- function(key, label, n) {
    if (is.null(n) || is.na(n) || n < 1L) return()
    out[[length(out) + 1L]] <<- list(key = key, label = label, n = as.integer(n))
  }

  add("fit", "Model fit indices", length(rz$fit_indices %||% list()))
  add("parameters", "Parameter estimates", length(params))
  add("loadings", "Factor loadings", length(payload_filter_parameters(params, "=~")))
  add("paths", "Regression paths", length(payload_filter_parameters(params, "~")))
  add("defined", "Defined and indirect effects", length(payload_filter_parameters(params, ":=")))
  add("covariances", "Covariances", length(payload_filter_parameters(params, "~~", same = FALSE)))
  add("variances", "Variances", length(payload_filter_parameters(params, "~~", same = TRUE)))
  add("r_squared", "R-squared values", length(fields$structural_r_squared %||% list()))
  add("reliability", "Reliability estimates", length(iz$reliability %||% list()))
  add("modification_indices", "Modification indices", length(iz$modification_indices %||% list()))
  add("observed_variables", "Observed variable summaries", length(fields$observed_variables %||% list()))

  out
}

payload_normalize_section <- function(section) {
  section <- as.character(section[[1L]] %||% "")
  tolower(gsub("[^a-z0-9]+", "_", trimws(section)))
}

payload_structural_section_df <- function(payload, key) {
  fields <- payload$fields %||% list()
  rz <- fields$report_zone %||% list()
  iz <- fields$inspection_zone %||% list()
  params <- iz$parameters %||% list()

  switch(key,
    fit = payload_fit_indices_to_df(rz$fit_indices %||% list()),
    parameters = payload_parameters_to_df(params),
    loadings = payload_parameters_to_df(payload_filter_parameters(params, "=~")),
    paths = payload_parameters_to_df(payload_filter_parameters(params, "~")),
    defined = payload_parameters_to_df(payload_filter_parameters(params, ":=")),
    covariances = payload_parameters_to_df(payload_filter_parameters(params, "~~", same = FALSE)),
    variances = payload_parameters_to_df(payload_filter_parameters(params, "~~", same = TRUE)),
    r_squared = payload_named_rows_to_df(fields$structural_r_squared %||% list(),
      c(variable = "Variable", r_squared = "R\u00B2")),
    reliability = payload_named_rows_to_df(iz$reliability %||% list(),
      c(factor = "Factor", omega = "\u03C9", ave = "AVE", n_indicators = "Indicators")),
    modification_indices = payload_named_rows_to_df(iz$modification_indices %||% list(),
      c(lhs = "lhs", op = "op", rhs = "rhs", mi = "MI", epc = "EPC",
        sepc.all = "Std. EPC")),
    observed_variables = payload_named_rows_to_df(fields$observed_variables %||% list(),
      c(variable = "Variable", mean = "M", sd = "SD", n = "N")),
    cli::cli_abort("Unknown structural table section {.val {key}}.")
  )
}

payload_filter_parameters <- function(params, op, same = NULL) {
  rows <- Filter(function(row) is.list(row) && identical(as.character(row$op %||% ""), op), params)
  if (identical(op, "~~") && !is.null(same)) {
    rows <- Filter(function(row) {
      lhs <- as.character(row$lhs %||% "")
      rhs <- as.character(row$rhs %||% "")
      identical(lhs == rhs, same)
    }, rows)
  }
  rows
}

payload_fit_indices_to_df <- function(rows) {
  if (!is.list(rows) || length(rows) == 0L) {
    cli::cli_abort("This payload has no fit-index table.")
  }
  out <- lapply(rows, function(row) {
    ci <- row$ci
    list(
      "Fit index" = payload_scalar(row$name),
      Value = payload_scalar(row$value),
      df = payload_scalar(row$df),
      p = payload_scalar(row$p_value),
      "CI level" = payload_scalar(row$ci_level),
      "Lower CI" = if (!is.null(ci) && length(ci) >= 2L) payload_scalar(ci[[1L]]) else NA,
      "Upper CI" = if (!is.null(ci) && length(ci) >= 2L) payload_scalar(ci[[2L]]) else NA,
      Scaled = if (isTRUE(row$scaled)) "yes" else NA
    )
  })
  payload_row_list_to_df(out)
}

payload_parameters_to_df <- function(rows) {
  if (!is.list(rows) || length(rows) == 0L) {
    cli::cli_abort("This payload has no rows for the requested parameter section.")
  }
  out <- lapply(rows, function(row) {
    list(
      Parameter = payload_parameter_label(row),
      lhs = payload_scalar(row$lhs),
      op = payload_scalar(row$op),
      rhs = payload_scalar(row$rhs),
      Estimate = payload_scalar(row$estimate),
      SE = payload_scalar(row$std_error),
      z = payload_scalar(row$statistic),
      p = payload_scalar(row$p_value),
      "Lower CI" = payload_scalar(row$ci_lower),
      "Upper CI" = payload_scalar(row$ci_upper),
      "Std. Estimate" = payload_scalar(row$std_estimate),
      Group = payload_scalar(row$group),
      Label = payload_scalar(row$label)
    )
  })
  payload_row_list_to_df(out)
}

payload_parameter_label <- function(row) {
  op <- as.character(row$op %||% "")
  lhs <- as.character(row$lhs %||% "")
  rhs <- as.character(row$rhs %||% "")
  if (identical(op, ":=")) return(lhs)
  sep <- switch(op,
    "=~" = " -> ",
    "~" = " ~ ",
    "~~" = " ~~ ",
    paste0(" ", op, " ")
  )
  paste0(lhs, sep, rhs)
}

payload_named_rows_to_df <- function(rows, columns) {
  if (!is.list(rows) || length(rows) == 0L) {
    cli::cli_abort("This payload has no rows for the requested section.")
  }
  out <- lapply(rows, function(row) {
    values <- lapply(names(columns), function(key) payload_scalar(row[[key]]))
    names(values) <- unname(columns)
    values
  })
  payload_row_list_to_df(out)
}

payload_row_list_to_df <- function(rows) {
  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  compacted <- lapply(rows, function(row) {
    row <- row[all_names]
    row[vapply(row, is.null, logical(1))] <- NA
    row
  })
  df <- as.data.frame(do.call(rbind, lapply(compacted, as.data.frame,
    check.names = FALSE, stringsAsFactors = FALSE)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  keep <- vapply(df, function(col) !all(is.na(col) | !nzchar(as.character(col))), logical(1))
  df[, keep, drop = FALSE]
}

payload_columns_from_rows <- function(rows) {
  keys <- unique(unlist(lapply(rows, names), use.names = FALSE))
  keys <- keys[nzchar(keys %||% "")]
  lapply(keys, function(key) list(key = key, label = key, format = "text"))
}

payload_cell_value <- function(row, key, column = NULL) {
  fmt <- as.character(column$format %||% "")
  if (identical(fmt, "ci")) {
    lo <- row$ci_lower
    hi <- row$ci_upper
    if (is.null(lo) || is.null(hi) || is.na(lo) || is.na(hi)) return(NA_character_)
    return(paste0("[", lo, ", ", hi, "]"))
  }

  value <- row[[key]]
  payload_scalar(value)
}

payload_scalar <- function(value) {
  if (is.null(value) || length(value) == 0L) return(NA)
  if (inherits(value, "AsIs")) value <- unclass(value)
  if (is.factor(value)) value <- as.character(value)
  if (inherits(value, c("Date", "POSIXct", "POSIXlt"))) value <- as.character(value)
  if (is.atomic(value) && length(value) == 1L) return(value[[1L]])
  if (is.atomic(value)) return(paste(as.character(value), collapse = ", "))
  paste(as.character(unlist(value, use.names = FALSE)), collapse = ", ")
}

payload_simplify_column <- function(values) {
  non_na <- values[!vapply(values, function(value) {
    length(value) == 1L && is.na(value)
  }, logical(1))]

  if (length(non_na) > 0L &&
      all(vapply(non_na, function(value) is.numeric(value) || is.integer(value), logical(1)))) {
    return(as.numeric(vapply(values, function(value) {
      if (length(value) == 1L && is.na(value)) return(NA_real_)
      as.numeric(value)
    }, numeric(1))))
  }

  vapply(values, function(value) {
    if (length(value) == 1L && is.na(value)) return(NA_character_)
    as.character(value)
  }, character(1))
}

payload_add_scalar <- function(row, label, value) {
  value <- payload_scalar(value)
  if (length(value) == 1L && is.na(value)) return(row)
  row[[label]] <- value
  row
}

payload_add_df <- function(row, df) {
  if (is.null(df) || length(df) == 0L) return(row)
  df <- payload_scalar(df)
  if (length(df) == 1L && is.na(df)) return(row)
  row[["df"]] <- df
  row
}

payload_add_ci <- function(row, ci, conf_level = NULL, prefix = "") {
  if (is.null(ci) || length(ci) < 2L) return(row)
  ci <- as.numeric(ci)
  if (any(is.na(ci[1:2]))) return(row)
  level <- conf_level %||% 0.95
  label <- paste0(prefix, as.integer(round(level * 100)), "% CI")
  row[[label]] <- paste0("[", ci[[1L]], ", ", ci[[2L]], "]")
  row
}

payload_inline_label <- function(x) {
  x <- as.character(x %||% "")
  switch(x,
    mean_diff = "Mean difference",
    odds_ratio = "Odds ratio",
    cramers_v = "Cramer's V",
    eta_sq_partial = "partial eta squared",
    std_alpha = "standardized alpha",
    average_r = "average r",
    logLik = "logLik",
    x
  )
}
