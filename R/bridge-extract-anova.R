# R bridge -- anova() extractor.
#
# Two common shapes:
#
# 1. Single-model anova(model)
#      colnames: Df, Sum Sq, Mean Sq, F value, Pr(>F)
#      One row per term + Residuals row at the bottom.
#      -> emit the LAST non-Residuals row (the focal term in sequential
#        coding) as the inline headline.
#
# 2. Model comparison anova(m1, m2)
#      colnames: Res.Df, RSS, Df, Sum of Sq, F, Pr(>F)
#      Last row contains the focal comparison F.
#      -> emit the LAST row.
#
# In both cases the full table goes into raw_output for reference.
# v1.5 will revisit this as part of the structural archetype where
# anova tables can also become Tables assets.

#' @rdname mellio_payload
#' @name mellio_payload
#' @export
mellio_payload.aov <- function(x, focal = NULL, controls = NULL, ..., .call = NULL) {
  call_str <- if (!is.null(.call)) {
    .call
  } else {
    formula_txt <- tryCatch(
      paste(deparse(stats::formula(x), width.cutoff = 500L), collapse = " "),
      error = function(e) NA_character_
    )
    if (!is.na(formula_txt) && nzchar(formula_txt)) formula_txt else ms_model_call_string(x)
  }
  ancova_spec <- ms_ancova_model_spec(x)
  if (is.null(focal) && !is.null(ancova_spec)) focal <- ancova_spec$factor$term
  if (is.null(controls) && !is.null(ancova_spec)) {
    controls <- vapply(ancova_spec$covariates, function(covariate) {
      covariate$term %||% covariate$variable
    }, character(1))
  }
  payload <- mellio_payload(
    stats::anova(x),
    focal = focal,
    controls = controls,
    ...,
    .call = call_str
  )
  balance_note <- ms_anova_balanced_design_note(x)
  if (!is.null(balance_note)) payload$fields$design_balance_note <- balance_note
  figure_focal <- focal
  if (is.null(figure_focal)) figure_focal <- payload$fields$term %||% NULL

  interaction_figure <- ms_anova_interaction_figure_data(x)
  if (!is.null(interaction_figure)) {
    payload$figure_data <- payload$figure_data %||% list()
    payload$figure_data$interaction_plot <- interaction_figure
    payload <- ms_add_available_figure(
      payload,
      type = "interaction_plot",
      label = "Interaction plot",
      default = length(payload$metadata$available_figures %||% list()) == 0L
    )
  }

  if (is.null(interaction_figure)) {
    means_figure <- ms_anova_means_figure_data(x)
    if (is.null(means_figure)) {
      means_figure <- ms_ancova_adjusted_means_figure_data(x, spec = ancova_spec)
    }
    if (is.null(means_figure)) {
      means_figure <- ms_factorial_main_effect_means_figure_data(
        x,
        focal = figure_focal,
        controls = controls
      )
    }
    if (!is.null(means_figure)) {
      payload$figure_data <- payload$figure_data %||% list()
      payload$figure_data$adjusted_means <- means_figure
      if (identical(means_figure$ci_method, "emmeans") ||
          grepl("emmeans", means_figure$source %||% "", fixed = TRUE)) {
        payload$packages <- ms_packages_basic(extras = "emmeans")
      }
      payload <- ms_add_available_figure(
        payload,
        type = "adjusted_means",
        label = "Means plot",
        default = length(payload$metadata$available_figures %||% list()) == 0L
      )
      if (identical(means_figure$source %||% "", "aov_one_way")) {
        payload <- ms_anova_attach_oneway_posthoc(payload, x, means_figure)
      }
    }
  }
  payload
}

#' @rdname mellio_payload
#' @export
mellio_payload.summary.aov <- function(x, focal = NULL, controls = NULL, ..., .call = NULL,
                                       .env = parent.frame()) {
  tables <- Filter(function(tbl) inherits(tbl, "anova") && is.data.frame(tbl), unclass(x))
  if (!length(tables)) {
    stop("summary.aov object does not contain an ANOVA table.", call. = FALSE)
  }
  if (length(tables) > 1L) {
    stop(
      "summary.aov objects with multiple ANOVA tables are not supported yet; pass the original aov or aovlist object.",
      call. = FALSE
    )
  }

  call_str <- if (!is.null(.call)) .call else "summary(aov(...))"
  table <- tables[[1L]]
  row_names <- rownames(table)
  if (!is.null(row_names)) rownames(table) <- trimws(row_names)

  payload <- mellio_payload(
    table,
    focal = focal,
    controls = controls,
    ...,
    .call = call_str
  )
  if (is.null(payload$fields$ss_type_label) || !nzchar(payload$fields$ss_type_label %||% "")) {
    payload$fields$ss_type <- "type_i_sequential"
    payload$fields$ss_type_label <- "Type I (sequential)"
    payload$fields$ss_type_note <- "Sequential tests depend on the order of terms in the model formula."
    payload$type_label <- ms_anova_type_label_with_ss(payload$type_label, "Type I (sequential)")
  }

  source_model <- ms_summary_aov_source_model(call_str, .env)
  if (!is.null(source_model)) {
    vars <- ms_anova_vars_from_call(ms_model_formula_string(source_model))
    source_terms <- tryCatch(
      attr(stats::terms(source_model), "term.labels") %||% character(0),
      error = function(e) character(0)
    )
    if (!is.null(vars$outcome) && !nzchar(payload$fields$outcome %||% "")) {
      payload$fields$outcome <- vars$outcome
    }
    if (!is.null(vars$predictor) &&
        (!nzchar(payload$fields$predictor %||% "") || length(source_terms) == 1L)) {
      payload$fields$predictor <- vars$predictor
    }
    balance_note <- ms_anova_balanced_design_note(source_model)
    if (!is.null(balance_note)) payload$fields$design_balance_note <- balance_note
  }
  payload$raw_output <- ms_capture_output(x)
  payload
}

ms_summary_aov_source_model <- function(call_str, env) {
  if (!is.environment(env)) return(NULL)
  target <- ms_summary_call_target_symbol(call_str)
  if (!nzchar(target)) return(NULL)
  model <- tryCatch(get(target, envir = env, inherits = TRUE), error = function(e) NULL)
  if (inherits(model, c("aov", "aovlist", "lm"))) model else NULL
}

ms_summary_call_target_symbol <- function(call_str) {
  call_str <- trimws(as.character(call_str %||% ""))
  if (!length(call_str) || is.na(call_str)) return("")
  m <- regexec("^summary\\(([^,\\)]+)(?:[,\\)].*)$", call_str, perl = TRUE)
  hit <- regmatches(call_str, m)[[1]]
  if (length(hit) < 2L) return("")
  target <- trimws(hit[[2]])
  if (grepl("^[A-Za-z.][A-Za-z0-9._]*$", target)) target else ""
}

ms_model_formula_string <- function(model) {
  formula_txt <- tryCatch(
    paste(deparse(stats::formula(model), width.cutoff = 500L), collapse = " "),
    error = function(e) NA_character_
  )
  if (!is.na(formula_txt) && nzchar(formula_txt)) formula_txt else ""
}

ms_anova_balanced_design_note <- function(model) {
  mf <- tryCatch(stats::model.frame(model), error = function(e) NULL)
  trm <- tryCatch(stats::terms(model), error = function(e) NULL)
  if (is.null(mf) || is.null(trm)) return(NULL)

  term_labels <- attr(trm, "term.labels") %||% character(0)
  if (!length(term_labels) || !any(grepl(":", term_labels, fixed = TRUE))) {
    return(NULL)
  }

  response <- as.character(attr(trm, "variables")[[2L]] %||% "")
  factor_vars <- setdiff(all.vars(stats::delete.response(trm)), response)
  factor_vars <- factor_vars[factor_vars %in% names(mf)]
  if (length(factor_vars) < 2L) return(NULL)

  is_discrete <- vapply(mf[factor_vars], function(value) {
    is.factor(value) || is.character(value) || is.logical(value)
  }, logical(1))
  if (!all(is_discrete)) return(NULL)

  complete <- stats::complete.cases(mf[factor_vars])
  if (!any(complete)) return(NULL)
  tab <- table(mf[complete, factor_vars, drop = FALSE])
  counts <- as.numeric(tab)
  if (!length(counts) || any(counts == 0L) || length(unique(counts)) != 1L) {
    return(NULL)
  }
  "Because the design is balanced, Type I, II, and III sums of squares coincide."
}

ms_anova_interaction_figure_data <- function(x) {
  figure <- ms_interaction_plot_data(x)
  if (is.null(figure)) return(NULL)
  if (ms_anova_blank_text(figure$source)) {
    figure$source <- "aov_interaction"
  }
  categorical_kinds <- c(
    "categorical_by_categorical",
    "categorical_by_categorical_by_categorical"
  )
  if (figure$interaction_kind %in% categorical_kinds) {
    if (ms_anova_blank_text(figure$mean_kind)) {
      figure$mean_kind <- "estimated_marginal"
    }
    if (!is.null(figure$outcome) && nzchar(figure$outcome %||% "")) {
      figure$y_label <- paste("Estimated marginal mean", figure$outcome)
    }
  }
  figure
}

ms_anova_blank_text <- function(value) {
  if (is.null(value) || length(value) == 0L) return(TRUE)
  text <- trimws(as.character(value[[1]]))
  is.na(text) || !nzchar(text)
}

#' @rdname mellio_payload
#' @export
mellio_payload.aovlist <- function(x, ..., .call = NULL) {
  rlang::check_installed("broom", reason = "to extract repeated-measures ANOVA results")
  call_str <- if (!is.null(.call)) .call else ms_model_call_string(x)
  tidy_df <- tryCatch(broom::tidy(x), error = function(e) NULL)
  if (is.null(tidy_df) || !nrow(tidy_df)) {
    stop("broom could not tidy this aovlist object.", call. = FALSE)
  }
  tidy_df <- as.data.frame(tidy_df)
  columns <- list()
  add_col <- function(key, label, format) {
    columns[[length(columns) + 1L]] <<- list(key = key, label = label, format = format)
  }
  if ("stratum" %in% names(tidy_df)) add_col("stratum", "Error stratum", "text")
  if ("term" %in% names(tidy_df)) add_col("term", "Term", "text")
  if ("df" %in% names(tidy_df)) add_col("df", "df", "integer")
  if ("sumsq" %in% names(tidy_df)) add_col("sumsq", "Sum Sq", "number")
  if ("meansq" %in% names(tidy_df)) add_col("meansq", "Mean Sq", "number")
  if ("statistic" %in% names(tidy_df)) add_col("statistic", "F", "number")
  if ("p.value" %in% names(tidy_df)) add_col("p_value", "p", "pvalue")

  rows <- lapply(seq_len(nrow(tidy_df)), function(i) {
    row <- list()
    if ("stratum" %in% names(tidy_df)) row$stratum <- as.character(tidy_df$stratum[i])
    if ("term" %in% names(tidy_df)) row$term <- as.character(tidy_df$term[i])
    if ("df" %in% names(tidy_df)) row$df <- ms_safe_numeric(tidy_df$df[i])
    if ("sumsq" %in% names(tidy_df)) row$sumsq <- ms_safe_numeric(tidy_df$sumsq[i])
    if ("meansq" %in% names(tidy_df)) row$meansq <- ms_safe_numeric(tidy_df$meansq[i])
    if ("statistic" %in% names(tidy_df)) row$statistic <- ms_safe_numeric(tidy_df$statistic[i])
    if ("p.value" %in% names(tidy_df)) row$p_value <- ms_safe_numeric(tidy_df$p.value[i])
    if (!is.null(row$statistic) && !is.na(row$statistic)) {
      row$f <- row$statistic
      row$df1 <- row$df
    }
    row
  })

  rows <- ms_aovlist_annotate_error_df(rows)
  repeated_metadata <- ms_aovlist_metadata(x)
  if (any(vapply(rows, function(row) {
    is.list(row) && !is.null(row$eta_sq_partial) && !is.na(row$eta_sq_partial)
  }, logical(1)))) {
    add_col("eta_sq_partial", "partial \u03b7\u00b2", "bounded")
  }
  rows <- ms_aovlist_label_strata(
    rows,
    repeated_metadata$subject_variable %||% repeated_metadata$subject %||% ""
  )
  rows <- ms_aovlist_publishable_rows(rows)
  means_figure <- ms_repeated_measures_means_figure_data(x, tidy_df)
  interaction_figure <- if (is.null(means_figure)) {
    ms_repeated_measures_interaction_figure_data(x, tidy_df)
  } else NULL
  figure_data <- list()
  if (!is.null(means_figure)) figure_data$adjusted_means <- means_figure
  if (!is.null(interaction_figure)) figure_data$interaction_plot <- interaction_figure

  ms_build_envelope(
    type = "custom_anova_table",
    type_label = "Repeated-measures ANOVA table",
    call = trimws(gsub("\\s+", " ", call_str)),
    fields = c(list(
      table_type = "anova",
      columns = columns,
      rows = rows,
      source = "R aovlist"
    ), repeated_metadata),
    raw_output = ms_capture_output(summary(x)),
    figure_data = if (length(figure_data)) figure_data else NULL,
    packages = ms_packages_basic(extras = "broom"),
    card_kind = "table"
  )
}

#' @rdname mellio_payload
#' @export
mellio_payload.afex_aov <- function(x, ..., .call = NULL) {
  anova_table <- tryCatch(as.data.frame(x$anova_table), error = function(e) NULL)
  if (is.null(anova_table) || !nrow(anova_table)) {
    stop("afex_aov object does not contain an ANOVA table.", call. = FALSE)
  }

  call_str <- if (!is.null(.call)) .call else ms_model_call_string(x)
  type_num <- attr(x, "type") %||% attr(anova_table, "type")
  correction <- attr(x$anova_table, "correction") %||% ""
  es <- attr(x$anova_table, "es") %||% ""
  outcome <- ms_model_clean_term(attr(x, "dv") %||% ms_afex_response_from_heading(x$anova_table))
  subject <- ms_model_clean_term(attr(x, "id") %||% "")
  within <- attr(x, "within") %||% list()
  between <- attr(x, "between") %||% list()
  within_terms <- names(within)
  between_terms <- names(between)

  rows <- ms_afex_anova_rows(anova_table)
  rows <- ms_afex_attach_sphericity(rows, x)
  sphericity <- ms_afex_sphericity_tests(x)
  corrections <- ms_afex_sphericity_corrections(x)
  n_subjects <- ms_afex_n_subjects(x)
  interaction_figure <- ms_afex_mixed_interaction_plot_data(x, rows)
  means_figure <- if (is.null(interaction_figure)) {
    ms_afex_repeated_measures_means_figure_data(x, rows)
  } else NULL
  figure_data <- list()
  if (!is.null(interaction_figure)) figure_data$interaction_plot <- interaction_figure
  if (!is.null(means_figure)) figure_data$adjusted_means <- means_figure
  if (!length(figure_data)) figure_data <- NULL

  type_label <- if (length(within_terms) && length(between_terms)) {
    "Mixed ANOVA table"
  } else if (length(within_terms)) {
    "Repeated-measures ANOVA table"
  } else {
    "ANOVA table"
  }
  if (!is.null(type_num) && length(type_num) && !is.na(type_num)) {
    type_label <- paste0(type_label, " (Type ", as.character(type_num[[1L]]), ")")
  }

  columns <- list(
    list(key = "term", label = "Effect", format = "text"),
    list(key = "df1", label = "df1", format = "number"),
    list(key = "df2", label = "df2", format = "number"),
    list(key = "mse", label = "MSE", format = "number"),
    list(key = "f", label = "F", format = "number")
  )
  if (any(vapply(rows, function(row) !is.null(row$ges), logical(1)))) {
    columns <- c(columns, list(list(key = "ges", label = "generalized \u03b7\u00b2", format = "bounded")))
  }
  if (any(vapply(rows, function(row) !is.null(row$eta_sq_partial), logical(1)))) {
    columns <- c(columns, list(list(key = "eta_sq_partial", label = "partial \u03b7\u00b2", format = "bounded")))
  }
  columns <- c(columns, list(list(key = "p_value", label = "p", format = "pvalue")))

  fields <- list(
    table_type = "anova",
    model_kind = if (length(within_terms) && length(between_terms)) {
      "mixed_anova"
    } else if (length(within_terms)) {
      "repeated_measures_anova"
    } else {
      "anova"
    },
    columns = columns,
    rows = rows,
    source = "R afex_aov",
    outcome = outcome,
    subject = subject,
    n = n_subjects,
    within_terms = unname(within_terms),
    between_terms = unname(between_terms),
    factors = ms_afex_factor_metadata(within, between),
    ss_type = if (!is.null(type_num) && length(type_num)) paste0("type_", as.character(type_num[[1L]])) else NULL,
    ss_type_label = if (!is.null(type_num) && length(type_num)) paste0("Type ", as.character(type_num[[1L]])) else NULL,
    effect_size = es,
    correction_method = correction,
    correction_label = ms_afex_correction_label(correction),
    sphericity_tests = sphericity,
    sphericity_corrections = corrections
  )

  ms_build_envelope(
    type = "custom_anova_table",
    type_label = type_label,
    call = trimws(gsub("\\s+", " ", call_str)),
    fields = fields,
    raw_output = ms_capture_output(summary(x)),
    figure_data = figure_data,
    packages = ms_packages_basic(extras = "afex"),
    card_kind = "table"
  )
}

ms_aovlist_annotate_error_df <- function(rows) {
  if (!length(rows)) return(rows)
  for (i in seq_along(rows)) {
    row <- rows[[i]]
    if (!is.list(row) || is.null(row$f) || is.na(row$f)) next
    stratum <- row$stratum %||% ""
    resid_idx <- which(vapply(rows, function(candidate) {
      is.list(candidate) &&
        identical(candidate$stratum %||% "", stratum) &&
        grepl("^Residuals?$", candidate$term %||% "", ignore.case = TRUE)
    }, logical(1)))
    if (length(resid_idx) > 0L) {
      rows[[i]]$df2 <- rows[[resid_idx[[1L]]]]$df
      rows[[i]]$residual_df <- rows[[resid_idx[[1L]]]]$df
      if (!is.null(rows[[resid_idx[[1L]]]]$sumsq)) {
        rows[[i]]$residual_sum_sq <- rows[[resid_idx[[1L]]]]$sumsq
      }
      if (!is.null(rows[[resid_idx[[1L]]]]$meansq)) {
        rows[[i]]$residual_mean_sq <- rows[[resid_idx[[1L]]]]$meansq
      }
      effect_ss <- ms_safe_numeric(row$sumsq %||% NA_real_)
      error_ss <- ms_safe_numeric(rows[[resid_idx[[1L]]]]$sumsq %||% NA_real_)
      if (!is.na(effect_ss) && !is.na(error_ss) && (effect_ss + error_ss) > 0) {
        eta_sq_partial <- ms_safe_numeric(effect_ss / (effect_ss + error_ss))
        rows[[i]]$eta_sq_partial <- eta_sq_partial
        rows[[i]]$eta_sq <- eta_sq_partial
        rows[[i]]$eta_sq_is_partial <- TRUE
        rows[[i]]$effect <- list(
          name = "eta_sq_partial",
          value = eta_sq_partial
        )
      }
    }
  }
  rows
}

ms_aovlist_label_strata <- function(rows, subject_var = "") {
  if (!length(rows)) return(rows)
  subject_var <- trimws(as.character(subject_var %||% ""))
  for (i in seq_along(rows)) {
    stratum <- trimws(as.character(rows[[i]]$stratum %||% ""))
    if (!nzchar(stratum)) next
    rows[[i]]$stratum_raw <- stratum
    rows[[i]]$stratum <- ms_aovlist_stratum_label(stratum, subject_var)
  }
  rows
}

ms_aovlist_publishable_rows <- function(rows) {
  if (!length(rows)) return(rows)
  Filter(function(row) {
    if (!is.list(row)) return(TRUE)
    stratum <- trimws(as.character(row$stratum %||% ""))
    term <- trimws(as.character(row$term %||% ""))
    !(identical(stratum, "Between subjects") &&
        grepl("^Residuals?$", term, ignore.case = TRUE))
  }, rows)
}

ms_aovlist_stratum_label <- function(stratum, subject_var = "") {
  stratum <- trimws(as.character(stratum %||% ""))
  if (!nzchar(stratum)) return(stratum)
  subject_var <- trimws(as.character(subject_var %||% ""))
  parts <- strsplit(stratum, ":", fixed = TRUE)[[1L]]
  parts <- trimws(parts[nzchar(parts)])
  if (nzchar(subject_var)) {
    parts <- parts[parts != subject_var]
  } else if (length(parts) > 1L) {
    parts <- parts[-1L]
  }
  if (!length(parts)) return("Between subjects")
  paste0(
    paste(vapply(parts, ms_model_clean_term, character(1)), collapse = " \u00d7 "),
    " \u00d7 subjects"
  )
}

ms_afex_response_from_heading <- function(x) {
  heading <- attr(x, "heading") %||% character(0)
  hit <- grep("^Response\\s*:", heading, value = TRUE)
  if (!length(hit)) return("")
  trimws(sub("^Response\\s*:\\s*", "", hit[[1L]]))
}

ms_afex_col <- function(df, candidates) {
  if (is.null(df) || !length(names(df))) return(NA_integer_)
  norm <- function(value) tolower(gsub("[^A-Za-z0-9]+", "", value))
  names_norm <- norm(names(df))
  cand_norm <- norm(candidates)
  hit <- match(cand_norm, names_norm, nomatch = 0L)
  hit <- hit[hit > 0L]
  if (length(hit)) hit[[1L]] else NA_integer_
}

ms_afex_value <- function(df, i, candidates) {
  col <- ms_afex_col(df, candidates)
  if (is.na(col)) return(NA_real_)
  ms_safe_numeric(df[[col]][[i]])
}

ms_afex_anova_rows <- function(anova_table) {
  term_names <- rownames(anova_table) %||% paste0("Effect ", seq_len(nrow(anova_table)))
  lapply(seq_len(nrow(anova_table)), function(i) {
    term <- as.character(term_names[[i]])
    df1 <- ms_afex_value(anova_table, i, c("num Df", "num.Df", "df1", "Df"))
    df2 <- ms_afex_value(anova_table, i, c("den Df", "den.Df", "df2", "denDF"))
    f_value <- ms_afex_value(anova_table, i, c("F", "F value", "F.value"))
    p_value <- ms_afex_value(anova_table, i, c("Pr(>F)", "p.value", "p", "Pr..F."))
    mse <- ms_afex_value(anova_table, i, c("MSE", "Mean Sq", "meansq"))
    ges <- ms_afex_value(anova_table, i, c("ges", "generalized eta squared"))
    pes <- ms_afex_value(anova_table, i, c("pes", "eta_sq_partial", "eta2_partial", "partial eta squared"))

    row <- list(
      term = term,
      label = ms_model_clean_term(term),
      df1 = df1,
      df2 = df2,
      df = df1,
      statistic = f_value,
      f = f_value,
      p_value = p_value,
      term_type = if (grepl(":", term, fixed = TRUE)) "interaction" else "main"
    )
    if (!is.na(mse)) row$mse <- mse
    if (!is.na(ges)) {
      row$ges <- ges
      row$effect <- list(name = "eta_sq_generalized", value = ges)
    }
    if (!is.na(pes)) {
      row$eta_sq_partial <- pes
      row$effect <- list(name = "eta_sq_partial", value = pes)
    }
    row
  })
}

ms_afex_summary <- function(x) {
  tryCatch(
    suppressMessages(suppressWarnings(summary(x))),
    error = function(e) NULL
  )
}

ms_afex_summary_table_df <- function(x) {
  if (is.null(x)) return(NULL)
  mat <- tryCatch(as.matrix(unclass(x)), error = function(e) NULL)
  if (is.null(mat) || !length(mat)) return(NULL)
  df <- as.data.frame(mat, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(df) <- rownames(mat)
  df
}

ms_afex_sphericity_tests <- function(x) {
  s <- ms_afex_summary(x)
  tab <- ms_afex_summary_table_df(s$sphericity.tests)
  if (is.null(tab) || !nrow(tab)) return(list())
  terms <- rownames(tab) %||% paste0("Effect ", seq_len(nrow(tab)))
  lapply(seq_len(nrow(tab)), function(i) {
    p_value <- ms_afex_value(tab, i, c("p-value", "p.value", "p"))
    row <- list(
      term = as.character(terms[[i]]),
      label = ms_model_clean_term(terms[[i]]),
      statistic = "Mauchly's W",
      w = ms_afex_value(tab, i, c("Test statistic", "Mauchly's W", "W")),
      p_value = p_value
    )
    if (!is.na(p_value)) row$violated <- isTRUE(p_value < 0.05)
    row
  })
}

ms_afex_sphericity_corrections <- function(x) {
  s <- ms_afex_summary(x)
  tab <- ms_afex_summary_table_df(s$pval.adjustments)
  if (is.null(tab) || !nrow(tab)) return(list())
  terms <- rownames(tab) %||% paste0("Effect ", seq_len(nrow(tab)))
  lapply(seq_len(nrow(tab)), function(i) {
    list(
      term = as.character(terms[[i]]),
      label = ms_model_clean_term(terms[[i]]),
      gg_epsilon = ms_afex_value(tab, i, c("GG eps", "GGe", "GG.eps")),
      gg_p_value = ms_afex_value(tab, i, c("Pr(>F[GG])", "Pr..F.GG..", "p[GG]")),
      hf_epsilon = ms_afex_value(tab, i, c("HF eps", "HFe", "HF.eps")),
      hf_p_value = ms_afex_value(tab, i, c("Pr(>F[HF])", "Pr..F.HF..", "p[HF]"))
    )
  })
}

ms_afex_attach_sphericity <- function(rows, x) {
  if (!length(rows)) return(rows)
  correction <- attr(x$anova_table, "correction") %||% ""
  correction_label <- ms_afex_correction_label(correction)
  sphericity <- ms_afex_sphericity_tests(x)
  corrections <- ms_afex_sphericity_corrections(x)

  for (i in seq_along(rows)) {
    term <- rows[[i]]$term %||% ""
    spher <- sphericity[vapply(sphericity, function(row) identical(row$term, term), logical(1))]
    if (length(spher)) {
      rows[[i]]$sphericity <- spher[[1L]]
    }
    corr <- corrections[vapply(corrections, function(row) identical(row$term, term), logical(1))]
    if (length(corr)) {
      rows[[i]]$sphericity_correction <- corr[[1L]]
    }
    if (nzchar(correction) && !identical(tolower(correction), "none") && length(corr)) {
      rows[[i]]$correction_method <- correction
      rows[[i]]$correction_label <- correction_label
      rows[[i]]$corrected <- TRUE
    }
  }
  rows
}

ms_afex_correction_label <- function(correction) {
  key <- toupper(trimws(as.character(correction %||% "")))
  if (!nzchar(key) || identical(key, "NONE")) return("")
  labels <- c(
    GG = "Greenhouse-Geisser",
    HF = "Huynh-Feldt"
  )
  labels[[key]] %||% correction
}

ms_afex_factor_metadata <- function(within, between) {
  out <- list()
  append_terms <- function(terms, role) {
    if (!length(terms)) return()
    for (term in names(terms)) {
      levels <- as.character(terms[[term]] %||% character(0))
      out[[length(out) + 1L]] <<- list(
        variable = term,
        label = ms_model_clean_term(term),
        role = role,
        levels = lapply(levels, function(level) {
          list(value = level, label = ms_model_clean_term(level))
        })
      )
    }
  }
  append_terms(within, "within")
  append_terms(between, "between")
  out
}

ms_afex_n_subjects <- function(x) {
  id <- attr(x, "id") %||% ""
  long <- tryCatch(x$data$long, error = function(e) NULL)
  if (is.null(long) || !is.data.frame(long) || !nzchar(id) || !id %in% names(long)) {
    return(NA_real_)
  }
  length(unique(long[[id]][!is.na(long[[id]])]))
}

ms_afex_repeated_measures_means_figure_data <- function(x, rows = NULL,
                                                        conf.level = 0.95,
                                                        max_levels = 12L) {
  within <- attr(x, "within") %||% list()
  between <- attr(x, "between") %||% list()
  if (length(within) != 1L || length(between) > 0L) return(NULL)

  outcome <- attr(x, "dv") %||% ""
  subject_var <- attr(x, "id") %||% ""
  within_var <- names(within)[[1L]]
  long <- tryCatch(x$data$long, error = function(e) NULL)
  if (is.null(long) || !is.data.frame(long) ||
      !outcome %in% names(long) ||
      !subject_var %in% names(long) ||
      !within_var %in% names(long)) {
    return(NULL)
  }

  y <- long[[outcome]]
  within_values <- long[[within_var]]
  subject_values <- long[[subject_var]]
  if (!is.numeric(y)) return(NULL)

  ok <- is.finite(y) & !is.na(within_values) & !is.na(subject_values)
  y <- y[ok]
  within_values <- within_values[ok]
  subject_values <- subject_values[ok]
  if (length(y) < 2L) return(NULL)

  levels <- as.character(within[[within_var]] %||% ms_anova_observed_levels(within_values))
  observed <- ms_anova_observed_levels(within_values)
  levels <- levels[levels %in% observed]
  if (!length(levels)) levels <- observed
  if (length(levels) < 2L || length(levels) > max_levels) return(NULL)

  subjects <- unique(as.character(subject_values))
  if (length(subjects) < 2L) return(NULL)

  cell_frame <- data.frame(
    subject = as.character(subject_values),
    level = factor(as.character(within_values), levels = levels),
    y = y,
    stringsAsFactors = FALSE
  )
  cells <- stats::aggregate(
    y ~ subject + level,
    data = cell_frame,
    FUN = function(value) mean(value, na.rm = TRUE),
    drop = FALSE
  )

  mat <- matrix(NA_real_, nrow = length(subjects), ncol = length(levels),
                dimnames = list(subjects, levels))
  for (i in seq_len(nrow(cells))) {
    subject <- as.character(cells$subject[[i]])
    level <- as.character(cells$level[[i]])
    if (!subject %in% rownames(mat) || !level %in% colnames(mat)) next
    mat[subject, level] <- ms_safe_numeric(cells$y[[i]])
  }
  mat <- mat[stats::complete.cases(mat), , drop = FALSE]
  if (nrow(mat) < 2L || ncol(mat) < 2L) return(NULL)

  means <- colMeans(mat, na.rm = TRUE)
  subject_means <- rowMeans(mat, na.rm = TRUE)
  grand_mean <- mean(mat, na.rm = TRUE)
  normalized <- sweep(mat, 1L, subject_means, FUN = "-") + grand_mean
  morey <- sqrt(ncol(mat) / (ncol(mat) - 1L))
  df_ci <- nrow(mat) - 1L
  critical <- stats::qt(1 - ((1 - conf.level) / 2), df = df_ci)
  se <- apply(normalized, 2L, stats::sd, na.rm = TRUE) / sqrt(nrow(mat)) * morey

  groups <- lapply(levels, function(level) {
    mean_value <- ms_safe_numeric(means[[level]])
    se_value <- ms_safe_numeric(se[[level]])
    raw_sd <- ms_safe_numeric(stats::sd(mat[, level], na.rm = TRUE))
    row <- list(
      level = level,
      label = ms_model_clean_term(level),
      n = nrow(mat),
      mean = mean_value,
      df = df_ci
    )
    if (!is.na(raw_sd)) row$sd <- raw_sd
    if (!is.na(se_value)) {
      row$se <- se_value
      row$ci_lower <- ms_safe_numeric(mean_value - critical * se_value)
      row$ci_upper <- ms_safe_numeric(mean_value + critical * se_value)
    }
    row
  })
  groups <- Filter(function(row) {
    is.list(row) &&
      !is.null(row$level) &&
      !is.null(row$mean) &&
      !is.na(ms_safe_numeric(row$mean))
  }, groups)
  if (length(groups) < 2L) return(NULL)

  focal_row <- NULL
  if (length(rows)) {
    idx <- which(vapply(rows, function(row) identical(row$term, within_var), logical(1)))
    if (length(idx)) focal_row <- rows[[idx[[1L]]]]
  }
  focal_test <- if (!is.null(focal_row)) {
    test <- list(
      statistic = list(name = "F", value = focal_row$f, df = I(c(focal_row$df1, focal_row$df2))),
      p_value = focal_row$p_value
    )
    if (!is.null(focal_row$effect)) test$effect <- focal_row$effect
    if (!is.null(focal_row$correction_method)) test$correction_method <- focal_row$correction_method
    if (!is.null(focal_row$correction_label)) test$correction_label <- focal_row$correction_label
    test
  } else NULL

  out <- list(
    mean_kind = "within_subject",
    source = "afex_aov",
    factor = list(
      variable = within_var,
      term = within_var,
      label = ms_model_clean_term(within_var),
      levels = lapply(levels, function(level) {
        list(value = level, label = ms_model_clean_term(level))
      })
    ),
    groups = groups,
    subject = list(
      variable = subject_var,
      label = ms_model_clean_term(subject_var),
      n = nrow(mat)
    ),
    outcome = ms_model_clean_term(outcome),
    y_label = if (nzchar(outcome)) paste("Mean", ms_model_clean_term(outcome)) else "Mean outcome",
    ci_level = conf.level,
    ci_method = "within_subject_morey"
  )
  if (!is.null(focal_test)) out$focal_test <- focal_test
  out
}

ms_afex_mixed_interaction_plot_data <- function(x, rows = NULL,
                                                conf.level = 0.95,
                                                max_levels = 8L) {
  within <- attr(x, "within") %||% list()
  between <- attr(x, "between") %||% list()
  if (length(within) != 1L || length(between) != 1L) return(NULL)

  outcome <- attr(x, "dv") %||% ""
  subject_var <- attr(x, "id") %||% ""
  within_var <- names(within)[[1L]]
  between_var <- names(between)[[1L]]
  long <- tryCatch(x$data$long, error = function(e) NULL)
  if (is.null(long) || !is.data.frame(long) ||
      !outcome %in% names(long) ||
      !subject_var %in% names(long) ||
      !within_var %in% names(long) ||
      !between_var %in% names(long)) {
    return(NULL)
  }

  y <- long[[outcome]]
  subject_values <- long[[subject_var]]
  within_values <- long[[within_var]]
  between_values <- long[[between_var]]
  if (!is.numeric(y)) return(NULL)

  ok <- is.finite(y) & !is.na(subject_values) & !is.na(within_values) & !is.na(between_values)
  y <- y[ok]
  subject_values <- subject_values[ok]
  within_values <- within_values[ok]
  between_values <- between_values[ok]
  if (length(y) < 2L) return(NULL)

  within_levels <- as.character(within[[within_var]] %||% ms_anova_observed_levels(within_values))
  between_levels <- as.character(between[[between_var]] %||% ms_anova_observed_levels(between_values))
  observed_within <- ms_anova_observed_levels(within_values)
  observed_between <- ms_anova_observed_levels(between_values)
  within_levels <- within_levels[within_levels %in% observed_within]
  between_levels <- between_levels[between_levels %in% observed_between]
  if (!length(within_levels)) within_levels <- observed_within
  if (!length(between_levels)) between_levels <- observed_between
  if (length(within_levels) < 2L || length(between_levels) < 2L) return(NULL)
  if (length(within_levels) > max_levels || length(between_levels) > max_levels) return(NULL)

  cell_frame <- data.frame(
    subject = as.character(subject_values),
    x_level = factor(as.character(within_values), levels = within_levels),
    moderator_level = factor(as.character(between_values), levels = between_levels),
    y = y,
    stringsAsFactors = FALSE
  )
  subject_cells <- stats::aggregate(
    y ~ subject + moderator_level + x_level,
    data = cell_frame,
    FUN = function(value) mean(value, na.rm = TRUE)
  )
  if (!nrow(subject_cells)) return(NULL)

  critical <- function(df) {
    if (is.na(df) || df <= 0) return(NA_real_)
    stats::qt(1 - ((1 - conf.level) / 2), df = df)
  }
  grid <- list()
  for (mod_level in between_levels) {
    for (x_index in seq_along(within_levels)) {
      x_level <- within_levels[[x_index]]
      values <- subject_cells$y[
        as.character(subject_cells$moderator_level) == mod_level &
          as.character(subject_cells$x_level) == x_level
      ]
      values <- values[is.finite(values)]
      n <- length(values)
      if (!n) next
      estimate <- ms_safe_numeric(mean(values, na.rm = TRUE))
      df_ci <- n - 1L
      sd_value <- if (n > 1L) ms_safe_numeric(stats::sd(values, na.rm = TRUE)) else NA_real_
      se_value <- if (!is.na(sd_value) && n > 1L) ms_safe_numeric(sd_value / sqrt(n)) else NA_real_
      row <- list(
        x = x_index,
        x_value = x_level,
        x_label = ms_model_clean_term(x_level),
        moderator_value = mod_level,
        moderator_label = ms_model_clean_term(mod_level),
        estimate = estimate,
        n = n,
        df = df_ci
      )
      if (!is.na(sd_value)) row$sd <- sd_value
      if (!is.na(se_value)) {
        row$se <- se_value
        ci <- critical(df_ci)
        if (!is.na(ci)) {
          row$ci_lower <- ms_safe_numeric(estimate - ci * se_value)
          row$ci_upper <- ms_safe_numeric(estimate + ci * se_value)
        }
      }
      grid[[length(grid) + 1L]] <- row
    }
  }
  if (length(grid) < length(within_levels) * length(between_levels)) return(NULL)

  interaction_term <- paste(between_var, within_var, sep = ":")
  interaction_row <- NULL
  if (length(rows)) {
    wanted <- sort(c(within_var, between_var))
    idx <- which(vapply(rows, function(row) {
      term <- row$term %||% ""
      identical(sort(strsplit(term, ":", fixed = TRUE)[[1L]]), wanted)
    }, logical(1)))
    if (length(idx)) {
      interaction_row <- rows[[idx[[1L]]]]
      interaction_term <- interaction_row$term %||% interaction_term
    }
  }
  focal_test <- if (!is.null(interaction_row)) {
    test <- list(
      statistic = list(name = "F", value = interaction_row$f, df = I(c(interaction_row$df1, interaction_row$df2))),
      p_value = interaction_row$p_value
    )
    if (!is.null(interaction_row$effect)) test$effect <- interaction_row$effect
    if (!is.null(interaction_row$correction_method)) test$correction_method <- interaction_row$correction_method
    if (!is.null(interaction_row$correction_label)) test$correction_label <- interaction_row$correction_label
    test
  } else NULL

  out <- list(
    interaction_term = interaction_term,
    interaction_kind = "categorical_by_categorical",
    source = "afex_aov",
    mean_kind = "observed_cell_means",
    variables = c(within_var, between_var),
    x = list(
      variable = within_var,
      term = within_var,
      label = ms_model_clean_term(within_var),
      type = "categorical",
      range = I(c(1, length(within_levels))),
      levels = lapply(within_levels, function(level) {
        list(value = level, label = ms_model_clean_term(level))
      })
    ),
    moderator = list(
      variable = between_var,
      term = between_var,
      label = ms_model_clean_term(between_var),
      type = "categorical",
      levels = lapply(between_levels, function(level) {
        list(value = level, label = ms_model_clean_term(level))
      })
    ),
    grid = grid,
    subject = list(
      variable = subject_var,
      label = ms_model_clean_term(subject_var),
      n = length(unique(as.character(subject_values)))
    ),
    outcome = ms_model_clean_term(outcome),
    y_label = if (nzchar(outcome)) paste("Mean", ms_model_clean_term(outcome)) else "Mean outcome",
    scale = "response",
    ci_level = conf.level,
    ci_method = "cell_standard_error",
    bounded_response = FALSE
  )
  if (!is.null(focal_test)) out$focal_test <- focal_test
  out
}

#' @rdname mellio_payload
#' @export
mellio_payload.anova <- function(x, focal = NULL, controls = NULL, ..., .call = NULL) {
  call_str <- if (!is.null(.call)) {
    .call
  } else {
    user_call <- match.call()$x
    if (!is.null(user_call) && !identical(user_call, as.name("x"))) {
      paste(deparse(user_call, width.cutoff = 500L), collapse = " ")
    } else NA_character_
  }

  is_comparison <- "Sum of Sq" %in% colnames(x) || "Res.Df" %in% colnames(x)

  ss_info <- if (!is_comparison) ms_anova_sum_squares_info(x) else NULL
  role_defaults <- if (!is_comparison &&
                       !is.null(ss_info) &&
                       ss_info$type %in% c("type_ii", "type_iii")) {
    ms_anova_infer_role_defaults(x)
  } else {
    list()
  }
  focal_for_roles <- focal %||% role_defaults$focal
  controls_for_roles <- controls %||% role_defaults$controls

  if (is_comparison) {
    info <- ms_anova_comparison(x)
  } else {
    info <- ms_anova_single(x, focal = focal_for_roles, controls = controls_for_roles)
    if (!is.null(ss_info)) {
      info$fields$ss_type <- ss_info$type
      info$fields$ss_type_label <- ss_info$label
      if (!is.null(ss_info$note)) info$fields$ss_type_note <- ss_info$note
      info$type_label <- ms_anova_type_label_with_ss(info$type_label, ss_info$label)
    }
  }

  vars <- ms_anova_vars_from_call(call_str)
  if (is.null(vars$outcome)) {
    vars$outcome <- ms_anova_outcome_from_heading(x)
  }
  if (!is.null(vars$outcome) && !nzchar(info$fields$outcome %||% "")) {
    info$fields$outcome <- vars$outcome
  }
  if (!is.null(vars$predictor) && !nzchar(info$fields$predictor %||% "")) {
    info$fields$predictor <- vars$predictor
  }

  if (!is_comparison) {
    role_info <- ms_anova_term_roles(x, focal = focal_for_roles, controls = controls_for_roles)
    if (length(role_info$terms) > 0L) info$fields$terms <- role_info$terms
    if (length(role_info$focal_terms) > 0L) info$fields$focal_terms <- role_info$focal_terms
    if (length(role_info$control_terms) > 0L) info$fields$control_terms <- role_info$control_terms
    if (!is.null(role_info$predictor) && !nzchar(info$fields$predictor %||% "")) {
      info$fields$predictor <- role_info$predictor
    }
    if (!is.null(role_info$model_kind)) info$fields$model_kind <- role_info$model_kind
    if (identical(info$fields$model_kind, "ancova")) {
      info$type_label <- sub("^ANOVA", "ANCOVA", info$type_label)
    }
  }

  ms_build_envelope(
    type       = info$type,
    type_label = info$type_label,
    call       = trimws(gsub("\\s+", " ", call_str)),
    fields     = info$fields,
    raw_output = ms_capture_output(x)
  )
}

ms_anova_infer_role_defaults <- function(x) {
  rn <- rownames(x)
  if (is.null(rn) || length(rn) == 0L) return(list())

  is_resid <- grepl("^Residuals?$", rn, ignore.case = TRUE)
  is_intercept <- grepl("^\\(?Intercept\\)?$", rn, ignore.case = TRUE)
  term_names <- rn[!is_resid & !is_intercept]
  if (!length(term_names)) return(list())

  is_interaction <- grepl(":", term_names, fixed = TRUE)
  factor_wrapped <- vapply(term_names, function(term) {
    variable <- ms_ancova_single_variable(term)
    !is.null(variable) &&
      !identical(trimws(term), trimws(variable)) &&
      ms_ancova_supported_factor_term(term, variable)
  }, logical(1))

  if (sum(factor_wrapped & !is_interaction) != 1L) return(list())

  focal <- term_names[which(factor_wrapped & !is_interaction)[[1L]]]
  controls <- term_names[term_names != focal & !is_interaction]
  list(focal = focal, controls = controls)
}

# Model-comparison case: anova(m1, m2). Headline is the last row's F.
ms_anova_comparison <- function(x) {
  # Last row has the focal comparison
  i <- nrow(x)
  f_col <- if ("F" %in% colnames(x)) "F" else "F value"
  p_col <- "Pr(>F)"
  comparison_rows <- ms_anova_comparison_rows(x)
  comparison_columns <- ms_anova_comparison_columns(comparison_rows)

  f_val  <- ms_safe_numeric(x[[f_col]][i])
  df1    <- ms_safe_numeric(x[["Df"]][i])
  df2    <- ms_safe_numeric(x[["Res.Df"]][i])
  pval   <- ms_safe_numeric(x[[p_col]][i])
  sum_col <- if ("Sum of Sq" %in% colnames(x)) "Sum of Sq" else NA_character_
  rss_col <- if ("RSS" %in% colnames(x)) "RSS" else NA_character_

  fields <- list(
    statistic = list(
      name  = "F",
      value = f_val,
      df    = I(c(df1, df2))
    ),
    p_value = pval,
    comparison = paste0("model ", max(1L, i - 1L), " vs. model ", i),
    model_count = nrow(x),
    df_change = df1,
    residual_df = df2,
    table_type = "model_comparison",
    rows = comparison_rows,
    columns = comparison_columns,
    table_note = "Each row compares the listed model with the previous nested model. F tests use residual sums of squares from the ANOVA comparison."
  )
  if (!is.na(sum_col)) {
    sum_sq <- ms_safe_numeric(x[[sum_col]][i])
    if (!is.na(sum_sq)) fields$sum_of_sq <- sum_sq
  }
  if (!is.na(rss_col)) {
    rss <- ms_safe_numeric(x[[rss_col]][i])
    if (!is.na(rss)) fields$residual_sum_sq <- rss
  }

  list(
    type       = "anova_model_comparison",
    type_label = paste0("Model comparison (", nrow(x), " models)"),
    fields     = fields
  )
}

ms_anova_comparison_rows <- function(x) {
  if (nrow(x) < 2L) return(list())
  f_col <- if ("F" %in% colnames(x)) {
    "F"
  } else if ("F value" %in% colnames(x)) {
    "F value"
  } else {
    NA_character_
  }
  p_col <- if ("Pr(>F)" %in% colnames(x)) "Pr(>F)" else NA_character_
  sum_col <- if ("Sum of Sq" %in% colnames(x)) "Sum of Sq" else NA_character_
  rss_col <- if ("RSS" %in% colnames(x)) "RSS" else NA_character_

  lapply(seq.int(2L, nrow(x)), function(i) {
    df1 <- if ("Df" %in% colnames(x)) ms_safe_numeric(x[["Df"]][i]) else NA_real_
    df2 <- if ("Res.Df" %in% colnames(x)) ms_safe_numeric(x[["Res.Df"]][i]) else NA_real_
    df <- if (!is.na(df1) && !is.na(df2)) {
      paste0("(", payload_format_df(df1), ", ", payload_format_df(df2), ")")
    } else {
      NA_character_
    }

    row <- list(
      comparison = paste0("model ", i - 1L, " vs. model ", i),
      df = df,
      model = paste("Model", i)
    )
    if (!is.na(f_col)) {
      f_val <- ms_safe_numeric(x[[f_col]][i])
      if (!is.na(f_val)) row$f <- f_val
    }
    if (!is.na(p_col)) {
      p_val <- ms_safe_numeric(x[[p_col]][i])
      if (!is.na(p_val)) row$p <- p_val
    }
    if (!is.na(sum_col)) {
      sum_sq <- ms_safe_numeric(x[[sum_col]][i])
      if (!is.na(sum_sq)) row$sum_of_sq <- sum_sq
    }
    if (!is.na(rss_col)) {
      rss <- ms_safe_numeric(x[[rss_col]][i])
      if (!is.na(rss)) row$residual_sum_sq <- rss
    }
    row
  })
}

ms_anova_comparison_columns <- function(rows) {
  has_key <- function(key) {
    any(vapply(rows, function(row) {
      !is.null(row[[key]]) && length(row[[key]]) > 0L &&
        !(length(row[[key]]) == 1L && is.na(row[[key]]))
    }, logical(1)))
  }

  columns <- list(
    list(key = "comparison", label = "Comparison", format = "text"),
    list(key = "df", label = "df", format = "df")
  )
  if (has_key("f")) {
    columns <- c(columns, list(list(key = "f", label = "F", format = "statistic")))
  }
  if (has_key("p")) {
    columns <- c(columns, list(list(key = "p", label = "p", format = "pvalue")))
  }
  if (has_key("sum_of_sq")) {
    columns <- c(columns, list(list(key = "sum_of_sq", label = "Sum of Sq", format = "number")))
  }
  if (has_key("residual_sum_sq")) {
    columns <- c(columns, list(list(key = "residual_sum_sq", label = "Residual SS", format = "number")))
  }
  columns <- c(columns, list(list(key = "model", label = "Model", format = "text")))
  columns
}

# Single-model case: pick the last non-Residuals row unless focal/controls
# identify a more specific headline term.
ms_anova_single <- function(x, focal = NULL, controls = NULL) {
  rn <- rownames(x)
  is_resid <- grepl("^Residuals?$", rn, ignore.case = TRUE)
  is_intercept <- grepl("^\\(?Intercept\\)?$", rn, ignore.case = TRUE)
  candidate <- which(!is_resid & !is_intercept)
  if (length(candidate) == 0L) {
    stop("anova table has no non-Residuals rows to extract", call. = FALSE)
  }
  i <- ms_anova_select_term_index(rn, candidate, focal = focal, controls = controls)
  term <- rn[i]

  f_col <- if ("F" %in% colnames(x)) "F" else "F value"
  p_col <- "Pr(>F)"

  f_val <- ms_safe_numeric(x[[f_col]][i])
  df1   <- ms_safe_numeric(x[["Df"]][i])
  sum_sq <- if ("Sum Sq" %in% colnames(x)) {
    ms_safe_numeric(x[["Sum Sq"]][i])
  } else NA_real_
  mean_sq <- if ("Mean Sq" %in% colnames(x)) {
    ms_safe_numeric(x[["Mean Sq"]][i])
  } else NA_real_

  # Residual df is on the Residuals row when present
  resid_idx <- which(is_resid)
  df2 <- if (length(resid_idx) > 0L) {
    ms_safe_numeric(x[["Df"]][resid_idx[1]])
  } else {
    NA_real_
  }
  resid_sum_sq <- if (length(resid_idx) > 0L && "Sum Sq" %in% colnames(x)) {
    ms_safe_numeric(x[["Sum Sq"]][resid_idx[1]])
  } else NA_real_
  resid_mean_sq <- if (length(resid_idx) > 0L && "Mean Sq" %in% colnames(x)) {
    ms_safe_numeric(x[["Mean Sq"]][resid_idx[1]])
  } else NA_real_

  pval <- ms_safe_numeric(x[[p_col]][i])

  fields <- list(
    statistic = list(
      name  = "F",
      value = f_val,
      df    = I(c(df1, df2))
    ),
    p_value = pval,
    term = term,
    df_effect = df1
  )
  if (!is.na(df2)) fields$df_error <- df2
  if (!is.na(sum_sq)) fields$sum_sq <- sum_sq
  if (!is.na(mean_sq)) fields$mean_sq <- mean_sq
  if (!is.na(resid_sum_sq)) fields$residual_sum_sq <- resid_sum_sq
  if (!is.na(resid_mean_sq)) fields$residual_mean_sq <- resid_mean_sq

  # P2: partial eta-squared for the focal term. SS_effect / (SS_effect +
  # SS_residual). For a single-term anova() table this equals
  # effectsize::eta_squared(..., partial = TRUE); for multi-term tables
  # the focal-term value matches the same row. Requires both sum_sq
  # columns; omit when either is missing (e.g. anova() on a list of
  # contrasts that drops Sum Sq).
  if (!is.na(sum_sq) && !is.na(resid_sum_sq) && (sum_sq + resid_sum_sq) > 0) {
    eta_p <- ms_safe_numeric(sum_sq / (sum_sq + resid_sum_sq))
    if (!is.na(eta_p)) {
      fields$effect <- list(name = "eta_sq_partial", value = eta_p)
    }
  }

  # Additive: per-term stats for two-way / factorial ANOVA narration and
  # the multi-row Tables projection. The card stays card_kind = "inline"
  # and every existing field-level scalar above is untouched, so the
  # EXTRACTED FIELDS panel, Result tab, and partial-\u03b7\u00b2 subscript all keep
  # working. fields$all_terms is consumed by the web bridge's inline anova
  # archetype to render one sentence per term; if it's absent (or a
  # consumer doesn't know about it), behavior falls back to the original
  # single-focal-term path.
  has_sumsq <- "Sum Sq" %in% colnames(x)
  has_meansq <- "Mean Sq" %in% colnames(x)
  all_terms <- lapply(candidate, function(j) {
    term_j <- rn[j]
    df_j   <- ms_safe_numeric(x[["Df"]][j])
    f_j    <- ms_safe_numeric(x[[f_col]][j])
    p_j    <- ms_safe_numeric(x[[p_col]][j])
    ss_j   <- if (has_sumsq) ms_safe_numeric(x[["Sum Sq"]][j]) else NA_real_
    ms_j   <- if (has_meansq) ms_safe_numeric(x[["Mean Sq"]][j]) else NA_real_
    eta_j  <- if (!is.na(ss_j) && !is.na(resid_sum_sq) &&
                  (ss_j + resid_sum_sq) > 0) {
      ms_safe_numeric(ss_j / (ss_j + resid_sum_sq))
    } else NA_real_

    row <- list(
      term      = term_j,
      term_type = if (grepl(":", term_j, fixed = TRUE)) "interaction" else "main",
      is_focal  = identical(j, i),
      df1       = df_j,
      df2       = df2,
      f         = f_j,
      p_value   = p_j
    )
    if (!is.na(ss_j))  row$sum_sq  <- ss_j
    if (!is.na(ms_j))  row$mean_sq <- ms_j
    if (!is.na(eta_j)) {
      # Keep eta_sq as a compatibility alias for older browser/table code,
      # but expose the measure explicitly for report wording and new tables.
      row$eta_sq_partial <- eta_j
      row$eta_sq <- eta_j
      row$effect <- list(name = "eta_sq_partial", value = eta_j)
    }
    row
  })
  fields$all_terms <- all_terms

  list(
    type       = "anova_single_model",
    type_label = paste0(
      "ANOVA \u00b7 ",
      ms_anova_single_label_suffix(
        rn[candidate],
        selected_term = term,
        focal = focal,
        controls = controls
      )
    ),
    fields     = fields
  )
}

ms_anova_single_label_suffix <- function(terms, selected_term = NULL,
                                         focal = NULL, controls = NULL) {
  terms <- trimws(as.character(terms %||% character(0)))
  terms <- terms[nzchar(terms)]
  focal_keys <- ms_role_keys(focal)
  control_keys <- ms_role_keys(controls)

  if (length(terms) <= 1L || length(focal_keys) > 0L || length(control_keys) > 0L) {
    return(paste0("focal term: ", selected_term %||% terms[[length(terms)]] %||% "term"))
  }

  if (length(terms) <= 3L) {
    return(paste0("terms: ", paste(terms, collapse = ", ")))
  }

  paste0(
    length(terms),
    " terms: ",
    paste(utils::head(terms, 3L), collapse = ", "),
    ", ..."
  )
}

ms_anova_select_term_index <- function(terms, candidate, focal = NULL, controls = NULL) {
  focal_keys <- ms_role_keys(focal)
  if (length(focal_keys) > 0L) {
    keys <- vapply(terms[candidate], ms_term_key, character(1))
    hit <- which(keys %in% focal_keys)
    if (length(hit) > 0L) return(candidate[hit[[1]]])
  }

  control_keys <- ms_role_keys(controls)
  if (length(control_keys) > 0L) {
    keys <- vapply(terms[candidate], ms_term_key, character(1))
    non_control <- which(!keys %in% control_keys)
    if (length(non_control) > 0L) return(candidate[non_control[[length(non_control)]]])
  }

  candidate[length(candidate)]
}

ms_anova_sum_squares_info <- function(x) {
  heading <- paste(attr(x, "heading") %||% character(0), collapse = " ")
  heading <- trimws(gsub("\\s+", " ", heading))
  if (!nzchar(heading)) return(NULL)

  if (grepl("type\\s+iii", heading, ignore.case = TRUE, perl = TRUE)) {
    return(list(
      type = "type_iii",
      label = "Type III",
      note = "Type III tests depend on the contrast coding used when the model was fitted."
    ))
  }
  if (grepl("type\\s+ii", heading, ignore.case = TRUE, perl = TRUE)) {
    return(list(type = "type_ii", label = "Type II"))
  }
  if (grepl("type\\s+i", heading, ignore.case = TRUE, perl = TRUE) ||
      grepl("^analysis of variance table", heading, ignore.case = TRUE, perl = TRUE)) {
    return(list(
      type = "type_i_sequential",
      label = "Type I (sequential)",
      note = "Sequential tests depend on the order of terms in the model formula."
    ))
  }
  NULL
}

ms_anova_type_label_with_ss <- function(type_label, ss_label) {
  type_label <- as.character(type_label %||% "ANOVA")
  ss_label <- as.character(ss_label %||% "")
  if (!nzchar(ss_label) || grepl(ss_label, type_label, fixed = TRUE)) {
    return(type_label)
  }
  sub("^ANOVA", paste0("ANOVA \u00b7 ", ss_label), type_label)
}

ms_anova_term_roles <- function(x, focal = NULL, controls = NULL) {
  rn <- rownames(x)
  is_resid <- grepl("^Residuals?$", rn, ignore.case = TRUE)
  is_intercept <- grepl("^\\(?Intercept\\)?$", rn, ignore.case = TRUE)
  term_names <- rn[!is_resid & !is_intercept]
  terms <- lapply(term_names, function(term) {
    list(
      name = term,
      label = ms_model_clean_term(term),
      role = "term",
      type = if (grepl(":", term, fixed = TRUE)) "interaction" else "other"
    )
  })

  roles <- ms_assign_term_roles(
    terms,
    focal = focal,
    controls = controls,
    infer_covariate_roles = TRUE
  )
  terms <- roles$terms

  focal_terms <- vapply(Filter(function(t) identical(t$role, "focal"), terms),
                        function(t) t$name, character(1))
  control_terms <- vapply(Filter(function(t) identical(t$role, "control"), terms),
                          function(t) t$name, character(1))
  predictor_terms <- if (length(focal_terms) > 0L) {
    Filter(function(t) identical(t$role, "focal"), terms)
  } else {
    Filter(function(t) !identical(t$role, "control"), terms)
  }
  # Use $name (raw identifier) so the field value lists variables verbatim.
  predictor <- ms_model_term_phrase(vapply(predictor_terms, function(t) {
    t$name %||% t$label
  }, character(1)))

  # Explicit controls on an ANOVA table usually mean ANCOVA-style
  # reporting, even though the anova object alone does not preserve
  # variable classes.
  model_kind <- roles$model_kind
  if (length(control_terms) > 0L && length(focal_terms) > 0L) {
    model_kind <- "ancova"
  }

  list(
    terms = terms,
    focal_terms = focal_terms,
    control_terms = control_terms,
    predictor = if (nzchar(predictor)) predictor else NULL,
    model_kind = model_kind
  )
}

ms_anova_vars_from_call <- function(call_str) {
  if (is.null(call_str) || length(call_str) == 0L || is.na(call_str[1])) {
    return(list())
  }

  call_str <- trimws(gsub("\\s+", " ", paste(call_str, collapse = " ")))
  m <- regexec("(?:^|[(,])\\s*([^,()~]+?)\\s*~\\s*([^,]+)", call_str, perl = TRUE)
  hit <- regmatches(call_str, m)[[1]]
  if (length(hit) < 3L) return(list())

  lhs <- ms_anova_clean_formula_side(hit[2])
  rhs <- ms_anova_clean_formula_side(hit[3])
  predictor <- ms_anova_primary_rhs_term(rhs)

  out <- list()
  if (nzchar(lhs)) out$outcome <- lhs
  if (nzchar(predictor)) out$predictor <- predictor
  out
}

ms_anova_outcome_from_heading <- function(x) {
  heading <- attr(x, "heading")
  if (is.null(heading) || length(heading) == 0L) return(NULL)

  heading <- paste(as.character(heading), collapse = "\n")
  m <- regexec("(?im)^\\s*Response:\\s*(.+?)\\s*$", heading, perl = TRUE)
  hit <- regmatches(heading, m)[[1]]
  if (length(hit) < 2L) return(NULL)

  outcome <- ms_anova_clean_formula_side(hit[2])
  if (nzchar(outcome)) outcome else NULL
}

ms_anova_primary_rhs_term <- function(rhs) {
  rhs <- strsplit(rhs %||% "", "[+*:|]", perl = TRUE)[[1]][1] %||% ""
  ms_anova_clean_formula_side(rhs)
}

ms_anova_clean_formula_side <- function(value) {
  value <- trimws(as.character(value %||% ""))
  value <- sub("^formula\\s*=\\s*", "", value)
  value <- sub("^`(.+)`$", "\\1", value)
  value <- sub("^.*\\$", "", value)
  value <- gsub("\\[[^]]*\\]", "", value)
  value <- gsub(
    "^(as\\.factor|factor|scale|I)\\(([^()]+)\\)$",
    "\\2",
    value,
    perl = TRUE
  )
  trimws(value)
}

ms_anova_means_figure_data <- function(x, conf.level = 0.95, max_levels = 12L) {
  if (!inherits(x, "aov") || inherits(x, c("aovlist", "maov", "manova"))) return(NULL)

  tf <- tryCatch(stats::terms(x), error = function(e) NULL)
  mf <- tryCatch(stats::model.frame(x), error = function(e) NULL)
  if (is.null(tf) || is.null(mf) || nrow(mf) < 2L) return(NULL)

  term_labels <- attr(tf, "term.labels") %||% character(0)
  orders <- attr(tf, "order") %||% integer(0)
  response_idx <- attr(tf, "response") %||% 0L
  if (length(term_labels) != 1L || length(orders) != 1L || orders[[1L]] != 1L) return(NULL)
  if (response_idx < 1L || response_idx > ncol(mf)) return(NULL)

  term <- term_labels[[1L]]
  group_values <- if (term %in% names(mf)) {
    mf[[term]]
  } else {
    vars <- tryCatch(all.vars(stats::as.formula(paste0("~", term))),
                     error = function(e) character(0))
    if (length(vars) == 1L && vars[[1L]] %in% names(mf)) {
      mf[[vars[[1L]]]]
    } else {
      NULL
    }
  }
  if (is.null(group_values)) return(NULL)
  if (!(is.factor(group_values) || is.character(group_values) || is.logical(group_values))) {
    return(NULL)
  }

  y <- mf[[response_idx]]
  if (!is.numeric(y) || length(y) != length(group_values)) return(NULL)
  ok <- is.finite(y) & !is.na(group_values)
  y <- y[ok]
  group_values <- group_values[ok]
  if (length(y) < 2L) return(NULL)

  levels <- ms_anova_observed_levels(group_values)
  if (length(levels) < 2L || length(levels) > max_levels) return(NULL)

  residual_df <- ms_safe_numeric(stats::df.residual(x))
  residual_mean_sq <- NA_real_
  if (!is.na(residual_df) && residual_df > 0) {
    rss <- sum(stats::residuals(x)^2, na.rm = TRUE)
    residual_mean_sq <- ms_safe_numeric(rss / residual_df)
  }
  critical <- if (!is.na(residual_df) && residual_df > 0) {
    stats::qt(1 - ((1 - conf.level) / 2), df = residual_df)
  } else {
    stats::qnorm(1 - ((1 - conf.level) / 2))
  }

  groups <- lapply(levels, function(level) {
    idx <- as.character(group_values) == as.character(level)
    vals <- y[idx]
    n <- length(vals)
    mean_value <- ms_safe_numeric(mean(vals, na.rm = TRUE))
    sd_value <- ms_safe_numeric(stats::sd(vals, na.rm = TRUE))
    se_value <- if (!is.na(residual_mean_sq) && n > 0) {
      ms_safe_numeric(sqrt(residual_mean_sq / n))
    } else if (!is.na(sd_value) && n > 0) {
      ms_safe_numeric(sd_value / sqrt(n))
    } else {
      NA_real_
    }

    row <- list(
      level = as.character(level),
      label = as.character(level),
      n = n,
      mean = mean_value
    )
    if (!is.na(sd_value)) row$sd <- sd_value
    if (!is.na(se_value)) {
      row$se <- se_value
      row$ci_lower <- ms_safe_numeric(mean_value - critical * se_value)
      row$ci_upper <- ms_safe_numeric(mean_value + critical * se_value)
    }
    row
  })
  groups <- Filter(function(row) {
    is.list(row) &&
      !is.null(row$level) &&
      !is.null(row$mean) &&
      !is.na(ms_safe_numeric(row$mean))
  }, groups)
  if (length(groups) < 2L) return(NULL)

  outcome <- ms_model_clean_term(names(mf)[[response_idx]])
  factor_variable <- ms_anova_clean_formula_side(term)
  if (!nzchar(factor_variable)) factor_variable <- ms_model_clean_term(term)

  out <- list(
    mean_kind = "estimated",
    source = "aov_one_way",
    factor = list(
      variable = factor_variable,
      term = term,
      label = ms_model_clean_term(term),
      levels = lapply(levels, function(level) {
        list(value = as.character(level), label = as.character(level))
      })
    ),
    groups = groups,
    outcome = outcome,
    y_label = if (nzchar(outcome)) paste("Mean", outcome) else "Mean outcome",
    ci_level = conf.level,
    ci_method = if (!is.na(residual_mean_sq)) "residual_mean_square" else "group_standard_error"
  )
  if (!is.na(residual_df)) out$residual_df <- residual_df
  if (!is.na(residual_mean_sq)) out$residual_mean_sq <- residual_mean_sq
  out
}

ms_anova_observed_levels <- function(values) {
  levels <- if (is.factor(values)) {
    observed <- unique(as.character(values[!is.na(values)]))
    levels(values)[levels(values) %in% observed]
  } else if (is.logical(values)) {
    as.character(sort(unique(values[!is.na(values)])))
  } else {
    unique(as.character(values[!is.na(values)]))
  }
  ms_anova_natural_order_levels(levels)
}

ms_anova_natural_order_levels <- function(levels) {
  levels <- as.character(levels)
  if (length(levels) < 2L) return(levels)

  numeric_levels <- suppressWarnings(as.numeric(levels))
  if (all(!is.na(numeric_levels)) && !anyDuplicated(numeric_levels)) {
    return(levels[order(numeric_levels)])
  }

  ranks <- vapply(levels, ms_anova_level_rank, numeric(1))
  if (all(!is.na(ranks)) && !anyDuplicated(ranks)) {
    return(levels[order(ranks)])
  }

  levels
}

ms_anova_level_rank <- function(level) {
  key <- tolower(trimws(as.character(level %||% "")))
  key <- gsub("[^a-z0-9]+", "", key)
  if (!nzchar(key)) return(NA_real_)

  temporal <- c(
    baseline = 0, base = 0,
    before = 1, prior = 1, pre = 1, pretest = 1,
    after = 2, post = 2, posttest = 2,
    followup = 3, followup1 = 3.01, followup2 = 3.02, followup3 = 3.03
  )
  intensity <- c(
    none = 0, zero = 0,
    low = 1, lower = 1,
    medium = 2, mid = 2, moderate = 2,
    high = 3, higher = 3
  )
  if (key %in% names(temporal)) return(unname(temporal[[key]]))
  if (key %in% names(intensity)) return(unname(intensity[[key]]))

  t_match <- regexec("^t([0-9]+)$", key)
  t_hit <- regmatches(key, t_match)[[1L]]
  if (length(t_hit) == 2L) return(100 + as.numeric(t_hit[[2L]]))

  follow_match <- regexec("^followup([0-9]+)$", key)
  follow_hit <- regmatches(key, follow_match)[[1L]]
  if (length(follow_hit) == 2L) return(3 + as.numeric(follow_hit[[2L]]) / 100)

  NA_real_
}

# --- One-way post-hoc ---------------------------------------------------------
#
# Two additions for a one-way aov, both reusing existing machinery so they stay
# in lock-step with the nonparametric tests and the emmeans pairwise path:
#
#   * a raw-data distribution figure (`nonparametric_group_plot`) built from the
#     model frame via ms_htest_group_descriptives() -- identical shape to the
#     Kruskal-Wallis boxplot, plus `mean_overlay = TRUE` so the renderer marks
#     the group mean (the quantity the ANOVA actually tests).
#   * pairwise contrasts via emmeans, run straight through mellio_payload.emmGrid
#     so the APA pairwise table + `pairwise_forest` figure come for free and
#     respect whatever multiple-comparison `adjust` is requested.
#
# Significance brackets are composed Lab-side from the pairwise rows' level_1 /
# level_2, so nothing bracket-specific is emitted here.

ms_anova_oneway_emm_contrasts <- function(x, focal_term, adjust = "tukey") {
  if (!requireNamespace("emmeans", quietly = TRUE)) return(NULL)
  if (is.null(focal_term) || !nzchar(focal_term)) return(NULL)
  emm <- tryCatch(
    suppressMessages(emmeans::emmeans(x, specs = focal_term)),
    error = function(e) NULL
  )
  if (is.null(emm)) return(NULL)
  tryCatch(
    suppressMessages(emmeans::contrast(emm, method = "pairwise", adjust = adjust)),
    error = function(e) NULL
  )
}

ms_anova_oneway_pairwise_payload <- function(x, focal_term, adjust = "tukey",
                                             call = NULL) {
  grid <- ms_anova_oneway_emm_contrasts(x, focal_term, adjust = adjust)
  if (is.null(grid)) return(NULL)
  call <- call %||% sprintf(
    "emmeans::contrast(emmeans::emmeans(model, ~%s), \"pairwise\", adjust = \"%s\")",
    focal_term, adjust
  )
  tryCatch(mellio_payload(grid, .call = call), error = function(e) NULL)
}

ms_anova_raw_group_figure <- function(x, means_figure, max_points = 500L) {
  if (is.null(means_figure)) return(NULL)
  tf <- tryCatch(stats::terms(x), error = function(e) NULL)
  mf <- tryCatch(stats::model.frame(x), error = function(e) NULL)
  if (is.null(tf) || is.null(mf)) return(NULL)
  response_idx <- attr(tf, "response") %||% 0L
  if (response_idx < 1L || response_idx > ncol(mf)) return(NULL)

  term <- means_figure$factor$term
  group_values <- if (!is.null(term) && term %in% names(mf)) {
    mf[[term]]
  } else {
    vars <- tryCatch(all.vars(stats::as.formula(paste0("~", term))),
                     error = function(e) character(0))
    if (length(vars) == 1L && vars[[1L]] %in% names(mf)) mf[[vars[[1L]]]] else NULL
  }
  y <- mf[[response_idx]]
  if (is.null(group_values) || !is.numeric(y) || length(y) != length(group_values)) {
    return(NULL)
  }
  keep <- is.finite(y) & !is.na(group_values)
  y <- y[keep]
  group_chr <- as.character(group_values)[keep]
  if (length(y) < 2L) return(NULL)

  levels_chr <- vapply(means_figure$factor$levels,
                       function(level) as.character(level$value), character(1))
  levels_chr <- levels_chr[levels_chr %in% unique(group_chr)]
  if (length(levels_chr) < 2L) return(NULL)
  g_factor <- factor(group_chr, levels = levels_chr)

  groups <- ms_htest_group_descriptives(y, g_factor, levels = levels_chr)
  if (is.null(groups) || length(groups) < 2L) return(NULL)

  observations <- lapply(seq_along(y), function(i) {
    list(id = i, group = as.character(g_factor[[i]]), value = ms_safe_numeric(y[[i]]))
  })
  total_n <- length(observations)
  truncated <- total_n > max_points
  if (truncated) observations <- observations[seq_len(max_points)]

  outcome <- as.character(means_figure$outcome %||% "")

  list(
    source = "aov_one_way",
    plot_kind = "independent_groups",
    factor = means_figure$factor,
    groups = groups,
    observations = observations,
    point_display = list(
      total_n = total_n,
      included_n = length(observations),
      truncated = truncated
    ),
    outcome = if (nzchar(outcome)) outcome else NULL,
    y_label = if (nzchar(outcome)) outcome else "Value",
    mean_overlay = TRUE
  )
}

ms_anova_attach_oneway_posthoc <- function(payload, x, means_figure) {
  payload$figure_data <- payload$figure_data %||% list()
  focal_term <- means_figure$factor$term

  raw_group <- tryCatch(ms_anova_raw_group_figure(x, means_figure),
                        error = function(e) NULL)
  if (!is.null(raw_group)) {
    payload$figure_data$nonparametric_group_plot <- raw_group
    payload <- ms_add_available_figure(
      payload, type = "nonparametric_group_plot", label = "Raw data", default = FALSE
    )
  }

  ph <- tryCatch(ms_anova_oneway_pairwise_payload(x, focal_term, adjust = "tukey"),
                 error = function(e) NULL)
  if (!is.null(ph)) {
    if (!is.null(ph$figure_data$pairwise_forest)) {
      payload$figure_data$pairwise_forest <- ph$figure_data$pairwise_forest
      payload <- ms_add_available_figure(
        payload, type = "pairwise_forest", label = "Pairwise contrasts", default = FALSE
      )
    }
    if (!is.null(ph$fields)) payload$fields$pairwise <- ph$fields
    payload$packages <- ms_packages_basic(extras = "emmeans")
  }
  payload
}

#' Recompute one-way ANOVA pairwise contrasts under a different correction
#'
#' Lets the Lab switch the multiple-comparison adjustment without refitting the
#' model: it reuses an aov object already held in the (WebR) session and runs the
#' contrasts straight through [mellio_payload()]'s emmGrid method.
#'
#' @param fit A fitted one-way `aov`/`lm` object.
#' @param by Optional grouping factor name; inferred from the model when `NULL`.
#' @param adjust Multiple-comparison adjustment passed to `emmeans::contrast()`
#'   (e.g. `"tukey"`, `"bonferroni"`, `"holm"`, `"sidak"`, `"fdr"`, `"none"`).
#' @return A pairwise-comparisons payload (table fields + `pairwise_forest`).
#' @export
mellio_anova_pairwise <- function(fit, by = NULL, adjust = "tukey") {
  focal <- by
  if (is.null(focal)) {
    mfig <- ms_anova_means_figure_data(fit)
    focal <- if (!is.null(mfig)) mfig$factor$term else NULL
  }
  if (is.null(focal) || !nzchar(focal)) {
    stop("mellio_anova_pairwise(): could not determine the grouping factor.", call. = FALSE)
  }
  payload <- ms_anova_oneway_pairwise_payload(fit, focal, adjust = adjust)
  if (is.null(payload)) {
    stop("mellio_anova_pairwise(): pairwise contrasts unavailable (is emmeans installed?).",
         call. = FALSE)
  }
  payload
}

ms_aovlist_metadata <- function(x, max_levels = 12L) {
  if (!inherits(x, "aovlist")) return(list())

  tf <- tryCatch(attr(x, "terms"), error = function(e) NULL)
  mf <- tryCatch(stats::model.frame(x), error = function(e) NULL)
  if (is.null(tf) || is.null(mf) || !is.data.frame(mf) || nrow(mf) < 1L) {
    return(list())
  }

  response_idx <- attr(tf, "response") %||% 0L
  if (response_idx < 1L || response_idx > ncol(mf)) return(list())

  term_labels <- attr(tf, "term.labels") %||% character(0)
  error_terms <- term_labels[grepl("^Error\\(", term_labels)]
  within_terms <- setdiff(term_labels, error_terms)
  within_terms <- within_terms[nzchar(within_terms)]
  if (!length(within_terms) || length(error_terms) != 1L) return(list())

  within_variables <- ms_aovlist_within_variables(within_terms, mf)
  subject_var <- ms_aovlist_subject_variable_multi(
    error_terms[[1L]],
    within_terms,
    within_variables,
    mf
  )

  out <- list(
    model_kind = "repeated_measures_anova",
    outcome = ms_model_clean_term(names(mf)[[response_idx]]),
    within_terms = within_terms
  )

  if (!is.null(subject_var) && subject_var %in% names(mf)) {
    subject_values <- mf[[subject_var]]
    subject_values <- subject_values[!is.na(subject_values)]
    out$subject <- ms_model_clean_term(subject_var)
    out$subject_variable <- subject_var
    out$n <- as.integer(length(unique(as.character(subject_values))))
  }

  factors <- ms_aovlist_factor_metadata(within_terms, within_variables, mf, max_levels = max_levels)
  if (length(factors)) out$factors <- factors
  sphericity_note <- ms_aovlist_sphericity_note(factors)
  if (nzchar(sphericity_note)) out$sphericity_note <- sphericity_note
  out$table_note <- ms_aovlist_table_note(out, sphericity_note)

  out
}

ms_aovlist_within_variables <- function(within_terms, mf) {
  vars <- character(0)
  for (term in within_terms) {
    term <- trimws(as.character(term %||% ""))
    if (!nzchar(term)) next
    if (term %in% names(mf)) {
      vars <- c(vars, term)
      next
    }
    cleaned <- ms_anova_clean_formula_side(term)
    if (nzchar(cleaned) && cleaned %in% names(mf)) {
      vars <- c(vars, cleaned)
      next
    }
    raw_vars <- tryCatch(all.vars(stats::as.formula(paste0("~", term))),
                         error = function(e) character(0))
    raw_vars <- raw_vars[raw_vars %in% names(mf)]
    if (length(raw_vars)) vars <- c(vars, raw_vars)
  }
  unique(vars)
}

ms_aovlist_subject_variable_multi <- function(error_term, within_terms, within_variables, mf) {
  inner <- sub("^Error\\((.*)\\)$", "\\1", error_term)
  vars <- tryCatch(all.vars(stats::as.formula(paste0("~", inner))),
                   error = function(e) character(0))
  vars <- vars[vars %in% names(mf)]
  excludes <- unique(c(
    within_variables,
    unlist(lapply(within_terms, function(term) {
      tryCatch(all.vars(stats::as.formula(paste0("~", term))),
               error = function(e) character(0))
    }), use.names = FALSE),
    vapply(within_terms, ms_anova_clean_formula_side, character(1))
  ))
  vars <- vars[!vars %in% excludes]
  if (!length(vars)) return(NULL)
  vars[[1L]]
}

ms_aovlist_factor_metadata <- function(within_terms, within_variables, mf, max_levels = 12L) {
  out <- list()
  if (!length(within_variables)) return(out)

  for (variable in within_variables) {
    if (!variable %in% names(mf)) next
    values <- mf[[variable]]
    if (!(is.factor(values) || is.character(values) || is.logical(values))) next
    levels <- ms_anova_observed_levels(values)
    if (length(levels) < 2L || length(levels) > max_levels) next
    term_idx <- match(variable, within_terms)
    term <- if (!is.na(term_idx)) within_terms[[term_idx]] else variable
    out[[length(out) + 1L]] <- list(
      variable = variable,
      term = term,
      label = ms_model_clean_term(term),
      role = "within",
      levels = lapply(levels, function(level) {
        list(value = as.character(level), label = ms_model_clean_term(level))
      })
    )
  }

  out
}

ms_aovlist_sphericity_note <- function(factors) {
  if (!length(factors)) return("")
  level_counts <- vapply(factors, function(factor) {
    length(factor$levels %||% list())
  }, integer(1))
  if (length(level_counts) && all(level_counts <= 2L)) {
    return("Sphericity is not applicable because all within-subject factors have two levels.")
  }
  paste(
    "Sphericity tests and Greenhouse-Geisser/Huynh-Feldt corrections",
    "are not available for this result; reported within-subject tests are uncorrected."
  )
}

ms_aovlist_table_note <- function(metadata, sphericity_note = "") {
  parts <- c(
    "All listed factors are within-subjects.",
    "Degrees of freedom for each effect and its matching error term are shown; F tests use the matching error stratum.",
    "Partial \u03b7\u00b2 is computed from each effect sum of squares and its matching error sum of squares."
  )
  if (nzchar(sphericity_note)) parts <- c(parts, sphericity_note)
  paste(parts, collapse = " ")
}

ms_repeated_measures_means_figure_data <- function(x, tidy_df = NULL,
                                                   conf.level = 0.95,
                                                   max_levels = 12L) {
  if (!inherits(x, "aovlist")) return(NULL)

  tf <- tryCatch(attr(x, "terms"), error = function(e) NULL)
  mf <- tryCatch(stats::model.frame(x), error = function(e) NULL)
  if (is.null(tf) || is.null(mf) || nrow(mf) < 2L) return(NULL)

  response_idx <- attr(tf, "response") %||% 0L
  if (response_idx < 1L || response_idx > ncol(mf)) return(NULL)
  y <- mf[[response_idx]]
  if (!is.numeric(y)) return(NULL)

  terms <- attr(tf, "term.labels") %||% character(0)
  error_terms <- terms[grepl("^Error\\(", terms)]
  within_terms <- setdiff(terms, error_terms)
  if (length(within_terms) != 1L || length(error_terms) != 1L) return(NULL)

  within_term <- within_terms[[1L]]
  within_var <- ms_anova_clean_formula_side(within_term)
  if (!nzchar(within_var) || !within_var %in% names(mf)) return(NULL)
  within_values <- mf[[within_var]]
  if (!(is.factor(within_values) || is.character(within_values) || is.logical(within_values))) {
    return(NULL)
  }

  subject_var <- ms_aovlist_subject_variable(error_terms[[1L]], within_var, mf)
  if (is.null(subject_var) || !subject_var %in% names(mf)) return(NULL)
  subject_values <- mf[[subject_var]]

  ok <- is.finite(y) & !is.na(within_values) & !is.na(subject_values)
  y <- y[ok]
  within_values <- within_values[ok]
  subject_values <- subject_values[ok]
  if (length(y) < 2L) return(NULL)

  levels <- ms_anova_observed_levels(within_values)
  if (length(levels) < 2L || length(levels) > max_levels) return(NULL)
  subjects <- unique(as.character(subject_values))
  if (length(subjects) < 2L) return(NULL)

  cell_frame <- data.frame(
    subject = as.character(subject_values),
    level = factor(as.character(within_values), levels = levels),
    y = y,
    stringsAsFactors = FALSE
  )
  cells <- stats::aggregate(
    y ~ subject + level,
    data = cell_frame,
    FUN = function(value) mean(value, na.rm = TRUE),
    drop = FALSE
  )

  mat <- matrix(NA_real_, nrow = length(subjects), ncol = length(levels),
                dimnames = list(subjects, levels))
  for (i in seq_len(nrow(cells))) {
    subject <- as.character(cells$subject[[i]])
    level <- as.character(cells$level[[i]])
    if (!subject %in% rownames(mat) || !level %in% colnames(mat)) next
    mat[subject, level] <- ms_safe_numeric(cells$y[[i]])
  }
  complete <- stats::complete.cases(mat)
  mat <- mat[complete, , drop = FALSE]
  if (nrow(mat) < 2L || ncol(mat) < 2L) return(NULL)

  means <- colMeans(mat, na.rm = TRUE)
  subject_means <- rowMeans(mat, na.rm = TRUE)
  grand_mean <- mean(mat, na.rm = TRUE)
  normalized <- sweep(mat, 1L, subject_means, FUN = "-") + grand_mean
  morey <- sqrt(ncol(mat) / (ncol(mat) - 1L))
  df_ci <- nrow(mat) - 1L
  critical <- stats::qt(1 - ((1 - conf.level) / 2), df = df_ci)
  se <- apply(normalized, 2L, stats::sd, na.rm = TRUE) / sqrt(nrow(mat)) * morey

  groups <- lapply(seq_along(levels), function(i) {
    level <- as.character(levels[[i]])
    mean_value <- ms_safe_numeric(means[[level]])
    se_value <- ms_safe_numeric(se[[level]])
    raw_sd <- ms_safe_numeric(stats::sd(mat[, level], na.rm = TRUE))
    row <- list(
      level = level,
      label = level,
      n = nrow(mat),
      mean = mean_value,
      df = df_ci
    )
    if (!is.na(raw_sd)) row$sd <- raw_sd
    if (!is.na(se_value)) {
      row$se <- se_value
      row$ci_lower <- ms_safe_numeric(mean_value - critical * se_value)
      row$ci_upper <- ms_safe_numeric(mean_value + critical * se_value)
    }
    row
  })
  groups <- Filter(function(row) {
    is.list(row) &&
      !is.null(row$level) &&
      !is.null(row$mean) &&
      !is.na(ms_safe_numeric(row$mean))
  }, groups)
  if (length(groups) < 2L) return(NULL)

  focal_test <- ms_aovlist_focal_test(tidy_df, within_term)
  outcome <- ms_model_clean_term(names(mf)[[response_idx]])
  out <- list(
    mean_kind = "within_subject",
    source = "repeated_measures_aov",
    factor = list(
      variable = within_var,
      term = within_term,
      label = ms_model_clean_term(within_term),
      levels = lapply(levels, function(level) {
        list(value = as.character(level), label = as.character(level))
      })
    ),
    groups = groups,
    subject = list(
      variable = subject_var,
      label = ms_model_clean_term(subject_var),
      n = nrow(mat)
    ),
    outcome = outcome,
    y_label = if (nzchar(outcome)) paste("Mean", outcome) else "Mean outcome",
    ci_level = conf.level,
    ci_method = "within_subject_morey"
  )
  if (!is.null(focal_test)) {
    out$focal_test <- focal_test
    stat <- focal_test$statistic %||% list()
    df <- stat$df %||% numeric(0)
    if (length(df) >= 2L) out$residual_df <- ms_safe_numeric(df[[2L]])
    if (!is.null(focal_test$residual_mean_sq)) {
      out$residual_mean_sq <- focal_test$residual_mean_sq
    }
  }
  out
}

ms_repeated_measures_interaction_figure_data <- function(x, tidy_df = NULL,
                                                         conf.level = 0.95,
                                                         max_levels = 8L) {
  if (!inherits(x, "aovlist")) return(NULL)

  tf <- tryCatch(attr(x, "terms"), error = function(e) NULL)
  mf <- tryCatch(stats::model.frame(x), error = function(e) NULL)
  if (is.null(tf) || is.null(mf) || nrow(mf) < 2L) return(NULL)

  response_idx <- attr(tf, "response") %||% 0L
  if (response_idx < 1L || response_idx > ncol(mf)) return(NULL)
  y <- mf[[response_idx]]
  if (!is.numeric(y)) return(NULL)

  terms <- attr(tf, "term.labels") %||% character(0)
  error_terms <- terms[grepl("^Error\\(", terms)]
  within_terms <- setdiff(terms, error_terms)
  within_terms <- within_terms[nzchar(within_terms)]
  within_variables <- ms_aovlist_within_variables(within_terms, mf)
  if (length(error_terms) != 1L || length(within_variables) != 2L) return(NULL)

  interaction_terms <- within_terms[vapply(within_terms, function(term) {
    parts <- strsplit(term, ":", fixed = TRUE)[[1L]]
    identical(sort(parts), sort(within_variables))
  }, logical(1))]
  if (!length(interaction_terms)) return(NULL)
  interaction_term <- interaction_terms[[1L]]

  x_var <- within_variables[[1L]]
  moderator_var <- within_variables[[2L]]
  if (!all(c(x_var, moderator_var) %in% names(mf))) return(NULL)

  subject_var <- ms_aovlist_subject_variable_multi(
    error_terms[[1L]],
    within_terms,
    within_variables,
    mf
  )
  if (is.null(subject_var) || !subject_var %in% names(mf)) return(NULL)

  x_values <- mf[[x_var]]
  moderator_values <- mf[[moderator_var]]
  subject_values <- mf[[subject_var]]
  if (!(is.factor(x_values) || is.character(x_values) || is.logical(x_values))) return(NULL)
  if (!(is.factor(moderator_values) || is.character(moderator_values) || is.logical(moderator_values))) return(NULL)

  ok <- is.finite(y) & !is.na(subject_values) & !is.na(x_values) & !is.na(moderator_values)
  y <- y[ok]
  subject_values <- subject_values[ok]
  x_values <- x_values[ok]
  moderator_values <- moderator_values[ok]
  if (length(y) < 2L) return(NULL)

  x_levels <- ms_anova_observed_levels(x_values)
  moderator_levels <- ms_anova_observed_levels(moderator_values)
  if (length(x_levels) < 2L || length(moderator_levels) < 2L) return(NULL)
  if (length(x_levels) > max_levels || length(moderator_levels) > max_levels) return(NULL)

  cell_frame <- data.frame(
    subject = as.character(subject_values),
    x_level = factor(as.character(x_values), levels = x_levels),
    moderator_level = factor(as.character(moderator_values), levels = moderator_levels),
    y = y,
    stringsAsFactors = FALSE
  )
  subject_cells <- stats::aggregate(
    y ~ subject + x_level + moderator_level,
    data = cell_frame,
    FUN = function(value) mean(value, na.rm = TRUE),
    drop = FALSE
  )
  if (!nrow(subject_cells)) return(NULL)

  subjects <- unique(as.character(subject_cells$subject))
  cell_ids <- as.vector(outer(x_levels, moderator_levels, paste, sep = "|||"))
  mat <- matrix(NA_real_, nrow = length(subjects), ncol = length(cell_ids),
                dimnames = list(subjects, cell_ids))
  for (i in seq_len(nrow(subject_cells))) {
    subject <- as.character(subject_cells$subject[[i]])
    x_level <- as.character(subject_cells$x_level[[i]])
    moderator_level <- as.character(subject_cells$moderator_level[[i]])
    cell_id <- paste(x_level, moderator_level, sep = "|||")
    if (!subject %in% rownames(mat) || !cell_id %in% colnames(mat)) next
    mat[subject, cell_id] <- ms_safe_numeric(subject_cells$y[[i]])
  }
  mat <- mat[stats::complete.cases(mat), , drop = FALSE]
  if (nrow(mat) < 2L || ncol(mat) < 2L) return(NULL)

  means <- colMeans(mat, na.rm = TRUE)
  subject_means <- rowMeans(mat, na.rm = TRUE)
  grand_mean <- mean(mat, na.rm = TRUE)
  normalized <- sweep(mat, 1L, subject_means, FUN = "-") + grand_mean
  morey <- sqrt(ncol(mat) / (ncol(mat) - 1L))
  df_ci <- nrow(mat) - 1L
  critical <- stats::qt(1 - ((1 - conf.level) / 2), df = df_ci)
  se <- apply(normalized, 2L, stats::sd, na.rm = TRUE) / sqrt(nrow(mat)) * morey

  grid <- list()
  for (moderator_level in moderator_levels) {
    for (x_index in seq_along(x_levels)) {
      x_level <- x_levels[[x_index]]
      cell_id <- paste(x_level, moderator_level, sep = "|||")
      mean_value <- ms_safe_numeric(means[[cell_id]])
      se_value <- ms_safe_numeric(se[[cell_id]])
      raw_sd <- ms_safe_numeric(stats::sd(mat[, cell_id], na.rm = TRUE))
      if (is.na(mean_value)) next
      row <- list(
        x = x_index,
        x_value = x_level,
        x_label = ms_model_clean_term(x_level),
        moderator_value = moderator_level,
        moderator_label = ms_model_clean_term(moderator_level),
        estimate = mean_value,
        n = nrow(mat),
        df = df_ci
      )
      if (!is.na(raw_sd)) row$sd <- raw_sd
      if (!is.na(se_value)) {
        row$se <- se_value
        row$ci_lower <- ms_safe_numeric(mean_value - critical * se_value)
        row$ci_upper <- ms_safe_numeric(mean_value + critical * se_value)
      }
      grid[[length(grid) + 1L]] <- row
    }
  }
  if (length(grid) < length(x_levels) * length(moderator_levels)) return(NULL)

  focal_test <- ms_aovlist_focal_test(tidy_df, interaction_term)

  out <- list(
    interaction_term = interaction_term,
    interaction_kind = "categorical_by_categorical",
    source = "repeated_measures_aov",
    mean_kind = "observed_cell_means",
    variables = c(x_var, moderator_var),
    x = list(
      variable = x_var,
      term = x_var,
      label = ms_model_clean_term(x_var),
      type = "categorical",
      range = I(c(1, length(x_levels))),
      levels = lapply(x_levels, function(level) {
        list(value = level, label = ms_model_clean_term(level))
      })
    ),
    moderator = list(
      variable = moderator_var,
      term = moderator_var,
      label = ms_model_clean_term(moderator_var),
      type = "categorical",
      levels = lapply(moderator_levels, function(level) {
        list(value = level, label = ms_model_clean_term(level))
      })
    ),
    grid = grid,
    subject = list(
      variable = subject_var,
      label = ms_model_clean_term(subject_var),
      n = nrow(mat)
    ),
    outcome = ms_model_clean_term(names(mf)[[response_idx]]),
    y_label = paste("Mean", ms_model_clean_term(names(mf)[[response_idx]])),
    scale = "response",
    ci_level = conf.level,
    ci_method = "within_subject_morey",
    bounded_response = FALSE
  )
  if (!is.null(focal_test)) out$focal_test <- focal_test
  out
}

ms_aovlist_subject_variable <- function(error_term, within_var, mf) {
  inner <- sub("^Error\\((.*)\\)$", "\\1", error_term)
  vars <- tryCatch(all.vars(stats::as.formula(paste0("~", inner))),
                   error = function(e) character(0))
  vars <- vars[vars %in% names(mf)]
  vars <- vars[vars != within_var]
  if (!length(vars)) return(NULL)
  vars[[1L]]
}

ms_aovlist_focal_test <- function(tidy_df, term) {
  if (is.null(tidy_df) || !nrow(tidy_df) || !"term" %in% names(tidy_df)) return(NULL)
  term_values <- as.character(tidy_df$term)
  hit <- which(term_values == term)
  if (!length(hit)) return(NULL)
  i <- hit[[1L]]
  f_value <- if ("statistic" %in% names(tidy_df)) {
    ms_safe_numeric(tidy_df$statistic[[i]])
  } else NA_real_
  df1 <- if ("df" %in% names(tidy_df)) ms_safe_numeric(tidy_df$df[[i]]) else NA_real_
  p_value <- if ("p.value" %in% names(tidy_df)) {
    ms_safe_numeric(tidy_df$p.value[[i]])
  } else NA_real_
  if (is.na(f_value) || is.na(df1)) return(NULL)

  stratum <- if ("stratum" %in% names(tidy_df)) as.character(tidy_df$stratum[[i]]) else ""
  resid_mask <- term_values %in% c("Residuals", "Residual")
  if ("stratum" %in% names(tidy_df)) {
    resid_mask <- resid_mask & as.character(tidy_df$stratum) == stratum
  }
  resid_idx <- which(resid_mask)
  df2 <- if (length(resid_idx) && "df" %in% names(tidy_df)) {
    ms_safe_numeric(tidy_df$df[[resid_idx[[1L]]]])
  } else NA_real_
  ss <- if ("sumsq" %in% names(tidy_df)) ms_safe_numeric(tidy_df$sumsq[[i]]) else NA_real_
  resid_ss <- if (length(resid_idx) && "sumsq" %in% names(tidy_df)) {
    ms_safe_numeric(tidy_df$sumsq[[resid_idx[[1L]]]])
  } else NA_real_
  resid_ms <- if (length(resid_idx) && "meansq" %in% names(tidy_df)) {
    ms_safe_numeric(tidy_df$meansq[[resid_idx[[1L]]]])
  } else NA_real_

  out <- list(
    statistic = list(name = "F", value = f_value),
    p_value = p_value
  )
  if (!is.na(df2)) out$statistic$df <- I(c(df1, df2))
  if (!is.na(ss)) out$sum_sq <- ss
  if (!is.na(resid_ss)) out$residual_sum_sq <- resid_ss
  if (!is.na(resid_ms)) out$residual_mean_sq <- resid_ms
  if (!is.na(ss) && !is.na(resid_ss) && (ss + resid_ss) > 0) {
    out$effect <- list(
      name = "eta_sq_partial",
      value = ms_safe_numeric(ss / (ss + resid_ss))
    )
  }
  out
}

ms_ancova_adjusted_means_figure_data <- function(x, spec = NULL,
                                                 conf.level = 0.95,
                                                 max_levels = 12L) {
  if (!inherits(x, "lm") || inherits(x, c("glm", "mlm", "maov", "manova", "aovlist"))) {
    return(NULL)
  }
  if (!requireNamespace("emmeans", quietly = TRUE)) return(NULL)

  spec <- spec %||% ms_ancova_model_spec(x, max_levels = max_levels)
  if (is.null(spec)) return(NULL)

  at_values <- setNames(
    lapply(spec$covariates, function(covariate) covariate$value),
    vapply(spec$covariates, function(covariate) covariate$variable, character(1))
  )
  emm_formula <- tryCatch(stats::as.formula(paste("~", spec$factor$term)),
                          error = function(e) NULL)
  if (is.null(emm_formula)) return(NULL)
  emm <- tryCatch(
    emmeans::emmeans(x, specs = emm_formula, at = at_values),
    error = function(e) NULL
  )
  if (is.null(emm)) return(NULL)
  emm_summary <- tryCatch(
    as.data.frame(summary(emm, infer = c(TRUE, FALSE), level = conf.level)),
    error = function(e) NULL
  )
  if (is.null(emm_summary) || nrow(emm_summary) < 2L) return(NULL)

  measure_cols <- c("emmean", "SE", "df", "lower.CL", "upper.CL",
                    "asymp.LCL", "asymp.UCL")
  group_cols <- setdiff(names(emm_summary), measure_cols)
  group_col <- if (length(group_cols)) group_cols[[1L]] else NA_character_
  if (is.na(group_col) || !group_col %in% names(emm_summary)) return(NULL)

  lower_col <- if ("lower.CL" %in% names(emm_summary)) "lower.CL" else "asymp.LCL"
  upper_col <- if ("upper.CL" %in% names(emm_summary)) "upper.CL" else "asymp.UCL"
  if (!all(c("emmean", "SE", lower_col, upper_col) %in% names(emm_summary))) return(NULL)

  counts <- spec$factor$counts
  groups <- lapply(seq_len(nrow(emm_summary)), function(i) {
    level <- as.character(emm_summary[[group_col]][[i]])
    n_value <- if (level %in% names(counts)) {
      as.integer(counts[[level]])
    } else {
      NA_integer_
    }
    row <- list(
      level = level,
      label = level,
      n = n_value,
      mean = ms_safe_numeric(emm_summary[["emmean"]][[i]]),
      se = ms_safe_numeric(emm_summary[["SE"]][[i]]),
      ci_lower = ms_safe_numeric(emm_summary[[lower_col]][[i]]),
      ci_upper = ms_safe_numeric(emm_summary[[upper_col]][[i]])
    )
    df_value <- if ("df" %in% names(emm_summary)) {
      ms_safe_numeric(emm_summary[["df"]][[i]])
    } else NA_real_
    if (!is.na(df_value)) row$df <- df_value
    row
  })
  groups <- Filter(function(row) {
    is.list(row) &&
      !is.null(row$level) &&
      !is.null(row$mean) &&
      !is.na(ms_safe_numeric(row$mean))
  }, groups)
  if (length(groups) < 2L) return(NULL)

  out <- list(
    mean_kind = "estimated_marginal",
    source = "ancova_emmeans",
    factor = spec$factor[c("variable", "term", "label", "levels")],
    groups = groups,
    covariates = lapply(spec$covariates, function(covariate) {
      list(
        variable = covariate$variable,
        term = covariate$term,
        label = covariate$label,
        value = covariate$value,
        value_label = ms_interaction_format_number(covariate$value),
        rule = "sample_mean"
      )
    }),
    adjustment = list(
      rule = "sample_mean",
      label = "sample means"
    ),
    outcome = spec$outcome,
    y_label = if (nzchar(spec$outcome)) {
      paste("Estimated marginal mean", spec$outcome)
    } else {
      "Estimated marginal mean"
    },
    ci_level = conf.level,
    ci_method = "emmeans"
  )
  focal_test <- ms_ancova_focal_test(x, spec$factor$term)
  if (!is.null(focal_test)) out$focal_test <- focal_test
  out
}

ms_factorial_main_effect_means_figure_data <- function(x, focal = NULL,
                                                       controls = NULL,
                                                       conf.level = 0.95,
                                                       max_levels = 12L) {
  if (!inherits(x, "lm") || inherits(x, c("glm", "mlm", "maov", "manova", "aovlist"))) {
    return(NULL)
  }
  if (!requireNamespace("emmeans", quietly = TRUE)) return(NULL)

  tf <- tryCatch(stats::terms(x), error = function(e) NULL)
  mf <- tryCatch(stats::model.frame(x), error = function(e) NULL)
  if (is.null(tf) || is.null(mf) || nrow(mf) < 2L) return(NULL)

  term_labels <- attr(tf, "term.labels") %||% character(0)
  orders <- attr(tf, "order") %||% integer(0)
  response_idx <- attr(tf, "response") %||% 0L
  if (!length(term_labels) || length(term_labels) != length(orders)) return(NULL)
  if (any(orders != 1L)) return(NULL)
  if (response_idx < 1L || response_idx > ncol(mf)) return(NULL)

  term_info <- lapply(term_labels, function(term) {
    values <- ms_ancova_term_values(term, mf)
    type <- if (is.null(values)) "other" else ms_model_term_type(term, mf)
    variable <- ms_ancova_single_variable(term)
    if (identical(type, "numeric") && !identical(term, variable)) type <- "other"
    if (identical(type, "factor") && !ms_ancova_supported_factor_term(term, variable)) {
      type <- "other"
    }
    list(
      term = term,
      variable = variable %||% term,
      label = ms_model_clean_term(term),
      type = type,
      values = values
    )
  })

  is_factor <- vapply(term_info, function(info) identical(info$type, "factor"), logical(1))
  is_numeric <- vapply(term_info, function(info) identical(info$type, "numeric"), logical(1))
  is_other <- !(is_factor | is_numeric)
  if (sum(is_factor) < 2L || any(is_other)) return(NULL)

  factor_idx <- which(is_factor)
  focal_keys <- ms_role_keys(focal)
  control_keys <- ms_role_keys(controls)
  if (length(control_keys) > 0L) {
    factor_idx <- factor_idx[!vapply(term_info[factor_idx], function(info) {
      ms_factorial_term_matches_role(info, control_keys)
    }, logical(1))]
  }
  if (!length(factor_idx)) return(NULL)

  if (length(focal_keys) > 0L) {
    focal_hits <- factor_idx[vapply(term_info[factor_idx], function(info) {
      ms_factorial_term_matches_role(info, focal_keys)
    }, logical(1))]
    if (!length(focal_hits)) return(NULL)
    selected_idx <- focal_hits[[1L]]
  } else {
    selected_idx <- factor_idx[[1L]]
  }

  factor_info <- term_info[[selected_idx]]
  factor_values <- factor_info$values
  levels <- ms_anova_observed_levels(factor_values)
  if (length(levels) < 2L || length(levels) > max_levels) return(NULL)
  counts <- table(as.character(factor_values[!is.na(factor_values)]))
  level_rows <- lapply(levels, function(level) {
    list(value = as.character(level), label = as.character(level))
  })

  covariates <- lapply(term_info[is_numeric], function(info) {
    values <- suppressWarnings(as.numeric(info$values))
    values <- values[is.finite(values)]
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
  emm_args <- list(object = x, specs = emm_formula)
  if (length(at_values) > 0L) emm_args$at <- at_values
  emm <- tryCatch(
    do.call(emmeans::emmeans, emm_args),
    error = function(e) NULL
  )
  if (is.null(emm)) return(NULL)
  emm_summary <- tryCatch(
    as.data.frame(summary(emm, infer = c(TRUE, FALSE), level = conf.level)),
    error = function(e) NULL
  )
  if (is.null(emm_summary) || nrow(emm_summary) < 2L) return(NULL)

  measure_cols <- c("emmean", "SE", "df", "lower.CL", "upper.CL",
                    "asymp.LCL", "asymp.UCL")
  group_cols <- setdiff(names(emm_summary), measure_cols)
  group_col <- if (length(group_cols)) group_cols[[1L]] else NA_character_
  if (is.na(group_col) || !group_col %in% names(emm_summary)) return(NULL)

  lower_col <- if ("lower.CL" %in% names(emm_summary)) "lower.CL" else "asymp.LCL"
  upper_col <- if ("upper.CL" %in% names(emm_summary)) "upper.CL" else "asymp.UCL"
  if (!all(c("emmean", "SE", lower_col, upper_col) %in% names(emm_summary))) return(NULL)

  groups <- lapply(seq_len(nrow(emm_summary)), function(i) {
    level <- as.character(emm_summary[[group_col]][[i]])
    n_value <- if (level %in% names(counts)) {
      as.integer(counts[[level]])
    } else {
      NA_integer_
    }
    row <- list(
      level = level,
      label = level,
      n = n_value,
      mean = ms_safe_numeric(emm_summary[["emmean"]][[i]]),
      se = ms_safe_numeric(emm_summary[["SE"]][[i]]),
      ci_lower = ms_safe_numeric(emm_summary[[lower_col]][[i]]),
      ci_upper = ms_safe_numeric(emm_summary[[upper_col]][[i]])
    )
    df_value <- if ("df" %in% names(emm_summary)) {
      ms_safe_numeric(emm_summary[["df"]][[i]])
    } else NA_real_
    if (!is.na(df_value)) row$df <- df_value
    row
  })
  groups <- Filter(function(row) {
    is.list(row) &&
      !is.null(row$level) &&
      !is.null(row$mean) &&
      !is.na(ms_safe_numeric(row$mean))
  }, groups)
  if (length(groups) < 2L) return(NULL)

  other_factors <- term_info[is_factor]
  other_factors <- other_factors[vapply(other_factors, function(info) {
    !identical(info$term, factor_info$term)
  }, logical(1))]

  out <- list(
    mean_kind = "estimated_marginal",
    source = "factorial_emmeans",
    factor = list(
      variable = factor_info$variable,
      term = factor_info$term,
      label = factor_info$label,
      levels = level_rows
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
    outcome = ms_model_clean_term(names(mf)[[response_idx]]),
    ci_level = conf.level,
    ci_method = "emmeans"
  )
  out$y_label <- if (nzchar(out$outcome)) {
    paste("Estimated marginal mean", out$outcome)
  } else {
    "Estimated marginal mean"
  }
  if (length(covariates) > 0L) {
    out$covariates <- lapply(covariates, function(covariate) {
      list(
        variable = covariate$variable,
        term = covariate$term,
        label = covariate$label,
        value = covariate$value,
        value_label = ms_interaction_format_number(covariate$value),
        rule = "sample_mean"
      )
    })
    out$adjustment <- list(
      rule = "sample_mean",
      label = "sample means"
    )
  }
  out
}

ms_factorial_term_matches_role <- function(info, keys) {
  if (!length(keys)) return(FALSE)
  term_keys <- unique(ms_term_key(c(
    info$term %||% "",
    info$variable %||% "",
    info$label %||% ""
  )))
  length(intersect(term_keys, keys)) > 0L
}

ms_ancova_model_spec <- function(x, max_levels = 12L) {
  if (!inherits(x, "lm") || inherits(x, c("glm", "mlm", "maov", "manova", "aovlist"))) {
    return(NULL)
  }
  tf <- tryCatch(stats::terms(x), error = function(e) NULL)
  mf <- tryCatch(stats::model.frame(x), error = function(e) NULL)
  if (is.null(tf) || is.null(mf) || nrow(mf) < 2L) return(NULL)

  term_labels <- attr(tf, "term.labels") %||% character(0)
  orders <- attr(tf, "order") %||% integer(0)
  response_idx <- attr(tf, "response") %||% 0L
  if (!length(term_labels) || length(term_labels) != length(orders)) return(NULL)
  if (any(orders != 1L)) return(NULL)
  if (response_idx < 1L || response_idx > ncol(mf)) return(NULL)

  term_info <- lapply(term_labels, function(term) {
    values <- ms_ancova_term_values(term, mf)
    type <- if (is.null(values)) "other" else ms_model_term_type(term, mf)
    variable <- ms_ancova_single_variable(term)
    if (identical(type, "numeric") && !identical(term, variable)) type <- "other"
    if (identical(type, "factor") && !ms_ancova_supported_factor_term(term, variable)) {
      type <- "other"
    }
    list(
      term = term,
      variable = variable %||% term,
      label = ms_model_clean_term(term),
      type = type,
      values = values
    )
  })
  is_factor <- vapply(term_info, function(info) identical(info$type, "factor"), logical(1))
  is_numeric <- vapply(term_info, function(info) identical(info$type, "numeric"), logical(1))
  is_other <- !(is_factor | is_numeric)
  if (sum(is_factor) != 1L || sum(is_numeric) < 1L || any(is_other)) return(NULL)

  factor_info <- term_info[[which(is_factor)[[1L]]]]
  factor_values <- factor_info$values
  levels <- ms_anova_observed_levels(factor_values)
  if (length(levels) < 2L || length(levels) > max_levels) return(NULL)
  counts <- table(as.character(factor_values[!is.na(factor_values)]))
  level_rows <- lapply(levels, function(level) {
    list(value = as.character(level), label = as.character(level))
  })

  covariates <- lapply(term_info[is_numeric], function(info) {
    values <- suppressWarnings(as.numeric(info$values))
    values <- values[is.finite(values)]
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
  if (!length(covariates)) return(NULL)

  outcome <- ms_model_clean_term(names(mf)[[response_idx]])
  list(
    outcome = outcome,
    factor = list(
      variable = factor_info$variable,
      term = factor_info$term,
      label = factor_info$label,
      levels = level_rows,
      counts = counts
    ),
    covariates = covariates
  )
}

ms_ancova_term_values <- function(term, mf) {
  if (is.null(mf) || !nzchar(term)) return(NULL)
  if (term %in% names(mf)) return(mf[[term]])
  variable <- ms_ancova_single_variable(term)
  if (!is.null(variable) && variable %in% names(mf)) return(mf[[variable]])
  NULL
}

ms_ancova_single_variable <- function(term) {
  vars <- tryCatch(all.vars(stats::as.formula(paste0("~", term))),
                   error = function(e) character(0))
  if (length(vars) != 1L || !nzchar(vars[[1L]])) return(NULL)
  vars[[1L]]
}

ms_ancova_supported_factor_term <- function(term, variable) {
  term <- trimws(as.character(term %||% ""))
  variable <- trimws(as.character(variable %||% ""))
  if (!nzchar(term) || !nzchar(variable)) return(FALSE)
  if (identical(term, variable)) return(TRUE)

  parsed <- tryCatch(parse(text = term)[[1L]], error = function(e) NULL)
  if (!is.call(parsed) || length(parsed) != 2L) return(FALSE)
  fun <- as.character(parsed[[1L]])
  arg <- parsed[[2L]]
  fun %in% c("factor", "as.factor") && identical(as.character(arg), variable)
}

ms_ancova_focal_test <- function(x, factor_term) {
  if (!nzchar(factor_term)) return(NULL)
  tbl <- tryCatch(stats::drop1(x, test = "F"), error = function(e) NULL)
  if (is.null(tbl) || !nrow(tbl)) return(NULL)
  rn <- rownames(tbl)
  hit <- match(factor_term, rn, nomatch = 0L)
  if (hit <= 0L) {
    key <- ms_term_key(factor_term)
    keys <- vapply(rn, ms_term_key, character(1))
    hit <- match(key, keys, nomatch = 0L)
  }
  if (hit <= 0L) return(NULL)
  df1 <- if ("Df" %in% names(tbl)) ms_safe_numeric(tbl[["Df"]][[hit]]) else NA_real_
  f_value <- if ("F value" %in% names(tbl)) {
    ms_safe_numeric(tbl[["F value"]][[hit]])
  } else NA_real_
  p_value <- if ("Pr(>F)" %in% names(tbl)) {
    ms_safe_numeric(tbl[["Pr(>F)"]][[hit]])
  } else NA_real_
  sum_sq <- if ("Sum of Sq" %in% names(tbl)) {
    ms_safe_numeric(tbl[["Sum of Sq"]][[hit]])
  } else NA_real_
  full_idx <- match("<none>", rn, nomatch = 0L)
  residual_sum_sq <- if (full_idx > 0L && "RSS" %in% names(tbl)) {
    ms_safe_numeric(tbl[["RSS"]][[full_idx]])
  } else {
    rss <- tryCatch(sum(stats::residuals(x)^2, na.rm = TRUE), error = function(e) NA_real_)
    ms_safe_numeric(rss)
  }
  df2 <- ms_safe_numeric(stats::df.residual(x))
  if (is.na(df1) || is.na(f_value)) return(NULL)
  out <- list(
    statistic = list(name = "F", value = f_value, df = I(c(df1, df2))),
    p_value = p_value,
    df_effect = df1,
    df_error = df2,
    term = factor_term
  )
  if (!is.na(sum_sq)) out$sum_sq <- sum_sq
  if (!is.na(residual_sum_sq)) out$residual_sum_sq <- residual_sum_sq
  if (!is.na(sum_sq) && !is.na(residual_sum_sq) && (sum_sq + residual_sum_sq) > 0) {
    out$effect <- list(
      name = "eta_sq_partial",
      value = ms_safe_numeric(sum_sq / (sum_sq + residual_sum_sq))
    )
  }
  out
}
