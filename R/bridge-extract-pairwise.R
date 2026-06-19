# R bridge -- pairwise comparisons / post-hoc tests.

#' @rdname mellio_payload
#' @export
mellio_payload.TukeyHSD <- function(x, ..., .call = NULL) {
  rows <- list()
  term_names <- names(x) %||% rep("", length(x))
  for (term in term_names) {
    mat <- as.matrix(x[[term]])
    if (!nrow(mat)) next
    rn <- rownames(mat) %||% paste0("contrast_", seq_len(nrow(mat)))
    for (i in seq_len(nrow(mat))) {
      row <- ms_pairwise_add_levels(list(
        family = as.character(term),
        contrast = ms_pairwise_clean_contrast(rn[[i]]),
        estimate = ms_pairwise_cell(mat, i, "diff"),
        ci_lower = ms_pairwise_cell(mat, i, "lwr"),
        ci_upper = ms_pairwise_cell(mat, i, "upr"),
        p_value = ms_pairwise_cell(mat, i, "p adj")
      ), rn[[i]])
      rows[[length(rows) + 1L]] <- row
    }
  }

  conf_level <- ms_safe_numeric(attr(x, "conf.level") %||% NA_real_)
  ms_pairwise_payload(
    rows = rows,
    method = "Tukey HSD",
    adjustment_method = "tukey",
    call = ms_pairwise_call(.call, match.call()$x),
    raw_output = ms_capture_output(x),
    conf_level = conf_level,
    source = "stats::TukeyHSD"
  )
}

#' @rdname mellio_payload
#' @export
mellio_payload.pairwise.htest <- function(x, ..., .call = NULL, .env = parent.frame()) {
  pmat <- x$p.value
  if (is.null(pmat)) {
    stop("pairwise.htest object has no p-value matrix.", call. = FALSE)
  }
  pmat <- as.matrix(pmat)
  rn <- rownames(pmat) %||% paste0("group_", seq_len(nrow(pmat)))
  cn <- colnames(pmat) %||% paste0("group_", seq_len(ncol(pmat)))
  context <- ms_pairwise_htest_context(x, .call = .call, .env = .env)
  level_order <- context$group_levels %||% ms_pairwise_matrix_level_order(rn, cn)

  rows <- list()
  for (i in seq_len(nrow(pmat))) {
    for (j in seq_len(ncol(pmat))) {
      p <- ms_safe_numeric(pmat[i, j])
      if (is.na(p)) next
      row <- list(
        contrast = ms_pairwise_clean_contrast(paste(rn[[i]], cn[[j]], sep = " - ")),
        level_1 = as.character(rn[[i]]),
        level_2 = as.character(cn[[j]]),
        p_value = p
      )
      row <- ms_pairwise_canonicalize_row(row, level_order = level_order, flip_values = FALSE)
      row <- ms_pairwise_enrich_htest_row(row, context)
      rows[[length(rows) + 1L]] <- row
    }
  }

  ms_pairwise_payload(
    rows = rows,
    method = as.character(x$method %||% "Pairwise comparisons"),
    adjustment_method = as.character(x$p.adjust.method %||% ""),
    call = ms_pairwise_call(.call, match.call()$x),
    raw_output = ms_capture_output(x),
    grouping_note = context$grouping_note %||% NULL,
    level_order = level_order,
    source = "stats::pairwise.htest"
  )
}

#' @rdname mellio_payload
#' @export
mellio_payload.dunn_test <- function(x, ..., .call = NULL) {
  df <- as.data.frame(x, stringsAsFactors = FALSE)
  args <- attr(x, "args") %||% list()
  rows <- ms_dunn_rows_from_rstatix(df, args)
  adjustment <- as.character(args$p.adjust.method %||% ms_dunn_adjustment_from_columns(df))

  ms_pairwise_payload(
    rows = rows,
    method = "Dunn's Kruskal-Wallis post hoc test",
    adjustment_method = adjustment,
    call = ms_pairwise_call(.call, match.call()$x, fallback = "rstatix::dunn_test(...)"),
    raw_output = ms_capture_output(x),
    packages = ms_packages_basic(extras = "rstatix"),
    grouping_note = ms_rstatix_grouping_note(args),
    level_order = ms_rstatix_group_levels(args),
    source = "rstatix::dunn_test"
  )
}

#' @rdname mellio_payload
#' @export
mellio_payload.wilcox_test <- function(x, ..., .call = NULL) {
  df <- as.data.frame(x, stringsAsFactors = FALSE)
  args <- attr(x, "args") %||% list()
  rows <- ms_wilcox_rows_from_rstatix(df, args, .call = .call)
  paired <- vapply(rows, function(row) isTRUE(row$paired), logical(1))
  method <- if (length(paired) && all(paired)) {
    "Wilcoxon signed-rank test"
  } else if (length(paired) && any(paired)) {
    "Wilcoxon pairwise comparisons"
  } else {
    "Wilcoxon rank-sum test"
  }
  adjustment <- as.character(args$p.adjust.method %||% ms_rstatix_adjustment_from_columns(df))

  ms_pairwise_payload(
    rows = rows,
    method = method,
    adjustment_method = adjustment,
    call = ms_pairwise_call(.call, match.call()$x, fallback = "rstatix::wilcox_test(...)"),
    raw_output = ms_capture_output(x),
    packages = ms_packages_basic(extras = "rstatix"),
    grouping_note = ms_rstatix_grouping_note(args),
    level_order = ms_rstatix_group_levels(args),
    source = "rstatix::wilcox_test"
  )
}

#' @rdname mellio_payload
#' @export
mellio_payload.t_test <- function(x, ..., .call = NULL) {
  df <- as.data.frame(x, stringsAsFactors = FALSE)
  args <- attr(x, "args") %||% list()
  rows <- ms_ttest_rows_from_rstatix(df, args, .call = .call)
  paired <- vapply(rows, function(row) isTRUE(row$paired), logical(1))
  method <- ms_rstatix_ttest_method(args, paired)
  adjustment <- as.character(args$p.adjust.method %||% ms_rstatix_adjustment_from_columns(df))

  ms_pairwise_payload(
    rows = rows,
    method = method,
    adjustment_method = adjustment,
    call = ms_pairwise_call(.call, match.call()$x, fallback = "rstatix::t_test(...)"),
    raw_output = ms_capture_output(x),
    conf_level = ms_safe_numeric(args$conf.level %||% NA_real_),
    packages = ms_packages_basic(extras = "rstatix"),
    grouping_note = ms_rstatix_grouping_note(args),
    level_order = ms_rstatix_group_levels(args),
    source = "rstatix::t_test"
  )
}

#' @rdname mellio_payload
#' @export
mellio_payload.dunnTest <- function(x, ..., .call = NULL) {
  res <- x$res %||% x$dtres %||% NULL
  if (is.null(res)) {
    stop("FSA::dunnTest object does not contain a comparison table.", call. = FALSE)
  }
  df <- as.data.frame(res, stringsAsFactors = FALSE)
  rows <- ms_dunn_rows_from_comparison_df(df)
  method <- as.character(x$method %||% "Dunn's Kruskal-Wallis post hoc test")
  adjustment <- ms_dunn_adjustment_from_text(method)

  ms_pairwise_payload(
    rows = rows,
    method = "Dunn's Kruskal-Wallis post hoc test",
    adjustment_method = adjustment,
    call = ms_pairwise_call(.call, match.call()$x, fallback = "FSA::dunnTest(...)"),
    raw_output = ms_capture_output(x),
    packages = ms_packages_basic(extras = "FSA"),
    source = "FSA::dunnTest"
  )
}

#' @rdname mellio_payload
#' @export
mellio_payload.dunn.test <- function(x, ..., .call = NULL) {
  rows <- ms_dunn_rows_from_dunn_test_object(x)
  adjustment <- ms_dunn_adjustment_from_text(as.character(x$method %||% ""))

  ms_pairwise_payload(
    rows = rows,
    method = "Dunn's Kruskal-Wallis post hoc test",
    adjustment_method = adjustment,
    call = ms_pairwise_call(.call, match.call()$x, fallback = "dunn.test::dunn.test(...)"),
    raw_output = ms_capture_output(x),
    packages = ms_packages_basic(extras = "dunn.test"),
    source = "dunn.test::dunn.test"
  )
}

#' @rdname mellio_payload
#' @export
mellio_payload.emmGrid <- function(x, ..., .call = NULL) {
  rlang::check_installed("emmeans", reason = "to extract pairwise emmeans contrasts")

  info <- ms_emm_pairwise_summary(x)
  s <- info$summary
  contrast_meta <- ms_emm_pairwise_contrast_meta(x, nrow(s))
  est_col <- attr(s, "estName") %||% "estimate"
  low_col <- (attr(s, "clNames") %||% c("lower.CL", "upper.CL"))[[1]]
  high_col <- (attr(s, "clNames") %||% c("lower.CL", "upper.CL"))[[2]]
  stat_col <- ms_pairwise_first_col(names(s), c("t.ratio", "z.ratio", "statistic"))
  stat_label <- if (!is.na(stat_col)) sub("\\.ratio$", "", stat_col) else NULL
  by_vars <- attr(s, "by.vars") %||% character(0)

  rows <- lapply(seq_len(nrow(s)), function(i) {
    meta <- contrast_meta[[i]] %||% NULL
    estimate <- ms_pairwise_df_cell(s, i, est_col)
    ci_lower <- ms_pairwise_df_cell(s, i, low_col)
    ci_upper <- ms_pairwise_df_cell(s, i, high_col)
    statistic <- if (!is.na(stat_col)) ms_pairwise_df_cell(s, i, stat_col) else NA_real_
    if (!is.null(meta) && isTRUE(meta$flip)) {
      estimate <- -estimate
      interval <- ms_pairwise_flip_interval(ci_lower, ci_upper)
      ci_lower <- interval[[1L]]
      ci_upper <- interval[[2L]]
      statistic <- -statistic
    }

    row <- list(
      contrast = (meta$contrast %||% ms_pairwise_clean_contrast(as.character(s$contrast[[i]]))),
      estimate = estimate,
      std_error = ms_pairwise_df_cell(s, i, "SE"),
      df = ms_pairwise_df_cell(s, i, "df"),
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      p_value = ms_pairwise_df_cell(s, i, "p.value")
    )
    if (!is.null(meta) && !is.null(meta$level_1) && !is.null(meta$level_2)) {
      row$level_1 <- meta$level_1
      row$level_2 <- meta$level_2
    } else {
      row <- ms_pairwise_add_levels(row, as.character(s$contrast[[i]]))
    }
    if (!is.na(stat_col)) {
      row$statistic <- statistic
      row$statistic_label <- stat_label
    }
    if (length(by_vars)) {
      parts <- vapply(by_vars, function(v) {
        if (v %in% names(s)) paste0(v, "=", as.character(s[[v]][[i]])) else ""
      }, character(1))
      parts <- parts[nzchar(parts)]
      if (length(parts)) row$by <- paste(parts, collapse = ", ")
    }
    row
  })

  ms_pairwise_payload(
    rows = rows,
    method = info$method,
    adjustment_method = info$adjustment,
    call = ms_pairwise_call(.call, match.call()$x, fallback = "emmeans::pairs(...)"),
    raw_output = ms_capture_output(s),
    conf_level = ms_safe_numeric((info$misc$level %||% attr(s, "level")) %||% NA_real_),
    packages = ms_packages_basic(extras = "emmeans"),
    source = "emmeans::emmGrid"
  )
}

#' @rdname mellio_payload
#' @export
mellio_payload.emm_list <- function(x, ..., .call = NULL) {
  rlang::check_installed("emmeans", reason = "to extract pairwise emmeans contrasts")

  grid <- NULL
  if (!is.null(x$contrasts) && inherits(x$contrasts, "emmGrid")) {
    grid <- x$contrasts
  } else if (!is.null(x$emmeans) && inherits(x$emmeans, "emmGrid")) {
    grid <- x$emmeans
  }
  if (is.null(grid)) {
    stop("emm_list object does not contain emmeans contrasts.", call. = FALSE)
  }

  call_str <- ms_pairwise_call(
    .call,
    match.call()$x,
    fallback = "emmeans::emmeans(..., pairwise ~ ...)"
  )
  payload <- mellio_payload(grid, ..., .call = call_str)
  payload$raw_output <- ms_capture_output(x)
  payload$fields$source <- "emmeans::emm_list"
  if (!is.null(x$emmeans) && inherits(x$emmeans, "emmGrid")) {
    emm_table <- ms_emm_means_table(x$emmeans)
    if (!is.null(emm_table)) {
      payload$fields$emmeans_rows <- emm_table$rows
      payload$fields$emmeans_columns <- emm_table$columns
      payload$fields$emmeans_note <- emm_table$note
      payload$fields$emmeans_title <- emm_table$title
      if (!is.null(emm_table$conf_level)) payload$fields$emmeans_conf_level <- emm_table$conf_level
    }
  }
  payload
}

#' @rdname mellio_payload
#' @export
mellio_payload.glht <- function(x, ..., .call = NULL) {
  ms_glht_pairwise_payload(x, .call = .call, expr = match.call()$x)
}

#' @rdname mellio_payload
#' @export
mellio_payload.summary.glht <- function(x, ..., .call = NULL) {
  ms_glht_pairwise_payload(
    x,
    summary_obj = x,
    .call = .call,
    expr = match.call()$x,
    summary_input = TRUE
  )
}

ms_pairwise_payload <- function(rows, method, adjustment_method, call,
                                raw_output, conf_level = NA_real_,
                                packages = NULL, grouping_note = NULL,
                                level_order = NULL,
                                source = "R") {
  rows <- Filter(function(row) length(row) > 0L, rows)
  if (!length(rows)) {
    stop("No pairwise comparison rows could be extracted.", call. = FALSE)
  }
  rows <- ms_pairwise_drop_single_family(rows)
  rows <- lapply(rows, ms_pairwise_canonicalize_row, level_order = level_order)

  fields <- list(
    table_type = "pairwise_comparisons",
    method = method,
    adjustment_method = adjustment_method,
    adjustment_label = ms_pairwise_adjustment_label(adjustment_method),
    n_comparisons = length(rows),
    columns = ms_pairwise_columns(rows, conf_level = conf_level),
    rows = rows,
    source = source,
    note = ms_pairwise_note(rows, method, adjustment_method, conf_level)
  )
  if (!is.na(conf_level)) fields$conf_level <- conf_level
  if (!is.null(grouping_note) && nzchar(trimws(as.character(grouping_note)))) {
    fields$grouping_note <- trimws(as.character(grouping_note))
  }

  pairwise_forest <- ms_pairwise_forest_data(
    rows = rows,
    method = method,
    adjustment_method = adjustment_method,
    conf_level = conf_level,
    source = source
  )
  figure_data <- NULL
  available_figures <- NULL
  if (!is.null(pairwise_forest)) {
    figure_data <- list(pairwise_forest = pairwise_forest)
    available_figures <- list(list(
      type = "pairwise_forest",
      label = "Pairwise forest plot",
      default = TRUE
    ))
  }

  ms_build_envelope(
    type = "pairwise_comparisons",
    type_label = "Pairwise comparisons",
    call = trimws(gsub("\\s+", " ", call %||% NA_character_)),
    fields = fields,
    raw_output = raw_output,
    packages = packages,
    card_kind = "table",
    figure_data = figure_data,
    available_figures = available_figures
  )
}

ms_pairwise_columns <- function(rows, conf_level = NA_real_) {
  has_key <- function(key) {
    any(vapply(rows, function(row) !is.null(row[[key]]) && !is.na(row[[key]]), logical(1)))
  }

  # Derive the statistic column label from the first row that carries
  # a statistic_label (e.g., "t" or "z"). Every row in a pairwise table
  # shares the same stat type, so reading once is safe. Falls back to
  # "Statistic" only when the extractor didn't set a label.
  stat_label <- NULL
  if (has_key("statistic")) {
    for (row in rows) {
      if (!is.null(row$statistic_label) && nzchar(row$statistic_label)) {
        stat_label <- row$statistic_label
        break
      }
    }
    if (is.null(stat_label)) stat_label <- "Statistic"
  }

  cols <- list()
  if (has_key("family")) {
    cols <- c(cols, list(list(key = "family", label = "Family", format = "text")))
  }
  if (has_key("by")) {
    cols <- c(cols, list(list(key = "by", label = "By", format = "text")))
  }
  cols <- c(cols, list(list(key = "contrast", label = "Contrast", format = "text")))
  if (has_key("estimate")) {
    cols <- c(cols, list(list(key = "estimate", label = "Difference", format = "number")))
  }
  if (has_key("n_pairs")) {
    cols <- c(cols, list(list(key = "n_pairs", label = "n pairs", format = "integer")))
  }
  if (has_key("n_1")) {
    cols <- c(cols, list(list(key = "n_1", label = "n 1", format = "integer")))
  }
  if (has_key("mean_1")) {
    cols <- c(cols, list(list(key = "mean_1", label = "M 1", format = "number")))
  }
  if (has_key("median_1")) {
    cols <- c(cols, list(list(key = "median_1", label = "Mdn 1", format = "number")))
  }
  if (has_key("n_2")) {
    cols <- c(cols, list(list(key = "n_2", label = "n 2", format = "integer")))
  }
  if (has_key("mean_2")) {
    cols <- c(cols, list(list(key = "mean_2", label = "M 2", format = "number")))
  }
  if (has_key("median_2")) {
    cols <- c(cols, list(list(key = "median_2", label = "Mdn 2", format = "number")))
  }
  if (has_key("median_difference")) {
    cols <- c(cols, list(list(key = "median_difference", label = "Median difference", format = "number")))
  }
  if (has_key("rank_mean_difference")) {
    cols <- c(cols, list(list(key = "rank_mean_difference", label = "Mean rank difference", format = "number")))
  }
  if (has_key("std_error")) {
    cols <- c(cols, list(list(key = "std_error", label = "SE", format = "number")))
  }
  # df: use format "df" so the JS renderer's formatDfCell() picks integer
  # rendering for residual df (945) and decimal rendering for fractional
  # df (e.g., Satterthwaite 24.65 from lmer). Was "number" which forced
  # 2 decimals on integer df.
  if (has_key("df")) {
    cols <- c(cols, list(list(key = "df", label = "df", format = "df")))
  }
  # Statistic: column header reads the row's statistic_label ("t"/"z")
  # instead of the generic "Statistic" so the header italicizes
  # correctly and matches reporter conventions.
  if (!is.null(stat_label)) {
    cols <- c(cols, list(list(key = "statistic", label = stat_label, format = "statistic")))
  }
  if (has_key("ci_lower")) {
    ci_label <- ms_pairwise_ci_column_label(conf_level)
    cols <- c(cols, list(
      list(key = "ci_lower", label = paste(ci_label, "lower"), format = "number"),
      list(key = "ci_upper", label = paste(ci_label, "upper"), format = "number")
    ))
  }
  if (has_key("effect_size")) {
    cols <- c(cols, list(list(
      key = "effect_size",
      label = ms_pairwise_effect_column_label(rows),
      format = "bounded"
    )))
  }
  if (has_key("p_raw")) {
    cols <- c(cols, list(list(key = "p_raw", label = "p (raw)", format = "pvalue")))
  }
  cols <- c(cols, list(list(key = "p_value", label = "p (adjusted)", format = "pvalue")))
  cols
}

ms_pairwise_note <- function(rows, method, adjustment_method, conf_level = NA_real_) {
  bits <- character(0)
  method <- trimws(as.character(method %||% ""))
  adjustment <- ms_pairwise_adjustment_label(adjustment_method)
  has_estimate <- any(vapply(rows, function(row) {
    !is.null(row$estimate) && !is.na(row$estimate)
  }, logical(1)))
  has_ci <- any(vapply(rows, function(row) {
    !is.null(row$ci_lower) && !is.na(row$ci_lower) &&
      !is.null(row$ci_upper) && !is.na(row$ci_upper)
  }, logical(1)))
  has_median_difference <- any(vapply(rows, function(row) {
    !is.null(row$median_difference) && !is.na(row$median_difference)
  }, logical(1)))
  has_rank_mean_difference <- any(vapply(rows, function(row) {
    !is.null(row$rank_mean_difference) && !is.na(row$rank_mean_difference)
  }, logical(1)))
  effect_names <- unique(vapply(rows, function(row) {
    trimws(as.character(row$effect_size_name %||% ""))
  }, character(1)))
  effect_names <- effect_names[nzchar(effect_names)]
  has_paired <- any(vapply(rows, function(row) isTRUE(row$paired), logical(1)))

  if (identical(tolower(method), "estimated marginal means contrasts")) {
    bits <- c(bits, "Pairwise comparisons are based on estimated marginal means (EMMs) from the fitted model.")
  } else if (nzchar(method)) {
    bits <- c(bits, paste0("Pairwise comparisons use ", method, "."))
  }
  if (nzchar(adjustment) && !tolower(adjustment) %in% c("none", "no adjustment")) {
    if (identical(tolower(adjustment), "single-step") &&
        grepl("Tukey contrasts", method, ignore.case = TRUE)) {
      bits <- c(bits, "p values are adjusted using the single-step method for Tukey contrasts.")
    } else {
      bits <- c(bits, paste0("p values are adjusted using the ", adjustment, " method."))
    }
  } else {
    bits <- c(bits, "p values are unadjusted.")
  }
  if (has_estimate) {
    bits <- c(bits, "Mean differences are reported as the first contrast level minus the second.")
  }
  if (has_median_difference) {
    median_note <- "Median differences are reported as the first contrast level minus the second."
    if (has_paired) median_note <- paste(median_note, "For paired tests, this is the median of paired differences.")
    bits <- c(bits, median_note)
  }
  if (has_rank_mean_difference) {
    bits <- c(bits, "Mean rank differences are reported as the first contrast level minus the second.")
  }
  if (has_ci) {
    ci_text <- ms_pairwise_ci_text_label(conf_level)
    if (grepl("tukey", paste(method, adjustment), ignore.case = TRUE)) {
      ci_text <- paste("Family-wise", ci_text)
    }
    bits <- c(bits, paste0(ci_text, " are reported where available."))
  }
  if (grepl("dunn", method, ignore.case = TRUE)) {
    bits <- c(bits, "Dunn z statistics are oriented to match the displayed contrast.")
  }
  if ("cliffs_delta" %in% effect_names) {
    bits <- c(bits, "Cliff's delta is positive when the first contrast level tends to have larger values.")
  }
  if ("rank_biserial" %in% effect_names) {
    bits <- c(bits, "The rank-biserial effect size is positive when paired differences favor the first contrast level.")
  }
  if ("cohens_d" %in% effect_names) {
    d_methods <- unique(vapply(rows, function(row) {
      if (!identical(trimws(as.character(row$effect_size_name %||% "")), "cohens_d")) return("")
      trimws(as.character(row$effect_size_method %||% ""))
    }, character(1)))
    d_methods <- d_methods[nzchar(d_methods)]
    if ("all_groups_pooled_sd" %in% d_methods) {
      bits <- c(bits, "Cohen's d uses the all-groups pooled SD from the pooled pairwise t test.")
    }
    if (!length(d_methods) || "pairwise_pooled_sd" %in% d_methods) {
      bits <- c(bits, "Cohen's d uses the pooled SD for each independent contrast.")
    }
  }
  if ("cohens_dz" %in% effect_names) {
    bits <- c(bits, "Cohen's dz uses the SD of the paired differences.")
  }
  paste(bits, collapse = " ")
}

ms_pairwise_effect_column_label <- function(rows) {
  for (row in rows) {
    name <- trimws(as.character(row$effect_size_name %||% ""))
    if (nzchar(name)) return(ms_pairwise_effect_label(name))
    label <- trimws(as.character(row$effect_size_label %||% ""))
    if (nzchar(label)) return(label)
  }
  "Effect size"
}

ms_pairwise_effect_label <- function(name) {
  lower <- tolower(trimws(as.character(name %||% "")))
  if (identical(lower, "cliffs_delta")) return("Cliff's delta")
  if (identical(lower, "rank_biserial")) return("Rank-biserial r")
  if (identical(lower, "cohens_d")) return("Cohen's d")
  if (identical(lower, "cohens_dz")) return("Cohen's dz")
  if (identical(lower, "eta_sq_h")) return("eta^2_H")
  if (identical(lower, "epsilon_sq")) return("epsilon^2")
  if (nzchar(name)) gsub("_", " ", name)
  else "Effect size"
}

ms_pairwise_forest_data <- function(rows, method, adjustment_method,
                                    conf_level = NA_real_, source = "R") {
  plot_rows <- Filter(Negate(is.null), lapply(rows, ms_pairwise_forest_row))
  if (!length(plot_rows)) return(NULL)

  out <- list(
    rows = plot_rows,
    estimate_label = "Mean difference",
    ci_label = ms_pairwise_ci_column_label(conf_level),
    method = as.character(method %||% ""),
    adjustment_method = as.character(adjustment_method %||% ""),
    adjustment_label = ms_pairwise_adjustment_label(adjustment_method),
    source = as.character(source %||% "R")
  )
  if (!is.na(ms_safe_numeric(conf_level))) {
    out$ci_level <- ms_safe_numeric(conf_level)
  }
  out
}

ms_pairwise_forest_row <- function(row) {
  estimate <- ms_safe_numeric(row$estimate %||% NA_real_)
  ci_lower <- ms_safe_numeric(row$ci_lower %||% NA_real_)
  ci_upper <- ms_safe_numeric(row$ci_upper %||% NA_real_)
  if (any(is.na(c(estimate, ci_lower, ci_upper)))) return(NULL)

  out <- list(
    label = ms_pairwise_forest_label(row),
    contrast = as.character(row$contrast %||% ""),
    estimate = estimate,
    ci_lower = min(ci_lower, ci_upper),
    ci_upper = max(ci_lower, ci_upper)
  )
  p_value <- ms_safe_numeric(row$p_value %||% NA_real_)
  if (!is.na(p_value)) out$p_value <- p_value
  family <- trimws(as.character(row$family %||% ""))
  if (nzchar(family)) out$family <- family
  by <- trimws(as.character(row$by %||% ""))
  if (nzchar(by)) out$by <- by
  out
}

ms_pairwise_forest_label <- function(row) {
  contrast <- trimws(as.character(row$contrast %||% ""))
  parts <- trimws(as.character(c(row$family %||% "", row$by %||% "")))
  parts <- parts[nzchar(parts)]
  if (length(parts)) paste(paste(parts, collapse = ", "), contrast, sep = ": ") else contrast
}

ms_pairwise_adjustment_label <- function(adjustment_method) {
  adjustment <- trimws(as.character(adjustment_method %||% ""))
  if (!nzchar(adjustment)) return("")
  lower <- tolower(adjustment)
  labels <- c(
    tukey = "Tukey",
    holm = "Holm",
    bonferroni = "Bonferroni",
    hochberg = "Hochberg",
    hommel = "Hommel",
    fdr = "FDR",
    bh = "BH",
    by = "BY",
    sidak = "Sidak",
    none = "none",
    `single-step` = "single-step"
  )
  if (lower %in% names(labels)) return(unname(labels[[lower]]))
  adjustment
}

ms_pairwise_ci_column_label <- function(conf_level = NA_real_) {
  level <- ms_safe_numeric(conf_level)
  if (is.na(level)) return("CI")
  paste0(as.integer(round(level * 100)), "% CI")
}

ms_pairwise_ci_text_label <- function(conf_level = NA_real_) {
  level <- ms_safe_numeric(conf_level)
  if (is.na(level)) return("Confidence intervals")
  paste0(as.integer(round(level * 100)), "% confidence intervals")
}

ms_pairwise_clean_contrast <- function(x) {
  x <- trimws(as.character(x %||% ""))
  x <- gsub("\\s*-\\s*", " - ", x)
  x <- gsub("\\s+", " ", x)
  x
}

ms_pairwise_matrix_level_order <- function(rn, cn) {
  values <- unique(as.character(c(cn, rn)))
  values[nzchar(values)]
}

ms_pairwise_canonicalize_row <- function(row, level_order = NULL, flip_values = TRUE) {
  if (is.null(row$level_1) || is.null(row$level_2)) {
    row <- ms_pairwise_add_levels(row, row$contrast %||% "")
  }
  a <- trimws(as.character(row$level_1 %||% ""))
  b <- trimws(as.character(row$level_2 %||% ""))
  if (!nzchar(a) || !nzchar(b)) return(row)
  if (!ms_pairwise_should_flip_levels(a, b, level_order)) return(row)

  stat_label <- toupper(trimws(as.character(row$statistic_label %||% "")))
  if (isTRUE(flip_values) && stat_label %in% c("W", "V")) return(row)

  row$level_1 <- b
  row$level_2 <- a
  row$contrast <- ms_pairwise_clean_contrast(paste(b, a, sep = " - "))
  if (!isTRUE(flip_values)) return(row)

  row <- ms_pairwise_swap(row, "n_1", "n_2")
  row <- ms_pairwise_swap(row, "mean_1", "mean_2")
  row <- ms_pairwise_swap(row, "median_1", "median_2")
  row <- ms_pairwise_flip_numeric(row, "estimate")
  row <- ms_pairwise_flip_numeric(row, "median_difference")
  row <- ms_pairwise_flip_numeric(row, "rank_mean_difference")
  if (stat_label %in% c("T", "Z")) {
    row <- ms_pairwise_flip_numeric(row, "statistic")
  }
  effect_name <- tolower(trimws(as.character(row$effect_size_name %||% "")))
  if (effect_name %in% c("cohens_d", "cohens_dz", "cliffs_delta", "rank_biserial")) {
    row <- ms_pairwise_flip_numeric(row, "effect_size")
  }
  ci_lower <- ms_safe_numeric(row$ci_lower %||% NA_real_)
  ci_upper <- ms_safe_numeric(row$ci_upper %||% NA_real_)
  if (!is.na(ci_lower) && !is.na(ci_upper)) {
    row$ci_lower <- ms_safe_numeric(-ci_upper)
    row$ci_upper <- ms_safe_numeric(-ci_lower)
  }
  row
}

ms_pairwise_should_flip_levels <- function(a, b, level_order = NULL) {
  order <- trimws(as.character(level_order %||% character(0)))
  order <- order[nzchar(order)]
  if (length(order)) {
    ia <- match(a, order)
    ib <- match(b, order)
    if (!is.na(ia) && !is.na(ib)) return(ia > ib)
  }
  nums <- suppressWarnings(as.numeric(c(a, b)))
  if (all(!is.na(nums))) return(nums[[1L]] > nums[[2L]])
  tolower(a) > tolower(b)
}

ms_pairwise_swap <- function(row, a, b) {
  tmp <- row[[a]]
  row[[a]] <- row[[b]]
  row[[b]] <- tmp
  row
}

ms_pairwise_flip_numeric <- function(row, key) {
  value <- ms_safe_numeric(row[[key]] %||% NA_real_)
  if (!is.na(value)) row[[key]] <- ms_safe_numeric(-value)
  row
}

ms_pairwise_cell <- function(mat, i, col) {
  if (!col %in% colnames(mat)) return(NA_real_)
  ms_safe_numeric(mat[i, col])
}

ms_pairwise_df_cell <- function(df, i, col) {
  if (is.na(col) || !col %in% names(df)) return(NA_real_)
  ms_safe_numeric(df[[col]][[i]])
}

ms_pairwise_first_col <- function(names, candidates) {
  hit <- candidates[candidates %in% names]
  if (length(hit)) hit[[1L]] else NA_character_
}

ms_pairwise_call <- function(.call, expr, fallback = NA_character_) {
  call_text <- ms_pairwise_call_text(.call)
  if (nzchar(call_text) && !ms_pairwise_is_symbol_call(call_text)) return(call_text)

  expr_text <- ms_pairwise_expr_text(expr)
  if (nzchar(expr_text) && !ms_pairwise_is_symbol_call(expr_text)) return(expr_text)

  fallback <- ms_pairwise_call_text(fallback)
  if (nzchar(fallback)) return(fallback)
  if (nzchar(call_text)) return(call_text)
  if (nzchar(expr_text)) return(expr_text)
  NA_character_
}

ms_pairwise_call_text <- function(x) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) return("")
  trimws(gsub("\\s+", " ", paste(x, collapse = " ")))
}

ms_pairwise_expr_text <- function(expr) {
  if (is.null(expr) || identical(expr, as.name("x"))) return("")
  trimws(gsub("\\s+", " ", paste(deparse(expr, width.cutoff = 500L), collapse = " ")))
}

ms_pairwise_is_symbol_call <- function(x) {
  x <- ms_pairwise_call_text(x)
  if (!nzchar(x)) return(FALSE)
  parsed <- tryCatch(parse(text = x)[[1L]], error = function(e) NULL)
  is.symbol(parsed)
}

ms_pairwise_add_levels <- function(row, contrast) {
  parts <- strsplit(trimws(as.character(contrast %||% "")), "\\s*-\\s*")[[1]]
  if (length(parts) == 2L && nzchar(parts[[1]]) && nzchar(parts[[2]])) {
    row$level_1 <- parts[[1]]
    row$level_2 <- parts[[2]]
  }
  row
}

ms_pairwise_drop_single_family <- function(rows) {
  families <- unique(vapply(rows, function(row) {
    value <- row$family %||% ""
    trimws(as.character(value))
  }, character(1)))
  families <- families[nzchar(families)]
  if (length(families) != 1L) return(rows)
  lapply(rows, function(row) {
    row$family <- NULL
    row
  })
}

ms_pairwise_htest_context <- function(x, .call = NULL, .env = parent.frame()) {
  method <- tolower(as.character(x$method %||% ""))
  is_t_test <- grepl("t tests?", method)
  is_wilcoxon <- grepl("wilcoxon", method)
  if (!is_t_test && !is_wilcoxon) return(NULL)

  args <- ms_pairwise_htest_call_args(.call)
  x_expr <- ms_pairwise_arg(args, "x", 1L)
  g_expr <- ms_pairwise_arg(args, "g", 2L)
  x_values <- ms_pairwise_eval_or_null(x_expr, .env)
  g_values <- ms_pairwise_eval_or_null(g_expr, .env)

  if (is.null(x_values) || is.null(g_values)) {
    data_exprs <- ms_pairwise_data_name_exprs(x$data.name %||% "")
    if (is.null(x_values)) x_values <- ms_pairwise_eval_or_null(data_exprs$x, .env)
    if (is.null(g_values)) g_values <- ms_pairwise_eval_or_null(data_exprs$g, .env)
  }

  if (is.null(x_values) || is.null(g_values)) return(NULL)
  if (length(x_values) != length(g_values)) return(NULL)
  grouping_note <- ms_pairwise_numeric_grouping_note(g_expr, g_values)

  alternative <- ms_pairwise_eval_character(
    ms_pairwise_arg(args, "alternative", NA_integer_), .env, "two.sided"
  )

  if (is_wilcoxon) {
    paired <- ms_pairwise_eval_logical(
      ms_pairwise_arg(args, "paired", NA_integer_),
      .env,
      grepl("paired|signed rank", method)
    )
    exact <- ms_pairwise_eval_nullable_logical(ms_pairwise_arg(args, "exact", NA_integer_), .env)
    correct <- ms_pairwise_eval_logical(ms_pairwise_arg(args, "correct", NA_integer_), .env, TRUE)
    return(list(
      kind = "wilcoxon",
      x = ms_safe_numeric(x_values),
      g = factor(g_values),
      group_levels = levels(factor(g_values)),
      paired = isTRUE(paired),
      alternative = alternative,
      exact = exact,
      correct = isTRUE(correct),
      grouping_note = grouping_note,
      has_call_args = length(args) > 0L
    ))
  }

  paired <- ms_pairwise_eval_logical(ms_pairwise_arg(args, "paired", NA_integer_), .env, FALSE)
  pool_default <- !isTRUE(paired)
  pool_sd <- ms_pairwise_eval_logical(ms_pairwise_arg(args, "pool.sd", NA_integer_), .env, pool_default)
  if (grepl("paired", method)) paired <- TRUE
  if (grepl("non-pooled", method)) pool_sd <- FALSE
  if (grepl("pooled sd", method) && !grepl("non-pooled", method)) pool_sd <- TRUE
  var_equal <- ms_pairwise_eval_logical(ms_pairwise_arg(args, "var.equal", NA_integer_), .env, FALSE)

  list(
    kind = "t_test",
    x = ms_safe_numeric(x_values),
    g = factor(g_values),
    group_levels = levels(factor(g_values)),
    paired = isTRUE(paired),
    pool_sd = isTRUE(pool_sd),
    alternative = alternative,
    var_equal = isTRUE(var_equal),
    grouping_note = grouping_note,
    has_call_args = length(args) > 0L
  )
}

ms_pairwise_htest_call_args <- function(.call) {
  call_text <- trimws(as.character(.call %||% ""))
  if (!nzchar(call_text)) return(list())
  expr <- tryCatch(parse(text = call_text)[[1]], error = function(e) NULL)
  if (!is.call(expr)) return(list())
  head <- paste(deparse(expr[[1]], width.cutoff = 500L), collapse = "")
  if (!grepl("pairwise\\.(t|wilcox)\\.test$", head)) return(list())
  as.list(expr)[-1L]
}

ms_pairwise_arg <- function(args, name, position) {
  if (!length(args)) return(NULL)
  nms <- names(args) %||% rep("", length(args))
  hit <- which(nms == name)
  if (length(hit)) return(args[[hit[[1L]]]])
  if (!is.na(position) && length(args) >= position && !nzchar(nms[[position]])) {
    return(args[[position]])
  }
  NULL
}

ms_pairwise_data_name_exprs <- function(data_name) {
  parts <- strsplit(as.character(data_name %||% ""), "\\s+and\\s+")[[1]]
  if (length(parts) != 2L) return(list(x = NULL, g = NULL))
  list(
    x = tryCatch(parse(text = parts[[1]])[[1]], error = function(e) NULL),
    g = tryCatch(parse(text = parts[[2]])[[1]], error = function(e) NULL)
  )
}

ms_pairwise_eval_or_null <- function(expr, env) {
  if (is.null(expr) || is.null(env)) return(NULL)
  tryCatch(eval(expr, envir = env), error = function(e) NULL)
}

ms_pairwise_eval_logical <- function(expr, env, default) {
  value <- ms_pairwise_eval_or_null(expr, env)
  if (is.null(value) || length(value) == 0L || is.na(value[[1]])) return(default)
  isTRUE(value[[1]])
}

ms_pairwise_eval_nullable_logical <- function(expr, env) {
  value <- ms_pairwise_eval_or_null(expr, env)
  if (is.null(value) || length(value) == 0L || is.na(value[[1]])) return(NULL)
  isTRUE(value[[1]])
}

ms_pairwise_eval_character <- function(expr, env, default) {
  value <- ms_pairwise_eval_or_null(expr, env)
  if (is.null(value) || length(value) == 0L || is.na(value[[1]])) return(default)
  as.character(value[[1]])
}

ms_pairwise_numeric_grouping_note <- function(group_expr, group_values) {
  if (!is.numeric(group_values)) return(NULL)
  values <- unique(group_values[!is.na(group_values)])
  if (length(values) < 2L) return(NULL)
  label <- ms_pairwise_expr_label(group_expr)
  if (!nzchar(label)) label <- "The grouping variable"
  paste0(
    label,
    " is numeric in the supplied data; Mellio treats its distinct values as group levels to match the R test."
  )
}

ms_pairwise_expr_label <- function(expr) {
  if (is.null(expr)) return("")
  text <- trimws(gsub("\\s+", " ", paste(deparse(expr, width.cutoff = 500L), collapse = " ")))
  if (!nzchar(text)) return("")
  if (exists("ms_htest_reference_name", mode = "function")) {
    label <- tryCatch(ms_htest_reference_name(text), error = function(e) "")
    if (nzchar(label)) return(label)
  }
  text
}

ms_pairwise_enrich_htest_row <- function(row, context) {
  if (is.null(context) || is.null(row$level_1) || is.null(row$level_2)) return(row)
  keep <- !is.na(context$x) & !is.na(context$g)
  x <- context$x[keep]
  g <- droplevels(context$g[keep])
  i <- row$level_1
  j <- row$level_2
  xi <- x[as.character(g) == i]
  xj <- x[as.character(g) == j]
  if (!length(xi) || !length(xj)) return(row)

  if (identical(context$kind, "wilcoxon")) {
    row$median_1 <- ms_safe_numeric(stats::median(xi, na.rm = TRUE))
    row$median_2 <- ms_safe_numeric(stats::median(xj, na.rm = TRUE))
    row$paired <- isTRUE(context$paired)

    if (isTRUE(context$paired)) {
      if (length(xi) != length(xj)) return(row)
      row$n_pairs <- as.integer(length(xi))
      diffs <- xi - xj
      row$median_difference <- ms_safe_numeric(stats::median(diffs, na.rm = TRUE))
      effect <- ms_safe_numeric(ms_signed_rank_biserial(diffs))
      if (!is.na(effect)) {
        row$effect_size <- effect
        row$effect_size_name <- "rank_biserial"
        row$effect_size_label <- "Rank-biserial r"
      }
    } else {
      row$n_1 <- as.integer(length(xi))
      row$n_2 <- as.integer(length(xj))
      row$median_difference <- ms_safe_numeric(row$median_1 - row$median_2)
      effect <- ms_safe_numeric(ms_cliffs_delta(xi, xj))
      if (!is.na(effect)) {
        row$effect_size <- effect
        row$effect_size_name <- "cliffs_delta"
        row$effect_size_label <- "Cliff's delta"
      }
    }

    wt_args <- list(
      x = xi,
      y = xj,
      paired = isTRUE(context$paired),
      alternative = context$alternative %||% "two.sided",
      correct = isTRUE(context$correct)
    )
    if (!is.null(context$exact)) wt_args$exact <- context$exact
    wt <- tryCatch(
      suppressWarnings(do.call(stats::wilcox.test, wt_args)),
      error = function(e) NULL
    )
    if (!is.null(wt)) {
      row$statistic <- ms_safe_numeric(unname(wt$statistic[[1]] %||% NA_real_))
      row$p_raw <- ms_safe_numeric(wt$p.value %||% NA_real_)
      stat_names <- names(wt$statistic)
      row$statistic_label <- if (length(stat_names) && nzchar(stat_names[[1]])) {
        stat_names[[1]]
      } else if (isTRUE(context$paired)) {
        "V"
      } else {
        "W"
      }
    }
    return(row)
  }

  row <- ms_pairwise_add_t_descriptives(row, xi, xj, paired = isTRUE(context$paired))
  row$estimate <- ms_safe_numeric(mean(xi, na.rm = TRUE) - mean(xj, na.rm = TRUE))
  if (isTRUE(context$pool_sd) && !isTRUE(context$paired)) {
    pooled <- ms_pairwise_all_groups_pooled_sd(x, g)
    if (!is.null(pooled) && all(c(i, j) %in% pooled$levels)) {
      se_diff <- pooled$sd * sqrt(1 / pooled$n[[i]] + 1 / pooled$n[[j]])
      if (is.finite(se_diff) && se_diff > 0) {
        row$std_error <- ms_safe_numeric(se_diff)
        row$df <- ms_safe_numeric(pooled$df)
        row$statistic <- ms_safe_numeric(row$estimate / se_diff)
        row$statistic_label <- "t"
        row$p_raw <- ms_safe_numeric(
          2 * stats::pt(abs(row$statistic), df = pooled$df, lower.tail = FALSE)
        )
        row <- ms_pairwise_set_cohens_d(row, pooled$sd, "all_groups_pooled_sd")
      }
    }
  } else {
    tt <- tryCatch(
      stats::t.test(
        xi, xj,
        paired = isTRUE(context$paired),
        alternative = context$alternative %||% "two.sided",
        var.equal = isTRUE(context$var_equal)
      ),
      error = function(e) NULL
    )
    if (!is.null(tt)) {
      row$df <- ms_safe_numeric(unname(tt$parameter[[1]] %||% NA_real_))
      row$statistic <- ms_safe_numeric(unname(tt$statistic[[1]] %||% NA_real_))
      row$statistic_label <- "t"
      row$p_raw <- ms_safe_numeric(tt$p.value %||% NA_real_)
      if (!is.null(tt$stderr)) row$std_error <- ms_safe_numeric(tt$stderr)
    }
  }
  row
}

ms_pairwise_add_t_descriptives <- function(row, xi, xj, paired = FALSE) {
  xi <- xi[is.finite(xi)]
  xj <- xj[is.finite(xj)]
  if (!length(xi) || !length(xj)) return(row)

  row$mean_1 <- ms_safe_numeric(mean(xi))
  row$mean_2 <- ms_safe_numeric(mean(xj))

  if (isTRUE(paired) && length(xi) == length(xj)) {
    row$paired <- TRUE
    row$n_pairs <- as.integer(length(xi))
    row$n_1 <- NULL
    row$n_2 <- NULL
    effect <- ms_pairwise_cohens_dz(xi - xj)
    if (!is.na(effect)) {
      row$effect_size <- effect
      row$effect_size_name <- "cohens_dz"
      row$effect_size_label <- "Cohen's dz"
    }
    return(row)
  }

  row$paired <- FALSE
  row$n_1 <- as.integer(length(xi))
  row$n_2 <- as.integer(length(xj))
  effect <- ms_pairwise_cohens_d_independent(xi, xj)
  if (!is.na(effect)) {
    row$effect_size <- effect
    row$effect_size_name <- "cohens_d"
    row$effect_size_label <- "Cohen's d"
    row$effect_size_method <- "pairwise_pooled_sd"
  }
  row
}

ms_pairwise_cohens_d_independent <- function(a, b) {
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  n1 <- length(a)
  n2 <- length(b)
  if (n1 < 2L || n2 < 2L) return(NA_real_)
  s1 <- stats::sd(a)
  s2 <- stats::sd(b)
  denom <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  value <- (mean(a) - mean(b)) / denom
  if (is.finite(value)) ms_safe_numeric(value) else NA_real_
}

ms_pairwise_all_groups_pooled_sd <- function(x, g) {
  keep <- is.finite(x) & !is.na(g)
  x <- x[keep]
  g <- droplevels(factor(g[keep]))
  if (!length(x) || length(levels(g)) < 2L) return(NULL)
  by_group <- split(x, g)
  n <- vapply(by_group, length, integer(1))
  s <- vapply(by_group, stats::sd, numeric(1))
  degf <- n - 1L
  total_degf <- sum(degf)
  if (!is.finite(total_degf) || total_degf <= 0) return(NULL)
  pooled_sd <- sqrt(sum(ifelse(degf, s^2, 0) * degf) / total_degf)
  if (!is.finite(pooled_sd) || pooled_sd <= 0) return(NULL)
  list(sd = pooled_sd, df = total_degf, n = n, levels = levels(g))
}

ms_pairwise_set_cohens_d <- function(row, denominator, method) {
  estimate <- ms_safe_numeric(row$estimate %||% NA_real_)
  denominator <- ms_safe_numeric(denominator %||% NA_real_)
  if (is.na(estimate) || is.na(denominator) || denominator <= 0) return(row)
  value <- estimate / denominator
  if (!is.finite(value)) return(row)
  row$effect_size <- ms_safe_numeric(value)
  row$effect_size_name <- "cohens_d"
  row$effect_size_label <- "Cohen's d"
  row$effect_size_method <- method
  row
}

ms_pairwise_cohens_dz <- function(diffs) {
  diffs <- diffs[is.finite(diffs)]
  if (length(diffs) < 2L) return(NA_real_)
  denom <- stats::sd(diffs)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  value <- mean(diffs) / denom
  if (is.finite(value)) ms_safe_numeric(value) else NA_real_
}

ms_ttest_rows_from_rstatix <- function(df, args = list(), .call = NULL) {
  required <- c("group1", "group2", "p")
  if (!all(required %in% names(df))) {
    stop("rstatix::t_test output is missing required comparison columns.", call. = FALSE)
  }

  context <- ms_rstatix_ttest_context(args, .call = .call)
  by_cols <- setdiff(
    names(df),
    c("estimate", "estimate1", "estimate2", ".y.", "group1", "group2",
      "n1", "n2", "statistic", "df", "p", "p.adj", "p.adj.signif",
      "p.signif", "method", "alternative", "conf.low", "conf.high")
  )

  lapply(seq_len(nrow(df)), function(i) {
    group1 <- as.character(df$group1[[i]])
    group2 <- as.character(df$group2[[i]])
    row <- list(
      contrast = ms_pairwise_clean_contrast(paste(group1, group2, sep = " - ")),
      level_1 = group1,
      level_2 = group2,
      estimate = ms_dunn_cell(df, i, c("estimate", "mean.diff", "difference")),
      mean_1 = ms_dunn_cell(df, i, c("estimate1", "mean1")),
      mean_2 = ms_dunn_cell(df, i, c("estimate2", "mean2")),
      n_1 = ms_rstatix_count_cell(df, i, "n1"),
      n_2 = ms_rstatix_count_cell(df, i, "n2"),
      statistic = ms_dunn_cell(df, i, c("statistic", "t")),
      statistic_label = "t",
      df = ms_dunn_cell(df, i, "df"),
      std_error = ms_dunn_cell(df, i, c("se", "SE", "stderr", "std.error")),
      ci_lower = ms_dunn_cell(df, i, c("conf.low", "conf.low.", "ci.low")),
      ci_upper = ms_dunn_cell(df, i, c("conf.high", "conf.high.", "ci.high")),
      p_value = ms_dunn_cell(df, i, c("p.adj", "p.adjusted", "p"))
    )
    raw_p <- ms_dunn_cell(df, i, c("p", "p.unadj", "p.raw"))
    if (!is.na(raw_p) && any(c("p.adj", "p.adjusted") %in% names(df))) {
      row$p_raw <- raw_p
    }
    if (length(by_cols)) {
      by <- ms_dunn_by_label(df[i, , drop = FALSE], by_cols)
      if (nzchar(by)) row$by <- by
    }
    row <- ms_pairwise_canonicalize_row(row, level_order = context$group_levels)
    row <- ms_rstatix_ttest_enrich(row, context, by_cols, df[i, , drop = FALSE])
    row <- ms_rstatix_ttest_count_fallback(row, df, i, context$paired)
    row
  })
}

ms_rstatix_ttest_context <- function(args = list(), .call = NULL) {
  data <- args$data %||% NULL
  formula <- args$formula %||% NULL
  vars <- if (!is.null(formula)) all.vars(formula) else character(0)
  paired <- ms_rstatix_ttest_paired_arg(args, .call = .call)
  list(
    data = data,
    outcome = vars[[1L]] %||% NULL,
    group = vars[[2L]] %||% NULL,
    group_levels = ms_rstatix_group_levels(args),
    paired = paired,
    pool_sd = isTRUE(args$pool.sd %||% FALSE),
    var_equal = isTRUE(args$var.equal %||% FALSE),
    alternative = as.character(args$alternative %||% "two.sided"),
    conf_level = ms_safe_numeric(args$conf.level %||% 0.95)
  )
}

ms_rstatix_ttest_paired_arg <- function(args = list(), .call = NULL) {
  if (!is.null(args$paired) && length(args$paired) && !is.na(args$paired[[1L]])) {
    return(isTRUE(args$paired[[1L]]))
  }
  ms_rstatix_ttest_call_paired(.call)
}

ms_rstatix_ttest_call_paired <- function(.call = NULL) {
  call_text <- ms_pairwise_call_text(.call)
  if (!nzchar(call_text)) return(NA)
  expr <- tryCatch(parse(text = call_text)[[1L]], error = function(e) NULL)
  if (!is.call(expr)) return(NA)
  head <- paste(deparse(expr[[1L]], width.cutoff = 500L), collapse = "")
  if (!grepl("(?:^|::)(?:pairwise_)?t_test$", head)) return(NA)
  args <- as.list(expr)[-1L]
  nms <- names(args) %||% rep("", length(args))
  hit <- which(nms == "paired")
  if (!length(hit)) return(NA)
  value <- args[[hit[[1L]]]]
  if (identical(value, TRUE)) return(TRUE)
  if (identical(value, FALSE)) return(FALSE)
  if (is.logical(value) && length(value) == 1L && !is.na(value)) return(isTRUE(value))
  NA
}

ms_rstatix_ttest_enrich <- function(row, context, by_cols = character(0), by_row = NULL) {
  if (is.null(context$data) || is.null(context$outcome) || is.null(context$group)) {
    return(row)
  }
  data <- tryCatch(as.data.frame(context$data, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(data) || !all(c(context$outcome, context$group) %in% names(data))) return(row)

  if (length(by_cols) && !is.null(by_row)) {
    for (by_col in by_cols) {
      if (!by_col %in% names(data) || !by_col %in% names(by_row)) next
      value <- as.character(by_row[[by_col]][[1L]])
      data <- data[as.character(data[[by_col]]) == value, , drop = FALSE]
    }
  }

  x <- ms_safe_numeric(data[[context$outcome]])
  g <- as.character(data[[context$group]])
  keep <- is.finite(x) & !is.na(g)
  x <- x[keep]
  g <- g[keep]
  xi <- x[g == row$level_1]
  xj <- x[g == row$level_2]
  if (!length(xi) || !length(xj)) return(row)

  row <- ms_pairwise_add_t_descriptives(row, xi, xj, paired = isTRUE(context$paired))
  row$estimate <- ms_safe_numeric(row$mean_1 - row$mean_2)
  context$data <- data
  fit <- ms_rstatix_ttest_try_test(row, context, xi, xj)
  if (!is.null(fit)) {
    row$df <- ms_safe_numeric(unname(fit$parameter[[1L]] %||% row$df %||% NA_real_))
    row$statistic <- ms_safe_numeric(unname(fit$statistic[[1L]] %||% row$statistic %||% NA_real_))
    row$statistic_label <- "t"
    if (!is.null(fit$stderr)) row$std_error <- ms_safe_numeric(fit$stderr)
    if (!is.null(fit$conf.int) && length(fit$conf.int) >= 2L) {
      row$ci_lower <- ms_safe_numeric(unname(fit$conf.int[[1L]]))
      row$ci_upper <- ms_safe_numeric(unname(fit$conf.int[[2L]]))
    }
  }
  if (isTRUE(context$pool_sd) && !isTRUE(context$paired)) {
    pooled <- ms_pairwise_all_groups_pooled_sd(x, g)
    if (!is.null(pooled)) {
      row <- ms_pairwise_set_cohens_d(row, pooled$sd, "all_groups_pooled_sd")
    }
  }
  row
}

ms_rstatix_ttest_try_test <- function(row, context, xi, xj) {
  if (isTRUE(context$pool_sd) && !isTRUE(context$paired)) {
    return(ms_rstatix_pooled_ttest(row, context, xi, xj))
  }
  args <- list(
    x = xi,
    y = xj,
    paired = isTRUE(context$paired),
    alternative = context$alternative %||% "two.sided",
    var.equal = isTRUE(context$var_equal),
    conf.level = ms_safe_numeric(context$conf_level %||% 0.95)
  )
  args$conf.level <- if (is.na(args$conf.level)) 0.95 else args$conf.level
  tryCatch(
    suppressWarnings(do.call(stats::t.test, args)),
    error = function(e) NULL
  )
}

ms_rstatix_pooled_ttest <- function(row, context, xi, xj) {
  data <- tryCatch(as.data.frame(context$data, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(data) || !all(c(context$outcome, context$group) %in% names(data))) return(NULL)
  x <- ms_safe_numeric(data[[context$outcome]])
  g <- as.character(data[[context$group]])
  keep <- is.finite(x) & !is.na(g)
  x <- x[keep]
  g <- factor(g[keep])
  pooled <- ms_pairwise_all_groups_pooled_sd(x, g)
  if (is.null(pooled)) return(NULL)
  n1 <- length(xi)
  n2 <- length(xj)
  se <- pooled$sd * sqrt(1 / n1 + 1 / n2)
  if (!is.finite(se) || se <= 0) return(NULL)
  estimate <- mean(xi) - mean(xj)
  t_value <- estimate / se
  p_value <- 2 * stats::pt(abs(t_value), df = pooled$df, lower.tail = FALSE)
  structure(
    list(
      statistic = stats::setNames(t_value, "t"),
      parameter = stats::setNames(pooled$df, "df"),
      p.value = p_value,
      stderr = se
    ),
    class = "htest"
  )
}

ms_rstatix_ttest_count_fallback <- function(row, df, i, paired = NA) {
  n1 <- ms_rstatix_count_cell(df, i, "n1")
  n2 <- ms_rstatix_count_cell(df, i, "n2")
  if (isTRUE(row$paired) || isTRUE(paired)) {
    if ((is.null(row$n_pairs) || is.na(row$n_pairs)) && !is.na(n1) && !is.na(n2)) {
      row$n_pairs <- min(n1, n2)
    }
    row$n_1 <- NULL
    row$n_2 <- NULL
    row$paired <- TRUE
  } else {
    if (is.null(row$n_1) || is.na(row$n_1)) row$n_1 <- n1
    if (is.null(row$n_2) || is.na(row$n_2)) row$n_2 <- n2
    row$paired <- FALSE
  }
  row
}

ms_rstatix_count_cell <- function(df, i, col) {
  if (!col %in% names(df)) return(NA_integer_)
  value <- as.integer(ms_safe_numeric(df[[col]][[i]]))
  if (is.na(value)) NA_integer_ else value
}

ms_rstatix_ttest_method <- function(args = list(), paired = logical(0)) {
  if (length(paired) && all(paired)) return("paired t test")
  if (isTRUE(args$pool.sd %||% FALSE)) return("t tests with pooled SD")
  if (isTRUE(args$var.equal %||% FALSE)) return("Student's t test")
  "Welch t test"
}

ms_rstatix_grouping_note <- function(args = list()) {
  data <- args$data %||% NULL
  formula <- args$formula %||% NULL
  if (is.null(data) || is.null(formula)) return(NULL)
  vars <- all.vars(formula)
  group <- vars[[2L]] %||% NULL
  if (is.null(group)) return(NULL)
  data <- tryCatch(as.data.frame(data, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(data) || !group %in% names(data)) return(NULL)
  ms_pairwise_numeric_grouping_note(as.name(group), data[[group]])
}

ms_rstatix_group_levels <- function(args = list()) {
  data <- args$data %||% NULL
  formula <- args$formula %||% NULL
  if (is.null(data) || is.null(formula)) return(NULL)
  vars <- all.vars(formula)
  group <- vars[[2L]] %||% NULL
  if (is.null(group)) return(NULL)
  data <- tryCatch(as.data.frame(data, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(data) || !group %in% names(data)) return(NULL)
  values <- data[[group]]
  if (is.factor(values)) return(levels(droplevels(values)))
  values <- values[!is.na(values)]
  if (!length(values)) return(NULL)
  if (is.numeric(values)) return(as.character(sort(unique(values))))
  levels(factor(values))
}

ms_wilcox_rows_from_rstatix <- function(df, args = list(), .call = NULL) {
  required <- c("group1", "group2", "statistic", "p")
  if (!all(required %in% names(df))) {
    stop("rstatix::wilcox_test output is missing required comparison columns.", call. = FALSE)
  }

  context <- ms_rstatix_wilcox_context(args, .call = .call)
  by_cols <- setdiff(
    names(df),
    c("estimate", ".y.", "group1", "group2", "n1", "n2", "statistic",
      "p", "p.adj", "p.adj.signif", "p.signif", "method", "alternative",
      "conf.low", "conf.high")
  )

  lapply(seq_len(nrow(df)), function(i) {
    group1 <- as.character(df$group1[[i]])
    group2 <- as.character(df$group2[[i]])
    row <- list(
      contrast = ms_pairwise_clean_contrast(paste(group1, group2, sep = " - ")),
      level_1 = group1,
      level_2 = group2,
      statistic = ms_safe_numeric(df$statistic[[i]]),
      p_value = ms_dunn_cell(df, i, c("p.adj", "p.adjusted", "p"))
    )
    raw_p <- ms_dunn_cell(df, i, c("p", "p.unadj", "p.raw"))
    if (!is.na(raw_p) && any(c("p.adj", "p.adjusted") %in% names(df))) {
      row$p_raw <- raw_p
    }
    if (length(by_cols)) {
      by <- ms_dunn_by_label(df[i, , drop = FALSE], by_cols)
      if (nzchar(by)) row$by <- by
    }
    row <- ms_pairwise_canonicalize_row(row, level_order = context$group_levels, flip_values = FALSE)
    row <- ms_rstatix_wilcox_enrich(row, context, by_cols, df[i, , drop = FALSE])
    row <- ms_rstatix_wilcox_count_fallback(row, df, i, context$paired)
    if (isTRUE(row$paired)) {
      row$statistic_label <- "V"
    } else {
      row$statistic_label <- "W"
    }
    row
  })
}

ms_rstatix_wilcox_context <- function(args = list(), .call = NULL) {
  data <- args$data %||% NULL
  formula <- args$formula %||% NULL
  vars <- if (!is.null(formula)) all.vars(formula) else character(0)
  paired <- ms_rstatix_wilcox_paired_arg(args, .call = .call)
  list(
    data = data,
    outcome = vars[[1L]] %||% NULL,
    group = vars[[2L]] %||% NULL,
    group_levels = ms_rstatix_group_levels(args),
    paired = paired,
    alternative = as.character(args$alternative %||% "two.sided"),
    exact = args$exact %||% NULL
  )
}

ms_rstatix_wilcox_paired_arg <- function(args = list(), .call = NULL) {
  if (!is.null(args$paired) && length(args$paired) && !is.na(args$paired[[1L]])) {
    return(isTRUE(args$paired[[1L]]))
  }
  ms_rstatix_wilcox_call_paired(.call)
}

ms_rstatix_wilcox_call_paired <- function(.call = NULL) {
  call_text <- ms_pairwise_call_text(.call)
  if (!nzchar(call_text)) return(NA)
  expr <- tryCatch(parse(text = call_text)[[1L]], error = function(e) NULL)
  if (!is.call(expr)) return(NA)
  head <- paste(deparse(expr[[1L]], width.cutoff = 500L), collapse = "")
  if (!grepl("(?:^|::)(?:pairwise_)?wilcox_test$", head)) return(NA)
  args <- as.list(expr)[-1L]
  nms <- names(args) %||% rep("", length(args))
  hit <- which(nms == "paired")
  if (!length(hit)) return(NA)
  value <- args[[hit[[1L]]]]
  if (identical(value, TRUE)) return(TRUE)
  if (identical(value, FALSE)) return(FALSE)
  if (is.logical(value) && length(value) == 1L && !is.na(value)) return(isTRUE(value))
  NA
}

ms_rstatix_wilcox_enrich <- function(row, context, by_cols = character(0), by_row = NULL) {
  if (is.null(context$data) || is.null(context$outcome) || is.null(context$group)) {
    return(row)
  }
  data <- tryCatch(as.data.frame(context$data, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(data) || !all(c(context$outcome, context$group) %in% names(data))) return(row)

  if (length(by_cols) && !is.null(by_row)) {
    for (by_col in by_cols) {
      if (!by_col %in% names(data) || !by_col %in% names(by_row)) next
      value <- as.character(by_row[[by_col]][[1L]])
      data <- data[as.character(data[[by_col]]) == value, , drop = FALSE]
    }
  }

  x <- ms_safe_numeric(data[[context$outcome]])
  g <- as.character(data[[context$group]])
  keep <- is.finite(x) & !is.na(g)
  x <- x[keep]
  g <- g[keep]
  xi <- x[g == row$level_1]
  xj <- x[g == row$level_2]
  if (!length(xi) || !length(xj)) return(row)

  row$median_1 <- ms_safe_numeric(stats::median(xi, na.rm = TRUE))
  row$median_2 <- ms_safe_numeric(stats::median(xj, na.rm = TRUE))
  paired <- ms_rstatix_wilcox_row_paired(row, context, xi, xj)
  row$paired <- isTRUE(paired)

  if (isTRUE(paired) && length(xi) == length(xj)) {
    row$n_pairs <- as.integer(length(xi))
    row$n_1 <- NULL
    row$n_2 <- NULL
    diffs <- xi - xj
    row$median_difference <- ms_safe_numeric(stats::median(diffs, na.rm = TRUE))
    effect <- ms_safe_numeric(ms_signed_rank_biserial(diffs))
    if (!is.na(effect)) {
      row$effect_size <- effect
      row$effect_size_name <- "rank_biserial"
      row$effect_size_label <- "Rank-biserial r"
    }
  } else {
    row$paired <- FALSE
    row$n_1 <- as.integer(length(xi))
    row$n_2 <- as.integer(length(xj))
    row$median_difference <- ms_safe_numeric(row$median_1 - row$median_2)
    effect <- ms_safe_numeric(ms_cliffs_delta(xi, xj))
    if (!is.na(effect)) {
      row$effect_size <- effect
      row$effect_size_name <- "cliffs_delta"
      row$effect_size_label <- "Cliff's delta"
    }
  }
  fit <- ms_rstatix_wilcox_try_test(row, context, xi, xj, paired = isTRUE(row$paired))
  if (!is.null(fit)) {
    row$statistic <- ms_safe_numeric(unname(fit$statistic[[1L]] %||% row$statistic %||% NA_real_))
    if (is.null(row$p_raw) || is.na(row$p_raw)) {
      row$p_raw <- ms_safe_numeric(fit$p.value %||% NA_real_)
    }
  }
  row
}

ms_rstatix_wilcox_row_paired <- function(row, context, xi, xj) {
  if (isTRUE(context$paired)) return(TRUE)
  if (identical(context$paired, FALSE)) return(FALSE)
  if (length(xi) != length(xj)) return(FALSE)

  paired_fit <- ms_rstatix_wilcox_try_test(row, context, xi, xj, paired = TRUE)
  unpaired_fit <- ms_rstatix_wilcox_try_test(row, context, xi, xj, paired = FALSE)
  paired_score <- ms_rstatix_wilcox_match_score(row, paired_fit)
  unpaired_score <- ms_rstatix_wilcox_match_score(row, unpaired_fit)
  if (is.na(paired_score) || is.na(unpaired_score)) return(FALSE)
  paired_score + sqrt(.Machine$double.eps) < unpaired_score
}

ms_rstatix_wilcox_try_test <- function(row, context, xi, xj, paired = FALSE) {
  args <- list(
    x = xi,
    y = xj,
    paired = isTRUE(paired),
    alternative = context$alternative %||% "two.sided"
  )
  if (!is.null(context$exact)) args$exact <- context$exact
  tryCatch(
    suppressWarnings(do.call(stats::wilcox.test, args)),
    error = function(e) NULL
  )
}

ms_rstatix_wilcox_match_score <- function(row, fit) {
  if (is.null(fit)) return(NA_real_)
  score <- 0
  used <- FALSE
  row_stat <- ms_safe_numeric(row$statistic %||% NA_real_)
  fit_stat <- ms_safe_numeric(unname(fit$statistic[[1L]] %||% NA_real_))
  if (!is.na(row_stat) && !is.na(fit_stat)) {
    score <- score + abs(row_stat - fit_stat)
    used <- TRUE
  }
  row_p <- ms_safe_numeric(row$p_raw %||% row$p_value %||% NA_real_)
  fit_p <- ms_safe_numeric(fit$p.value %||% NA_real_)
  if (!is.na(row_p) && !is.na(fit_p)) {
    score <- score + abs(row_p - fit_p)
    used <- TRUE
  }
  if (used) score else NA_real_
}

ms_rstatix_wilcox_count_fallback <- function(row, df, i, paired = NA) {
  n1 <- if ("n1" %in% names(df)) as.integer(ms_safe_numeric(df$n1[[i]])) else NA_integer_
  n2 <- if ("n2" %in% names(df)) as.integer(ms_safe_numeric(df$n2[[i]])) else NA_integer_
  missing_n_pairs <- is.null(row$n_pairs) || is.na(row$n_pairs)
  if (isTRUE(row$paired) || (isTRUE(paired) && missing_n_pairs)) {
    if (missing_n_pairs && !is.na(n1) && !is.na(n2)) row$n_pairs <- min(n1, n2)
    row$n_1 <- NULL
    row$n_2 <- NULL
    row$paired <- TRUE
  } else {
    if (is.null(row$n_1) || is.na(row$n_1)) row$n_1 <- n1
    if (is.null(row$n_2) || is.na(row$n_2)) row$n_2 <- n2
    row$paired <- FALSE
  }
  row
}

ms_rstatix_adjustment_from_columns <- function(df) {
  if (any(c("p.adj", "p.adjusted") %in% names(df))) return("adjusted")
  "none"
}

ms_dunn_rows_from_rstatix <- function(df, args = list()) {
  required <- c("group1", "group2", "statistic")
  if (!all(required %in% names(df))) {
    stop("rstatix::dunn_test output is missing required comparison columns.", call. = FALSE)
  }

  data <- args$data %||% NULL
  formula <- args$formula %||% NULL
  vars <- if (!is.null(formula)) all.vars(formula) else character(0)
  outcome <- vars[[1L]] %||% NULL
  group <- vars[[2L]] %||% NULL
  level_order <- ms_rstatix_group_levels(args)
  by_cols <- setdiff(
    names(df),
    c(".y.", "group1", "group2", "n1", "n2", "statistic", "p", "p.adj",
      "p.adj.signif", "p.signif", "method", "estimate", "estimate1", "estimate2")
  )

  lapply(seq_len(nrow(df)), function(i) {
    group1 <- as.character(df$group1[[i]])
    group2 <- as.character(df$group2[[i]])
    row <- list(
      contrast = ms_pairwise_clean_contrast(paste(group2, group1, sep = " - ")),
      level_1 = group2,
      level_2 = group1,
      statistic = ms_safe_numeric(df$statistic[[i]]),
      statistic_label = "z",
      p_value = ms_dunn_cell(df, i, c("p.adj", "p.adjusted", "p"))
    )
    raw_p <- ms_dunn_cell(df, i, c("p", "p.unadj", "p.raw"))
    if (!is.na(raw_p) && "p.adj" %in% names(df)) row$p_raw <- raw_p
    if ("n2" %in% names(df)) row$n_1 <- as.integer(ms_safe_numeric(df$n2[[i]]))
    if ("n1" %in% names(df)) row$n_2 <- as.integer(ms_safe_numeric(df$n1[[i]]))
    if ("estimate" %in% names(df)) {
      row$rank_mean_difference <- ms_dunn_cell(df, i, "estimate")
    }
    if (length(by_cols)) {
      by <- ms_dunn_by_label(df[i, , drop = FALSE], by_cols)
      if (nzchar(by)) row$by <- by
    }
    row <- ms_pairwise_canonicalize_row(row, level_order = level_order)
    ms_dunn_enrich_from_data(row, data, outcome, group, by_cols, df[i, , drop = FALSE])
  })
}

ms_dunn_rows_from_comparison_df <- function(df) {
  comparison_col <- ms_dunn_first_col(names(df), c("Comparison", "comparison", "comparisons", "contrast"))
  if (is.na(comparison_col)) {
    stop("Dunn test comparison table is missing a comparison column.", call. = FALSE)
  }

  lapply(seq_len(nrow(df)), function(i) {
    comparison <- as.character(df[[comparison_col]][[i]])
    levels <- ms_dunn_comparison_levels(comparison)
    list(
      contrast = ms_pairwise_clean_contrast(paste(levels[[1L]], levels[[2L]], sep = " - ")),
      level_1 = levels[[1L]],
      level_2 = levels[[2L]],
      statistic = ms_dunn_cell(df, i, c("Z", "z", "statistic")),
      statistic_label = "z",
      p_raw = ms_dunn_cell(df, i, c("P.unadj", "P.unadjusted", "P.raw", "P", "p")),
      p_value = ms_dunn_cell(df, i, c("P.adj", "P.adjusted", "P.adjust", "p.adj", "p.adjusted", "p"))
    )
  })
}

ms_dunn_rows_from_dunn_test_object <- function(x) {
  comparisons <- x$comparisons %||% x$Comparison %||% x$comparison %||% names(x$Z)
  if (is.null(comparisons) || !length(comparisons)) {
    stop("dunn.test object is missing comparison labels.", call. = FALSE)
  }
  n <- length(comparisons)
  df <- data.frame(
    Comparison = as.character(comparisons),
    Z = ms_dunn_recycle(x$Z %||% x$z %||% x$statistic, n),
    P = ms_dunn_recycle(x$P %||% x$p %||% x$P.unadj, n),
    P.adj = ms_dunn_recycle(x$P.adjusted %||% x$P.adj %||% x$p.adjusted %||% x$P, n),
    stringsAsFactors = FALSE
  )
  ms_dunn_rows_from_comparison_df(df)
}

ms_dunn_enrich_from_data <- function(row, data, outcome, group, by_cols = character(0),
                                     by_row = NULL) {
  if (is.null(data) || is.null(outcome) || is.null(group)) return(row)
  data <- tryCatch(as.data.frame(data, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(data) || !all(c(outcome, group) %in% names(data))) return(row)

  if (length(by_cols) && !is.null(by_row)) {
    for (by_col in by_cols) {
      if (!by_col %in% names(data) || !by_col %in% names(by_row)) next
      value <- as.character(by_row[[by_col]][[1L]])
      data <- data[as.character(data[[by_col]]) == value, , drop = FALSE]
    }
  }

  x <- ms_safe_numeric(data[[outcome]])
  g <- as.character(data[[group]])
  keep <- is.finite(x) & !is.na(g)
  x <- x[keep]
  g <- g[keep]
  xi <- x[g == row$level_1]
  xj <- x[g == row$level_2]
  if (!length(xi) || !length(xj)) return(row)

  if (is.null(row$n_1) || is.na(row$n_1)) row$n_1 <- as.integer(length(xi))
  if (is.null(row$n_2) || is.na(row$n_2)) row$n_2 <- as.integer(length(xj))
  row$median_1 <- ms_safe_numeric(stats::median(xi, na.rm = TRUE))
  row$median_2 <- ms_safe_numeric(stats::median(xj, na.rm = TRUE))
  row$median_difference <- ms_safe_numeric(row$median_1 - row$median_2)
  effect <- ms_safe_numeric(ms_cliffs_delta(xi, xj))
  if (!is.na(effect)) {
    row$effect_size <- effect
    row$effect_size_name <- "cliffs_delta"
    row$effect_size_label <- "Cliff's delta"
  }
  row
}

ms_dunn_cell <- function(df, i, cols) {
  col <- ms_dunn_first_col(names(df), cols)
  if (is.na(col)) return(NA_real_)
  ms_safe_numeric(df[[col]][[i]])
}

ms_dunn_first_col <- function(names, candidates) {
  hit <- candidates[candidates %in% names]
  if (length(hit)) hit[[1L]] else NA_character_
}

ms_dunn_comparison_levels <- function(comparison) {
  comparison <- trimws(as.character(comparison %||% ""))
  parts <- strsplit(comparison, "(?i)\\s+-\\s+|\\s+vs\\.?\\s+|\\s+versus\\s+", perl = TRUE)[[1L]]
  if (length(parts) < 2L) parts <- strsplit(comparison, "\\s*-\\s*")[[1L]]
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]
  if (length(parts) >= 2L) return(parts[1:2])
  c(comparison, "")
}

ms_dunn_by_label <- function(row, by_cols) {
  parts <- vapply(by_cols, function(col) {
    if (!col %in% names(row)) return("")
    paste0(col, "=", as.character(row[[col]][[1L]]))
  }, character(1))
  paste(parts[nzchar(parts)], collapse = ", ")
}

ms_dunn_adjustment_from_columns <- function(df) {
  if ("p.adj" %in% names(df) || "P.adj" %in% names(df)) return("adjusted")
  ""
}

ms_dunn_adjustment_from_text <- function(text) {
  lower <- tolower(paste(as.character(text %||% ""), collapse = " "))
  known <- c("holm", "bonferroni", "hochberg", "hommel", "fdr", "bh", "by", "sidak", "none")
  hit <- known[vapply(known, function(method) grepl(method, lower, fixed = TRUE), logical(1))]
  if (length(hit)) hit[[1L]] else ""
}

ms_dunn_recycle <- function(x, n) {
  if (is.null(x) || !length(x)) return(rep(NA_real_, n))
  x <- ms_safe_numeric(x)
  if (length(x) >= n) return(x[seq_len(n)])
  rep(x, length.out = n)
}

ms_emm_pairwise_summary <- function(x) {
  misc <- tryCatch(x@misc, error = function(e) list())
  s <- tryCatch(
    suppressMessages(suppressWarnings(summary(x, infer = c(TRUE, TRUE)))),
    error = function(e) NULL
  )
  if (is.null(s) || !"contrast" %in% names(s)) {
    paired <- tryCatch(graphics::pairs(x), error = function(e) NULL)
    if (is.null(paired)) {
      stop("emmGrid object does not contain pairwise contrasts.", call. = FALSE)
    }
    misc <- tryCatch(paired@misc, error = function(e) misc)
    s <- suppressMessages(suppressWarnings(summary(paired, infer = c(TRUE, TRUE))))
  }
  s <- as.data.frame(s)
  method <- as.character(misc$methDesc %||% "Estimated marginal means contrasts")
  if (tolower(trimws(method)) %in% c("pairwise differences", "pairwise contrasts")) {
    method <- "estimated marginal means contrasts"
  }
  adjustment <- as.character((attr(s, "adjust") %||% misc$adjust) %||% "")
  list(summary = s, method = method, adjustment = adjustment, misc = misc)
}

ms_emm_means_table <- function(x) {
  s <- tryCatch(
    suppressMessages(suppressWarnings(summary(x, infer = c(TRUE, TRUE)))),
    error = function(e) NULL
  )
  if (is.null(s)) return(NULL)
  est_col <- attr(s, "estName") %||% "emmean"
  cl_names <- attr(s, "clNames") %||% c("lower.CL", "upper.CL")
  low_col <- cl_names[[1L]] %||% "lower.CL"
  high_col <- cl_names[[2L]] %||% "upper.CL"
  pri_vars <- attr(s, "pri.vars") %||% character(0)
  by_vars <- attr(s, "by.vars") %||% character(0)
  conf_level <- ms_safe_numeric(attr(s, "level") %||% NA_real_)
  messages <- attr(s, "mesg") %||% character(0)
  df <- as.data.frame(s, stringsAsFactors = FALSE)
  if (!nrow(df)) return(NULL)

  level_cols <- unique(c(by_vars, pri_vars))
  level_cols <- level_cols[level_cols %in% names(df)]
  if (!length(level_cols)) {
    skip_cols <- c(est_col, "SE", "df", low_col, high_col, "t.ratio", "z.ratio", "p.value")
    level_cols <- setdiff(names(df), skip_cols)
  }
  level_cols <- level_cols[level_cols %in% names(df)]

  rows <- lapply(seq_len(nrow(df)), function(i) {
    row <- list()
    for (col in level_cols) row[[col]] <- as.character(df[[col]][[i]])
    if (est_col %in% names(df)) row$emmean <- ms_safe_numeric(df[[est_col]][[i]])
    if ("SE" %in% names(df)) row$std_error <- ms_safe_numeric(df$SE[[i]])
    if ("df" %in% names(df)) row$df <- ms_safe_numeric(df$df[[i]])
    if (low_col %in% names(df)) row$ci_lower <- ms_safe_numeric(df[[low_col]][[i]])
    if (high_col %in% names(df)) row$ci_upper <- ms_safe_numeric(df[[high_col]][[i]])
    row
  })

  columns <- lapply(level_cols, function(col) {
    list(key = col, label = col, format = "text")
  })
  columns <- c(columns, list(list(key = "emmean", label = "EMM", format = "number")))
  if ("SE" %in% names(df)) {
    columns <- c(columns, list(list(key = "std_error", label = "SE", format = "number")))
  }
  if ("df" %in% names(df)) {
    columns <- c(columns, list(list(key = "df", label = "df", format = "df")))
  }
  ci_label <- ms_pairwise_ci_column_label(conf_level)
  if (low_col %in% names(df) && high_col %in% names(df)) {
    columns <- c(columns, list(
      list(key = "ci_lower", label = paste(ci_label, "lower"), format = "number"),
      list(key = "ci_upper", label = paste(ci_label, "upper"), format = "number")
    ))
  }

  message_text <- trimws(as.character(messages))
  message_text <- message_text[nzchar(message_text)]
  message_text <- ifelse(grepl("[.!?]$", message_text), message_text, paste0(message_text, "."))
  note <- paste(message_text, collapse = " ")
  note <- trimws(gsub("\\s+", " ", note))
  if (!nzchar(note)) {
    note <- "Estimated marginal means are shown with confidence intervals where available."
  }
  interaction_note <- ms_emm_interaction_note(x, s, messages)
  if (nzchar(interaction_note) && !grepl("misleading due to involvement in interactions", note, ignore.case = TRUE)) {
    note <- trimws(paste(note, interaction_note))
  }

  list(
    title = "Estimated marginal means",
    rows = rows,
    columns = columns,
    note = note,
    conf_level = conf_level
  )
}

ms_emm_interaction_note <- function(x, summary_obj = NULL, messages = character(0)) {
  misc <- tryCatch(x@misc, error = function(e) list())
  pri_vars <- attr(summary_obj, "pri.vars") %||% misc$pri.vars %||% character(0)
  by_vars <- attr(summary_obj, "by.vars") %||% misc$by.vars %||% character(0)
  avgd_vars <- misc$avgd.over %||% character(0)
  if (!length(avgd_vars)) {
    avgd_line <- grep("averaged over the levels of:", messages, value = TRUE, ignore.case = TRUE)
    if (length(avgd_line)) {
      avgd_text <- sub("^.*averaged over the levels of:\\s*", "", avgd_line[[1L]], ignore.case = TRUE)
      avgd_vars <- unlist(strsplit(avgd_text, "\\s*,\\s*|\\s+and\\s+"), use.names = FALSE)
      avgd_vars <- gsub("[.;]\\s*$", "", trimws(avgd_vars))
    }
  }
  focal_vars <- unique(trimws(as.character(c(pri_vars, by_vars))))
  avgd_vars <- unique(trimws(as.character(avgd_vars)))
  focal_vars <- focal_vars[nzchar(focal_vars)]
  avgd_vars <- avgd_vars[nzchar(avgd_vars)]
  if (!length(focal_vars) || !length(avgd_vars)) return("")

  terms_obj <- tryCatch(x@model.info$terms, error = function(e) NULL)
  term_labels <- tryCatch(attr(stats::terms(terms_obj), "term.labels"), error = function(e) character(0))
  if (!length(term_labels)) {
    term_labels <- tryCatch(attr(terms_obj, "term.labels"), error = function(e) character(0))
  }
  interaction_terms <- term_labels[grepl(":", term_labels, fixed = TRUE)]
  if (!length(interaction_terms)) return("")

  for (term in interaction_terms) {
    pieces <- trimws(unlist(strsplit(term, ":", fixed = TRUE), use.names = FALSE))
    if (any(focal_vars %in% pieces) && any(avgd_vars %in% pieces)) {
      return("emmeans note: results may be misleading due to involvement in interactions.")
    }
  }
  ""
}

ms_emm_pairwise_contrast_meta <- function(x, n_rows) {
  empty <- vector("list", n_rows)
  misc <- tryCatch(x@misc, error = function(e) list())
  coefs <- misc$con.coef %||% NULL
  grid <- misc$orig.grid %||% NULL
  if (is.null(coefs) || is.null(grid)) return(empty)
  coefs <- tryCatch(as.matrix(coefs), error = function(e) NULL)
  grid <- tryCatch(as.data.frame(grid, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(coefs) || is.null(grid)) return(empty)
  if (nrow(coefs) < n_rows || ncol(coefs) != nrow(grid)) return(empty)

  lapply(seq_len(n_rows), function(i) {
    ms_emm_pairwise_contrast_row_meta(coefs[i, ], grid)
  })
}

ms_emm_pairwise_contrast_row_meta <- function(coefs, grid) {
  nz <- which(abs(coefs) > sqrt(.Machine$double.eps))
  if (length(nz) != 2L) return(NULL)
  plus <- nz[coefs[nz] > 0]
  minus <- nz[coefs[nz] < 0]
  if (length(plus) != 1L || length(minus) != 1L) return(NULL)
  if (!isTRUE(all.equal(abs(coefs[[plus]]), abs(coefs[[minus]]), tolerance = 1e-8))) {
    return(NULL)
  }

  plus_row <- grid[plus, , drop = FALSE]
  minus_row <- grid[minus, , drop = FALSE]
  differing <- names(grid)[vapply(names(grid), function(name) {
    !identical(as.character(plus_row[[name]][[1L]]), as.character(minus_row[[name]][[1L]]))
  }, logical(1))]
  if (length(differing) != 1L) return(NULL)

  focal <- differing[[1L]]
  plus_level <- as.character(plus_row[[focal]][[1L]])
  minus_level <- as.character(minus_row[[focal]][[1L]])
  level_order <- unique(as.character(grid[[focal]]))
  plus_pos <- match(plus_level, level_order)
  minus_pos <- match(minus_level, level_order)
  flip <- !is.na(plus_pos) && !is.na(minus_pos) && plus_pos < minus_pos
  level_1 <- if (flip) minus_level else plus_level
  level_2 <- if (flip) plus_level else minus_level

  list(
    contrast = ms_pairwise_clean_contrast(paste(level_1, level_2, sep = " - ")),
    level_1 = level_1,
    level_2 = level_2,
    flip = flip,
    factor = focal
  )
}

ms_pairwise_flip_interval <- function(ci_lower, ci_upper) {
  if (is.na(ci_lower) || is.na(ci_upper)) return(list(-ci_lower, -ci_upper))
  list(-ci_upper, -ci_lower)
}

ms_glht_pairwise_payload <- function(x, summary_obj = NULL, .call = NULL,
                                     expr = NULL, summary_input = FALSE) {
  if (is.null(summary_obj)) {
    summary_obj <- tryCatch(summary(x), error = function(e) NULL)
  }
  if (is.null(summary_obj) || is.null(summary_obj$test)) {
    stop("Could not extract multcomp pairwise comparisons.", call. = FALSE)
  }

  test <- summary_obj$test
  coefs <- test$coefficients %||% numeric(0)
  sigma <- test$sigma %||% rep(NA_real_, length(coefs))
  stat <- test$tstat %||% rep(NA_real_, length(coefs))
  pvals <- test$pvalues %||% rep(NA_real_, length(coefs))
  names <- names(coefs) %||% paste0("contrast_", seq_along(coefs))
  ci_obj <- tryCatch(stats::confint(x), error = function(e) NULL)
  ci <- ci_obj$confint %||% NULL
  conf_level <- ms_safe_numeric(attr(ci, "conf.level") %||% NA_real_)
  df <- ms_glht_df(x, summary_obj, test)
  stat_label <- if (!is.na(df)) "t" else "z"

  rows <- lapply(seq_along(coefs), function(i) {
    row <- ms_pairwise_add_levels(list(
      contrast = ms_pairwise_clean_contrast(names[[i]]),
      estimate = ms_safe_numeric(coefs[[i]]),
      std_error = ms_safe_numeric(sigma[[i]]),
      statistic = ms_safe_numeric(stat[[i]]),
      statistic_label = stat_label,
      p_value = ms_safe_numeric(pvals[[i]])
    ), names[[i]])
    if (!is.na(df)) row$df <- df
    if (!is.null(ci) && names[[i]] %in% rownames(ci)) {
      row$ci_lower <- ms_safe_numeric(ci[names[[i]], "lwr"])
      row$ci_upper <- ms_safe_numeric(ci[names[[i]], "upr"])
    }
    row
  })

  ms_pairwise_payload(
    rows = rows,
    method = ms_glht_pairwise_method(x, summary_obj),
    adjustment_method = as.character(test$type %||% ""),
    call = ms_pairwise_call(
      .call,
      expr,
      fallback = if (isTRUE(summary_input)) {
        "summary(multcomp::glht(...))"
      } else {
        "multcomp::glht(...)"
      }
    ),
    raw_output = ms_capture_output(summary_obj),
    conf_level = conf_level,
    packages = ms_packages_basic(extras = "multcomp"),
    source = "multcomp::glht"
  )
}

ms_glht_df <- function(x, summary_obj, test) {
  candidates <- list(test$df, summary_obj$df, x$df)
  for (candidate in candidates) {
    value <- ms_safe_numeric(candidate %||% NA_real_)
    if (!is.na(value) && is.finite(value) && value > 0) return(value)
  }
  NA_real_
}

ms_glht_pairwise_method <- function(x, summary_obj) {
  type <- trimws(as.character(
    (x$type %||% summary_obj$type %||% attr(x$linfct, "type")) %||% ""
  ))
  if (nzchar(type) && !tolower(type) %in% c("none", "general")) {
    return(paste(type, "contrasts"))
  }
  "General linear hypotheses"
}
