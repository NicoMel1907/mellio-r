# R bridge -- htest extractor.
#
# Covers: t.test (one-sample, two-sample, Welch, paired), cor.test
# (Pearson, Spearman, Kendall), wilcox.test (rank-sum and signed-rank),
# chisq.test, fisher.test, kruskal.test, friedman.test, and mcnemar.test.
#
# All return objects of class "htest" with a consistent shape:
# statistic, parameter (df), p.value, estimate, conf.int, alternative,
# method, data.name. Dispatch on $method to set the schema TypeKey
# and handle per-test quirks.
#
# Schema: docs/STATS-R-BRIDGE-SCHEMA.md

#' @rdname mellio_payload
#' @param .data Optional data.frame to enrich supported htest results with
#'   observed descriptives, effect sizes, and figure data. When omitted,
#'   Mellio tries to recover simple source variables from the calling
#'   environment before falling back to test-only output.
#' @param .env Optional environment used when recovering source variables for
#'   supported test objects.
#' @export
mellio_payload.htest <- function(x, .data = NULL, ..., .call = NULL,
                                 .env = parent.frame()) {
  # Capture the user's call. Priority: explicit .call passed by mellio_open
  # (since mellio_open's x is local-named "x" by that point); otherwise
  # use match.call()$x for direct mellio_payload(t.test(...)) callers.
  call_str <- if (!is.null(.call)) {
    .call
  } else {
    user_call <- match.call()$x
    if (!is.null(user_call) && !identical(user_call, as.name("x"))) {
      paste(deparse(user_call, width.cutoff = 500L), collapse = " ")
    } else {
      x$data.name %||% NA_character_
    }
  }

  type_info <- ms_htest_type(x$method)
  fields    <- ms_htest_fields(x, type_info$type, data = .data, env = .env)
  figure_data <- ms_htest_figure_data(x, type_info$type, fields, data = .data, env = .env)
  raw       <- ms_capture_output(x)

  ms_build_envelope(
    type       = type_info$type,
    type_label = type_info$type_label,
    call       = trimws(gsub("\\s+", " ", call_str)),
    fields     = fields,
    raw_output = raw,
    figure_data = figure_data
  )
}

# ── Type mapping ──────────────────────────────────────────────────────

# Map x$method to schema TypeKey + human label.
# x$method strings come from base R and are stable across versions.
ms_htest_type <- function(method) {
  m <- as.character(method)
  type <- if      (grepl("Welch.*Two Sample t",          m)) "welch_t_test"
          else if (grepl("Two Sample t",                 m)) "students_t_test"
          else if (grepl("Paired t",                     m)) "paired_t_test"
          else if (grepl("One Sample t",                 m)) "one_sample_t_test"
          else if (grepl("Pearson.*correlation",         m)) "pearson_correlation"
          else if (grepl("Spearman.*correlation",        m)) "spearman_correlation"
          else if (grepl("Kendall.*correlation",         m)) "kendall_correlation"
          else if (grepl("Wilcoxon rank sum",            m)) "wilcoxon_rank_sum"
          else if (grepl("Wilcoxon signed rank",         m)) "wilcoxon_signed_rank"
          else if (grepl("Kruskal-Wallis",               m)) "kruskal_wallis_test"
          else if (grepl("Friedman rank sum",             m)) "friedman_test"
          else if (grepl("McNemar",                      m)) "mcnemar_test"
          else if (grepl("Pearson.*Chi-squared",           m) ||
                   grepl("Chi-squared test for given probabilities", m)) "chi_squared_test"
          else if (grepl("Fisher.*Exact Test",           m)) "fisher_exact_test"
          else                                                "htest_other"

  list(type = type, type_label = m)
}

# ── Field extraction ──────────────────────────────────────────────────

# Extract the schema's InlineFields object from an htest result.
# `data` (P4) is an optional data.frame. `env` is a conservative fallback
# for the common script flow where the user stores a t.test() object and then
# calls mellio_open(result): htest objects retain only a data.name string, so
# we recover simple names like score, d$score, or d[["score"]] when they are
# still available in the caller's environment.
ms_htest_fields <- function(x, type, data = NULL, env = NULL) {
  fields <- list(p_value = ms_safe_numeric(x$p.value))

  # Statistic is optional: Fisher's exact has no test statistic, so
  # ms_htest_statistic returns NULL and the field is omitted entirely
  # rather than emitted as {name, value: NA}.
  stat <- ms_htest_statistic(x, type)
  if (!is.null(stat)) fields$statistic <- stat

  est <- ms_htest_estimate(x, type)
  if (!is.null(est)) fields$estimate <- est

  # P2: group level labels for two-sample t-tests. R's $estimate is
  # c("mean in group X" = ..., "mean in group Y" = ...), so the level
  # names live in names(est). Labels alone (no n/M/SD) ship without
  # .data; the P4 enrichment below adds n/mean/sd when data is supplied.
  groups <- ms_htest_groups(x, type)
  if (!is.null(groups)) fields$groups <- groups

  context <- ms_htest_context_fields(x, type, fields$groups)
  if (length(context)) {
    for (key in names(context)) fields[[key]] <- context[[key]]
  }

  # One-sample t-test: surface mu (htest$null.value) so the paragraph
  # builder can render "...differed from <mu>." cleanly.
  if (identical(type, "one_sample_t_test") && !is.null(x$null.value)) {
    fields$null_value <- ms_safe_numeric(unname(x$null.value)[[1]])
  }

  # ── P4 enrichment from .data / simple caller-env references ───────
  if (!is.null(data) && type %in% ms_htest_data_enrichment_types()) {
    enriched <- tryCatch(
      ms_htest_enrich_from_data(x, type, data, existing_groups = fields$groups),
      error = function(e) NULL
    )
    if (!is.null(enriched)) {
      if (!is.null(enriched$groups))          fields$groups          <- enriched$groups
      if (!is.null(enriched$paired_measures)) fields$paired_measures <- enriched$paired_measures
      if (!is.null(enriched$effect_size))     fields$effect_size     <- enriched$effect_size
      if (!is.null(enriched$grouping_note))   fields$grouping_note   <- enriched$grouping_note
    }
  }
  if (is.null(fields$effect_size) &&
      type %in% ms_htest_data_enrichment_types() &&
      is.environment(env)) {
    enriched <- tryCatch(
      ms_htest_enrich_from_env(x, type, env, existing_groups = fields$groups),
      error = function(e) NULL
    )
    if (!is.null(enriched)) {
      if (!is.null(enriched$groups))          fields$groups          <- enriched$groups
      if (!is.null(enriched$paired_measures)) fields$paired_measures <- enriched$paired_measures
      if (!is.null(enriched$effect_size))     fields$effect_size     <- enriched$effect_size
      if (!is.null(enriched$grouping_note))   fields$grouping_note   <- enriched$grouping_note
    }
  }

  # P2: chi-square gets sample size + Cramer's V; Fisher gets the
  # odds-ratio + CI from $estimate above. Cramer's V is computed from
  # the contingency table dimensions in $observed.
  if (identical(type, "chi_squared_test") && !is.null(x$observed)) {
    n <- ms_safe_numeric(sum(x$observed))
    if (!is.na(n) && n > 0) fields$sample_size <- n
    context <- ms_htest_chi_square_context(x)
    if (nzchar(context)) fields$test_context <- context
    count_table <- ms_htest_chi_square_count_table(x, context)
    if (!is.null(count_table)) {
      for (key in names(count_table)) fields[[key]] <- count_table[[key]]
    }

    v <- ms_cramers_v(x)
    if (!is.null(v)) {
      fields$effect <- list(name = "cramers_v", value = v)
    }
  } else if (identical(type, "fisher_exact_test")) {
    fields$test_context <- "association"
  }

  if (!is.null(x$alternative)) {
    fields$alternative <- as.character(x$alternative)
  }

  if (!is.null(x$conf.int)) {
    cl <- attr(x$conf.int, "conf.level")
    if (!is.null(cl)) fields$conf_level <- ms_safe_numeric(cl)
  }

  method_meta <- ms_htest_method_metadata(x, type)
  if (!is.null(method_meta$estimate_ci_method)) {
    fields$estimate_ci_method <- method_meta$estimate_ci_method
  }
  if (!is.null(method_meta$p_value_method)) {
    fields$p_value_method <- method_meta$p_value_method
  }

  if (is.null(fields$sample_size)) {
    n <- ms_htest_sample_size(x, type, fields)
    if (!is.null(n) && !is.na(n) && is.finite(n) && n > 0) {
      fields$sample_size <- ms_safe_numeric(n)
    }
  }
  if (is.null(fields$sample_size) &&
      type %in% c("pearson_correlation", "spearman_correlation",
                  "kendall_correlation")) {
    n <- tryCatch(
      ms_htest_correlation_sample_size(x, data = data, env = env),
      error = function(e) NULL
    )
    if (!is.null(n) && !is.na(n) && is.finite(n) && n > 0) {
      fields$sample_size <- ms_safe_numeric(n)
    }
  }

  fields
}

ms_htest_context_fields <- function(x, type, groups = NULL) {
  parsed <- ms_htest_parse_data_name(x$data.name)
  if (is.null(parsed)) return(list())

  out <- list()
  clean <- function(value) {
    label <- ms_model_clean_term(value)
    if (nzchar(label)) label else NULL
  }

  if (type %in% c("welch_t_test", "students_t_test") &&
      identical(parsed$form, "formula")) {
    outcome <- clean(parsed$outcome)
    group <- clean(parsed$by)
    if (!is.null(outcome)) out$outcome <- outcome
    if (!is.null(group)) {
      out$predictor <- group
      out$group <- group
    }
    group_labels <- ms_htest_group_labels(groups)
    if (length(group_labels) == 2L) {
      out$comparison <- paste(group_labels[[1L]], "vs.", group_labels[[2L]])
    }
  } else if (identical(type, "paired_t_test") &&
             identical(parsed$form, "pair")) {
    labels <- Filter(Negate(is.null), lapply(c(parsed$left, parsed$right), clean))
    if (length(labels) == 2L) out$comparison <- paste(labels[[1L]], "vs.", labels[[2L]])
  } else if (identical(type, "one_sample_t_test") &&
             identical(parsed$form, "single")) {
    outcome <- clean(parsed$name)
    if (!is.null(outcome)) out$outcome <- outcome
  } else if (identical(type, "kruskal_wallis_test") &&
             identical(parsed$form, "formula")) {
    outcome <- clean(parsed$outcome)
    group <- clean(parsed$by)
    if (!is.null(outcome)) out$outcome <- outcome
    if (!is.null(group)) {
      out$predictor <- group
      out$group <- group
    }
  } else if (identical(type, "wilcoxon_rank_sum") &&
             identical(parsed$form, "formula")) {
    outcome <- clean(parsed$outcome)
    group <- clean(parsed$by)
    if (!is.null(outcome)) out$outcome <- outcome
    if (!is.null(group)) {
      out$predictor <- group
      out$group <- group
    }
  } else if (identical(type, "wilcoxon_signed_rank") &&
             identical(parsed$form, "pair")) {
    labels <- Filter(Negate(is.null), lapply(c(parsed$left, parsed$right), clean))
    if (length(labels) == 2L) out$comparison <- paste(labels[[1L]], "vs.", labels[[2L]])
  } else if (identical(type, "wilcoxon_signed_rank") &&
             identical(parsed$form, "single")) {
    outcome <- clean(parsed$name)
    if (!is.null(outcome)) out$outcome <- outcome
  } else if (type %in% c("pearson_correlation", "spearman_correlation",
                         "kendall_correlation") &&
             identical(parsed$form, "pair")) {
    labels <- Filter(Negate(is.null), lapply(c(parsed$left, parsed$right), clean))
    if (length(labels) == 2L) out$comparison <- paste(labels[[1L]], "and", labels[[2L]])
  }

  out
}

ms_htest_group_labels <- function(groups) {
  if (!is.list(groups) || !length(groups)) return(character(0))
  labels <- vapply(groups, function(g) {
    as.character(g$label %||% NA_character_)
  }, character(1))
  labels <- labels[!is.na(labels) & nzchar(labels)]
  labels
}

ms_htest_data_enrichment_types <- function() {
  c("welch_t_test", "students_t_test", "paired_t_test", "one_sample_t_test",
    "wilcoxon_rank_sum", "wilcoxon_signed_rank", "kruskal_wallis_test")
}

ms_htest_method_metadata <- function(x, type) {
  out <- list()
  has_ci <- !is.null(x$conf.int) && length(x$conf.int) >= 2L

  if (type %in% c("welch_t_test", "students_t_test",
                  "paired_t_test", "one_sample_t_test")) {
    out$p_value_method <- "t_distribution"
    if (has_ci) out$estimate_ci_method <- "t_interval"
  } else if (identical(type, "pearson_correlation")) {
    out$p_value_method <- "t_distribution"
    if (has_ci) out$estimate_ci_method <- "fisher_z"
  } else if (identical(type, "spearman_correlation")) {
    out$p_value_method <- "spearman_rank_test"
  } else if (identical(type, "kendall_correlation")) {
    out$p_value_method <- "kendall_rank_test"
  } else if (type %in% c("wilcoxon_rank_sum", "wilcoxon_signed_rank")) {
    out$p_value_method <- "wilcoxon_rank_test"
  } else if (type %in% c("kruskal_wallis_test", "friedman_test", "mcnemar_test")) {
    out$p_value_method <- "chi_square_approximation"
  } else if (identical(type, "chi_squared_test")) {
    method <- as.character(x$method %||% "")
    out$p_value_method <- if (grepl("simulated p-value", method, ignore.case = TRUE)) {
      "monte_carlo_simulation"
    } else {
      "chi_square_approximation"
    }
  } else if (identical(type, "fisher_exact_test")) {
    out$p_value_method <- "fisher_exact"
    if (has_ci) out$estimate_ci_method <- "exact"
  }

  out
}

ms_htest_chi_square_context <- function(x) {
  obs <- x$observed
  if (is.null(obs)) return("")
  dims <- dim(obs)
  if (!is.null(dims) && length(dims) == 2L && all(dims >= 2L)) {
    return("association")
  }
  "goodness_of_fit"
}

ms_htest_chi_square_count_table <- function(x, context) {
  obs <- x$observed
  if (is.null(obs)) return(NULL)
  context <- as.character(context %||% "")

  if (identical(context, "association")) {
    dims <- dim(obs)
    if (is.null(dims) || length(dims) != 2L || any(dims < 1L)) return(NULL)
    rn <- rownames(obs)
    cn <- colnames(obs)
    if (is.null(rn) || length(rn) != dims[[1L]] || any(!nzchar(rn))) {
      rn <- paste("Row", seq_len(dims[[1L]]))
    }
    if (is.null(cn) || length(cn) != dims[[2L]] || any(!nzchar(cn))) {
      cn <- paste("Column", seq_len(dims[[2L]]))
    }

    columns <- c(
      list(list(key = "row", label = "Row", format = "text")),
      lapply(seq_len(dims[[2L]]), function(j) {
        list(key = paste0("col_", j), label = cn[[j]], format = "integer")
      })
    )
    rows <- lapply(seq_len(dims[[1L]]), function(i) {
      row <- list(row = rn[[i]])
      for (j in seq_len(dims[[2L]])) {
        row[[paste0("col_", j)]] <- ms_safe_numeric(obs[i, j])
      }
      row
    })
    n <- ms_safe_numeric(sum(obs))
    return(list(
      table_type = "chi_square_contingency",
      columns = columns,
      rows = rows,
      table_note = paste0(
        "Cells show observed counts for the contingency table. ",
        if (!is.na(n) && is.finite(n) && n > 0) paste0("N = ", round(n), ". ") else "",
        "The chi-square p value is based on the chi-square approximation."
      )
    ))
  }

  values <- as.numeric(obs)
  if (!length(values)) return(NULL)
  expected <- if (!is.null(x$expected)) as.numeric(x$expected) else rep(NA_real_, length(values))
  labels <- names(obs)
  if (is.null(labels) || length(labels) != length(values) || any(!nzchar(labels))) {
    labels <- paste("Category", seq_along(values))
  }
  rows <- lapply(seq_along(values), function(i) {
    row <- list(
      category = labels[[i]],
      observed = ms_safe_numeric(values[[i]])
    )
    exp_value <- ms_safe_numeric(expected[[i]])
    if (!is.na(exp_value) && is.finite(exp_value)) row$expected <- exp_value
    row
  })
  n <- ms_safe_numeric(sum(values))
  list(
    table_type = "chi_square_goodness_of_fit",
    columns = list(
      list(key = "category", label = "Category", format = "text"),
      list(key = "observed", label = "Observed", format = "integer"),
      list(key = "expected", label = "Expected", format = "number")
    ),
    rows = rows,
    table_note = paste0(
      "Observed and expected counts are shown. ",
      if (!is.na(n) && is.finite(n) && n > 0) paste0("N = ", round(n), ". ") else "",
      "The chi-square p value is based on the chi-square approximation."
    )
  )
}

ms_htest_sample_size <- function(x, type, fields = list()) {
  groups <- fields$groups
  if (is.list(groups) && length(groups) > 0L) {
    ns <- vapply(groups, function(g) {
      ms_safe_numeric(g$n %||% NA_real_)
    }, numeric(1))
    ns <- ns[!is.na(ns) & is.finite(ns) & ns > 0]
    if (length(ns) > 0L) return(sum(ns))
  }

  if (identical(type, "chi_squared_test") && !is.null(x$observed)) {
    n <- ms_safe_numeric(sum(x$observed))
    if (!is.na(n) && is.finite(n) && n > 0) return(n)
  }

  df <- if (!is.null(x$parameter) && length(x$parameter) >= 1L) {
    ms_safe_numeric(unname(x$parameter)[[1L]])
  } else {
    NA_real_
  }
  if (is.na(df) || !is.finite(df)) return(NULL)

  if (type %in% c("paired_t_test", "one_sample_t_test")) {
    return(df + 1)
  }
  if (identical(type, "pearson_correlation")) {
    return(df + 2)
  }

  NULL
}

ms_htest_correlation_sample_size <- function(x, data = NULL, env = NULL) {
  parsed <- ms_htest_parse_data_name(x$data.name)
  if (is.null(parsed) || !identical(parsed$form, "pair")) return(NULL)

  left <- right <- NULL
  if (!is.null(data) && is.data.frame(data)) {
    left <- ms_htest_resolve_var(parsed$left, data)
    right <- ms_htest_resolve_var(parsed$right, data)
  }
  if ((is.null(left) || is.null(right)) && is.environment(env)) {
    left <- ms_htest_resolve_var_env(parsed$left, env)
    right <- ms_htest_resolve_var_env(parsed$right, env)
  }

  if (is.null(left) || is.null(right) || length(left) != length(right)) return(NULL)
  keep <- stats::complete.cases(left, right)
  n <- sum(keep)
  if (!is.na(n) && is.finite(n) && n > 0) n else NULL
}

ms_htest_figure_data <- function(x, type, fields, data = NULL, env = NULL) {
  parsed <- ms_htest_parse_data_name(x$data.name)
  if (identical(type, "paired_t_test")) {
    source_data <- data
    if (is.null(source_data) && is.environment(env)) {
      source_data <- tryCatch(ms_htest_env_data(parsed, type, env, x = x),
                              error = function(e) NULL)
    }
    paired_plot <- ms_htest_paired_figure_data(x, fields, source_data, parsed)
    if (!is.null(paired_plot)) return(list(paired_difference_plot = paired_plot))
    return(NULL)
  }
  if (identical(type, "one_sample_t_test")) {
    one_sample_plot <- ms_htest_one_sample_figure_data(x, fields, parsed)
    if (!is.null(one_sample_plot)) return(list(one_sample_mean_plot = one_sample_plot))
    return(NULL)
  }

  if (type %in% c("wilcoxon_rank_sum", "kruskal_wallis_test")) {
    source_data <- data
    if (is.null(source_data) && is.environment(env)) {
      source_data <- tryCatch(ms_htest_env_data(parsed, type, env, x = x),
                              error = function(e) NULL)
    }
    group_plot <- ms_htest_nonparametric_group_figure_data(
      x, type, fields, source_data, parsed
    )
    if (!is.null(group_plot)) return(list(nonparametric_group_plot = group_plot))
    return(NULL)
  }

  if (type %in% c("pearson_correlation", "spearman_correlation",
                  "kendall_correlation")) {
    scatter <- ms_htest_scatter_figure_data(x, type, fields, data, parsed, env = env)
    if (!is.null(scatter)) return(list(scatter_plot = scatter))
    return(NULL)
  }

  if (is.null(data) || !(type %in% c("welch_t_test", "students_t_test"))) return(NULL)
  groups <- fields$groups
  if (!is.list(groups) || length(groups) < 2L || length(groups) > 12L) return(NULL)
  usable <- lapply(groups, function(g) {
    n <- ms_safe_numeric(g$n %||% NA_real_)
    mean_value <- ms_safe_numeric(g$mean %||% NA_real_)
    sd_value <- ms_safe_numeric(g$sd %||% NA_real_)
    label <- as.character(g$label %||% "")
    if (!nzchar(label) || is.na(n) || is.na(mean_value) || is.na(sd_value) ||
        !is.finite(n) || !is.finite(mean_value) || !is.finite(sd_value) ||
        n < 2L || sd_value < 0) {
      return(NULL)
    }
    se_value <- sd_value / sqrt(n)
    list(
      level = label,
      label = label,
      n = n,
      mean = mean_value,
      sd = sd_value,
      se = ms_safe_numeric(se_value),
      df = n - 1
    )
  })
  usable <- Filter(Negate(is.null), usable)
  if (length(usable) < 2L) return(NULL)

  if (is.null(parsed) || !identical(parsed$form, "formula")) return(NULL)
  outcome <- as.character(fields$outcome %||% ms_model_clean_term(parsed$outcome))
  group <- as.character(fields$group %||% fields$predictor %||% ms_model_clean_term(parsed$by))
  if (!nzchar(outcome) || !nzchar(group)) return(NULL)

  level_rows <- lapply(usable, function(g) {
    list(value = g$level, label = g$label)
  })
  adjusted_means <- list(
    mean_kind = "observed",
    source = "htest_two_sample",
    factor = list(
      variable = group,
      term = parsed$by,
      label = group,
      levels = level_rows
    ),
    groups = usable,
    outcome = outcome,
    y_label = paste("Mean", outcome),
    ci_level = ms_safe_numeric(fields$conf_level %||% 0.95),
    ci_method = "group_standard_error"
  )

  figures <- list(adjusted_means = adjusted_means)
  estimation_plot <- ms_htest_two_sample_estimation_figure_data(
    fields, data, parsed, usable, outcome, group
  )
  if (!is.null(estimation_plot)) figures$estimation_plot <- estimation_plot
  figures
}

# Gardner-Altman estimation plot for two-group t-tests: raw observations +
# group means on a shared outcome axis, with the mean difference and its CI
# on a floating delta axis. Two-group only; needs the source data frame for
# the raw dots. Returns NULL (caller still emits the means-plot fallback) when
# data is unavailable or the difference can't be formed.
ms_htest_two_sample_estimation_figure_data <- function(fields, data, parsed,
                                                       usable, outcome, group) {
  if (length(usable) != 2L) return(NULL)
  if (is.null(data) || !is.data.frame(data)) return(NULL)
  if (is.null(parsed) || !identical(parsed$form, "formula")) return(NULL)

  y <- ms_htest_resolve_var(parsed$outcome, data)
  g <- ms_htest_resolve_var(parsed$by, data)
  if (is.null(y) || is.null(g) || length(y) != length(g)) return(NULL)
  keep <- is.finite(y) & !is.na(g)
  y <- y[keep]
  g <- g[keep]
  if (length(y) < 2L) return(NULL)
  g_chr <- as.character(g)

  ref <- usable[[1L]]
  comparison <- usable[[2L]]
  diff_value <- ms_safe_numeric(comparison$mean - ref$mean)
  if (is.na(diff_value) || !is.finite(diff_value)) return(NULL)

  ci <- fields$estimate$ci %||% NULL
  if (!is.null(ci)) ci <- ms_safe_numeric(unname(ci))
  difference <- list(
    estimate   = diff_value,
    reference  = ref$label,
    comparison = comparison$label,
    conf_level = ms_safe_numeric(fields$conf_level %||% 0.95)
  )
  # R's conf.int is the symmetric CI of (mean1 - mean2); reuse its half-width
  # around our reference->comparison difference so direction stays consistent
  # with the plotted groups regardless of how R ordered $estimate.
  if (!is.null(ci) && length(ci) == 2L && all(is.finite(ci))) {
    half_width <- abs(ci[[2L]] - ci[[1L]]) / 2
    difference$ci_lower <- ms_safe_numeric(diff_value - half_width)
    difference$ci_upper <- ms_safe_numeric(diff_value + half_width)
  }

  observations <- lapply(seq_along(y), function(i) {
    list(id = i, group = g_chr[[i]], value = ms_safe_numeric(y[[i]]))
  })
  max_points <- 500L
  truncated <- length(observations) > max_points
  if (truncated) observations <- observations[seq_len(max_points)]

  level_rows <- lapply(usable, function(gr) {
    list(value = gr$level, label = gr$label)
  })

  list(
    source = "htest_two_sample",
    plot_kind = "two_group_estimation",
    factor = list(
      variable = group,
      term = parsed$by,
      label = group,
      levels = level_rows
    ),
    groups = usable,
    observations = observations,
    point_display = list(
      total_n = length(y),
      included_n = length(observations),
      truncated = truncated
    ),
    difference = difference,
    outcome = outcome,
    y_label = outcome
  )
}

# Scatter plot for correlation tests: the raw (x, y) cloud with axis labels and
# the correlation estimate. Needs the source data frame for the points (cor.test
# stores no observations). For Pearson it adds an OLS fit line + CI band; rank
# correlations (Spearman/Kendall) get points only, since a least-squares line
# would misrepresent a monotonic-rank relationship.
ms_htest_scatter_figure_data <- function(x, type, fields, data, parsed = NULL, env = NULL) {
  if (is.null(parsed)) parsed <- ms_htest_parse_data_name(x$data.name)
  if (is.null(parsed) || !identical(parsed$form, "pair")) return(NULL)

  xv <- yv <- NULL
  if (!is.null(data) && is.data.frame(data)) {
    xv <- ms_htest_resolve_var(parsed$left, data)
    yv <- ms_htest_resolve_var(parsed$right, data)
  }
  if ((is.null(xv) || is.null(yv)) && is.environment(env)) {
    xv <- ms_htest_resolve_var_env(parsed$left, env)
    yv <- ms_htest_resolve_var_env(parsed$right, env)
  }
  if (is.null(xv) || is.null(yv)) return(NULL)
  xv <- suppressWarnings(as.numeric(xv))
  yv <- suppressWarnings(as.numeric(yv))
  if (length(xv) != length(yv)) return(NULL)

  keep <- is.finite(xv) & is.finite(yv)
  total_n <- sum(keep)
  xv <- xv[keep]
  yv <- yv[keep]
  if (length(xv) < 3L) return(NULL)

  x_label <- ms_htest_reference_name(parsed$left)
  y_label <- ms_htest_reference_name(parsed$right)
  if (!nzchar(x_label)) x_label <- "x"
  if (!nzchar(y_label)) y_label <- "y"

  max_points <- 1000L
  idx <- seq_along(xv)
  truncated <- length(idx) > max_points
  if (truncated) idx <- idx[seq_len(max_points)]
  observations <- lapply(idx, function(i) {
    list(x = ms_safe_numeric(xv[[i]]), y = ms_safe_numeric(yv[[i]]))
  })

  est <- fields$estimate %||% list()
  est_ci <- est$ci %||% NULL
  if (!is.null(est_ci)) est_ci <- ms_safe_numeric(unname(est_ci))

  out <- list(
    source = "htest_correlation",
    plot_kind = "scatter",
    method = sub("_correlation$", "", type),
    x_label = x_label,
    y_label = y_label,
    observations = observations,
    point_display = list(
      total_n = total_n,
      included_n = length(observations),
      truncated = truncated
    ),
    estimate = list(
      name = as.character(est$name %||% "cor"),
      value = ms_safe_numeric(est$value %||% NA_real_),
      ci = est_ci
    ),
    p_value = ms_safe_numeric(fields$p_value %||% NA_real_),
    conf_level = ms_safe_numeric(fields$conf_level %||% 0.95),
    n = total_n
  )

  if (identical(type, "pearson_correlation") &&
      stats::sd(xv) > 0 && length(unique(xv)) >= 2L) {
    fit <- ms_htest_scatter_fit(xv, yv, conf_level = out$conf_level)
    if (!is.null(fit)) out$fit <- fit
  }
  out
}

# OLS fit line + pointwise confidence band over the observed x-range, sampled at
# n_points for a smooth polyline. Returns NULL when the model can't be formed.
ms_htest_scatter_fit <- function(xv, yv, conf_level = 0.95, n_points = 64L) {
  df <- data.frame(fx = as.numeric(xv), fy = as.numeric(yv))
  fit <- tryCatch(stats::lm(fy ~ fx, data = df), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  rng <- range(df$fx, finite = TRUE)
  if (!all(is.finite(rng)) || rng[[1L]] == rng[[2L]]) return(NULL)
  xs <- seq(rng[[1L]], rng[[2L]], length.out = n_points)
  pr <- tryCatch(
    stats::predict(fit, newdata = data.frame(fx = xs),
                   interval = "confidence", level = conf_level),
    error = function(e) NULL
  )
  if (is.null(pr)) return(NULL)
  co <- stats::coef(fit)
  line <- lapply(seq_along(xs), function(i) {
    list(
      x = ms_safe_numeric(xs[[i]]),
      y = ms_safe_numeric(pr[i, "fit"]),
      lower = ms_safe_numeric(pr[i, "lwr"]),
      upper = ms_safe_numeric(pr[i, "upr"])
    )
  })
  list(
    kind = "ols",
    slope = ms_safe_numeric(unname(co[[2L]])),
    intercept = ms_safe_numeric(unname(co[[1L]])),
    conf_level = ms_safe_numeric(conf_level),
    line = line
  )
}

ms_htest_nonparametric_group_figure_data <- function(x, type, fields, data, parsed = NULL) {
  if (is.null(parsed)) parsed <- ms_htest_parse_data_name(x$data.name)
  if (is.null(parsed) || !identical(parsed$form, "formula")) return(NULL)
  if (is.null(data) || !is.data.frame(data)) return(NULL)

  y <- ms_htest_resolve_var(parsed$outcome, data)
  g <- ms_htest_resolve_var(parsed$by, data)
  if (is.null(y) || is.null(g) || length(y) != length(g)) return(NULL)

  keep <- is.finite(y) & !is.na(g)
  y <- y[keep]
  g <- g[keep]
  if (length(y) < 2L) return(NULL)

  g_factor <- factor(g)
  levels <- levels(g_factor)
  if (identical(type, "wilcoxon_rank_sum") && length(levels) != 2L) return(NULL)
  if (identical(type, "kruskal_wallis_test") && length(levels) < 2L) return(NULL)

  groups <- ms_htest_group_descriptives(y, g_factor, levels = levels)
  if (is.null(groups) || length(groups) < 2L) return(NULL)

  observations <- lapply(seq_along(y), function(i) {
    list(
      id = i,
      group = as.character(g_factor[[i]]),
      value = ms_safe_numeric(y[[i]])
    )
  })
  max_points <- 500L
  truncated <- length(observations) > max_points
  if (truncated) observations <- observations[seq_len(max_points)]

  level_rows <- lapply(groups, function(row) {
    list(value = row$label, label = row$label)
  })
  outcome <- as.character(fields$outcome %||% ms_model_clean_term(parsed$outcome))
  group <- as.character(fields$group %||% fields$predictor %||% ms_model_clean_term(parsed$by))
  if (!nzchar(outcome)) outcome <- ms_htest_reference_name(parsed$outcome)
  if (!nzchar(group)) group <- ms_htest_reference_name(parsed$by)
  if (!nzchar(outcome) || !nzchar(group)) return(NULL)

  list(
    source = type,
    plot_kind = "independent_groups",
    factor = list(
      variable = group,
      term = parsed$by,
      label = group,
      levels = level_rows
    ),
    groups = groups,
    observations = observations,
    point_display = list(
      total_n = length(y),
      included_n = length(observations),
      truncated = truncated
    ),
    outcome = outcome,
    y_label = outcome
  )
}

ms_htest_one_sample_figure_data <- function(x, fields, parsed = NULL) {
  estimate <- fields$estimate %||% list()
  mean_value <- ms_safe_numeric(estimate$value %||%
                                  unname(x$estimate %||% NA_real_)[[1L]])
  ci <- estimate$ci %||% x$conf.int %||% NULL
  if (!is.null(ci)) ci <- ms_safe_numeric(unname(ci))
  if (is.na(mean_value) || !is.finite(mean_value) ||
      is.null(ci) || length(ci) < 2L ||
      any(is.na(ci[1:2])) || any(!is.finite(ci[1:2]))) {
    return(NULL)
  }

  sample_row <- if (is.list(fields$groups) && length(fields$groups)) {
    fields$groups[[1L]]
  } else {
    list()
  }
  label <- as.character(
    fields$outcome %||%
      sample_row$label %||%
      ms_htest_reference_name(parsed$name %||% "") %||%
      ""
  )
  if (!nzchar(label)) label <- "Sample"

  n <- ms_safe_numeric(sample_row$n %||% fields$sample_size %||% NA_real_)
  sd_value <- ms_safe_numeric(sample_row$sd %||% NA_real_)
  se_value <- if (!is.na(sd_value) && is.finite(sd_value) &&
                  !is.na(n) && is.finite(n) && n > 0) {
    ms_safe_numeric(sd_value / sqrt(n))
  } else {
    NA_real_
  }
  null_value <- ms_safe_numeric(fields$null_value %||%
                                  unname(x$null.value %||% 0)[[1L]])
  if (is.na(null_value) || !is.finite(null_value)) null_value <- 0
  conf_level <- ms_safe_numeric(fields$conf_level %||%
                                  attr(x$conf.int, "conf.level") %||% 0.95)

  list(
    source = "one_sample_t_test",
    sample = list(label = label, variable = parsed$name %||% ""),
    n = if (!is.na(n) && is.finite(n)) n else NULL,
    mean = mean_value,
    sd = if (!is.na(sd_value) && is.finite(sd_value)) sd_value else NULL,
    se = if (!is.na(se_value) && is.finite(se_value)) se_value else NULL,
    null_value = null_value,
    ci_lower = ci[[1L]],
    ci_upper = ci[[2L]],
    conf_level = if (!is.na(conf_level) && is.finite(conf_level)) conf_level else 0.95,
    outcome = label,
    y_label = label
  )
}

ms_htest_paired_figure_data <- function(x, fields, data, parsed = NULL) {
  if (is.null(parsed)) parsed <- ms_htest_parse_data_name(x$data.name)
  if (is.null(parsed) || !identical(parsed$form, "pair")) return(NULL)
  if (is.null(data) || !is.data.frame(data)) return(NULL)

  left <- ms_htest_resolve_var(parsed$left, data)
  right <- ms_htest_resolve_var(parsed$right, data)
  if (is.null(left) || is.null(right) || length(left) != length(right)) return(NULL)
  keep <- is.finite(left) & is.finite(right)
  left <- left[keep]
  right <- right[keep]
  n <- length(left)
  if (n < 2L) return(NULL)

  left_label <- ms_htest_reference_name(parsed$left)
  right_label <- ms_htest_reference_name(parsed$right)
  if (!nzchar(left_label)) left_label <- "Measure 1"
  if (!nzchar(right_label)) right_label <- "Measure 2"

  diffs <- left - right
  pairs <- lapply(seq_len(n), function(i) {
    list(
      id = i,
      first = ms_safe_numeric(left[[i]]),
      second = ms_safe_numeric(right[[i]]),
      difference = ms_safe_numeric(diffs[[i]])
    )
  })

  stat <- fields$statistic %||% list()
  estimate <- fields$estimate %||% list()
  ci <- estimate$ci %||% NULL
  if (!is.null(ci)) ci <- ms_safe_numeric(unname(ci))
  conf_level <- ms_safe_numeric(fields$conf_level %||% attr(x$conf.int, "conf.level") %||% 0.95)

  list(
    source = "paired_t_test",
    mean_kind = "within_subject",
    first_measure = list(label = left_label, variable = parsed$left),
    second_measure = list(label = right_label, variable = parsed$right),
    pairs = pairs,
    n = n,
    outcome = paste(left_label, "and", right_label),
    y_label = "Observed value",
    direction = paste(left_label, "minus", right_label),
    difference = list(
      mean = ms_safe_numeric(mean(diffs)),
      sd = ms_safe_numeric(stats::sd(diffs)),
      se = ms_safe_numeric(stats::sd(diffs) / sqrt(n)),
      df = ms_safe_numeric(stat$df %||% (n - 1L)),
      ci_lower = if (!is.null(ci) && length(ci) == 2L) ci[[1L]] else NULL,
      ci_upper = if (!is.null(ci) && length(ci) == 2L) ci[[2L]] else NULL,
      conf_level = conf_level
    )
  )
}

# Extract per-group level labels and retained group means for two-sample
# t-tests. Returns NULL for any other test type. n/sd are added later
# (P4) when raw data is available.
ms_htest_groups <- function(x, type) {
  if (!type %in% c("welch_t_test", "students_t_test")) return(NULL)
  est <- x$estimate
  if (is.null(est) || length(est) != 2L) return(NULL)
  nms <- names(est)
  if (is.null(nms) || length(nms) != 2L) return(NULL)
  # Only formula-call t-tests have names of the form "mean in group <X>".
  # Vector-input calls (t.test(x, y)) get names like "mean of x"/"mean of y",
  # which aren't real group labels — return NULL so the paragraph falls
  # back to placeholders rather than leaking "mean of x".
  if (!all(grepl("^mean in group\\s+", nms))) return(NULL)
  labels <- sub("^mean in group\\s+", "", nms)
  if (any(!nzchar(labels))) return(NULL)
  means <- ms_safe_numeric(unname(est))
  unname(Map(function(lbl, mean_value) {
    row <- list(label = lbl)
    if (!is.na(mean_value) && is.finite(mean_value)) row$mean <- mean_value
    row
  }, labels, means))
}

# Cramer's V for a chisq.test() result. Returns NULL when the result
# doesn't have a usable observed table (e.g. one-dimensional GoF tests
# where V isn't defined).
ms_cramers_v <- function(x) {
  obs <- x$observed
  chi <- ms_safe_numeric(unname(x$statistic))
  n   <- if (!is.null(obs)) sum(obs) else NA_real_
  if (is.na(chi) || is.na(n) || n <= 0) return(NULL)

  dims <- dim(obs)
  if (is.null(dims) || length(dims) != 2L || any(dims < 2L)) return(NULL)

  k <- min(dims) - 1L
  if (k < 1L) return(NULL)
  v <- sqrt(chi / (n * k))
  if (is.na(v) || !is.finite(v)) return(NULL)
  ms_safe_numeric(v)
}

# Build the statistic block: name + value + optional df. Returns NULL
# when x$statistic is missing (e.g. fisher.test, which has no test
# statistic — only p-value, odds ratio, and CI). Callers must omit the
# field entirely rather than emit {value: NA}; the JS validator rejects
# inline cards with NA statistics.
ms_htest_statistic <- function(x, type = NULL) {
  if (is.null(x$statistic)) return(NULL)
  stat_val  <- ms_safe_numeric(unname(x$statistic))
  stat_name <- names(x$statistic)
  if (is.null(stat_name) || !nzchar(stat_name)) stat_name <- "stat"
  if (identical(type, "kruskal_wallis_test")) stat_name <- "H"

  out <- list(name = stat_name, value = stat_val)

  # df: x$parameter is named numeric. Single value common (t-test, chi-sq);
  # two values for F-tests (not produced by htest path here). We still
  # serialise as a vector when length == 2.
  if (!is.null(x$parameter) && length(x$parameter) > 0) {
    df_vals <- ms_safe_numeric(unname(x$parameter))
    if (length(df_vals) == 1L) {
      out$df <- df_vals
    } else {
      # Force JSON serialisation as an array, not a scalar.
      out$df <- I(df_vals)
    }
  }
  out
}

# Build the estimate block: name + value + optional CI.
# Test-type-specific because the meaning of $estimate varies.
ms_htest_estimate <- function(x, type) {
  if (is.null(x$estimate)) return(NULL)

  est <- x$estimate
  ci  <- if (!is.null(x$conf.int)) ms_safe_numeric(unname(x$conf.int)) else NULL

  # Two-sample t-tests: $estimate is c(mean of x, mean of y).
  # We emit the mean difference rather than two separate means.
  if (type %in% c("welch_t_test", "students_t_test") && length(est) == 2L) {
    return(list(
      name  = "mean_diff",
      value = ms_safe_numeric(unname(est[1] - est[2])),
      ci    = if (!is.null(ci) && length(ci) == 2L) I(ci) else NULL
    ))
  }

  # Paired and one-sample t-tests: $estimate is the mean of differences
  # or the sample mean.
  if (type %in% c("paired_t_test", "one_sample_t_test")) {
    nm <- names(est)[1] %||% "mean"
    return(list(
      name  = nm,
      value = ms_safe_numeric(unname(est[1])),
      ci    = if (!is.null(ci) && length(ci) == 2L) I(ci) else NULL
    ))
  }

  # Correlation: $estimate is c(cor = r) or similar.
  if (type %in% c("pearson_correlation", "spearman_correlation",
                  "kendall_correlation")) {
    nm <- names(est)[1] %||% "estimate"
    return(list(
      name  = nm,
      value = ms_safe_numeric(unname(est[1])),
      ci    = if (!is.null(ci) && length(ci) == 2L) I(ci) else NULL
    ))
  }

  # Wilcoxon with conf.int=TRUE produces a pseudo-median estimate.
  # chisq/fisher have an $estimate for the odds ratio (fisher) or no estimate.
  if (identical(type, "fisher_exact_test")) {
    return(list(
      name  = "odds_ratio",
      value = ms_safe_numeric(unname(est[1])),
      ci    = if (!is.null(ci) && length(ci) == 2L) I(ci) else NULL
    ))
  }

  nm <- names(est)[1] %||% "estimate"
  list(
    name  = nm,
    value = ms_safe_numeric(unname(est[1])),
    ci    = if (!is.null(ci) && length(ci) == 2L) I(ci) else NULL
  )
}

# ── P4: .data enrichment for t-tests ────────────────────────────────
#
# Parse the htest$data.name string into a structured spec describing the
# call shape so we know how to look up variables in `.data`:
#   "y by g"   → formula two-sample / paired-from-formula
#   "x and y"  → two-vector call (paired or two-sample-from-vectors)
#   "x"        → one-sample
# Returns NULL when the pattern is unrecognised.
ms_htest_parse_data_name <- function(data_name) {
  s <- trimws(as.character(data_name %||% ""))
  if (!nzchar(s)) return(NULL)
  if (grepl("\\s+by\\s+", s)) {
    parts <- strsplit(s, "\\s+by\\s+", perl = TRUE)[[1]]
    if (length(parts) == 2L) {
      return(list(form = "formula",
                  outcome = trimws(parts[[1]]),
                  by      = trimws(parts[[2]])))
    }
  }
  if (grepl("\\s+and\\s+", s)) {
    parts <- strsplit(s, "\\s+and\\s+", perl = TRUE)[[1]]
    if (length(parts) == 2L) {
      return(list(form = "pair",
                  left  = trimws(parts[[1]]),
                  right = trimws(parts[[2]])))
    }
  }
  list(form = "single", name = s)
}

ms_htest_reference_name <- function(name) {
  s <- trimws(as.character(name %||% ""))
  if (!nzchar(s)) return("")

  simple <- regexec("^`?([A-Za-z.][A-Za-z0-9._]*)`?$", s, perl = TRUE)
  m <- regmatches(s, simple)[[1]]
  if (length(m) == 2L) return(m[[2L]])

  dollar <- regexec(
    "^[A-Za-z.][A-Za-z0-9._]*\\$`?([A-Za-z.][A-Za-z0-9._]*)`?$",
    s,
    perl = TRUE
  )
  m <- regmatches(s, dollar)[[1]]
  if (length(m) == 2L) return(m[[2L]])

  bracket <- regexec(
    "^[A-Za-z.][A-Za-z0-9._]*\\[\\[\\s*['\"]([^'\"]+)['\"]\\s*\\]\\]$",
    s,
    perl = TRUE
  )
  m <- regmatches(s, bracket)[[1]]
  if (length(m) == 2L) return(m[[2L]])

  comma <- regexec(
    "^[A-Za-z.][A-Za-z0-9._]*\\[\\s*,\\s*['\"]([^'\"]+)['\"]\\s*\\]$",
    s,
    perl = TRUE
  )
  m <- regmatches(s, comma)[[1]]
  if (length(m) == 2L) return(m[[2L]])

  ""
}

# Resolve a variable reference (bare name or "df$col") in `.data`. Returns
# NULL when the column isn't present.
ms_htest_resolve_var <- function(name, data) {
  if (is.null(data) || !is.data.frame(data)) return(NULL)
  nm <- ms_htest_reference_name(name)
  if (!nzchar(nm) || !(nm %in% names(data))) return(NULL)
  data[[nm]]
}

ms_htest_resolve_var_env <- function(name, env) {
  if (!is.environment(env)) return(NULL)
  s <- trimws(as.character(name %||% ""))
  if (!nzchar(s)) return(NULL)

  simple <- regexec("^`?([A-Za-z.][A-Za-z0-9._]*)`?$", s, perl = TRUE)
  m <- regmatches(s, simple)[[1]]
  if (length(m) == 2L) {
    nm <- m[[2L]]
    if (exists(nm, envir = env, inherits = TRUE)) {
      return(get(nm, envir = env, inherits = TRUE))
    }
    return(NULL)
  }

  dollar <- regexec(
    "^([A-Za-z.][A-Za-z0-9._]*)\\$`?([A-Za-z.][A-Za-z0-9._]*)`?$",
    s,
    perl = TRUE
  )
  m <- regmatches(s, dollar)[[1]]
  if (length(m) == 3L && exists(m[[2L]], envir = env, inherits = TRUE)) {
    obj <- get(m[[2L]], envir = env, inherits = TRUE)
    key <- m[[3L]]
    if ((is.data.frame(obj) || is.list(obj)) && key %in% names(obj)) {
      return(obj[[key]])
    }
    return(NULL)
  }

  bracket <- regexec(
    "^([A-Za-z.][A-Za-z0-9._]*)\\[\\[\\s*['\"]([^'\"]+)['\"]\\s*\\]\\]$",
    s,
    perl = TRUE
  )
  m <- regmatches(s, bracket)[[1]]
  if (length(m) == 3L && exists(m[[2L]], envir = env, inherits = TRUE)) {
    obj <- get(m[[2L]], envir = env, inherits = TRUE)
    key <- m[[3L]]
    if ((is.data.frame(obj) || is.list(obj)) && key %in% names(obj)) {
      return(obj[[key]])
    }
    return(NULL)
  }

  comma <- regexec(
    "^([A-Za-z.][A-Za-z0-9._]*)\\[\\s*,\\s*['\"]([^'\"]+)['\"]\\s*\\]$",
    s,
    perl = TRUE
  )
  m <- regmatches(s, comma)[[1]]
  if (length(m) == 3L && exists(m[[2L]], envir = env, inherits = TRUE)) {
    obj <- get(m[[2L]], envir = env, inherits = TRUE)
    key <- m[[3L]]
    if (is.data.frame(obj) && key %in% names(obj)) {
      return(obj[[key]])
    }
  }

  NULL
}

# Per-group observed descriptives as a list of GroupSummary rows. Group
# ordering follows the order of `levels` when supplied (so labels in
# fields$groups line up with the test object's direction).
ms_htest_group_descriptives <- function(values, group, levels = NULL) {
  if (length(values) != length(group)) return(NULL)
  keep <- is.finite(values) & !is.na(group)
  values <- values[keep]
  group  <- group[keep]
  if (length(values) == 0L) return(NULL)
  g <- if (is.factor(group)) group else factor(group)
  lvls <- if (!is.null(levels)) levels else levels(g)
  out <- lapply(lvls, function(lvl) {
    idx <- which(as.character(g) == as.character(lvl))
    if (!length(idx)) return(NULL)
    vals <- values[idx]
    qs <- ms_safe_quantile(vals, probs = c(0.25, 0.5, 0.75))
    list(
      label = as.character(lvl),
      n     = length(idx),
      mean  = ms_safe_numeric(mean(vals)),
      sd    = ms_safe_numeric(stats::sd(vals)),
      median = ms_safe_numeric(qs[[2L]]),
      q1 = ms_safe_numeric(qs[[1L]]),
      q3 = ms_safe_numeric(qs[[3L]]),
      iqr = ms_safe_numeric(stats::IQR(vals)),
      min = ms_safe_numeric(min(vals)),
      max = ms_safe_numeric(max(vals))
    )
  })
  out <- Filter(Negate(is.null), out)
  if (!length(out)) return(NULL)
  out
}

ms_safe_quantile <- function(values, probs) {
  tryCatch(
    ms_safe_numeric(stats::quantile(values, probs = probs, names = FALSE, na.rm = TRUE)),
    error = function(e) rep(NA_real_, length(probs))
  )
}

ms_htest_env_data <- function(parsed, type, env, x = NULL) {
  if (is.null(parsed) || !is.environment(env)) return(NULL)

  make_data <- function(refs, values) {
    if (any(vapply(values, is.null, logical(1)))) return(NULL)
    lens <- vapply(values, length, integer(1))
    if (!length(lens) || any(lens != lens[[1L]])) return(NULL)
    nms <- vapply(refs, ms_htest_reference_name, character(1))
    if (any(!nzchar(nms)) || anyDuplicated(nms)) return(NULL)
    out <- as.data.frame(values, optional = TRUE, stringsAsFactors = FALSE)
    names(out) <- nms
    out
  }

  if (type %in% c("welch_t_test", "students_t_test",
                  "wilcoxon_rank_sum", "kruskal_wallis_test") &&
      identical(parsed$form, "formula")) {
    refs <- c(parsed$outcome, parsed$by)
    values <- lapply(refs, ms_htest_resolve_var_env, env = env)
    data <- make_data(refs, values)
    if (!is.null(data) && ms_htest_candidate_matches(data, x, type)) {
      return(data)
    }
    return(ms_htest_find_env_data_frame(refs, env, x = x, type = type))
  }

  if (type %in% c("paired_t_test", "wilcoxon_signed_rank") &&
      identical(parsed$form, "pair")) {
    refs <- c(parsed$left, parsed$right)
    values <- lapply(refs, ms_htest_resolve_var_env, env = env)
    return(make_data(refs, values))
  }

  if (type %in% c("one_sample_t_test", "wilcoxon_signed_rank") &&
      identical(parsed$form, "single")) {
    refs <- parsed$name
    values <- list(ms_htest_resolve_var_env(parsed$name, env))
    return(make_data(refs, values))
  }

  NULL
}

ms_htest_find_env_data_frame <- function(refs, env, x = NULL, type = NULL) {
  if (!is.environment(env)) return(NULL)
  cols <- vapply(refs, ms_htest_reference_name, character(1))
  if (any(!nzchar(cols)) || anyDuplicated(cols)) return(NULL)

  find_candidates <- function(search_env) {
    names_in_env <- ls(envir = search_env, all.names = TRUE)
    candidates <- list()
    for (nm in names_in_env) {
      obj <- tryCatch(get(nm, envir = search_env, inherits = FALSE),
                      error = function(e) NULL)
      if (!is.data.frame(obj)) next
      if (!all(cols %in% names(obj))) next
      candidates[[length(candidates) + 1L]] <- obj[cols]
    }
    candidates
  }

  resolve_candidates <- function(candidates) {
    if (!length(candidates)) return(NULL)
    matches <- Filter(function(candidate) {
      ms_htest_candidate_matches(candidate, x, type)
    }, candidates)
    if (length(matches) == 1L) return(matches[[1L]])
    if (length(matches) > 1L) {
      signatures <- vapply(matches, ms_htest_candidate_signature,
                           character(1), type = type)
      if (length(unique(signatures)) == 1L) return(matches[[1L]])
    }
    if (is.null(x) && length(candidates) == 1L) return(candidates[[1L]])
    NULL
  }

  # Prefer local data frames, but when there are several, only use one if
  # it reproduces the stored t-test exactly. This handles workspaces with
  # sleep_df/sleep_a/sleep_b copies without blindly picking by name.
  local_candidates <- find_candidates(env)
  local_match <- resolve_candidates(local_candidates)
  if (!is.null(local_match)) return(local_match)
  if (length(local_candidates) > 0L) return(NULL)

  # If no local data frame is present, walk parent environments. This
  # recovers common examples such as datasets::sleep while keeping the
  # same "exactly one match" guard.
  parent_candidates <- list()
  seen <- character(0)
  search_env <- parent.env(env)
  while (is.environment(search_env) && !identical(search_env, emptyenv())) {
    env_name <- environmentName(search_env)
    env_key <- if (nzchar(env_name)) env_name else utils::capture.output(print(search_env))[1]
    if (env_key %in% seen) break
    seen <- c(seen, env_key)
    found <- find_candidates(search_env)
    if (length(found)) parent_candidates <- c(parent_candidates, found)
    if (identical(search_env, baseenv())) break
    search_env <- parent.env(search_env)
  }
  resolve_candidates(parent_candidates)
}

ms_htest_candidate_matches <- function(data, x, type) {
  if (is.null(x) || is.null(type)) return(TRUE)
  if (!type %in% c("welch_t_test", "students_t_test",
                   "wilcoxon_rank_sum", "kruskal_wallis_test")) return(TRUE)
  if (!is.data.frame(data) || ncol(data) < 2L) return(FALSE)

  y <- data[[1L]]
  g <- data[[2L]]
  if (length(y) != length(g)) return(FALSE)
  keep <- !is.na(y) & !is.na(g)
  y <- y[keep]
  g <- g[keep]
  if (length(y) < 2L) return(FALSE)

  fit_data <- data.frame(.y = y, .g = g)
  fit <- tryCatch({
    if (type %in% c("welch_t_test", "students_t_test")) {
      stats::t.test(.y ~ .g, data = fit_data,
                    var.equal = identical(type, "students_t_test"))
    } else if (identical(type, "wilcoxon_rank_sum")) {
      suppressWarnings(stats::wilcox.test(.y ~ .g, data = fit_data))
    } else if (identical(type, "kruskal_wallis_test")) {
      stats::kruskal.test(.y ~ .g, data = fit_data)
    } else {
      NULL
    }
  }, error = function(e) NULL)
  if (is.null(fit)) return(FALSE)

  same_numeric <- function(a, b, tolerance = 1e-7) {
    a <- ms_safe_numeric(unname(a))
    b <- ms_safe_numeric(unname(b))
    if (length(a) != length(b)) return(FALSE)
    all(is.finite(a) == is.finite(b)) &&
      isTRUE(all.equal(a, b, tolerance = tolerance, check.attributes = FALSE))
  }

  same_numeric(fit$statistic, x$statistic) &&
    same_numeric(fit$parameter %||% numeric(0), x$parameter %||% numeric(0)) &&
    same_numeric(fit$p.value, x$p.value) &&
    same_numeric(fit$estimate %||% numeric(0), x$estimate %||% numeric(0))
}

ms_htest_candidate_signature <- function(data, type = NULL) {
  if (!is.data.frame(data) || ncol(data) < 2L) return("")
  desc <- ms_htest_group_descriptives(data[[1L]], data[[2L]])
  if (is.null(desc)) return("")
  paste(vapply(desc, function(row) {
    paste(
      row$label %||% "",
      row$n %||% "",
      signif(row$mean %||% NA_real_, 12),
      signif(row$sd %||% NA_real_, 12),
      sep = ":"
    )
  }, character(1)), collapse = "|")
}

ms_htest_enrich_from_env <- function(x, type, env, existing_groups = NULL) {
  parsed <- ms_htest_parse_data_name(x$data.name)
  data <- ms_htest_env_data(parsed, type, env, x = x)
  if (is.null(data)) return(NULL)
  ms_htest_enrich_from_data(x, type, data, existing_groups = existing_groups)
}

# Cohen's d (point estimate only — CI is out of scope per the plan).
# Method tag tells the renderer which denominator was used, so users
# downstream can distinguish pooled-SD / paired-dz / one-sample.
ms_htest_cohens_d <- function(x, type, data, parsed) {
  if (type %in% c("welch_t_test", "students_t_test")) {
    if (is.null(parsed) || !identical(parsed$form, "formula")) return(NULL)
    y <- ms_htest_resolve_var(parsed$outcome, data)
    g <- ms_htest_resolve_var(parsed$by, data)
    if (is.null(y) || is.null(g)) return(NULL)
    keep <- !is.na(y) & !is.na(g)
    y <- y[keep]; g <- factor(g[keep])
    lvls <- levels(g)
    if (length(lvls) != 2L) return(NULL)
    a <- y[as.character(g) == lvls[[1]]]
    b <- y[as.character(g) == lvls[[2]]]
    n1 <- length(a); n2 <- length(b)
    if (n1 < 2L || n2 < 2L) return(NULL)
    s1 <- stats::sd(a); s2 <- stats::sd(b)
    denominator <- if (identical(type, "welch_t_test")) {
      sqrt((s1^2 + s2^2) / 2)
    } else {
      sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
    }
    if (!is.finite(denominator) || denominator <= 0) return(NULL)
    d <- (mean(a) - mean(b)) / denominator
    return(list(
      name   = "cohens_d",
      value  = ms_safe_numeric(d),
      method = if (identical(type, "welch_t_test")) "average_group_sd" else "pooled_sd"
    ))
  }
  if (identical(type, "paired_t_test")) {
    # Build the difference vector from `.data`. Two-vector form
    # ("data: pre and post") looks up both columns; formula form
    # ("data: y by g") splits by the grouping factor — but paired
    # formulas are rare and order-fragile, so we focus on the
    # two-vector form which is the common paired t-test idiom.
    diffs <- NULL
    if (!is.null(parsed) && identical(parsed$form, "pair")) {
      l <- ms_htest_resolve_var(parsed$left, data)
      r <- ms_htest_resolve_var(parsed$right, data)
      if (!is.null(l) && !is.null(r) && length(l) == length(r)) {
        diffs <- l - r
      }
    }
    if (is.null(diffs)) return(NULL)
    diffs <- diffs[is.finite(diffs)]
    n <- length(diffs)
    if (n < 2L) return(NULL)
    s <- stats::sd(diffs)
    if (!is.finite(s) || s <= 0) return(NULL)
    return(list(
      name   = "cohens_dz",
      value  = ms_safe_numeric(mean(diffs) / s),
      method = "paired_dz"
    ))
  }
  if (identical(type, "one_sample_t_test")) {
    if (is.null(parsed) || !identical(parsed$form, "single")) return(NULL)
    y <- ms_htest_resolve_var(parsed$name, data)
    if (is.null(y)) return(NULL)
    y <- y[is.finite(y)]
    if (length(y) < 2L) return(NULL)
    s <- stats::sd(y)
    if (!is.finite(s) || s <= 0) return(NULL)
    mu0 <- ms_safe_numeric(unname(x$null.value %||% 0)[[1]])
    if (is.na(mu0)) mu0 <- 0
    return(list(
      name   = "cohens_d",
      value  = ms_safe_numeric((mean(y) - mu0) / s),
      method = "one_sample"
    ))
  }
  NULL
}

ms_htest_nonparametric_effect_size <- function(x, type, data, parsed) {
  if (identical(type, "wilcoxon_rank_sum")) {
    if (is.null(parsed) || !identical(parsed$form, "formula")) return(NULL)
    y <- ms_htest_resolve_var(parsed$outcome, data)
    g <- ms_htest_resolve_var(parsed$by, data)
    if (is.null(y) || is.null(g) || length(y) != length(g)) return(NULL)
    keep <- is.finite(y) & !is.na(g)
    y <- y[keep]
    g <- factor(g[keep])
    lvls <- levels(g)
    if (length(lvls) != 2L) return(NULL)
    a <- y[as.character(g) == lvls[[1L]]]
    b <- y[as.character(g) == lvls[[2L]]]
    delta <- ms_cliffs_delta(a, b)
    if (is.na(delta)) return(NULL)
    return(list(
      name = "cliffs_delta",
      value = ms_safe_numeric(delta),
      method = "pairwise_order_probability",
      group1 = as.character(lvls[[1L]]),
      group2 = as.character(lvls[[2L]])
    ))
  }

  if (identical(type, "wilcoxon_signed_rank")) {
    diffs <- NULL
    if (!is.null(parsed) && identical(parsed$form, "pair")) {
      l <- ms_htest_resolve_var(parsed$left, data)
      r <- ms_htest_resolve_var(parsed$right, data)
      if (!is.null(l) && !is.null(r) && length(l) == length(r)) diffs <- l - r
    } else if (!is.null(parsed) && identical(parsed$form, "single")) {
      y <- ms_htest_resolve_var(parsed$name, data)
      if (!is.null(y)) {
        mu0 <- ms_safe_numeric(unname(x$null.value %||% 0)[[1]])
        if (is.na(mu0)) mu0 <- 0
        diffs <- y - mu0
      }
    }
    rbc <- ms_signed_rank_biserial(diffs)
    if (is.na(rbc)) return(NULL)
    return(list(
      name = "rank_biserial",
      value = ms_safe_numeric(rbc),
      method = "matched_pairs_rank_biserial"
    ))
  }

  if (identical(type, "kruskal_wallis_test")) {
    groups <- NULL
    if (!is.null(parsed) && identical(parsed$form, "formula")) {
      y <- ms_htest_resolve_var(parsed$outcome, data)
      g <- ms_htest_resolve_var(parsed$by, data)
      if (!is.null(y) && !is.null(g) && length(y) == length(g)) {
        keep <- is.finite(y) & !is.na(g)
        groups <- factor(g[keep])
      }
    }
    if (is.null(groups)) return(NULL)
    n <- length(groups)
    k <- length(levels(groups))
    h <- ms_safe_numeric(unname(x$statistic %||% NA_real_)[[1L]])
    if (is.na(h) || !is.finite(h) || n <= k || k < 2L) return(NULL)
    eta_sq_h <- (h - k + 1) / (n - k)
    if (!is.finite(eta_sq_h)) return(NULL)
    eta_sq_h <- max(0, min(1, eta_sq_h))
    return(list(
      name = "eta_sq_h",
      value = ms_safe_numeric(eta_sq_h),
      method = "kruskal_wallis_eta_squared_h",
      formula = "(H - k + 1) / (N - k)"
    ))
  }

  NULL
}

ms_cliffs_delta <- function(a, b) {
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  if (length(a) == 0L || length(b) == 0L) return(NA_real_)
  comparisons <- outer(a, b, "-")
  delta <- mean(sign(comparisons))
  if (is.finite(delta)) delta else NA_real_
}

ms_signed_rank_biserial <- function(diffs) {
  if (is.null(diffs)) return(NA_real_)
  diffs <- diffs[is.finite(diffs) & diffs != 0]
  if (length(diffs) == 0L) return(NA_real_)
  ranks <- rank(abs(diffs), ties.method = "average")
  denom <- sum(ranks)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  value <- (sum(ranks[diffs > 0]) - sum(ranks[diffs < 0])) / denom
  if (is.finite(value)) value else NA_real_
}

# Compose the two enrichment outputs together. Group descriptives use
# the same label ordering as fields$groups (recovered from R's
# `mean in group <X>` names) so display lines up with the test object.
ms_htest_enrich_from_data <- function(x, type, data, existing_groups = NULL) {
  parsed <- ms_htest_parse_data_name(x$data.name)
  out <- list(
    groups = NULL,
    paired_measures = NULL,
    effect_size = NULL,
    grouping_note = NULL
  )

  # ── Group-wise n/M/SD ────────────────────────────────────────────
  if (type %in% c("welch_t_test", "students_t_test") &&
      !is.null(parsed) && identical(parsed$form, "formula")) {
    y <- ms_htest_resolve_var(parsed$outcome, data)
    g <- ms_htest_resolve_var(parsed$by, data)
    if (!is.null(y) && !is.null(g)) {
      # Use the existing labels order (from R's $estimate names) so
      # the descriptives align with the mean-diff direction the test
      # computed.
      lvls <- if (!is.null(existing_groups)) {
        vapply(existing_groups, function(grp) grp$label %||% NA_character_,
               character(1))
      } else NULL
      desc <- ms_htest_group_descriptives(y, g, levels = lvls)
      if (!is.null(desc)) out$groups <- desc
    }
  } else if (identical(type, "paired_t_test") &&
             !is.null(parsed) && identical(parsed$form, "pair")) {
    l <- ms_htest_resolve_var(parsed$left, data)
    r <- ms_htest_resolve_var(parsed$right, data)
    if (!is.null(l) && !is.null(r) && length(l) == length(r)) {
      keep <- is.finite(l) & is.finite(r)
      l_complete <- l[keep]
      r_complete <- r[keep]
      diffs <- l_complete - r_complete
      if (length(diffs) >= 2L) {
        left_label <- ms_htest_reference_name(parsed$left)
        right_label <- ms_htest_reference_name(parsed$right)
        if (!nzchar(left_label)) left_label <- "Measure 1"
        if (!nzchar(right_label)) right_label <- "Measure 2"
        out$paired_measures <- list(
          list(
            label = left_label,
            role  = "first_measure",
            n     = length(l_complete),
            mean  = ms_safe_numeric(mean(l_complete)),
            sd    = ms_safe_numeric(stats::sd(l_complete))
          ),
          list(
            label = right_label,
            role  = "second_measure",
            n     = length(r_complete),
            mean  = ms_safe_numeric(mean(r_complete)),
            sd    = ms_safe_numeric(stats::sd(r_complete))
          )
        )
        out$groups <- list(list(
          label = "Differences",
          n     = length(diffs),
          mean  = ms_safe_numeric(mean(diffs)),
          sd    = ms_safe_numeric(stats::sd(diffs))
        ))
      }
    }
  } else if (identical(type, "one_sample_t_test") &&
             !is.null(parsed) && identical(parsed$form, "single")) {
    y <- ms_htest_resolve_var(parsed$name, data)
    if (!is.null(y)) {
      y <- y[is.finite(y)]
      if (length(y) >= 2L) {
        out$groups <- list(list(
          label = ms_htest_reference_name(parsed$name),
          n     = length(y),
          mean  = ms_safe_numeric(mean(y)),
          sd    = ms_safe_numeric(stats::sd(y)),
          median = ms_safe_numeric(stats::median(y)),
          q1 = ms_safe_numeric(stats::quantile(y, 0.25, names = FALSE)),
          q3 = ms_safe_numeric(stats::quantile(y, 0.75, names = FALSE)),
          iqr = ms_safe_numeric(stats::IQR(y)),
          min = ms_safe_numeric(min(y)),
          max = ms_safe_numeric(max(y))
        ))
      }
    }
  } else if (type %in% c("wilcoxon_rank_sum", "kruskal_wallis_test") &&
             !is.null(parsed) && identical(parsed$form, "formula")) {
    y <- ms_htest_resolve_var(parsed$outcome, data)
    g <- ms_htest_resolve_var(parsed$by, data)
    if (!is.null(y) && !is.null(g)) {
      desc <- ms_htest_group_descriptives(y, g)
      if (!is.null(desc)) out$groups <- desc
      out$grouping_note <- ms_numeric_grouping_note(parsed$by, g)
    }
  } else if (identical(type, "wilcoxon_signed_rank") &&
             !is.null(parsed) && identical(parsed$form, "pair")) {
    l <- ms_htest_resolve_var(parsed$left, data)
    r <- ms_htest_resolve_var(parsed$right, data)
    if (!is.null(l) && !is.null(r) && length(l) == length(r)) {
      keep <- is.finite(l) & is.finite(r)
      l_complete <- l[keep]
      r_complete <- r[keep]
      diffs <- l_complete - r_complete
      if (length(diffs) >= 1L) {
        left_label <- ms_htest_reference_name(parsed$left)
        right_label <- ms_htest_reference_name(parsed$right)
        if (!nzchar(left_label)) left_label <- "Measure 1"
        if (!nzchar(right_label)) right_label <- "Measure 2"
        out$paired_measures <- list(
          list(
            label = left_label,
            role = "first_measure",
            n = length(l_complete),
            mean = ms_safe_numeric(mean(l_complete)),
            sd = ms_safe_numeric(stats::sd(l_complete)),
            median = ms_safe_numeric(stats::median(l_complete)),
            iqr = ms_safe_numeric(stats::IQR(l_complete))
          ),
          list(
            label = right_label,
            role = "second_measure",
            n = length(r_complete),
            mean = ms_safe_numeric(mean(r_complete)),
            sd = ms_safe_numeric(stats::sd(r_complete)),
            median = ms_safe_numeric(stats::median(r_complete)),
            iqr = ms_safe_numeric(stats::IQR(r_complete))
          )
        )
        qs <- ms_safe_quantile(diffs, probs = c(0.25, 0.5, 0.75))
        out$groups <- list(list(
          label = "Differences",
          n = length(diffs),
          mean = ms_safe_numeric(mean(diffs)),
          sd = ms_safe_numeric(stats::sd(diffs)),
          median = ms_safe_numeric(qs[[2L]]),
          q1 = ms_safe_numeric(qs[[1L]]),
          q3 = ms_safe_numeric(qs[[3L]]),
          iqr = ms_safe_numeric(stats::IQR(diffs)),
          min = ms_safe_numeric(min(diffs)),
          max = ms_safe_numeric(max(diffs))
        ))
      }
    }
  } else if (identical(type, "wilcoxon_signed_rank") &&
             !is.null(parsed) && identical(parsed$form, "single")) {
    y <- ms_htest_resolve_var(parsed$name, data)
    if (!is.null(y)) {
      mu0 <- ms_safe_numeric(unname(x$null.value %||% 0)[[1]])
      if (is.na(mu0)) mu0 <- 0
      diffs <- y - mu0
      diffs <- diffs[is.finite(diffs)]
      if (length(diffs) >= 1L) {
        qs <- ms_safe_quantile(diffs, probs = c(0.25, 0.5, 0.75))
        out$groups <- list(list(
          label = "Differences",
          n = length(diffs),
          mean = ms_safe_numeric(mean(diffs)),
          sd = ms_safe_numeric(stats::sd(diffs)),
          median = ms_safe_numeric(qs[[2L]]),
          q1 = ms_safe_numeric(qs[[1L]]),
          q3 = ms_safe_numeric(qs[[3L]]),
          iqr = ms_safe_numeric(stats::IQR(diffs)),
          min = ms_safe_numeric(min(diffs)),
          max = ms_safe_numeric(max(diffs))
        ))
      }
    }
  }

  # ── Cohen's d (no CI; out of scope per the plan) ─────────────────
  d <- tryCatch(ms_htest_cohens_d(x, type, data, parsed),
                error = function(e) NULL)
  if (!is.null(d)) out$effect_size <- d
  if (is.null(out$effect_size) &&
      type %in% c("wilcoxon_rank_sum", "wilcoxon_signed_rank", "kruskal_wallis_test")) {
    es <- tryCatch(ms_htest_nonparametric_effect_size(x, type, data, parsed),
                   error = function(e) NULL)
    if (!is.null(es)) out$effect_size <- es
  }

  if (is.null(out$groups) && is.null(out$paired_measures) &&
      is.null(out$effect_size) && is.null(out$grouping_note)) {
    return(NULL)
  }
  out
}

ms_numeric_grouping_note <- function(group_ref, group) {
  if (!is.numeric(group)) return(NULL)
  values <- unique(group[!is.na(group)])
  if (length(values) < 2L) return(NULL)
  label <- ms_htest_reference_name(group_ref)
  if (!nzchar(label)) label <- "The grouping variable"
  paste0(
    label,
    " is numeric in the supplied data; Mellio treats its distinct values as group levels to match the R test."
  )
}
