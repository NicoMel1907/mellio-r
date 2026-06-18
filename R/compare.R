# Side-by-side model comparison (stargazer-style)

#' Compare multiple models side-by-side
#'
#' Creates a stargazer-style comparison table with models as columns,
#' stacked coefficient/SE rows, significance stars, and model statistics.
#'
#' Accepts any model objects supported by `broom::tidy()` and
#' `broom::glance()`, including lm, glm, and lme4 (via broom.mixed).
#' You can call this helper directly, or pass multiple model objects to
#' [melliotab()] or [mellio_open()], which delegate here.
#'
#' @param ... Model objects to compare
#' @param style Citation style (default "apa7")
#' @param title Table title
#' @param number Table number
#' @param note Table note (NULL auto-generates significance note when stars = TRUE)
#' @param source Table source text
#' @param column.labels Character vector of model labels (e.g., c("Step 1", "Step 2")).
#'   Defaults to "Model 1", "Model 2", etc.
#' @param dep.var.labels Dependent variable label (rendered as spanner across
#'   all model columns)
#' @param se.type How to display standard errors: "parentheses" (default,
#'   separate indented row), "inline" (same cell), "none" (omit)
#' @param stars Logical, append significance stars (default TRUE)
#' @param stats Character vector of model statistics for bottom rows.
#'   Default: c("n", "r.squared", "adj.r.squared", "f.statistic").
#'   Available: "n", "r.squared", "adj.r.squared", "f.statistic",
#'   "aic", "bic", "loglik", "residual.se"
#' @param decimals Decimal places for coefficients (default 2)
#' @param p_decimals Decimal places for p-values (default 3)
#' @param omit Character vector of term names to exclude
#' @param intercept.bottom Place intercept at bottom of coefficient rows
#'   (default TRUE)
#' @return A melliotab object
#' @export
#'
#' @examples
#' m1 <- lm(Ozone ~ Temp, data = airquality)
#' m2 <- lm(Ozone ~ Temp + Wind, data = airquality)
#' m3 <- lm(Ozone ~ Temp + Wind + Solar.R, data = airquality)
#' mt_compare(
#'   m1, m2, m3,
#'   title = "Predictors of ozone concentration",
#'   column.labels = c("Step 1", "Step 2", "Step 3"),
#'   dep.var.labels = "Ozone concentration"
#' )
mt_compare <- function(...,
                       style = "apa7",
                       title = NULL,
                       number = NULL,
                       note = NULL,
                       source = NULL,
                       column.labels = NULL,
                       dep.var.labels = NULL,
                       se.type = c("parentheses", "inline", "none"),
                       stars = TRUE,
                       stats = c("n", "r.squared", "adj.r.squared",
                                 "f.statistic"),
                       decimals = 2L,
                       p_decimals = 3L,
                       omit = NULL,
                       intercept.bottom = TRUE) {

  rlang::check_installed("broom", reason = "to extract model coefficients")

  models <- list(...)
  n_models <- length(models)

  if (n_models == 0) {
    cli::cli_abort("mt_compare() requires at least one model object.")
  }

  se.type <- match.arg(se.type)
  style <- match.arg(style, list_styles())
  style_config <- get_style(style)
  remove_lz <- isTRUE(style_config$remove_leading_zeros)

  # --- Validate column.labels length ---
  if (!is.null(column.labels) && length(column.labels) != n_models) {
    cli::cli_abort(
      "{.arg column.labels} has {length(column.labels)} label{?s} but {n_models} model{?s} {?was/were} provided."
    )
  }
  model_labels <- column.labels %||% paste("Model", seq_len(n_models))

  # --- Extract tidy + glance from each model ---
  tidy_list <- vector("list", n_models)
  glance_list <- vector("list", n_models)

  for (k in seq_len(n_models)) {
    m <- models[[k]]
    cls <- paste(class(m), collapse = "/")

    # tidy
    tidy_list[[k]] <- tryCatch(
      broom::tidy(m),
      error = function(e) {
        cli::cli_abort(c(
          "Could not extract coefficients from model {k} ({.cls {cls}}).",
          "i" = "broom::tidy() failed: {conditionMessage(e)}",
          "i" = "Supported types include lm, glm, and other models with broom methods."
        ))
      }
    )

    # glance (non-fatal)
    glance_list[[k]] <- tryCatch(
      broom::glance(m),
      error = function(e) {
        cli::cli_warn(
          "Could not extract summary statistics from model {k} ({.cls {cls}}): {conditionMessage(e)}"
        )
        NULL
      }
    )
  }

  # --- Union all term names ---
  all_terms <- unique(unlist(lapply(tidy_list, function(td) td$term)))

  # Move intercept to bottom
  if (intercept.bottom) {
    int_idx <- which(all_terms == "(Intercept)")
    if (length(int_idx) > 0) {
      all_terms <- c(all_terms[-int_idx], "(Intercept)")
    }
  }

  # Apply omit filter
  if (!is.null(omit)) {
    all_terms <- all_terms[!all_terms %in% omit]
    if (length(all_terms) == 0) {
      cli::cli_abort("No terms remain after applying {.arg omit} filter.")
    }
  }

  # --- Format coefficients and SEs for each model ---
  coef_formatted <- vector("list", n_models)
  se_formatted <- vector("list", n_models)

  for (k in seq_len(n_models)) {
    td <- tidy_list[[k]]
    coefs <- setNames(rep("", length(all_terms)), all_terms)
    ses <- setNames(rep("", length(all_terms)), all_terms)

    for (term in all_terms) {
      row <- td[td$term == term, , drop = FALSE]
      if (nrow(row) == 0) next

      # Format coefficient
      coef_val <- format_apa_number(
        as.character(row$estimate[1]), "estimate",
        decimals = decimals, remove_lz = remove_lz
      )

      # Add significance stars
      if (stars && "p.value" %in% names(row)) {
        star_str <- get_sig_stars(row$p.value[1])
        if (nzchar(star_str)) {
          coef_val <- paste0(coef_val, star_str)
        }
      }
      coefs[term] <- coef_val

      # Format SE
      if ("std.error" %in% names(row)) {
        se_val <- format_apa_number(
          as.character(row$std.error[1]), "estimate",
          decimals = decimals, remove_lz = remove_lz
        )
        ses[term] <- paste0("(", se_val, ")")
      }
    }

    coef_formatted[[k]] <- coefs
    se_formatted[[k]] <- ses
  }

  # --- Clean term labels ---
  term_labels <- .compare_clean_terms(all_terms)

  # --- Build rows based on se.type ---
  if (se.type == "parentheses") {
    n_rows <- length(all_terms) * 2
    stub_col <- character(n_rows)
    model_cols <- matrix("", nrow = n_rows, ncol = n_models)
    se_row_indices <- integer(0)

    for (i in seq_along(all_terms)) {
      coef_row <- (i - 1) * 2 + 1
      se_row <- coef_row + 1

      stub_col[coef_row] <- term_labels[i]
      stub_col[se_row] <- ""

      se_row_indices <- c(se_row_indices, se_row)

      for (k in seq_len(n_models)) {
        model_cols[coef_row, k] <- coef_formatted[[k]][all_terms[i]]
        model_cols[se_row, k] <- se_formatted[[k]][all_terms[i]]
      }
    }
  } else if (se.type == "inline") {
    n_rows <- length(all_terms)
    stub_col <- term_labels
    model_cols <- matrix("", nrow = n_rows, ncol = n_models)
    se_row_indices <- integer(0)

    for (i in seq_along(all_terms)) {
      for (k in seq_len(n_models)) {
        coef <- coef_formatted[[k]][all_terms[i]]
        se <- se_formatted[[k]][all_terms[i]]
        if (nzchar(coef) && nzchar(se)) {
          model_cols[i, k] <- paste0(coef, " ", se)
        } else {
          model_cols[i, k] <- coef
        }
      }
    }
  } else {
    # se.type == "none"
    n_rows <- length(all_terms)
    stub_col <- term_labels
    model_cols <- matrix("", nrow = n_rows, ncol = n_models)
    se_row_indices <- integer(0)

    for (i in seq_along(all_terms)) {
      for (k in seq_len(n_models)) {
        model_cols[i, k] <- coef_formatted[[k]][all_terms[i]]
      }
    }
  }

  # --- Build model statistics rows ---
  stat_data <- .compare_build_stat_rows(
    glance_list, models, stats, decimals, remove_lz
  )
  stat_start <- length(stub_col) + 1L

  # --- Assemble wide data frame ---
  all_stub <- c(stub_col, stat_data$label)
  wide_df <- data.frame(` ` = all_stub, check.names = FALSE,
                         stringsAsFactors = FALSE)

  for (k in seq_len(n_models)) {
    col_vals <- c(model_cols[, k], stat_data[[k + 1]])
    wide_df[[model_labels[k]]] <- col_vals
  }

  # --- Build melliotab object ---
  col_types <- c("stub", rep("default", n_models))

  result <- structure(
    list(
      data = wide_df,
      raw_data = wide_df,
      column_types = col_types,
      style = style,
      style_config = style_config,
      title = title,
      number = number,
      note = note,
      source = source,
      spanners = list(),
      section_titles = list(),
      merged_regions = list(),
      indent_levels = integer(0),
      decimals = as.integer(decimals),
      p_decimals = as.integer(p_decimals),
      options = list(
        remove_leading_zeros = remove_lz,
        italic_stat_headers = isTRUE(style_config$italic_stat_headers),
        bold_section_titles = isTRUE(style_config$bold_section_titles),
        sig_stars = stars,
        format_ci = FALSE,
        diagonal_mode = "all",
        triangle = "all",
        is_correlation = FALSE,
        is_comparison = TRUE,
        comparison_stat_start = if (nrow(stat_data) > 0L) stat_start else NA_integer_
      ),
      model = models,
      model_summary = glance_list
    ),
    class = "melliotab"
  )

  # --- Add dep.var.labels spanner ---
  if (!is.null(dep.var.labels)) {
    result <- mt_spanner(result, label = "Dependent variable:",
                         columns = seq(2, n_models + 1), level = 2L)
    result <- mt_spanner(result, label = dep.var.labels,
                         columns = seq(2, n_models + 1), level = 1L)
    result$options$dependent_variable_label <- dep.var.labels
  }

  # --- Indent SE rows ---
  if (se.type == "parentheses" && length(se_row_indices) > 0) {
    indent <- result$indent_levels
    for (r in se_row_indices) {
      indent[as.character(r)] <- 1L
    }
    result$indent_levels <- indent
  }

  # --- Auto-generate significance note ---
  if (is.null(note) && stars) {
    result$note <- "*p < .05. **p < .01. ***p < .001."
  }

  result
}


# --- Internal helpers ---

#' Clean R term names to display labels
#' @keywords internal
.compare_clean_terms <- function(terms) {
  labels <- terms

  # (Intercept) -> Intercept
  labels <- sub("^\\(Intercept\\)$", "Intercept", labels)

  # factor(var)level -> var: level
  labels <- gsub("^factor\\(([^)]+)\\)(.+)$", "\\1: \\2", labels)

  # poly(var, n)k -> var^k
  labels <- gsub("^poly\\(([^,]+),\\s*\\d+\\)(\\d+)$", "\\1^\\2", labels)

  # I(var^2) -> var^2
  labels <- gsub("^I\\((.+)\\)$", "\\1", labels)

  labels
}


#' Build model statistics rows
#' @keywords internal
.compare_build_stat_rows <- function(glance_list, models, stats, decimals,
                                      remove_lz) {
  n_models <- length(models)

  # Stat metadata: key -> (glance column, display label, format type)
  stat_meta <- list(
    n = list(col = "nobs", label = "Observations", fmt = "integer"),
    r.squared = list(col = "r.squared", label = "R\u00B2", fmt = "bounded"),
    adj.r.squared = list(col = "adj.r.squared", label = "Adjusted R\u00B2",
                          fmt = "bounded"),
    f.statistic = list(col = "statistic", label = "F Statistic", fmt = "f"),
    aic = list(col = "AIC", label = "AIC", fmt = "decimal"),
    bic = list(col = "BIC", label = "BIC", fmt = "decimal"),
    loglik = list(col = "logLik", label = "Log Likelihood", fmt = "decimal"),
    residual.se = list(col = "sigma", label = "Residual Std. Error",
                        fmt = "decimal")
  )

  # Filter to requested stats
  stats <- stats[stats %in% names(stat_meta)]
  if (length(stats) == 0) {
    return(data.frame(label = character(0), stringsAsFactors = FALSE))
  }

  rows <- data.frame(label = character(length(stats)),
                      stringsAsFactors = FALSE)
  for (k in seq_len(n_models)) {
    rows[[k + 1]] <- character(length(stats))
  }

  for (s in seq_along(stats)) {
    meta <- stat_meta[[stats[s]]]
    rows$label[s] <- meta$label

    for (k in seq_len(n_models)) {
      gl <- glance_list[[k]]
      if (is.null(gl)) {
        rows[[k + 1]][s] <- ""
        next
      }

      val <- gl[[meta$col]]

      # For n, try nobs() as fallback
      if (stats[s] == "n" && (is.null(val) || is.na(val))) {
        val <- tryCatch(stats::nobs(models[[k]]), error = function(e) NA)
      }

      if (is.null(val) || is.na(val)) {
        rows[[k + 1]][s] <- ""
        next
      }

      # Format based on type
      formatted <- switch(meta$fmt,
        integer = as.character(as.integer(val)),
        bounded = {
          fv <- formatC(val, digits = decimals, format = "f")
          if (remove_lz) {
            fv <- sub("^0\\.", ".", fv)
            fv <- sub("^-0\\.", "-.", fv)
          }
          fv
        },
        decimal = formatC(val, digits = decimals, format = "f"),
        f = {
          f_val <- formatC(val, digits = decimals, format = "f")
          # Try to include df info
          num_df <- gl[["df"]]
          df2 <- gl[["df.residual"]]
          if (!is.null(num_df) && !is.null(df2) &&
              !is.na(num_df) && !is.na(df2)) {
            paste0(f_val, " (df = ", as.integer(num_df), "; ",
                   as.integer(df2), ")")
          } else {
            f_val
          }
        },
        as.character(val)
      )

      rows[[k + 1]][s] <- formatted
    }
  }

  rows
}
