# R bridge + table extractors for mediation::mediate objects.
#
# A mediate object is already a compact statistical result: ACME, ADE,
# total effect, and proportion mediated. Mellio keeps the model fitting
# in R and receives only those reportable rows.

#' @rdname mellio_payload
#' @export
mellio_payload.mediate <- function(x, ..., .call = NULL) {
  fields <- ms_mediate_payload_fields(x)

  call_str <- if (!is.null(.call)) {
    .call
  } else {
    ms_deparse_call(x$call)
  }

  ms_build_envelope(
    type       = "mediation_mediate",
    type_label = "Mediation analysis (mediation::mediate)",
    call       = trimws(gsub("\\s+", " ", call_str %||% "mediation::mediate(...)")),
    fields     = fields,
    raw_output = ms_mediate_raw_output(x),
    packages   = ms_packages_basic("mediation"),
    card_kind  = "table"
  )
}

#' @rdname mellio_payload
#' @export
mellio_payload.summary.mediate <- function(x, ..., .call = NULL) {
  mellio_payload.mediate(x, ..., .call = .call)
}

#' @rdname melliotab
#' @export
melliotab.mediate <- function(x, style = "apa7", title = NULL,
                              number = NULL, note = NULL,
                              source = NULL, decimals = 2L,
                              p_decimals = 3L, ...) {
  df <- ms_mediate_table_data(x, pretty_names = TRUE)
  meta <- ms_mediate_meta(x)

  if (is.null(note)) {
    bits <- c(
      if (!is.na(meta$treatment)) paste0("treatment = ", meta$treatment) else NULL,
      if (!is.na(meta$mediator)) paste0("mediator = ", meta$mediator) else NULL,
      if (!is.na(meta$outcome)) paste0("outcome = ", meta$outcome) else NULL,
      if (!is.na(meta$n)) paste0("N = ", meta$n) else NULL,
      if (!is.na(meta$sims)) paste0("simulations = ", meta$sims) else NULL
    )
    interval <- if (isTRUE(meta$boot)) "nonparametric bootstrap" else "quasi-Bayesian"
    note <- paste0(
      "Mediation estimates from mediation::mediate",
      if (length(bits) > 0L) paste0(" (", paste(bits, collapse = "; "), ")") else "",
      ". Intervals are ", interval, " confidence intervals."
    )
  }

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
#' @export
melliotab.summary.mediate <- function(x, style = "apa7", title = NULL,
                                      number = NULL, note = NULL,
                                      source = NULL, decimals = 2L,
                                      p_decimals = 3L, ...) {
  melliotab.mediate(
    x, style = style, title = title,
    number = number, note = note,
    source = source, decimals = decimals,
    p_decimals = p_decimals, ...
  )
}

# -- Mediation helpers --------------------------------------------------------

ms_is_mediate_result <- function(x) {
  inherits(x, "mediate") || inherits(x, "summary.mediate")
}

ms_mediate_summary <- function(x) {
  if (inherits(x, "summary.mediate")) return(x)
  rlang::check_installed("mediation", reason = "to extract mediate results")
  summary(x)
}

ms_mediate_meta <- function(x) {
  s <- ms_mediate_summary(x)
  labels <- ms_mediate_labels(x)
  list(
    treatment = ms_label_or(labels$treatment,
                            ms_chr_or_na(s$treat %||% x$treat)),
    mediator  = ms_label_or(labels$mediator,
                            ms_chr_or_na(s$mediator %||% x$mediator)),
    outcome   = ms_label_or(labels$outcome,
                            ms_mediate_outcome(s$model.y %||% x$model.y)),
    n         = as.integer(ms_safe_numeric(s$nobs %||% x$nobs)),
    sims      = as.integer(ms_safe_numeric(s$sims %||% x$sims)),
    boot      = isTRUE(s$boot %||% x$boot),
    ci_type   = ms_chr_or_na(s$boot.ci.type %||% x$boot.ci.type),
    conf_level = ms_safe_numeric(s$conf.level %||% x$conf.level),
    interaction = isTRUE(s$INT %||% x$INT)
  )
}

ms_mediate_labels <- function(x) {
  labels <- attr(x, "mellio_labels", exact = TRUE)
  if (!is.list(labels)) labels <- list()
  list(
    treatment = ms_chr_or_na(labels$treatment),
    mediator  = ms_chr_or_na(labels$mediator),
    outcome   = ms_chr_or_na(labels$outcome)
  )
}

ms_label_or <- function(label, fallback) {
  if (!is.na(label)) label else fallback
}

ms_mediate_outcome <- function(model) {
  f <- tryCatch(stats::formula(model), error = function(e) NULL)
  if (is.null(f) || length(f) < 2L) return(NA_character_)
  trimws(paste(deparse(f[[2]], width.cutoff = 500L), collapse = " "))
}

ms_mediate_payload_fields <- function(x) {
  rows <- ms_mediate_rows(x)
  has_contrast <- any(nzchar(vapply(rows, function(r) r$contrast %||% "", character(1))))

  # CI label tracks the actual confidence level the user passed to
  # mediate(). Hardcoded "95% CI" was wrong for users who set
  # conf.level = 0.99 etc -- the column label lied about the bounds.
  meta <- ms_mediate_meta(x)
  ci_label <- if (!is.na(meta$conf_level)) {
    sprintf("%d%% CI", round(meta$conf_level * 100))
  } else {
    "95% CI"
  }

  columns <- list(list(key = "effect", label = "Effect", format = "text"))
  if (has_contrast) {
    columns <- c(columns, list(list(key = "contrast", label = "Contrast", format = "text")))
  }
  columns <- c(columns, list(
    list(key = "estimate", label = "Estimate", format = "number"),
    list(key = "ci",       label = ci_label,   format = "ci"),
    list(key = "p_value",  label = "p",        format = "pvalue")
  ))

  fields <- list(
    table_type = "mediation",
    columns = columns,
    rows = rows,
    source = "mediation::mediate"
  )
  fields$treatment <- meta$treatment
  fields$mediator <- meta$mediator
  fields$outcome <- meta$outcome
  fields$n <- meta$n
  fields$sims <- meta$sims
  fields$boot <- meta$boot
  fields$ci_type <- meta$ci_type
  fields$conf_level <- meta$conf_level
  # Methodology note rendered below the table as an APA "Note." line.
  # Keeps provenance (bootstrap type, sim count, CI level) out of the
  # results prose -- it's methods, not findings. The JS renderer picks
  # this up automatically via tableProjectionNote().
  method_note <- ms_mediate_method_note(meta, rows = rows)
  if (nzchar(method_note)) fields$note <- method_note
  # Path coefficients (a, b, c') extracted from model.m and model.y on the
  # mediate object. mediate()'s own summary only reports the high-level
  # decomposition (ACME / ADE / Total / Prop mediated); a path diagram
  # needs the underlying regression coefficients. NULL when models or
  # variable names can't be recovered -- the figure renderer falls back to
  # showing the decomposition only.
  paths <- ms_mediate_paths(x, meta)
  if (!is.null(paths) && length(paths) > 0L) fields$paths <- paths

  # Controls / covariates -- any term in model.y / model.m that isn't on
  # the canonical triangle (treat, mediator, intercept, treat:mediator
  # interaction). Surfaces as an "Adjusted for: ..." caption on the
  # diagram and in the table note. Almost no published mediation runs
  # uncontrolled, so this is high-impact when present.
  controls <- ms_mediate_controls(x, meta)
  if (length(controls) > 0L) fields$controls <- controls
  fields
}

ms_mediate_rows <- function(x) {
  s <- ms_mediate_summary(x)
  meta <- ms_mediate_meta(x)
  ms_mediate_rows_one(s, meta)
}

ms_mediate_rows_one <- function(s, meta) {
  if (isTRUE(meta$interaction)) {
    specs <- list(
      list("ACME", "control", "d0", "d0.ci", "d0.p"),
      list("ACME", "treated", "d1", "d1.ci", "d1.p"),
      list("ADE", "control", "z0", "z0.ci", "z0.p"),
      list("ADE", "treated", "z1", "z1.ci", "z1.p"),
      list("Total effect", "", "tau.coef", "tau.ci", "tau.p"),
      list("Proportion mediated", "control", "n0", "n0.ci", "n0.p"),
      list("Proportion mediated", "treated", "n1", "n1.ci", "n1.p"),
      list("ACME", "average", "d.avg", "d.avg.ci", "d.avg.p"),
      list("ADE", "average", "z.avg", "z.avg.ci", "z.avg.p"),
      list("Proportion mediated", "average", "n.avg", "n.avg.ci", "n.avg.p")
    )
  } else {
    specs <- list(
      list("ACME", "", c("d.avg", "d0"), c("d.avg.ci", "d0.ci"), c("d.avg.p", "d0.p")),
      list("ADE", "", c("z.avg", "z0"), c("z.avg.ci", "z0.ci"), c("z.avg.p", "z0.p")),
      list("Total effect", "", "tau.coef", "tau.ci", "tau.p"),
      list("Proportion mediated", "", c("n.avg", "n0"), c("n.avg.ci", "n0.ci"), c("n.avg.p", "n0.p"))
    )
  }

  rows <- lapply(specs, function(sp) {
    ms_mediate_effect_row(
      s,
      effect = sp[[1]],
      contrast = sp[[2]],
      estimate_keys = sp[[3]],
      ci_keys = sp[[4]],
      p_keys = sp[[5]]
    )
  })
  Filter(function(row) !is.na(row$estimate), rows)
}

ms_mediate_effect_row <- function(s, effect, contrast,
                                  estimate_keys, ci_keys, p_keys) {
  ci <- ms_mediate_first_ci(s, ci_keys)
  row <- list(
    effect   = effect,
    contrast = contrast,
    estimate = ms_mediate_first_value(s, estimate_keys),
    ci_lower = ci[[1]],
    ci_upper = ci[[2]],
    p_value  = ms_mediate_first_value(s, p_keys)
  )
  row
}

ms_mediate_first_value <- function(x, keys) {
  for (key in keys) {
    value <- x[[key]]
    if (!is.null(value) && length(value) > 0L) {
      out <- ms_safe_numeric(value[[1]])
      if (!is.na(out)) return(out)
    }
  }
  NA_real_
}

ms_mediate_first_ci <- function(x, keys) {
  for (key in keys) {
    value <- x[[key]]
    if (!is.null(value) && length(value) >= 2L) {
      return(list(ms_safe_numeric(value[[1]]), ms_safe_numeric(value[[2]])))
    }
  }
  list(NA_real_, NA_real_)
}

ms_mediate_table_data <- function(x, pretty_names = FALSE) {
  rows <- ms_mediate_rows(x)
  has_contrast <- any(nzchar(vapply(rows, function(r) r$contrast %||% "", character(1))))

  cols <- "effect"
  if (has_contrast) cols <- c(cols, "contrast")
  cols <- c(cols, "estimate", "ci_lower", "ci_upper", "p_value")

  numeric_cols <- c("estimate", "ci_lower", "ci_upper", "p_value")
  df <- as.data.frame(
    setNames(lapply(cols, function(col) {
      if (col %in% numeric_cols) {
        return(vapply(rows, function(row) {
          value <- row[[col]]
          if (is.null(value)) NA_real_ else as.numeric(value)
        }, numeric(1)))
      }
      vapply(rows, function(row) {
        value <- row[[col]]
        if (is.null(value) || is.na(value)) NA_character_ else as.character(value)
      }, character(1))
    }), cols),
    stringsAsFactors = FALSE
  )

  if (pretty_names) {
    names(df) <- vapply(names(df), function(nm) {
      switch(nm,
        effect = "Effect",
        contrast = "Contrast",
        estimate = "Estimate",
        ci_lower = "Lower CI",
        ci_upper = "Upper CI",
        p_value = "p",
        nm
      )
    }, character(1), USE.NAMES = FALSE)
  }

  df
}

ms_mediate_raw_output <- function(x) {
  paste(utils::capture.output(print(ms_mediate_summary(x))), collapse = "\n")
}

# Compose the methodology "Note." line for the table footer. Surfaces:
#   - estimation method (bootstrap type or quasi-Bayesian Monte Carlo)
#   - simulation/resample count
#   - CI level
#   - N
#   - whether coefficients are standardized (currently always
#     unstandardized for mediate() output; left as a switch for when the
#     beta toggle ships in v1 figure work)
# Capitalisation note: the bootstrap-method label opens the sentence, so
# it's capitalised here ("Percentile bootstrap...", "Quasi-Bayesian...");
# the renderer prepends "Note. " upstream. Returns "" when meta is empty
# so the renderer can omit the note entirely.
ms_mediate_method_note <- function(meta, scale = "unstandardized", rows = NULL) {
  if (is.null(meta)) return("")
  bits <- character(0)
  sims_val <- ms_safe_numeric(meta$sims)
  if (isTRUE(meta$boot)) {
    ci_type <- tolower(as.character(meta$ci_type %||% ""))
    ci_desc <- switch(
      ci_type,
      "perc"  = "Percentile bootstrap",
      "bca"   = "BCa bootstrap",
      "norm"  = "Normal-approximation bootstrap",
      "basic" = "Basic bootstrap",
      "Nonparametric bootstrap"
    )
    if (!is.na(sims_val) && sims_val > 0) {
      bits <- c(bits, sprintf("%s with %d resamples", ci_desc,
                              as.integer(round(sims_val))))
    } else {
      bits <- c(bits, ci_desc)
    }
  } else if (!is.na(sims_val) && sims_val > 0) {
    bits <- c(bits,
              sprintf("Quasi-Bayesian Monte Carlo with %d simulations",
                      as.integer(round(sims_val))))
  }
  conf <- ms_safe_numeric(meta$conf_level)
  if (!is.na(conf)) {
    bits <- c(bits, sprintf("%d%% CIs", as.integer(round(conf * 100))))
  }
  n_val <- ms_safe_numeric(meta$n)
  if (!is.na(n_val)) {
    bits <- c(bits, sprintf("N = %d", as.integer(round(n_val))))
  }
  # Estimator type -- explicit, so a reader doesn't have to guess whether
  # the column is B or beta. Important once the beta toggle lands.
  has_prop <- ms_mediate_has_prop_mediated(rows)
  if (identical(scale, "standardized")) {
    bits <- c(bits, "Estimates are standardized (\u03b2)")
  } else if (isTRUE(has_prop)) {
    bits <- c(bits,
              "Effects are unstandardized (B) except proportion mediated, which is a ratio")
  } else {
    bits <- c(bits, "Estimates are unstandardized (B)")
  }
  if (isTRUE(has_prop) && ms_mediate_prop_caveat_needed(rows)) {
    bits <- c(bits,
              "proportion mediated can be unstable when the total effect is small and its CI may fall outside [0, 1]")
  }
  if (length(bits) == 0L) return("")
  paste0(paste(bits, collapse = "; "), ".")
}

ms_mediate_has_prop_mediated <- function(rows) {
  if (is.null(rows) || !length(rows)) return(FALSE)
  any(vapply(rows, function(row) {
    effect <- tolower(as.character(row$effect %||% ""))
    effect <- gsub("[._[:space:]-]+", ".", trimws(effect))
    effect %in% c("prop.mediated", "proportion.mediated")
  }, logical(1)))
}

ms_mediate_prop_caveat_needed <- function(rows) {
  if (is.null(rows) || !length(rows)) return(FALSE)
  norm_effect <- function(row) {
    effect <- tolower(as.character(row$effect %||% ""))
    gsub("[._[:space:]-]+", ".", trimws(effect))
  }
  total_rows <- Filter(function(row) {
    norm_effect(row) %in% c("total.effect", "total")
  }, rows)
  prop_rows <- Filter(function(row) {
    norm_effect(row) %in% c("prop.mediated", "proportion.mediated")
  }, rows)
  total_near_zero <- any(vapply(total_rows, function(row) {
    est <- ms_safe_numeric(row$estimate)
    lo <- ms_safe_numeric(row$ci_lower)
    hi <- ms_safe_numeric(row$ci_upper)
    if (!is.na(lo) && !is.na(hi) && min(lo, hi) <= 0 && max(lo, hi) >= 0) {
      return(TRUE)
    }
    !is.na(est) && abs(est) < 0.1
  }, logical(1)))
  prop_outside_unit <- any(vapply(prop_rows, function(row) {
    lo <- ms_safe_numeric(row$ci_lower)
    hi <- ms_safe_numeric(row$ci_upper)
    (!is.na(lo) && lo < 0) || (!is.na(hi) && hi > 1)
  }, logical(1)))
  total_near_zero || prop_outside_unit
}

ms_chr_or_na <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x[[1]])) return(NA_character_)
  as.character(x[[1]])
}

# ----------------------------------------------------------------------------
# Path coefficients for the mediation diagram
# ----------------------------------------------------------------------------
# A mediate object preserves both fitted models (model.m and model.y) but
# summary(mediate_obj) only reports the higher-level decomposition. The
# mediation diagram needs the per-arrow coefficients:
#   a  = effect of treat on M           (coef of treat in model.m)
#   b  = effect of M on Y, controlling for treat
#                                       (coef of mediator in model.y)
#   c' = direct effect of treat on Y, controlling for M
#                                       (coef of treat in model.y)
# We don't extract c (total) -- it's already the "Total effect" row from
# the decomposition. Returns a list with one entry per recovered path, or
# NULL when models / variable names aren't available.

ms_mediate_paths <- function(x, meta) {
  model_m <- tryCatch(x$model.m, error = function(e) NULL)
  model_y <- tryCatch(x$model.y, error = function(e) NULL)
  if (is.null(model_m) && is.null(model_y)) return(NULL)
  treat <- meta$treatment
  mediator <- meta$mediator
  if (is.na(treat) || is.na(mediator)) return(NULL)

  # Detect categorical treat (factor or binary numeric). beta = coef *
  # sd(predictor) / sd(outcome) is meaningless when "1 SD increase in X"
  # doesn't correspond to a real comparison -- for binary/factor X, the
  # coefficient already IS the across-level difference. Per Hayes &
  # Rockwood (2017) we omit beta for these cases rather than print a
  # number readers will misinterpret. The figure renderer / table note
  # surfaces "beta not shown for categorical predictor" so the omission is
  # explained instead of looking like missing data.
  treat_is_categorical <- ms_mediate_term_is_categorical(model_m, model_y, treat)

  # Pre-compute SDs from the model frame so each coef_row call doesn't
  # re-fetch. NA when the term isn't in that model's data.
  sd_treat <- ms_mediate_term_sd(model_m, treat)
  sd_mediator_in_y <- ms_mediate_term_sd(model_y, mediator)
  sd_outcome_m <- ms_mediate_response_sd(model_m)
  sd_outcome_y <- ms_mediate_response_sd(model_y)

  paths <- list()
  a <- ms_mediate_coef_row(model_m, treat, "a",
                           paste0(treat, " \u2192 ", mediator),
                           sd_predictor = sd_treat,
                           sd_outcome   = sd_outcome_m,
                           skip_std     = treat_is_categorical)
  if (!is.null(a)) paths[[length(paths) + 1L]] <- a
  b <- ms_mediate_coef_row(model_y, mediator, "b",
                           paste0(mediator, " \u2192 ", meta$outcome %||% "Y"),
                           sd_predictor = sd_mediator_in_y,
                           sd_outcome   = sd_outcome_y,
                           skip_std     = FALSE)
  if (!is.null(b)) paths[[length(paths) + 1L]] <- b
  cprime <- ms_mediate_coef_row(model_y, treat, "c_prime",
                                paste0(treat, " \u2192 ",
                                       (meta$outcome %||% "Y"),
                                       " (direct)"),
                                sd_predictor = sd_treat,
                                sd_outcome   = sd_outcome_y,
                                skip_std     = treat_is_categorical)
  if (!is.null(cprime)) paths[[length(paths) + 1L]] <- cprime

  if (length(paths) == 0L) NULL else paths
}

# Categorical detection used to decide whether to publish beta. Checks
# both model.m and model.y so a factor declared in either captures the
# limitation. Uses model$model (the model frame stored on lm/glm) when
# available; gracefully returns FALSE when introspection fails.
ms_mediate_term_is_categorical <- function(model_m, model_y, term) {
  for (model in list(model_m, model_y)) {
    if (is.null(model)) next
    mf <- tryCatch(stats::model.frame(model), error = function(e) NULL)
    if (is.null(mf) || !(term %in% names(mf))) next
    col <- mf[[term]]
    if (is.factor(col) || is.character(col) || is.logical(col)) return(TRUE)
    if (is.numeric(col)) {
      vals <- unique(col[is.finite(col)])
      if (length(vals) <= 2L) return(TRUE)  # binary numeric, e.g. 0/1
    }
  }
  FALSE
}

ms_mediate_term_sd <- function(model, term) {
  if (is.null(model) || is.null(term) || is.na(term) || !nzchar(term)) return(NA_real_)
  mf <- tryCatch(stats::model.frame(model), error = function(e) NULL)
  if (is.null(mf) || !(term %in% names(mf))) return(NA_real_)
  col <- mf[[term]]
  if (!is.numeric(col)) return(NA_real_)
  ms_safe_numeric(stats::sd(col, na.rm = TRUE))
}

ms_mediate_response_sd <- function(model) {
  if (is.null(model)) return(NA_real_)
  mf <- tryCatch(stats::model.frame(model), error = function(e) NULL)
  if (is.null(mf) || ncol(mf) < 1L) return(NA_real_)
  resp <- mf[[1L]]   # model.frame puts the response in column 1
  if (!is.numeric(resp)) return(NA_real_)
  ms_safe_numeric(stats::sd(resp, na.rm = TRUE))
}

# Extract a deduplicated list of control / covariate terms from the
# mediate object's underlying models. A term is a "control" when it
# appears in coef(model.m) or coef(model.y) and is NOT:
#   - "(Intercept)"
#   - the treatment variable
#   - the mediator variable
#   - the treatment:mediator interaction (or its reverse)
# Returns a character vector ordered by appearance. The JS renderer's
# mediationDiagramControlsList() strips I() wrappers cosmetically -- we
# return them as R wrote them so the audit trail is intact.
ms_mediate_controls <- function(x, meta) {
  treat <- meta$treatment
  mediator <- meta$mediator
  if (is.na(treat) || is.na(mediator)) return(character(0))

  # Drop the canonical triangle. For factor predictors R prefixes the
  # coef name (e.g. treat="X" with X a factor becomes Xtreated), so we
  # also drop anything that starts with the treat or mediator name --
  # otherwise factor-encoded treat would leak into the controls list.
  drop_set <- c("(Intercept)", treat, mediator,
                paste0(treat, ":", mediator),
                paste0(mediator, ":", treat))
  is_triangle <- function(name) {
    if (name %in% drop_set) return(TRUE)
    if (startsWith(name, treat))    return(TRUE)
    if (startsWith(name, mediator)) return(TRUE)
    FALSE
  }
  collect <- character(0)
  for (slot in c("model.m", "model.y")) {
    model <- tryCatch(x[[slot]], error = function(e) NULL)
    if (is.null(model)) next
    co <- tryCatch(names(stats::coef(model)), error = function(e) NULL)
    if (is.null(co)) next
    for (term in co) {
      if (is_triangle(term)) next
      if (!(term %in% collect)) collect <- c(collect, term)
    }
  }
  collect
}

# Pull a single coefficient row out of an lm/glm fit. Robust to both
# Pr(>|t|) (lm) and Pr(>|z|) (glm) column naming. Returns NULL when the
# requested term isn't in the model's coefficient table (e.g. factor
# treatment whose coef name doesn't match the bare variable name).
#
# Standardized estimate (beta) is the raw-SD form: beta = B * sd(X) / sd(Y).
# Matches the lavaan std.all convention for simple paths and is what
# psych/methods journals print for mediation diagrams. We omit it
# (std_estimate stays NA, std_estimate_skipped is TRUE) when the
# predictor is categorical -- see ms_mediate_term_is_categorical for
# why ("beta" for a 0/1 treat just rescales the binary contrast, which
# readers misinterpret as a real SD-units effect).
ms_mediate_coef_row <- function(model, term, path_id, label,
                                sd_predictor = NA_real_,
                                sd_outcome   = NA_real_,
                                skip_std     = FALSE) {
  if (is.null(model) || is.null(term) || is.na(term) || !nzchar(term)) {
    return(NULL)
  }
  co <- tryCatch(stats::coef(model), error = function(e) NULL)
  if (is.null(co)) return(NULL)

  # Resolve the coefficient name. Continuous / binary-numeric predictors
  # appear under the bare variable name (`X`). Factor predictors get
  # prefixed by R: a 2-level factor becomes `Xtreated`, a 3-level factor
  # becomes `Xmed` / `Xhigh` (one row per non-reference level). We try
  # bare name first, then prefix-match. For multi-level factors we pick
  # the first non-reference level -- multi-level mediation diagrams are
  # methodologically awkward (no convention surfaces it cleanly), but
  # showing one comparison is more useful than nothing as long as the
  # caption is honest about which level we drew.
  coef_name <- NULL
  level_used <- NA_character_
  multi_level <- FALSE
  if (term %in% names(co)) {
    coef_name <- term
  } else {
    matches <- names(co)[startsWith(names(co), term)]
    if (length(matches) == 0L) return(NULL)
    coef_name <- matches[1L]
    level_used <- substr(coef_name, nchar(term) + 1L, nchar(coef_name))
    if (length(matches) > 1L) multi_level <- TRUE
  }

  est <- ms_safe_numeric(co[[coef_name]])
  term <- coef_name  # downstream summary / confint lookups use the resolved name
  se <- NA_real_
  stat <- NA_real_
  pval <- NA_real_
  sm <- tryCatch(summary(model), error = function(e) NULL)
  if (!is.null(sm) && !is.null(sm$coefficients) &&
      term %in% rownames(sm$coefficients)) {
    sm_row <- sm$coefficients[term, , drop = TRUE]
    sm_cols <- names(sm_row)
    if ("Std. Error" %in% sm_cols) se <- ms_safe_numeric(sm_row[["Std. Error"]])
    if ("t value" %in% sm_cols) {
      stat <- ms_safe_numeric(sm_row[["t value"]])
    } else if ("z value" %in% sm_cols) {
      stat <- ms_safe_numeric(sm_row[["z value"]])
    }
    if ("Pr(>|t|)" %in% sm_cols) {
      pval <- ms_safe_numeric(sm_row[["Pr(>|t|)"]])
    } else if ("Pr(>|z|)" %in% sm_cols) {
      pval <- ms_safe_numeric(sm_row[["Pr(>|z|)"]])
    }
  }

  ci_bounds <- tryCatch(
    suppressMessages(stats::confint(model, parm = term)),
    error = function(e) NULL
  )
  ci_lower <- NA_real_
  ci_upper <- NA_real_
  if (!is.null(ci_bounds) && length(dim(ci_bounds)) == 2L &&
      nrow(ci_bounds) >= 1L && ncol(ci_bounds) >= 2L) {
    ci_lower <- ms_safe_numeric(ci_bounds[1L, 1L])
    ci_upper <- ms_safe_numeric(ci_bounds[1L, 2L])
  }

  std_estimate <- NA_real_
  std_skipped <- FALSE
  if (isTRUE(skip_std)) {
    std_skipped <- TRUE
  } else if (!is.na(est) && !is.na(sd_predictor) && !is.na(sd_outcome) &&
             sd_outcome > 0) {
    std_estimate <- ms_safe_numeric(est * sd_predictor / sd_outcome)
  }

  row <- list(
    path      = path_id,
    label     = label,
    term      = as.character(term),
    estimate  = est,
    std_error = se,
    statistic = stat,
    p_value   = pval,
    ci_lower  = ci_lower,
    ci_upper  = ci_upper
  )
  if (!is.na(std_estimate)) row$std_estimate <- std_estimate
  if (isTRUE(std_skipped))  row$std_estimate_skipped <- TRUE
  # Factor-treatment metadata so the figure caption can be honest about
  # which contrast it's showing. `level_used` is set when the bare
  # variable name didn't match coef() and we fell back to a prefix
  # search (R's standard factor dummy-coding naming). `multi_level` is
  # an additional flag when the factor has 3+ levels -- in that case
  # the diagram only shows ONE non-reference contrast and the caption
  # should say so.
  if (!is.na(level_used) && nzchar(level_used)) {
    row$level_used <- level_used
  }
  if (isTRUE(multi_level)) row$multi_level <- TRUE
  row
}
