# Extractors for the `interactions` package: johnson_neyman() region-of-
# significance results and sim_slopes() simple-slopes analyses.
#
# The Johnson-Neyman figure ships as five analytic parameters rather than
# the raw 1000-row cbands grid: the conditional slope is linear in the
# moderator (s0 + s1*m) and the confidence-band half-width squared is
# quadratic (w0 + w1*m + w2*m^2). Both are recovered from cbands itself,
# so they are exact under whatever vcov / FDR adjustment the user chose —
# we never recompute from the model. A thinned curve ships alongside as a
# fallback and parity check; the web app prefers the analytic form.

#' @rdname mellio_payload
#' @export
mellio_payload.johnson_neyman <- function(x, ..., .call = NULL) {
  call_str <- ms_interactions_call(.call, "interactions::johnson_neyman(...)")
  figure <- ms_jn_figure_data(x)
  ms_build_envelope(
    type = "johnson_neyman",
    type_label = "Johnson-Neyman interval",
    call = call_str,
    fields = ms_jn_fields(x),
    raw_output = ms_jn_raw_output(x),
    packages = ms_packages_basic(extras = "interactions"),
    card_kind = "table",
    figure_data = if (!is.null(figure)) list(johnson_neyman_plot = figure) else NULL
  )
}

#' @rdname mellio_payload
#' @export
mellio_payload.sim_slopes <- function(x, ..., .call = NULL) {
  call_str <- ms_interactions_call(.call, "interactions::sim_slopes(...)")

  pred <- as.character(attr(x, "pred") %||% NA_character_)
  modx <- as.character(attr(x, "modx") %||% NA_character_)
  mod2 <- attr(x, "mod2")
  resp <- attr(x, "resp")
  ci_width <- ms_safe_numeric(attr(x, "ci.width") %||% 0.95)

  parsed <- ms_sim_slopes_rows(x)
  fields <- list(
    table_type = "simple_slopes",
    source = "interactions",
    pred = pred,
    modx = modx,
    rows = parsed$rows,
    columns = parsed$columns,
    statistic_label = parsed$statistic_label,
    ci_width = ci_width
  )
  if (!is.null(mod2)) fields$mod2 <- as.character(mod2)
  if (!is.null(resp)) fields$outcome <- as.character(resp)

  note <- ms_sim_slopes_note(x, pred, modx, resp)

  # Embed the Johnson-Neyman figure when sim_slopes probed it (the default
  # for continuous moderators). Three-way probes carry one J-N object per
  # second-moderator level; a single figure would misrepresent them, so we
  # only attach the figure for two-way probes.
  figure <- NULL
  jn_list <- x$jn %||% attr(x, "jns")
  if (is.null(mod2) && is.list(jn_list) && length(jn_list) >= 1 &&
      inherits(jn_list[[1]], "johnson_neyman")) {
    figure <- ms_jn_figure_data(jn_list[[1]])
    jn_note <- ms_jn_note(jn_list[[1]])
    if (nzchar(jn_note)) note <- paste(note, jn_note)
  }
  if (nzchar(note)) fields$note <- trimws(note)

  ms_build_envelope(
    type = "simple_slopes",
    type_label = "Simple slopes analysis",
    call = call_str,
    fields = fields,
    raw_output = ms_capture_output(x),
    packages = ms_packages_basic(extras = "interactions"),
    card_kind = "table",
    figure_data = if (!is.null(figure)) list(johnson_neyman_plot = figure) else NULL
  )
}

#' @export
mellio_open_dispatch.johnson_neyman <- function(x, browse = TRUE, ..., .call = NULL) {
  send_payload_to_stats(mellio_payload(x, ..., .call = .call), browse = browse)
}

#' @export
mellio_open_dispatch.sim_slopes <- function(x, browse = TRUE, ..., .call = NULL) {
  send_payload_to_stats(mellio_payload(x, ..., .call = .call), browse = browse)
}

# ---------------------------------------------------------------------------
# johnson_neyman internals
# ---------------------------------------------------------------------------

ms_interactions_call <- function(.call, fallback) {
  out <- .call %||% fallback
  trimws(gsub("\\s+", " ", out))
}

# print.johnson_neyman draws the stored ggplot when attr(x, "plot") is TRUE,
# which opens a graphics device under R CMD check. Capture text only.
ms_jn_raw_output <- function(x) {
  attr(x, "plot") <- FALSE
  ms_capture_output(x)
}

ms_jn_columns <- function(x) {
  cb <- x$cbands
  if (!is.data.frame(cb) || nrow(cb) < 3) return(NULL)
  modx_col <- names(cb)[[1]]
  slope_col <- paste("Slope of", as.character(attr(x, "pred")))
  if (!slope_col %in% names(cb)) slope_col <- names(cb)[[2]]
  if (!all(c(slope_col, "Lower", "Upper") %in% names(cb))) return(NULL)
  list(cb = cb, modx_col = modx_col, slope_col = slope_col)
}

# Recover the analytic form from cbands: slope(m) = s0 + s1*m from the
# endpoints, half-width(m)^2 = w0 + w1*m + w2*m^2 from three points. Only
# trusted when it reproduces every cbands row to numerical precision.
ms_jn_analytic <- function(x) {
  cols <- ms_jn_columns(x)
  if (is.null(cols)) return(NULL)
  cb <- cols$cb
  m <- as.numeric(cb[[cols$modx_col]])
  slope <- as.numeric(cb[[cols$slope_col]])
  half <- (as.numeric(cb$Upper) - as.numeric(cb$Lower)) / 2
  keep <- is.finite(m) & is.finite(slope) & is.finite(half)
  m <- m[keep]; slope <- slope[keep]; half <- half[keep]
  n <- length(m)
  if (n < 3 || m[1] == m[n]) return(NULL)

  s1 <- (slope[n] - slope[1]) / (m[n] - m[1])
  s0 <- slope[1] - s1 * m[1]
  idx <- unique(c(1L, as.integer(ceiling(n / 2)), n))
  if (length(idx) < 3) return(NULL)
  w <- tryCatch(
    solve(cbind(1, m[idx], m[idx]^2), half[idx]^2),
    error = function(e) NULL
  )
  if (is.null(w) || any(!is.finite(w))) return(NULL)

  slope_hat <- s0 + s1 * m
  half_hat <- sqrt(pmax(w[1] + w[2] * m + w[3] * m^2, 0))
  tol <- 1e-6 * max(1, max(abs(slope), na.rm = TRUE), max(abs(half), na.rm = TRUE))
  if (max(abs(slope_hat - slope)) > tol || max(abs(half_hat - half)) > tol) {
    return(NULL)
  }
  list(s0 = s0, s1 = s1, w0 = w[1], w1 = w[2], w2 = w[3])
}

# Thin cbands to at most `points` rows: the client-side fallback when the
# analytic recovery is unavailable, and the parity reference in tests.
ms_jn_curve <- function(x, points = 81L) {
  cols <- ms_jn_columns(x)
  if (is.null(cols)) return(NULL)
  cb <- cols$cb
  n <- nrow(cb)
  idx <- unique(as.integer(round(seq(1L, n, length.out = min(points, n)))))
  lapply(idx, function(i) {
    list(
      m = ms_safe_numeric(as.numeric(cb[[cols$modx_col]][i])),
      slope = ms_safe_numeric(as.numeric(cb[[cols$slope_col]][i])),
      ci_lower = ms_safe_numeric(as.numeric(cb$Lower[i])),
      ci_upper = ms_safe_numeric(as.numeric(cb$Upper[i]))
    )
  })
}

ms_jn_bound <- function(x, which) {
  bounds <- x$bounds
  if (is.null(bounds) || length(bounds) < 2) return(NA_real_)
  ms_safe_numeric(as.numeric(bounds[[which]]))
}

ms_jn_format <- function(value, digits = 2) {
  if (is.null(value) || !is.finite(value)) return("NA")
  format(round(value, digits), trim = TRUE, scientific = FALSE, nsmall = 0)
}

ms_jn_note <- function(x) {
  pred <- as.character(attr(x, "pred") %||% "the predictor")
  modx <- as.character(attr(x, "modx") %||% "the moderator")
  alpha <- ms_safe_numeric(attr(x, "alpha") %||% 0.05)
  inside <- isTRUE(attr(x, "inside"))
  all_sig <- isTRUE(attr(x, "all_sig"))
  failed <- isTRUE(attr(x, "failed"))
  modrange <- suppressWarnings(as.numeric(attr(x, "modrange")))
  lo <- ms_jn_bound(x, 1L)
  hi <- ms_jn_bound(x, 2L)

  if (failed) {
    return(sprintf(
      "No Johnson-Neyman interval could be found: the significance of the slope of %s does not change across the range of %s examined.",
      pred, modx))
  }
  parts <- character(0)
  if (all_sig) {
    parts <- c(parts, sprintf(
      "The slope of %s is significant (p < %s) across the entire observed range of %s.",
      pred, ms_jn_format(alpha), modx))
  } else {
    parts <- c(parts, sprintf(
      "When %s is %s the interval [%s, %s], the slope of %s is significant at p < %s.",
      modx, if (inside) "inside" else "outside",
      ms_jn_format(lo), ms_jn_format(hi), pred, ms_jn_format(alpha)))
  }
  if (length(modrange) == 2 && all(is.finite(modrange))) {
    parts <- c(parts, sprintf(
      "Observed range of %s: [%s, %s].",
      modx, ms_jn_format(modrange[1]), ms_jn_format(modrange[2])))
    if (!all_sig && is.finite(lo) && is.finite(hi) &&
        (hi < modrange[1] || lo > modrange[2])) {
      parts <- c(parts, "Both interval bounds fall outside the observed range of the moderator.")
    }
  }
  if (isTRUE(attr(x, "control.fdr"))) {
    tcrit <- ms_safe_numeric(if (length(x$t_value)) as.numeric(x$t_value[[1]]) else NA_real_)
    parts <- c(parts, if (is.finite(tcrit)) {
      sprintf("Interval computed with a false discovery rate adjusted critical t of %s.",
              ms_jn_format(tcrit))
    } else {
      "Interval computed with a false discovery rate adjusted critical t."
    })
  }
  paste(parts, collapse = " ")
}

ms_jn_fields <- function(x) {
  pred <- as.character(attr(x, "pred") %||% NA_character_)
  modx <- as.character(attr(x, "modx") %||% NA_character_)
  alpha <- ms_safe_numeric(attr(x, "alpha") %||% 0.05)
  modrange <- suppressWarnings(as.numeric(attr(x, "modrange")))
  if (length(modrange) != 2) modrange <- c(NA_real_, NA_real_)

  fields <- list(
    table_type = "johnson_neyman",
    source = "interactions",
    pred = pred,
    modx = modx,
    alpha = alpha,
    control_fdr = isTRUE(attr(x, "control.fdr")),
    inside = isTRUE(attr(x, "inside")),
    all_sig = isTRUE(attr(x, "all_sig")),
    failed = isTRUE(attr(x, "failed")),
    bound_lower = ms_jn_bound(x, 1L),
    bound_upper = ms_jn_bound(x, 2L),
    modrange = I(ms_safe_numeric(modrange))
  )
  tcrit <- if (length(x$t_value)) ms_safe_numeric(as.numeric(x$t_value[[1]])) else NA_real_
  if (is.finite(tcrit)) fields$critical_t <- tcrit

  table_data <- ms_jn_key_rows(x)
  fields$rows <- table_data$rows
  fields$columns <- table_data$columns
  note <- ms_jn_note(x)
  if (nzchar(note)) fields$note <- note
  fields
}

# Key-point table: conditional slope + CI at the observed extremes of the
# moderator and at the J-N bounds, evaluated from the analytic form so the
# values match cbands exactly.
ms_jn_key_rows <- function(x) {
  pred <- as.character(attr(x, "pred") %||% "predictor")
  modx <- as.character(attr(x, "modx") %||% "moderator")
  alpha <- ms_safe_numeric(attr(x, "alpha") %||% 0.05)
  inside <- isTRUE(attr(x, "inside"))
  failed <- isTRUE(attr(x, "failed"))
  all_sig <- isTRUE(attr(x, "all_sig"))
  analytic <- ms_jn_analytic(x)
  modrange <- suppressWarnings(as.numeric(attr(x, "modrange")))
  lo <- ms_jn_bound(x, 1L)
  hi <- ms_jn_bound(x, 2L)

  eval_at <- function(m) {
    if (is.null(analytic) || !is.finite(m)) {
      return(list(slope = NA_real_, lower = NA_real_, upper = NA_real_))
    }
    slope <- analytic$s0 + analytic$s1 * m
    half <- sqrt(max(analytic$w0 + analytic$w1 * m + analytic$w2 * m^2, 0))
    list(slope = slope, lower = slope - half, upper = slope + half)
  }
  significance_at <- function(m) {
    if (failed || !is.finite(m)) return(NA_character_)
    if (all_sig) return("Significant")
    if (!is.finite(lo) || !is.finite(hi)) return(NA_character_)
    in_interval <- m >= lo & m <= hi
    if (identical(inside, in_interval)) "Significant" else "n.s."
  }

  points <- list()
  if (length(modrange) == 2 && is.finite(modrange[1])) {
    points <- c(points, list(list(label = "Min. observed", m = modrange[1], boundary = FALSE)))
  }
  if (!failed && is.finite(lo)) {
    points <- c(points, list(list(label = "J-N lower bound", m = lo, boundary = TRUE)))
  }
  if (!failed && is.finite(hi)) {
    points <- c(points, list(list(label = "J-N upper bound", m = hi, boundary = TRUE)))
  }
  if (length(modrange) == 2 && is.finite(modrange[2])) {
    points <- c(points, list(list(label = "Max. observed", m = modrange[2], boundary = FALSE)))
  }
  points <- points[order(vapply(points, function(p) p$m, numeric(1)))]

  rows <- lapply(points, function(p) {
    est <- eval_at(p$m)
    # Display rows only — the figure keeps full precision in figure_data.
    # The web app's "ci" cell format pastes values verbatim, so round here.
    list(
      point = p$label,
      modx_value = ms_safe_numeric(round(p$m, 2)),
      estimate = ms_safe_numeric(round(est$slope, 2)),
      ci_lower = ms_safe_numeric(round(est$lower, 2)),
      ci_upper = ms_safe_numeric(round(est$upper, 2)),
      significance = if (p$boundary) {
        paste0("boundary (p = ", ms_jn_format(alpha), ")")
      } else {
        significance_at(p$m)
      }
    )
  })
  rows <- Filter(function(row) is.finite(row$modx_value %||% NA_real_), rows)

  columns <- list(
    list(key = "point", label = "Point", format = "text"),
    list(key = "modx_value", label = modx, format = "number"),
    list(key = "estimate", label = paste("Slope of", pred), format = "number"),
    list(key = "ci", label = paste0(ms_jn_format((1 - alpha) * 100, 0), "% CI"), format = "ci"),
    list(key = "significance", label = "Significance", format = "text")
  )
  list(rows = rows, columns = columns)
}

ms_jn_figure_data <- function(x) {
  cols <- ms_jn_columns(x)
  if (is.null(cols)) return(NULL)
  pred <- as.character(attr(x, "pred") %||% NA_character_)
  modx <- as.character(attr(x, "modx") %||% NA_character_)
  modrange <- suppressWarnings(as.numeric(attr(x, "modrange")))
  if (length(modrange) != 2) modrange <- c(NA_real_, NA_real_)
  grid_m <- as.numeric(cols$cb[[cols$modx_col]])
  grid_range <- range(grid_m[is.finite(grid_m)])
  analytic <- ms_jn_analytic(x)
  curve <- ms_jn_curve(x)
  if (is.null(curve) || length(curve) < 2) return(NULL)

  out <- list(
    source = "interactions",
    pred = list(variable = pred, label = pred),
    moderator = list(variable = modx, label = modx),
    alpha = ms_safe_numeric(attr(x, "alpha") %||% 0.05),
    control_fdr = isTRUE(attr(x, "control.fdr")),
    inside = isTRUE(attr(x, "inside")),
    all_sig = isTRUE(attr(x, "all_sig")),
    failed = isTRUE(attr(x, "failed")),
    bounds = list(lower = ms_jn_bound(x, 1L), upper = ms_jn_bound(x, 2L)),
    modrange = I(ms_safe_numeric(modrange)),
    grid_range = I(ms_safe_numeric(grid_range)),
    curve = curve,
    x_label = modx,
    y_label = paste("Slope of", pred)
  )
  if (!is.null(analytic)) {
    out$analytic <- list(
      s0 = ms_safe_numeric(analytic$s0),
      s1 = ms_safe_numeric(analytic$s1),
      w0 = ms_safe_numeric(analytic$w0),
      w1 = ms_safe_numeric(analytic$w1),
      w2 = ms_safe_numeric(analytic$w2)
    )
  }
  tcrit <- if (length(x$t_value)) ms_safe_numeric(as.numeric(x$t_value[[1]])) else NA_real_
  if (is.finite(tcrit)) out$critical_t <- tcrit
  out
}

# ---------------------------------------------------------------------------
# sim_slopes internals
# ---------------------------------------------------------------------------

ms_sim_slopes_rows <- function(x) {
  modx <- as.character(attr(x, "modx") %||% "moderator")
  pred <- as.character(attr(x, "pred") %||% "predictor")
  mod2 <- attr(x, "mod2")
  ci_width <- ms_safe_numeric(attr(x, "ci.width") %||% 0.95)

  tables <- if (is.data.frame(x$slopes)) list(x$slopes) else as.list(x$slopes)
  mod2_labels <- ms_sim_slopes_level_labels(attr(x, "mod2.values"))
  modx_labels <- ms_sim_slopes_level_labels(attr(x, "modx.values"))

  statistic_label <- "t"
  rows <- list()
  for (ti in seq_along(tables)) {
    df <- tables[[ti]]
    if (!is.data.frame(df) || nrow(df) == 0) next
    value_col <- grep("^Value of ", names(df), value = TRUE)[1]
    ci_cols <- grep("%$", names(df), value = TRUE)
    stat_col <- grep("val\\.$", names(df), value = TRUE)[1]
    if (is.na(value_col) || length(ci_cols) < 2 || is.na(stat_col)) next
    statistic_label <- sub(" val\\.$", "", stat_col)

    for (ri in seq_len(nrow(df))) {
      value <- ms_safe_numeric(as.numeric(df[[value_col]][ri]))
      row <- list(
        modx_value = value,
        estimate = ms_safe_numeric(as.numeric(df[["Est."]][ri])),
        se = ms_safe_numeric(as.numeric(df[["S.E."]][ri])),
        # The web app's "ci" cell format pastes values verbatim — round for display.
        ci_lower = ms_safe_numeric(round(as.numeric(df[[ci_cols[1]]][ri]), 2)),
        ci_upper = ms_safe_numeric(round(as.numeric(df[[ci_cols[2]]][ri]), 2)),
        statistic = ms_safe_numeric(as.numeric(df[[stat_col]][ri])),
        p_value = ms_safe_numeric(as.numeric(df[["p"]][ri]))
      )
      level <- ms_sim_slopes_match_label(value, modx_labels)
      if (!is.na(level)) row$level <- level
      if (!is.null(mod2)) {
        row$mod2_label <- if (ti <= length(mod2_labels$labels)) {
          mod2_labels$labels[[ti]]
        } else {
          paste("Level", ti)
        }
      }
      rows <- c(rows, list(row))
    }
  }
  rows <- Filter(function(row) is.finite(row$estimate %||% NA_real_), rows)

  columns <- list()
  if (!is.null(mod2)) {
    columns <- c(columns, list(list(key = "mod2_label", label = as.character(mod2), format = "text")))
  }
  has_levels <- any(vapply(rows, function(row) !is.null(row$level), logical(1)))
  if (has_levels) {
    columns <- c(columns, list(list(key = "level", label = "Level", format = "text")))
  }
  columns <- c(columns, list(
    list(key = "modx_value", label = paste("Value of", modx), format = "number"),
    list(key = "estimate", label = paste("Slope of", pred), format = "number"),
    list(key = "se", label = "SE", format = "number"),
    list(key = "ci", label = paste0(ms_jn_format(ci_width * 100, 0), "% CI"), format = "ci"),
    list(key = "statistic", label = statistic_label, format = "statistic"),
    list(key = "p_value", label = "p", format = "pvalue")
  ))
  list(rows = rows, columns = columns, statistic_label = statistic_label)
}

# attr(x, "modx.values") is a named numeric when interactions picked the
# default probe values ("- 1 SD" / "Mean" / "+ 1 SD"); unnamed when the
# user supplied modx.values directly.
ms_sim_slopes_level_labels <- function(values) {
  if (is.null(values)) return(list(values = numeric(0), labels = character(0)))
  vals <- suppressWarnings(as.numeric(values))
  labels <- names(values)
  if (is.null(labels)) labels <- rep(NA_character_, length(vals))
  list(values = vals, labels = labels)
}

ms_sim_slopes_match_label <- function(value, level_labels) {
  if (!is.finite(value) || length(level_labels$values) == 0) return(NA_character_)
  diffs <- abs(level_labels$values - value)
  idx <- which.min(diffs)
  if (length(idx) != 1 || !is.finite(diffs[idx])) return(NA_character_)
  tol <- 1e-8 * max(1, abs(value))
  if (diffs[idx] > tol) return(NA_character_)
  label <- level_labels$labels[idx]
  if (is.na(label) || !nzchar(label)) NA_character_ else label
}

ms_sim_slopes_note <- function(x, pred, modx, resp) {
  parts <- sprintf("Conditional slopes of %s at selected values of %s.", pred, modx)
  if (!is.null(resp) && nzchar(resp)) {
    parts <- sprintf("Conditional slopes of %s on %s at selected values of %s.",
                     pred, as.character(resp), modx)
  }
  robust <- attr(x, "robust")
  if (is.character(robust) && nzchar(robust)) {
    parts <- paste(parts, sprintf("Robust standard errors (%s).", robust))
  } else if (isTRUE(robust)) {
    parts <- paste(parts, "Robust standard errors.")
  }
  parts
}
