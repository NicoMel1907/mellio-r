# R bridge -- effectsize package tables.

#' @rdname mellio_payload
#' @export
mellio_payload.effectsize_table <- function(x, ..., .call = NULL) {
  fields <- ms_effectsize_fields(x, source = "effectsize")

  call_str <- if (!is.null(.call)) .call else {
    paste(deparse(match.call()$x, width.cutoff = 500L), collapse = " ")
  }

  ms_build_envelope(
    type = "effect_size",
    type_label = "Effect size",
    call = trimws(gsub("\\s+", " ", call_str)),
    fields = fields,
    raw_output = ms_capture_output(x),
    packages = ms_packages_basic(extras = "effectsize"),
    card_kind = "table"
  )
}

ms_detect_effectsize_data_frame <- function(x) {
  if (!is.data.frame(x) || nrow(x) == 0L) return(NULL)
  lower <- ms_df_names(x)
  est_col <- ms_effectsize_estimate_col(x)
  if (is.na(est_col)) return(NULL)
  if (!lower[[est_col]] %in% ms_effectsize_estimate_names()) return(NULL)

  fields <- ms_effectsize_fields(x, source = "custom_data_frame")
  list(
    type = "effect_size",
    type_label = "Effect size",
    fields = fields
  )
}

ms_effectsize_fields <- function(x, source = "effectsize") {
  df <- as.data.frame(x)
  if (!nrow(df)) {
    stop("effectsize table has no rows to extract.", call. = FALSE)
  }

  est_col <- ms_effectsize_estimate_col(df)
  if (is.na(est_col)) {
    stop("Could not identify the effect-size estimate column.", call. = FALSE)
  }

  lower <- ms_df_names(df)
  variable_col <- ms_effectsize_variable_col(df, lower)
  low_col <- ms_df_col(lower, ms_df_low_cols())
  high_col <- ms_df_col(lower, ms_df_high_cols())
  ci_col <- ms_df_col(lower, c("ci", "ci_level", "confidence", "confidence_level"))
  mag_col <- ms_df_col(lower, c("magnitude", "interpretation", "effect_magnitude"))
  standardizer_col <- ms_df_col(lower, c("standardizer", "standardiser", "sd", "denominator"))

  effect_key <- names(df)[[est_col]]
  effect_label <- ms_effectsize_label(effect_key)

  rows <- lapply(seq_len(nrow(df)), function(i) {
    row <- list(
      variable = if (!is.na(variable_col)) {
        ms_df_text_value(df[[variable_col]][[i]])
      } else {
        "Effect size"
      },
      effect = effect_label,
      effect_key = effect_key,
      estimate = ms_safe_numeric(df[[est_col]][[i]])
    )
    if (!is.na(low_col) && !is.na(high_col)) {
      row$ci_lower <- ms_safe_numeric(df[[low_col]][[i]])
      row$ci_upper <- ms_safe_numeric(df[[high_col]][[i]])
    }
    if (!is.na(ci_col)) row$ci_level <- ms_safe_numeric(df[[ci_col]][[i]])
    if (!is.na(mag_col)) row$magnitude <- ms_df_text_value(df[[mag_col]][[i]])
    if (!is.na(standardizer_col)) row$standardizer <- ms_df_text_value(df[[standardizer_col]][[i]])
    row
  })

  columns <- list(
    list(key = "variable", label = "Term", format = "text"),
    list(key = "effect", label = "Effect", format = "text"),
    list(key = "estimate", label = "Estimate", format = "number")
  )
  if (!is.na(low_col) && !is.na(high_col)) {
    columns <- c(columns, list(
      list(key = "ci_lower", label = "95% CI lower", format = "number"),
      list(key = "ci_upper", label = "95% CI upper", format = "number")
    ))
  }
  if (!is.na(mag_col)) {
    columns <- c(columns, list(list(key = "magnitude", label = "Magnitude", format = "text")))
  }

  list(
    table_type = "effect_sizes",
    statistic_name = effect_label,
    method = if (identical(source, "effectsize")) "effectsize" else "effect size table",
    columns = columns,
    rows = rows,
    n_effects = length(rows),
    source = source
  )
}

#' @rdname mellio_payload
#' @export
mellio_payload.effectsize_difference <- function(x, ..., .call = NULL) {
  mellio_payload.effectsize_table(x, ..., .call = .call)
}

#' @rdname mellio_payload
#' @export
mellio_payload.effectsize_anova <- function(x, ..., .call = NULL) {
  mellio_payload.effectsize_table(x, ..., .call = .call)
}

#' @rdname mellio_payload
#' @export
mellio_payload.see_effectsize_table <- function(x, ..., .call = NULL) {
  mellio_payload.effectsize_table(x, ..., .call = .call)
}

ms_effectsize_estimate_col <- function(df) {
  lower <- ms_df_names(df)
  excluded <- c(
    "ci", "ci_level", "ci_low", "ci_lower", "ci_high", "ci_upper",
    "parameter", "term", "effect", "group", "group1", "group2",
    "magnitude", "interpretation", "standardizer", "standardiser"
  )
  candidates <- which(vapply(df, is.numeric, logical(1)) & !lower %in% excluded)
  if (!length(candidates)) return(NA_integer_)
  preferred <- ms_df_col(lower, ms_effectsize_estimate_names())
  if (!is.na(preferred) && preferred %in% candidates) return(preferred)
  candidates[[1]]
}

ms_effectsize_estimate_names <- function() {
  c(
    "cohens_d", "hedges_g", "glass_delta", "eta2", "eta2_partial",
    "omega2", "omega2_partial", "epsilon2", "r2", "r", "cramers_v",
    "odds_ratio"
  )
}

ms_effectsize_variable_col <- function(df, lower) {
  col <- ms_df_col(lower, c(
    "parameter", "term", "effect", "contrast", "comparison",
    "variable", "group", "group1"
  ))
  if (!is.na(col)) return(col)
  NA_integer_
}

ms_effectsize_label <- function(key) {
  normal <- tolower(gsub("[^A-Za-z0-9]+", "_", key))
  labels <- list(
    cohens_d = "Cohen's d",
    hedges_g = "Hedges' g",
    glass_delta = "Glass' delta",
    eta2 = "\u03b7\u00b2",
    eta2_partial = "\u03b7\u00b2\u209a",
    omega2 = "\u03c9\u00b2",
    omega2_partial = "\u03c9\u00b2\u209a",
    epsilon2 = "\u03b5\u00b2",
    cramers_v = "Cramer's V",
    odds_ratio = "OR",
    r2 = "R\u00b2"
  )
  labels[[normal]] %||% key
}
