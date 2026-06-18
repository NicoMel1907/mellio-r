# R bridge — descriptive statistics extractor.
#
# Supports the long-tail use case: "save M / SD / N for a variable so
# I can put it in a Method or Results paragraph." Two entry points:
#
#   mellio_payload(numeric_vector, name = "Age")
#       → computes M, SD, N, median, range, n_missing
#
#   mellio_payload(summary(numeric_vector))
#       → passes through what summary() gives (no SD)
#
# Both emit type = "descriptive_summary" inline cards. The JS renderer
# is descriptive-type aware so the headline becomes "M = X, SD = Y"
# instead of the test-statistic shape.
#
# Schema: docs/STATS-R-BRIDGE-SCHEMA.md

#' @rdname mellio_payload
#' @param name Human-readable variable name for descriptive payloads
#'   (e.g. `"Reaction time (ms)"`). Defaults to the deparsed argument
#'   when possible.
#' @export
mellio_payload.numeric <- function(x, name = NULL, ..., .call = NULL) {
  user_call <- match.call()$x
  call_str <- if (!is.null(.call)) {
    .call
  } else if (!is.null(user_call)) {
    paste(deparse(user_call, width.cutoff = 500L), collapse = " ")
  } else NA_character_

  var_name <- name %||% if (!is.null(user_call)) {
    paste(deparse(user_call, width.cutoff = 500L), collapse = " ")
  } else "variable"

  n_complete <- sum(!is.na(x))
  if (n_complete < 1L) {
    stop("No non-missing values to summarise.", call. = FALSE)
  }
  n_missing <- sum(is.na(x))

  m  <- ms_safe_numeric(mean(x,   na.rm = TRUE))
  sd_ <- ms_safe_numeric(stats::sd(x, na.rm = TRUE))
  md <- ms_safe_numeric(stats::median(x, na.rm = TRUE))
  mn <- ms_safe_numeric(min(x, na.rm = TRUE))
  mx <- ms_safe_numeric(max(x, na.rm = TRUE))

  fields <- list(
    statistic = list(name = "M", value = m),
    p_value   = NA_real_,
    estimate  = list(name = "SD", value = sd_),
    n         = n_complete,
    median    = md,
    range     = I(c(mn, mx))
  )
  if (n_missing > 0L) fields$n_missing <- n_missing

  ms_build_envelope(
    type       = "descriptive_summary",
    type_label = paste0("Descriptive statistics \u2014 ", trimws(var_name)),
    call       = trimws(gsub("\\s+", " ", call_str)),
    fields     = fields,
    raw_output = ms_capture_output(summary(x))
  )
}

#' @rdname mellio_payload
#' @export
mellio_payload.summaryDefault <- function(x, ..., .call = NULL) {
  # x is the named numeric returned by summary(numeric_vector).
  # Names: "Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", optional "NA's".
  user_call <- match.call()$x
  call_str <- if (!is.null(.call)) {
    .call
  } else if (!is.null(user_call)) {
    paste(deparse(user_call, width.cutoff = 500L), collapse = " ")
  } else NA_character_

  pick <- function(key) {
    if (key %in% names(x)) ms_safe_numeric(unname(x[[key]])) else NA_real_
  }

  fields <- list(
    statistic = list(name = "M", value = pick("Mean")),
    p_value   = NA_real_,
    # SD isn't available from summary() — surface this in the type
    # label so the user notices. The headline will show M, median,
    # range; users wanting SD should pass the raw vector to
    # mellio_payload() directly.
    median    = pick("Median"),
    range     = I(c(pick("Min."), pick("Max."))),
    quartiles = list(
      q1 = pick("1st Qu."),
      q3 = pick("3rd Qu.")
    )
  )
  if ("NA's" %in% names(x)) {
    fields$n_missing <- ms_safe_numeric(unname(x[["NA's"]]))
  }

  ms_build_envelope(
    type       = "descriptive_summary",
    type_label = "Descriptive statistics (from summary())",
    call       = trimws(gsub("\\s+", " ", call_str)),
    fields     = fields,
    raw_output = ms_capture_output(x)
  )
}
