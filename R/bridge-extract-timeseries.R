# R bridge -- time-series model and forecast extractors.

#' @rdname mellio_payload
#' @export
mellio_payload.Arima <- function(x, ..., .call = NULL) {
  ms_arima_payload(x, .call = .call)
}

#' @rdname mellio_payload
#' @export
mellio_payload.forecast_ARIMA <- function(x, ..., .call = NULL) {
  ms_arima_payload(x, .call = .call, package = "forecast")
}

#' @rdname mellio_payload
#' @export
mellio_payload.ARIMA <- function(x, ..., .call = NULL) {
  ms_arima_payload(x, .call = .call, package = "forecast")
}

#' @rdname mellio_payload
#' @export
mellio_payload.ets <- function(x, ..., .call = NULL) {
  coefs <- x$par %||% numeric(0)
  coefs <- coefs[is.finite(coefs)]
  rows <- lapply(seq_along(coefs), function(i) {
    list(
      term = ms_vector_name(coefs, i, "parameter"),
      estimate = ms_safe_numeric(coefs[[i]])
    )
  })
  if (!length(rows)) {
    rows <- list(list(term = "Model", estimate = NA_real_))
  }

  fields <- list(
    table_type = "time_series_coefficients",
    model_family = "ETS",
    method = as.character(x$method %||% "Exponential smoothing"),
    columns = list(
      list(key = "term", label = "Term", format = "text"),
      list(key = "estimate", label = "Estimate", format = "number")
    ),
    rows = rows,
    n_parameters = length(rows),
    source = "forecast::ets"
  )
  fields <- ms_timeseries_fit_fields(fields, x)

  ms_build_envelope(
    type = "time_series_model",
    type_label = "Time-series model",
    call = ms_timeseries_call(x, .call),
    fields = fields,
    raw_output = ms_capture_output(x),
    packages = ms_packages_basic(extras = "forecast"),
    card_kind = "table"
  )
}

#' @rdname mellio_payload
#' @export
mellio_payload.forecast <- function(x, ..., .call = NULL) {
  mean_values <- as.numeric(x$mean %||% numeric(0))
  if (!length(mean_values)) {
    stop("forecast object does not contain point forecasts.", call. = FALSE)
  }

  interval <- ms_forecast_interval(x)
  periods <- ms_forecast_periods(x$mean, length(mean_values))
  rows <- lapply(seq_along(mean_values), function(i) {
    row <- list(
      horizon = as.integer(i),
      period = periods[[i]],
      forecast = ms_safe_numeric(mean_values[[i]])
    )
    if (!is.null(interval)) {
      row$ci_lower <- ms_safe_numeric(interval$lower[[i]])
      row$ci_upper <- ms_safe_numeric(interval$upper[[i]])
    }
    row
  })

  columns <- list(
    list(key = "horizon", label = "Horizon", format = "integer"),
    list(key = "period", label = "Period", format = "text"),
    list(key = "forecast", label = "Forecast", format = "number")
  )
  if (!is.null(interval)) {
    columns <- c(columns, list(
      list(key = "ci_lower", label = paste0(interval$level, "% PI lower"), format = "number"),
      list(key = "ci_upper", label = paste0(interval$level, "% PI upper"), format = "number")
    ))
  }

  fields <- list(
    table_type = "forecasts",
    method = as.character(x$method %||% "Forecast"),
    columns = columns,
    rows = rows,
    n_forecasts = length(rows),
    interval_level = if (!is.null(interval)) interval$level else NA_real_,
    source = "forecast"
  )

  ms_build_envelope(
    type = "time_series_forecast",
    type_label = "Forecast",
    call = ms_timeseries_call(x, .call),
    fields = fields,
    raw_output = ms_capture_output(x),
    packages = ms_packages_basic(extras = "forecast"),
    card_kind = "table"
  )
}

ms_arima_payload <- function(x, .call = NULL, package = NULL) {
  coefs <- tryCatch(stats::coef(x), error = function(e) x$coef %||% numeric(0))
  coefs <- coefs[is.finite(coefs)]
  se <- ms_arima_se(x, names(coefs))

  rows <- lapply(seq_along(coefs), function(i) {
    estimate <- ms_safe_numeric(coefs[[i]])
    std_error <- ms_safe_numeric(se[[i]] %||% NA_real_)
    row <- list(
      term = ms_vector_name(coefs, i, "parameter"),
      estimate = estimate,
      std_error = std_error
    )
    if (!is.na(estimate) && !is.na(std_error) && std_error > 0) {
      z <- estimate / std_error
      row$statistic <- ms_safe_numeric(z)
      row$p_value <- ms_safe_numeric(2 * stats::pnorm(abs(z), lower.tail = FALSE))
    }
    row
  })
  if (!length(rows)) {
    rows <- list(list(term = "Model", estimate = NA_real_))
  }

  fields <- list(
    table_type = "time_series_coefficients",
    model_family = "ARIMA",
    method = as.character(x$method %||% "ARIMA"),
    columns = ms_arima_columns(rows),
    rows = rows,
    n_parameters = length(rows),
    source = if (is.null(package)) "stats::arima" else "forecast::Arima"
  )
  fields <- c(fields, ms_arima_order_fields(x))
  fields <- ms_timeseries_fit_fields(fields, x)

  ms_build_envelope(
    type = "time_series_model",
    type_label = "Time-series model",
    call = ms_timeseries_call(x, .call),
    fields = fields,
    raw_output = ms_capture_output(x),
    packages = ms_packages_basic(extras = package),
    card_kind = "table"
  )
}

ms_arima_se <- function(x, names) {
  vc <- x$var.coef %||% NULL
  if (is.null(vc) || !length(vc)) return(rep(NA_real_, length(names)))
  out <- sqrt(diag(as.matrix(vc)))
  if (!length(names)) return(out)
  out[names]
}

ms_arima_columns <- function(rows) {
  has <- function(key) {
    any(vapply(rows, function(row) !is.null(row[[key]]) && !is.na(row[[key]]), logical(1)))
  }
  cols <- list(
    list(key = "term", label = "Term", format = "text"),
    list(key = "estimate", label = "Estimate", format = "number")
  )
  if (has("std_error")) cols <- c(cols, list(list(key = "std_error", label = "SE", format = "number")))
  if (has("statistic")) cols <- c(cols, list(list(key = "statistic", label = "z", format = "statistic")))
  if (has("p_value")) cols <- c(cols, list(list(key = "p_value", label = "p", format = "pvalue")))
  cols
}

ms_arima_order_fields <- function(x) {
  arma <- x$arma %||% NULL
  if (is.null(arma) || length(arma) < 7L) return(list())
  p <- as.integer(arma[[1]])
  q <- as.integer(arma[[2]])
  p_seasonal <- as.integer(arma[[3]])
  q_seasonal <- as.integer(arma[[4]])
  period <- as.integer(arma[[5]])
  d <- as.integer(arma[[6]])
  d_seasonal <- as.integer(arma[[7]])

  fields <- list(
    model_order = paste0("(", p, ",", d, ",", q, ")"),
    p = p,
    d = d,
    q = q
  )
  if (p_seasonal != 0L || d_seasonal != 0L || q_seasonal != 0L) {
    fields$seasonal_order <- paste0(
      "(", p_seasonal, ",", d_seasonal, ",", q_seasonal, ")[", period, "]"
    )
    fields$seasonal_period <- period
  }
  fields
}

ms_timeseries_fit_fields <- function(fields, x) {
  add_num <- function(name, value) {
    value <- ms_safe_numeric(value)
    if (!is.na(value)) fields[[name]] <<- value
  }
  add_num("aic", x$aic %||% NA_real_)
  add_num("bic", x$bic %||% NA_real_)
  add_num("aicc", x$aicc %||% NA_real_)
  add_num("sigma2", x$sigma2 %||% NA_real_)
  add_num("log_likelihood", x$loglik %||% x$logLik %||% NA_real_)
  n <- x$nobs %||% x$n %||% NA_integer_
  if (!is.null(n) && length(n) > 0L && !is.na(n)) fields$n <- as.integer(n)
  fields
}

ms_forecast_interval <- function(x) {
  lower <- x$lower %||% NULL
  upper <- x$upper %||% NULL
  if (is.null(lower) || is.null(upper)) return(NULL)
  lower <- as.matrix(lower)
  upper <- as.matrix(upper)
  if (!nrow(lower) || !nrow(upper)) return(NULL)
  levels <- x$level %||% colnames(lower) %||% NA_real_
  level_nums <- suppressWarnings(as.numeric(levels))
  idx <- if (all(is.na(level_nums))) ncol(lower) else which.max(level_nums)
  list(
    level = level_nums[[idx]] %||% levels[[idx]],
    lower = as.numeric(lower[, idx]),
    upper = as.numeric(upper[, idx])
  )
}

ms_forecast_periods <- function(x, n) {
  times <- tryCatch(stats::time(x), error = function(e) NULL)
  if (!is.null(times) && length(times) >= n) {
    return(as.character(round(as.numeric(times)[seq_len(n)], 6)))
  }
  nm <- names(x)
  if (!is.null(nm) && length(nm) >= n) return(as.character(nm[seq_len(n)]))
  as.character(seq_len(n))
}

ms_timeseries_call <- function(x, .call = NULL) {
  if (!is.null(.call)) return(trimws(gsub("\\s+", " ", .call)))
  call_obj <- x$call %||% NULL
  if (!is.null(call_obj)) return(ms_deparse_call(call_obj))
  NA_character_
}

ms_vector_name <- function(x, i, prefix) {
  nm <- names(x)
  if (!is.null(nm) && length(nm) >= i && nzchar(nm[[i]])) return(nm[[i]])
  paste0(prefix, "_", i)
}
