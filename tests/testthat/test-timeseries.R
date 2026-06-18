test_that("base Arima payload produces time-series coefficient rows", {
  fit <- stats::arima(Nile, order = c(1, 0, 0))
  p <- mellio_payload(fit)

  expect_equal(p$card_kind, "table")
  expect_equal(p$type, "time_series_model")
  expect_equal(p$fields$table_type, "time_series_coefficients")
  expect_equal(p$fields$model_family, "ARIMA")
  expect_equal(p$fields$model_order, "(1,0,0)")
  expect_true(length(p$fields$rows) >= 1L)
  expect_true(is.numeric(p$fields$rows[[1]]$estimate))
})

test_that("forecast Arima and forecast objects produce structured rows", {
  skip_if_not_installed("forecast")

  fit <- forecast::Arima(Nile, order = c(1, 0, 0))
  p <- mellio_payload(fit)
  expect_equal(p$type, "time_series_model")
  expect_equal(p$fields$source, "forecast::Arima")

  fc <- forecast::forecast(fit, h = 3)
  pf <- mellio_payload(fc)
  expect_equal(pf$type, "time_series_forecast")
  expect_equal(pf$fields$table_type, "forecasts")
  expect_equal(length(pf$fields$rows), 3L)
  expect_true(is.numeric(pf$fields$rows[[1]]$forecast))
  expect_true(is.numeric(pf$fields$rows[[1]]$ci_lower))
})

test_that("forecast ets payload produces time-series parameter rows", {
  skip_if_not_installed("forecast")

  fit <- forecast::ets(Nile)
  p <- mellio_payload(fit)

  expect_equal(p$type, "time_series_model")
  expect_equal(p$fields$table_type, "time_series_coefficients")
  expect_equal(p$fields$model_family, "ETS")
  expect_true(length(p$fields$rows) >= 1L)
})
