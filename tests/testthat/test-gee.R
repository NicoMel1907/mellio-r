test_that("mellio_payload.geeglm reports Wald coefficients and GEE context", {
  skip_if_not_installed("geepack")

  geeglm <- getExportedValue("geepack", "geeglm")
  set.seed(2026)
  panel <- data.frame(
    participant = rep(seq_len(20), each = 4),
    time = rep(seq_len(4), 20),
    condition = factor(rep(c("A", "B"), each = 40))
  )
  panel$y <- 12 + 0.3 * panel$time + 1.7 * (panel$condition == "B") +
    stats::rnorm(nrow(panel), sd = 2)

  fit <- geeglm(
    y ~ time + condition,
    family = stats::gaussian,
    id = participant,
    data = panel,
    corstr = "exchangeable"
  )
  p <- mellio_payload(fit)

  expect_equal(p$type, "gee_model_summary")
  expect_equal(p$card_kind, "inline")
  expect_equal(p$fields$statistic_label, "Wald \u03c7\u00b2")
  expect_equal(p$fields$std_error_method, "robust_sandwich")
  expect_equal(p$fields$id_variable, "participant")
  expect_equal(p$fields$correlation_structure, "exchangeable")
  expect_equal(p$fields$working_correlation_parameter_name, "alpha")
  expect_true(is.numeric(p$fields$working_correlation_parameter))
  expect_true(is.numeric(p$fields$scale_parameter))

  coefs <- p$fields$coefficients
  expect_true(length(coefs) >= 2)
  expect_true(all(vapply(coefs, function(row) identical(row$estimate_name, "B"), logical(1))))
  expect_true(any(vapply(coefs, function(row) identical(row$statistic_label, "Wald \u03c7\u00b2"), logical(1))))
  expect_true(all(vapply(coefs, function(row) identical(row$df, 1), logical(1))))
})
