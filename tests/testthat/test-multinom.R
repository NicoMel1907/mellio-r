test_that("mellio_payload.multinom preserves comparison-by-term coefficients", {
  skip_if_not_installed("nnet")

  multinom <- getExportedValue("nnet", "multinom")
  set.seed(2026)
  survey <- data.frame(
    choice = factor(sample(c("car", "bus", "bike"), 180, replace = TRUE)),
    age = stats::rnorm(180, mean = 42, sd = 11),
    income = stats::rnorm(180, mean = 52000, sd = 12000)
  )

  fit <- multinom(choice ~ age + income, data = survey, trace = FALSE)
  p <- mellio_payload(fit)

  expect_equal(p$type, "multinomial_logistic_regression")
  expect_equal(p$card_kind, "table")
  expect_equal(p$fields$table_type, "multinomial_coefficients")
  expect_equal(p$fields$statistic_label, "z")
  expect_equal(p$fields$model_family, "multinomial logistic")
  expect_true(is.numeric(p$fields$aic))
  expect_true(is.numeric(p$fields$n))

  rows <- p$fields$rows
  expect_true(length(rows) >= 6)
  expect_true(all(vapply(rows, function(row) {
    !is.null(row$comparison) && nzchar(row$comparison)
  }, logical(1))))
  expect_true(any(vapply(rows, function(row) identical(row$term, "age"), logical(1))))
  expect_true(any(vapply(rows, function(row) !is.null(row$p_value), logical(1))))
})
