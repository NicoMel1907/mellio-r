test_that("mellio_payload.manova emits multivariate tests instead of univariate ANOVA", {
  fit <- manova(
    cbind(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width) ~ Species,
    data = iris
  )

  p <- mellio_payload(fit)

  expect_equal(p$card_kind, "table")
  expect_equal(p$type, "manova_multivariate_tests")
  expect_equal(p$fields$table_type, "manova_multivariate_tests")

  tests <- vapply(p$fields$rows, function(row) row$test, character(1))
  expect_true(all(c("Pillai's Trace", "Wilks' Lambda") %in% tests))
  expect_true(all(vapply(p$fields$rows, function(row) {
    is.numeric(row$f) && is.numeric(row$hypothesis_df) &&
      is.numeric(row$error_df) && is.numeric(row$p_value)
  }, logical(1))))
})

test_that("mellio_payload.summary.manova supports a selected multivariate criterion", {
  fit <- manova(
    cbind(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width) ~ Species,
    data = iris
  )
  pillai <- summary(fit, test = "Pillai")

  p <- mellio_payload(pillai)

  expect_equal(p$card_kind, "table")
  expect_equal(p$type, "manova_multivariate_tests")
  expect_equal(length(p$fields$rows), 1L)
  expect_equal(p$fields$rows[[1]]$test, "Pillai's Trace")
  expect_true(is.numeric(p$fields$rows[[1]]$eta_sq_partial))
})
