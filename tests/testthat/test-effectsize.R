test_that("effectsize::cohens_d payload produces an effect-size table", {
  skip_if_not_installed("effectsize")

  p <- mellio_payload(effectsize::cohens_d(mpg ~ am, data = mtcars))

  expect_equal(p$card_kind, "table")
  expect_equal(p$type, "effect_size")
  expect_equal(p$fields$table_type, "effect_sizes")
  expect_equal(p$fields$rows[[1]]$effect, "Cohen's d")
  expect_true(is.numeric(p$fields$rows[[1]]$estimate))
  expect_true(is.numeric(p$fields$rows[[1]]$ci_lower))
})

test_that("effectsize ANOVA measures preserve one row per term", {
  skip_if_not_installed("effectsize")

  fit <- aov(mpg ~ factor(cyl), data = mtcars)
  eta <- suppressMessages(effectsize::eta_squared(fit))
  omega <- suppressMessages(effectsize::omega_squared(fit))

  p_eta <- mellio_payload(eta)
  expect_equal(p_eta$fields$rows[[1]]$variable, "factor(cyl)")
  expect_equal(p_eta$fields$rows[[1]]$effect, "\u03b7\u00b2")
  expect_true(is.numeric(p_eta$fields$rows[[1]]$estimate))

  p_omega <- mellio_payload(omega)
  expect_equal(p_omega$fields$rows[[1]]$variable, "factor(cyl)")
  expect_equal(p_omega$fields$rows[[1]]$effect, "\u03c9\u00b2")
  expect_true(is.numeric(p_omega$fields$rows[[1]]$estimate))
})

test_that("effect-size-shaped data frames still produce effect-size cards", {
  df <- data.frame(
    Cohens_d = 0.49,
    CI = 0.95,
    CI_low = -0.02,
    CI_high = 1.01
  )

  p <- mellio_payload(df)

  expect_equal(p$type, "effect_size")
  expect_equal(p$fields$table_type, "effect_sizes")
  expect_equal(p$fields$source, "custom_data_frame")
  expect_equal(p$fields$rows[[1]]$effect, "Cohen's d")
  expect_equal(p$fields$rows[[1]]$estimate, 0.49)
  expect_equal(p$fields$rows[[1]]$ci_lower, -0.02)
})
