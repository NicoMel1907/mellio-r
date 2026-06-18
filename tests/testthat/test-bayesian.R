test_that("brmsfit payload extracts fixed-effect posterior summaries", {
  fixed <- matrix(
    c(0.40, 0.10, 0.20, 0.60, 1.00, 1200, 1100,
      1.25, 0.30, 0.65, 1.85, 1.01, 900, 850),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("Intercept", "x"),
      c("Estimate", "Est.Error", "l-95% CI", "u-95% CI",
        "Rhat", "Bulk_ESS", "Tail_ESS")
    )
  )
  fit <- structure(list(fixed = fixed), class = "brmsfit")
  p <- mellio_payload(fit)

  expect_equal(p$card_kind, "table")
  expect_equal(p$type, "bayesian_model_summary")
  expect_equal(p$fields$table_type, "bayesian_parameters")
  expect_equal(p$fields$method, "brms")
  expect_equal(length(p$fields$rows), 2L)
  expect_equal(p$fields$rows[[2]]$term, "x")
  expect_true(is.numeric(p$fields$rows[[2]]$ci_lower))
  expect_true(is.numeric(p$fields$rows[[2]]$rhat))
})

test_that("stanreg payload extracts posterior summaries", {
  posterior_summary <- matrix(
    c(0.40, 0.10, 0.20, 0.60, 1.00, 1200,
      1.25, 0.30, 0.65, 1.85, 1.01, 900),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("(Intercept)", "x"),
      c("mean", "sd", "2.5%", "97.5%", "Rhat", "n_eff")
    )
  )
  fit <- structure(list(posterior_summary = posterior_summary), class = "stanreg")
  p <- mellio_payload(fit)

  expect_equal(p$type, "bayesian_model_summary")
  expect_equal(p$fields$table_type, "bayesian_parameters")
  expect_equal(p$fields$method, "rstanarm")
  expect_equal(length(p$fields$rows), 2L)
  expect_equal(p$fields$rows[[2]]$term, "x")
  expect_true(is.numeric(p$fields$rows[[2]]$ci_upper))
  expect_true(is.numeric(p$fields$rows[[2]]$ess_bulk))
})
