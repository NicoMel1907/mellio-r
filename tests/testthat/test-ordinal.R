test_that("MASS::polr payload produces ordinal regression coefficients", {
  skip_if_not_installed("MASS")

  d <- mtcars
  d$gear_ord <- ordered(d$gear)
  fit <- suppressWarnings(MASS::polr(gear_ord ~ mpg, data = d, Hess = TRUE))

  p <- mellio_payload(fit)

  expect_equal(p$card_kind, "inline")
  expect_equal(p$type, "ordinal_regression")
  expect_equal(p$type_label, "Ordinal regression")
  expect_equal(p$fields$link_function, "logit")
  expect_equal(p$fields$coefficient_scale, "proportional_odds")
  expect_equal(p$fields$coefficients[[1]]$term, "mpg")
  expect_equal(p$fields$coefficients[[1]]$estimate_name, "OR")
  expect_true(p$fields$coefficients[[1]]$estimate > 0)
  expect_true(is.numeric(p$fields$coefficients[[1]]$p_value))

  plot <- p$figure_data$interaction_plot
  expect_equal(plot$source, "ordinal_predictions")
  expect_equal(plot$mean_kind, "predicted_probability")
  expect_equal(plot$interaction_kind, "continuous_by_categorical")
  expect_equal(plot$x$variable, "mpg")
  expect_equal(plot$moderator$variable, "outcome_level")
  expect_equal(length(plot$moderator$levels), length(levels(d$gear_ord)))
  expect_equal(length(plot$grid), 80 * length(levels(d$gear_ord)))
  expect_equal(p$metadata$available_figures[[1]]$type, "interaction_plot")
  expect_equal(p$metadata$available_figures[[1]]$label, "Predicted probability curves")
})

test_that("ordinal::clm payload produces ordinal regression coefficients", {
  skip_if_not_installed("ordinal")

  d <- mtcars
  d$gear_ord <- ordered(d$gear)
  fit <- ordinal::clm(gear_ord ~ mpg, data = d)

  p <- mellio_payload(fit)

  expect_equal(p$type, "ordinal_regression")
  expect_equal(p$fields$coefficient_scale, "proportional_odds")
  expect_equal(p$fields$coefficients[[1]]$term, "mpg")
  expect_equal(p$fields$coefficients[[1]]$estimate_name, "OR")
  expect_true(p$fields$coefficients[[1]]$estimate > 0)
  expect_equal(p$fields$coefficients[[1]]$ci_method, "profile_likelihood")
  expect_equal(p$fields$coefficient_ci_method, "profile_likelihood")
  expect_equal(round(p$fields$coefficients[[1]]$ci_lower, 2), 1.07)
  expect_equal(round(p$fields$coefficients[[1]]$ci_upper, 2), 1.40)
  expect_equal(p$figure_data$interaction_plot$source, "ordinal_predictions")
  expect_equal(p$figure_data$interaction_plot$interaction_kind, "continuous_by_categorical")
})

test_that("ordinal payload reports omnibus tests for multi-level factors", {
  skip_if_not_installed("ordinal")

  set.seed(2026)
  n <- 360
  d <- data.frame(
    group = factor(rep(c("A", "B", "C"), each = n / 3)),
    x = rnorm(n)
  )
  eta <- 0.45 * d$x + ifelse(d$group == "B", 0.7, ifelse(d$group == "C", -0.35, 0))
  p_low <- plogis(-0.8 - eta)
  p_mid <- plogis(0.9 - eta) - p_low
  probs <- cbind(p_low, p_mid, 1 - p_low - p_mid)
  d$y <- ordered(apply(probs, 1, function(pr) {
    sample(c("low", "mid", "high"), 1, prob = pr)
  }), levels = c("low", "mid", "high"))

  fit <- ordinal::clm(y ~ group + x, data = d)
  p <- mellio_payload(fit)
  tests <- p$fields$model_term_tests

  expect_equal(length(tests), 1)
  expect_equal(tests[[1]]$term, "group")
  expect_equal(tests[[1]]$method, "ordinal_factor_lrt")
  expect_equal(tests[[1]]$statistic$name, "chi2")
  expect_equal(as.numeric(tests[[1]]$statistic$df), 2)
  expect_true(is.numeric(tests[[1]]$p_value))
  expect_equal(p$fields$residual_df, stats::df.residual(fit))

  plot <- p$figure_data$interaction_plot
  expect_equal(plot$interaction_kind, "categorical_by_categorical")
  expect_equal(plot$x$variable, "group")
  expect_false(isTRUE(plot$connect_levels))
})

test_that("MASS::polr payload produces categorical ordinal probability figures", {
  skip_if_not_installed("MASS")

  d <- mtcars
  d$gear_ord <- ordered(d$gear)
  d$am_fac <- factor(d$am, labels = c("Automatic", "Manual"))
  fit <- suppressWarnings(MASS::polr(gear_ord ~ am_fac + mpg, data = d, Hess = TRUE))

  p <- mellio_payload(fit)
  plot <- p$figure_data$interaction_plot

  expect_equal(plot$source, "ordinal_predictions")
  expect_equal(plot$interaction_kind, "categorical_by_categorical")
  expect_equal(plot$x$variable, "am_fac")
  expect_equal(length(plot$x$levels), 2)
  expect_equal(length(plot$grid), 2 * length(levels(d$gear_ord)))
  expect_false(isTRUE(plot$connect_levels))
  expect_true(all(vapply(plot$held_constant, function(row) {
    isTRUE(row$variable == "mpg") && isTRUE(row$type == "mean")
  }, logical(1))))
  expect_equal(p$metadata$available_figures[[1]]$label, "Predicted probabilities")

  warnings <- vapply(p$fields$model_warnings, function(row) {
    if (is.null(row$message)) "" else row$message
  }, character(1))
  expect_true(any(grepl("Manual x 3", warnings, fixed = TRUE)))
  expect_true(any(grepl("Automatic x 5", warnings, fixed = TRUE)))
  am_row <- p$fields$coefficients[[which(vapply(p$fields$coefficients, function(row) {
    identical(row$term, "am_facManual")
  }, logical(1)))]]
  expect_true(isTRUE(am_row$estimate_unstable))
  expect_equal(am_row$estimate_note, "unstable/unbounded")
  expect_equal(am_row$statistic, "not interpreted")
  expect_equal(am_row$p_value, "not interpreted")
  expect_equal(am_row$ci_note, "not estimable: separation")

  mpg_row <- p$fields$coefficients[[which(vapply(p$fields$coefficients, function(row) {
    identical(row$term, "mpg")
  }, logical(1)))]]
  expect_false(isTRUE(mpg_row$estimate_unstable))
  expect_equal(mpg_row$ci_method, "wald")
  expect_true(isTRUE(mpg_row$ci_fallback))
  expect_true(is.numeric(mpg_row$ci_lower))
  expect_true(is.numeric(mpg_row$ci_upper))
})

test_that("ordinal::clm separated fits explain singular inference", {
  skip_if_not_installed("ordinal")

  d <- mtcars
  d$gear_ord <- ordered(d$gear)
  d$am_fac <- factor(d$am, labels = c("Automatic", "Manual"))
  fit <- suppressWarnings(ordinal::clm(gear_ord ~ am_fac + mpg, data = d))

  p <- mellio_payload(fit)
  messages <- vapply(p$fields$model_warnings, function(row) {
    if (is.null(row$message)) "" else row$message
  }, character(1))

  expect_true(any(grepl("information matrix is singular", messages, fixed = TRUE)))
  expect_true(any(grepl("Manual x 3", messages, fixed = TRUE)))
  expect_true(any(grepl("Automatic x 5", messages, fixed = TRUE)))

  mpg_row <- p$fields$coefficients[[which(vapply(p$fields$coefficients, function(row) {
    identical(row$term, "mpg")
  }, logical(1)))]]
  expect_false(isTRUE(mpg_row$estimate_unstable))
  expect_true(isTRUE(mpg_row$inference_unavailable))
  expect_equal(mpg_row$ci_note, "not available: singular information")
  expect_true(is.na(mpg_row$std_error))
})

test_that("MASS::polr interaction models warn without probability curves", {
  skip_if_not_installed("MASS")

  d <- mtcars
  d$gear_ord <- ordered(d$gear)
  d$am_fac <- factor(d$am, labels = c("Automatic", "Manual"))
  fit <- suppressWarnings(MASS::polr(gear_ord ~ mpg * am_fac, data = d, Hess = TRUE))

  p <- mellio_payload(fit)

  expect_null(p$figure_data$interaction_plot)
  expect_null(p$figure_data$coefficient_plot)
  messages <- vapply(p$fields$model_warnings, function(row) {
    if (is.null(row$message)) "" else row$message
  }, character(1))
  expect_true(any(grepl("main-effect proportional-odds models", messages, fixed = TRUE)))
  expect_true(any(grepl("Manual x 3", messages, fixed = TRUE)))
  expect_true(any(grepl("Automatic x 5", messages, fixed = TRUE)))
})
