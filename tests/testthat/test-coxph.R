test_that("mellio_payload.coxph reports hazard ratios and survival fit", {
  skip_if_not_installed("survival")

  fit <- survival::coxph(
    survival::Surv(time, status) ~ age + sex,
    data = survival::lung
  )
  p <- mellio_payload(fit)

  expect_equal(p$type, "cox_proportional_hazards")
  expect_equal(p$card_kind, "inline")
  expect_equal(p$fields$coefficient_scale, "hazard_ratio")
  expect_equal(p$fields$statistic$name, "\u03c7\u00b2")
  expect_true(is.numeric(p$fields$n))
  expect_true(is.numeric(p$fields$events))

  coefs <- p$fields$coefficients
  expect_true(length(coefs) >= 2)
  expect_true(all(vapply(coefs, function(row) identical(row$estimate_name, "HR"), logical(1))))
  expect_true(all(vapply(coefs, function(row) is.numeric(row$estimate), logical(1))))
  expect_true(any(vapply(coefs, function(row) !is.null(row$ci_lower), logical(1))))

  fig <- p$figure_data$coefficient_plot
  expect_equal(fig$coefficient_scale, "hazard_ratio")
  expect_equal(fig$estimate_label, "HR")
  expect_equal(fig$n, p$fields$n)
  expect_equal(fig$events, p$fields$events)
  age_fig <- fig$coefficients[[which(vapply(fig$coefficients, function(row) {
    identical(row$term, "age")
  }, logical(1)))[1]]]
  age_field <- coefs[[which(vapply(coefs, function(row) {
    identical(row$term, "age")
  }, logical(1)))[1]]]
  expect_equal(age_fig$estimate, age_field$estimate, tolerance = 1e-8)
  expect_equal(age_fig$log_estimate, age_field$log_estimate, tolerance = 1e-8)
})

test_that("mellio_payload.coxph warns and suppresses unstable hazard-ratio plots", {
  skip_if_not_installed("survival")

  d <- data.frame(
    time = c(1, 2, 3, 4, 5, 6, 7, 8),
    status = c(1, 1, 1, 1, 0, 0, 0, 0),
    exposed = c(1, 1, 1, 1, 0, 0, 0, 0)
  )
  fit <- suppressWarnings(survival::coxph(
    survival::Surv(time, status) ~ exposed,
    data = d
  ))
  p <- mellio_payload(fit)

  warning_types <- vapply(p$fields$model_warnings, function(row) row$type, character(1))
  expect_true("separation_or_boundary" %in% warning_types)
  expect_null(p$figure_data$coefficient_plot)
  figures <- p$metadata$available_figures %||% list()
  expect_false(any(vapply(figures, function(fig) {
    identical(fig$type, "coefficient_plot")
  }, logical(1))))
})

test_that("mellio_payload.coxph labels factor contrasts", {
  skip_if_not_installed("survival")

  lung2 <- transform(
    survival::lung,
    sex = factor(sex, levels = c(1, 2), labels = c("Male", "Female"))
  )
  fit <- survival::coxph(
    survival::Surv(time, status) ~ age + sex,
    data = lung2
  )
  p <- mellio_payload(fit)

  sex_row <- p$fields$coefficients[[which(vapply(p$fields$coefficients, function(row) {
    identical(row$term, "sexFemale")
  }, logical(1)))[1]]]
  expect_equal(sex_row$label, "sex: Female vs. Male")
  expect_equal(sex_row$term_source, "sex")
  expect_equal(sex_row$contrast_level, "Female")
  expect_equal(sex_row$contrast_reference, "Male")

  sex_fig <- p$figure_data$coefficient_plot$coefficients[[which(vapply(
    p$figure_data$coefficient_plot$coefficients,
    function(row) identical(row$term, "sexFemale"),
    logical(1)
  ))[1]]]
  expect_equal(sex_fig$label, "sex: Female vs. Male")
})
