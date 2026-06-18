test_that("mellio_payload.survfit returns single-arm KM payload with tracks and table", {
  skip_if_not_installed("survival")

  fit <- survival::survfit(
    survival::Surv(time, status) ~ 1,
    data = survival::lung
  )
  p <- mellio_payload(fit)

  expect_equal(p$type, "kaplan_meier_survival")
  expect_equal(p$card_kind, "inline")
  expect_true(is.numeric(p$fields$n))
  expect_true(is.numeric(p$fields$events))
  expect_false(p$fields$has_strata)
  expect_equal(p$fields$groups_count, 1L)
  expect_true(length(p$fields$rows) == 1L)
  expect_true(length(p$fields$columns) >= 4L)
  expect_equal(p$fields$table_type, "kaplan_meier_summary")
  expect_null(p$fields$statistic)

  km <- p$figure_data$km_curve
  expect_false(is.null(km))
  expect_true(length(km$tracks) == 1L)
  track <- km$tracks[[1L]]
  expect_true(length(track$time) > 1L)
  expect_equal(length(track$time), length(track$surv))
  expect_true(is.numeric(track$lower) && is.numeric(track$upper))
  expect_true(is.integer(track$n_risk) && is.integer(track$n_event) && is.integer(track$n_censor))
  expect_true(is.numeric(track$median))

  figures <- p$metadata$available_figures
  expect_true(is.list(figures) && length(figures) >= 1L)
  expect_true(any(vapply(figures, function(fig) identical(fig$type, "km_curve"), logical(1))))
})

test_that("mellio_payload.survfit produces per-stratum tracks and log-rank when grouped", {
  skip_if_not_installed("survival")

  fit <- survival::survfit(
    survival::Surv(time, status) ~ sex,
    data = survival::lung
  )
  p <- mellio_payload(fit, .env = environment())

  expect_equal(p$type, "kaplan_meier_survival")
  expect_true(p$fields$has_strata)
  expect_true(p$fields$groups_count >= 2L)
  expect_true(length(p$fields$rows) >= 2L)
  expect_true(!is.null(p$fields$predictor) && nzchar(p$fields$predictor))

  stat <- p$fields$statistic
  expect_false(is.null(stat))
  expect_equal(stat$name, "\u03c7\u00b2")
  expect_true(is.numeric(stat$value) && is.finite(stat$value))
  expect_true(is.numeric(stat$df) && stat$df >= 1)
  expect_true(is.numeric(p$fields$p_value))
  expect_equal(p$fields$test_name, "Log-rank test")

  km <- p$figure_data$km_curve
  expect_true(length(km$tracks) >= 2L)
  expect_true(km$has_strata)
  expect_false(is.null(km$log_rank))
  expect_equal(km$predictor, "sex")
  expect_equal(km$time_label, "time")
  expect_equal(km$time_unit, "days")
  expect_equal(km$event_label, "status")

  expect_true(all(vapply(km$tracks, function(tr) !is.null(tr$label) && nzchar(tr$label), logical(1))))
  expect_true(any(vapply(km$tracks, function(tr) {
    !is.null(tr$median) && !is.na(tr$median)
  }, logical(1))))
})

test_that("mellio_payload.survfit table rows carry n/events/median/CI per group", {
  skip_if_not_installed("survival")

  fit <- survival::survfit(
    survival::Surv(time, status) ~ sex,
    data = survival::lung
  )
  p <- mellio_payload(fit, .env = environment())

  rows <- p$fields$rows
  expect_true(length(rows) >= 2L)
  row <- rows[[1L]]
  expect_true(!is.null(row$group) && nzchar(as.character(row$group)))
  expect_true(is.numeric(row$n) || is.null(row$n))
  expect_true(is.numeric(row$events) || is.null(row$events))
  expect_true(!is.null(row$median))

  col_keys <- vapply(p$fields$columns, function(c) as.character(c$key), character(1))
  expect_true(all(c("group", "n", "events", "median", "ci") %in% col_keys))
})

test_that("mellio_payload.survfit preserves non-default confidence levels", {
  skip_if_not_installed("survival")

  fit <- survival::survfit(
    survival::Surv(time, status) ~ 1,
    data = survival::lung,
    conf.int = 0.9
  )
  p <- mellio_payload(fit)

  ci_col <- p$fields$columns[[which(vapply(p$fields$columns, function(col) {
    identical(col$key, "ci")
  }, logical(1)))[1L]]]
  expect_equal(p$fields$conf_level, 0.9)
  expect_equal(ci_col$label, "90% CI")
  expect_match(p$fields$note, "90% CI", fixed = TRUE)
  expect_equal(p$figure_data$km_curve$conf_level, 0.9)
  expect_true(!is.na(p$fields$rows[[1L]]$ci_lower))
  expect_true(!is.na(p$fields$rows[[1L]]$ci_upper))
})

test_that("mellio_payload.survfit labels counting-process Surv axes honestly", {
  skip_if_not_installed("survival")

  d <- data.frame(
    start = c(0, 0, 4, 4),
    stop = c(4, 6, 9, 12),
    status = c(0, 1, 1, 0),
    arm = factor(c("A", "A", "B", "B"))
  )
  fit <- survival::survfit(
    survival::Surv(start, stop, status) ~ arm,
    data = d
  )
  p <- mellio_payload(fit, .env = environment())

  expect_equal(p$fields$time_label, "stop")
  expect_equal(p$fields$event_label, "status")
  expect_equal(p$figure_data$km_curve$time_label, "stop")
  expect_equal(p$figure_data$km_curve$event_label, "status")
})
