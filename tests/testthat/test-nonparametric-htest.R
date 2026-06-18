test_that("Kruskal-Wallis htest maps to a named nonparametric payload", {
  p <- mellio_payload(kruskal.test(extra ~ group, data = datasets::sleep))

  expect_equal(p$card_kind, "inline")
  expect_equal(p$type, "kruskal_wallis_test")
  expect_match(p$type_label, "Kruskal-Wallis")
  expect_equal(p$fields$statistic$name, "H")
  expect_true(is.numeric(p$fields$statistic$df))
  expect_equal(p$fields$p_value_method, "chi_square_approximation")
  expect_equal(p$fields$outcome, "extra")
  expect_equal(p$fields$predictor, "group")
})

test_that("Friedman htest maps to a named nonparametric payload", {
  y <- c(10, 12, 11, 14, 13, 16, 9, 11, 12, 13, 14, 15)
  groups <- factor(rep(c("A", "B", "C"), times = 4))
  blocks <- factor(rep(seq_len(4), each = 3))

  p <- mellio_payload(friedman.test(y, groups, blocks))

  expect_equal(p$card_kind, "inline")
  expect_equal(p$type, "friedman_test")
  expect_match(p$type_label, "Friedman")
  expect_equal(p$fields$statistic$name, "Friedman chi-squared")
  expect_true(is.numeric(p$fields$statistic$df))
  expect_equal(p$fields$p_value_method, "chi_square_approximation")
})

test_that("McNemar htest maps to a named nonparametric payload", {
  tab <- matrix(c(32, 8, 19, 41), nrow = 2)
  p <- mellio_payload(mcnemar.test(tab))

  expect_equal(p$card_kind, "inline")
  expect_equal(p$type, "mcnemar_test")
  expect_match(p$type_label, "McNemar")
  expect_equal(p$fields$statistic$name, "McNemar's chi-squared")
  expect_true(is.numeric(p$fields$statistic$df))
  expect_equal(p$fields$p_value_method, "chi_square_approximation")
})

test_that("Wilcoxon htests record rank-test p-value provenance", {
  rank_sum <- suppressWarnings(mellio_payload(wilcox.test(extra ~ group, data = datasets::sleep)))
  signed_rank <- suppressWarnings(mellio_payload(wilcox.test(
    c(1, 2, 3, 4, 5),
    c(2, 3, 5, 5, 7),
    paired = TRUE
  )))

  expect_equal(rank_sum$type, "wilcoxon_rank_sum")
  expect_equal(rank_sum$fields$p_value_method, "wilcoxon_rank_test")
  expect_equal(signed_rank$type, "wilcoxon_signed_rank")
  expect_equal(signed_rank$fields$p_value_method, "wilcoxon_rank_test")
})

test_that("Kruskal-Wallis htest enriches medians, effect size, and figure data from .data", {
  df <- data.frame(
    score = c(4, 5, 6, 7, 8, 8, 9, 10, 12, 13, 15, 16),
    group = factor(rep(c("A", "B", "C"), each = 4))
  )
  fit <- kruskal.test(score ~ group, data = df)

  p <- mellio_payload(fit, .data = df)

  expect_equal(p$type, "kruskal_wallis_test")
  expect_length(p$fields$groups, 3L)
  expect_equal(p$fields$groups[[1]]$median, median(df$score[df$group == "A"]))
  expect_equal(p$fields$groups[[2]]$iqr, IQR(df$score[df$group == "B"]))
  expect_equal(p$fields$sample_size, nrow(df))
  expect_equal(p$fields$effect_size$name, "eta_sq_h")
  expect_equal(p$fields$effect_size$method, "kruskal_wallis_eta_squared_h")
  expect_equal(p$fields$effect_size$formula, "(H - k + 1) / (N - k)")
  expect_true(p$fields$effect_size$value >= 0)

  fig <- p$figure_data$nonparametric_group_plot
  expect_equal(fig$source, "kruskal_wallis_test")
  expect_equal(fig$outcome, "score")
  expect_equal(fig$factor$variable, "group")
  expect_length(fig$groups, 3L)
  expect_length(fig$observations, nrow(df))
  expect_true(any(vapply(p$metadata$available_figures, function(f) {
    identical(f$type, "nonparametric_group_plot") && isTRUE(f$default)
  }, logical(1))))
})

test_that("Rank tests flag numeric variables treated as grouping factors", {
  fit <- kruskal.test(len ~ dose, data = ToothGrowth)
  p <- mellio_payload(fit, .data = ToothGrowth)

  expect_equal(p$type, "kruskal_wallis_test")
  expect_match(p$fields$grouping_note, "dose is numeric")
  expect_match(p$fields$grouping_note, "distinct values as group levels")
})

test_that("Wilcoxon rank-sum htest enriches group medians and Cliff's delta", {
  sleep_df <- datasets::sleep
  fit <- suppressWarnings(wilcox.test(extra ~ group, data = sleep_df, exact = FALSE))

  p <- mellio_payload(fit, .data = sleep_df)

  expect_equal(p$type, "wilcoxon_rank_sum")
  expect_length(p$fields$groups, 2L)
  expect_equal(p$fields$groups[[1]]$median, median(sleep_df$extra[sleep_df$group == 1]))
  expect_equal(p$fields$groups[[2]]$q3, unname(quantile(sleep_df$extra[sleep_df$group == 2], 0.75)))
  expect_equal(p$fields$sample_size, nrow(sleep_df))
  expect_equal(p$fields$effect_size$name, "cliffs_delta")

  a <- sleep_df$extra[sleep_df$group == 1]
  b <- sleep_df$extra[sleep_df$group == 2]
  expect_equal(p$fields$effect_size$value, mean(sign(outer(a, b, "-"))), tolerance = 1e-9)
  expect_equal(p$figure_data$nonparametric_group_plot$source, "wilcoxon_rank_sum")
})

test_that("Nonparametric htests can recover simple source data from caller env", {
  p <- local({
    sleep_df <- datasets::sleep
    fit <- kruskal.test(extra ~ group, data = sleep_df)
    mellio_payload(fit)
  })

  expect_equal(p$type, "kruskal_wallis_test")
  expect_equal(p$fields$sample_size, nrow(datasets::sleep))
  expect_length(p$figure_data$nonparametric_group_plot$observations, nrow(datasets::sleep))
})
