test_that("TukeyHSD payload produces pairwise comparison rows", {
  fit <- aov(mpg ~ factor(cyl), data = mtcars)
  p <- mellio_payload(TukeyHSD(fit))

  expect_equal(p$card_kind, "table")
  expect_equal(p$type, "pairwise_comparisons")
  expect_equal(p$fields$table_type, "pairwise_comparisons")
  expect_equal(p$fields$method, "Tukey HSD")
  expect_equal(p$fields$adjustment_method, "tukey")
  expect_equal(p$fields$adjustment_label, "Tukey")
  expect_match(p$fields$note, "Tukey")
  expect_match(p$fields$note, "first contrast level minus the second")
  expect_equal(p$metadata$available_figures[[1]]$type, "pairwise_forest")
  expect_equal(p$metadata$available_figures[[1]]$label, "Pairwise forest plot")
  expect_true(isTRUE(p$metadata$available_figures[[1]]$default))
  expect_equal(p$figure_data$pairwise_forest$estimate_label, "Mean difference")
  expect_equal(p$figure_data$pairwise_forest$rows[[1]]$contrast, "4 - 6")
  expect_equal(p$figure_data$pairwise_forest$rows[[1]]$estimate, p$fields$rows[[1]]$estimate)
  expect_false("Family" %in% vapply(p$fields$columns, function(col) col$label, character(1)))
  expect_equal(length(p$fields$rows), 3L)
  expect_equal(p$fields$rows[[1]]$contrast, "4 - 6")
  expect_equal(p$fields$rows[[1]]$level_1, "4")
  expect_equal(p$fields$rows[[1]]$level_2, "6")
  expect_true(is.numeric(p$fields$rows[[1]]$estimate))
  expect_true(is.numeric(p$fields$rows[[1]]$p_value))

  p90 <- mellio_payload(TukeyHSD(fit, conf.level = 0.90))
  ci_cols <- vapply(p90$fields$columns, function(col) col$label, character(1))
  expect_true("90% CI lower" %in% ci_cols)
  expect_match(p90$fields$note, "90% confidence intervals")
  expect_equal(p90$figure_data$pairwise_forest$ci_label, "90% CI")

  fit2 <- aov(mpg ~ factor(cyl) + factor(am), data = mtcars)
  p2 <- mellio_payload(TukeyHSD(fit2))
  expect_true("Family" %in% vapply(p2$fields$columns, function(col) col$label, character(1)))
})

test_that("pairwise.htest payload preserves adjusted p-value matrix rows", {
  p <- mellio_payload(pairwise.t.test(mtcars$mpg, mtcars$cyl))

  expect_equal(p$type, "pairwise_comparisons")
  expect_equal(p$fields$adjustment_method, "holm")
  expect_equal(p$fields$adjustment_label, "Holm")
  expect_match(p$fields$note, "Holm")
  expect_null(p$figure_data)
  expect_null(p$metadata$available_figures)
  expect_equal(length(p$fields$rows), 3L)
  expect_equal(p$fields$rows[[1]]$contrast, "4 - 6")
  expect_match(p$fields$note, "first contrast level minus the second")
  expect_true("Difference" %in% vapply(p$fields$columns, function(col) col$label, character(1)))
  expect_true("t" %in% vapply(p$fields$columns, function(col) col$label, character(1)))
  expect_true(is.numeric(p$fields$rows[[1]]$p_value))
  expect_equal(p$fields$rows[[1]]$estimate, 6.920779, tolerance = 1e-6)
  expect_equal(p$fields$rows[[1]]$df, 29)
  expect_equal(p$fields$rows[[1]]$statistic, 4.441099, tolerance = 1e-6)
  group_6 <- mtcars$mpg[mtcars$cyl == 6]
  group_4 <- mtcars$mpg[mtcars$cyl == 4]
  by_cyl <- split(mtcars$mpg, mtcars$cyl)
  group_n <- vapply(by_cyl, length, integer(1))
  group_sd <- vapply(by_cyl, stats::sd, numeric(1))
  pooled_sd <- sqrt(sum((group_n - 1) * group_sd^2) / sum(group_n - 1))
  expect_equal(p$fields$rows[[1]]$n_1, length(group_4))
  expect_equal(p$fields$rows[[1]]$n_2, length(group_6))
  expect_equal(p$fields$rows[[1]]$mean_1, mean(group_4))
  expect_equal(p$fields$rows[[1]]$mean_2, mean(group_6))
  expect_true(is.numeric(p$fields$rows[[1]]$p_raw))
  expect_equal(p$fields$rows[[1]]$effect_size_name, "cohens_d")
  expect_equal(p$fields$rows[[1]]$effect_size_method, "all_groups_pooled_sd")
  expect_equal(
    p$fields$rows[[1]]$effect_size,
    (mean(group_4) - mean(group_6)) / pooled_sd,
    tolerance = 1e-9
  )
  expect_match(p$fields$note, "all-groups pooled SD")
})

test_that("mellio_open enriches inline pairwise.htest calls", {
  url <- withr::with_options(
    list(mellio.editor_url = "https://example.com"),
    mellio_open(pairwise.t.test(mtcars$mpg, mtcars$cyl), browse = FALSE)
  )
  b64 <- sub("&.*$", "", sub(".*payload=", "", url))
  std <- chartr("-_", "+/", b64)
  pad <- (4 - nchar(std) %% 4) %% 4
  if (pad > 0) std <- paste0(std, strrep("=", pad))
  parsed <- jsonlite::fromJSON(rawToChar(jsonlite::base64_dec(std)), simplifyVector = FALSE)

  expect_equal(parsed$call, "pairwise.t.test(mtcars$mpg, mtcars$cyl)")
  expect_equal(parsed$fields$rows[[1]]$estimate, 6.920779, tolerance = 1e-6)
  expect_equal(parsed$fields$rows[[1]]$df, 29)
})

test_that("pairwise.wilcox.test payload enriches nonparametric rows", {
  sleep_df <- datasets::sleep
  fit <- pairwise.wilcox.test(
    sleep_df$extra,
    sleep_df$group,
    p.adjust.method = "holm",
    exact = FALSE
  )
  p <- mellio_payload(fit)
  row <- p$fields$rows[[1]]
  group_2 <- sleep_df$extra[sleep_df$group == "2"]
  group_1 <- sleep_df$extra[sleep_df$group == "1"]
  wt <- suppressWarnings(stats::wilcox.test(group_2, group_1, exact = FALSE))
  delta <- mean(sign(outer(group_2, group_1, "-")))
  labels <- vapply(p$fields$columns, function(col) col$label, character(1))

  expect_equal(p$type, "pairwise_comparisons")
  expect_match(p$fields$method, "Wilcoxon")
  expect_equal(p$fields$adjustment_method, "holm")
  expect_equal(row$contrast, "1 - 2")
  expect_equal(row$n_1, 10L)
  expect_equal(row$n_2, 10L)
  expect_equal(row$median_1, stats::median(group_1))
  expect_equal(row$median_2, stats::median(group_2))
  expect_equal(row$median_difference, stats::median(group_1) - stats::median(group_2))
  expect_equal(row$statistic_label, "W")
  wt_canonical <- suppressWarnings(stats::wilcox.test(group_1, group_2, exact = FALSE))
  expect_equal(row$statistic, unname(wt_canonical$statistic[[1]]))
  expect_equal(row$effect_size_name, "cliffs_delta")
  expect_equal(row$effect_size, -delta)
  expect_true("Mdn 1" %in% labels)
  expect_true("Mdn 2" %in% labels)
  expect_true("Median difference" %in% labels)
  expect_true("Cliff's delta" %in% labels)
  expect_match(p$fields$note, "Median differences")
  expect_match(p$fields$note, "Cliff's delta")
})

test_that("mellio_open enriches inline pairwise.wilcox.test calls", {
  sleep_df <- datasets::sleep
  url <- withr::with_options(
    list(mellio.editor_url = "https://example.com"),
    mellio_open(
      pairwise.wilcox.test(
        sleep_df$extra,
        sleep_df$group,
        p.adjust.method = "holm",
        exact = FALSE
      ),
      browse = FALSE
    )
  )
  b64 <- sub("&.*$", "", sub(".*payload=", "", url))
  std <- chartr("-_", "+/", b64)
  pad <- (4 - nchar(std) %% 4) %% 4
  if (pad > 0) std <- paste0(std, strrep("=", pad))
  parsed <- jsonlite::fromJSON(rawToChar(jsonlite::base64_dec(std)), simplifyVector = FALSE)

  expect_match(parsed$call, "pairwise\\.wilcox\\.test")
  expect_equal(parsed$fields$rows[[1]]$statistic_label, "W")
  expect_equal(parsed$fields$rows[[1]]$effect_size_name, "cliffs_delta")
  expect_equal(parsed$fields$rows[[1]]$median_difference, -1.4)
})

test_that("paired pairwise.wilcox.test rows use signed-rank effect sizes", {
  paired_df <- data.frame(
    score = c(4, 5, 5, 6, 7, 7, 8, 9, 5, 6, 6, 7, 8, 8, 9, 10),
    condition = factor(rep(c("before", "after"), each = 8), levels = c("before", "after"))
  )
  fit <- pairwise.wilcox.test(
    paired_df$score,
    paired_df$condition,
    paired = TRUE,
    exact = FALSE
  )
  p <- mellio_payload(fit)
  row <- p$fields$rows[[1]]
  labels <- vapply(p$fields$columns, function(col) col$label, character(1))

  expect_equal(row$contrast, "before - after")
  expect_true(isTRUE(row$paired))
  expect_equal(row$n_pairs, 8L)
  expect_null(row$n_1)
  expect_null(row$n_2)
  expect_equal(row$statistic_label, "V")
  expect_equal(row$median_difference, -1)
  expect_equal(row$effect_size_name, "rank_biserial")
  expect_equal(row$effect_size, -1)
  expect_true("n pairs" %in% labels)
  expect_false("n 1" %in% labels)
  expect_false("n 2" %in% labels)
  expect_match(p$fields$note, "median of paired differences")
  expect_match(p$fields$note, "rank-biserial")
})

test_that("rstatix pairwise_wilcox_test payload enriches rank-sum rows", {
  skip_if_not_installed("rstatix")

  fit <- rstatix::pairwise_wilcox_test(
    ToothGrowth,
    len ~ dose,
    p.adjust.method = "holm",
    exact = FALSE
  )
  p <- mellio_payload(fit)
  row <- p$fields$rows[[1]]
  group_05 <- ToothGrowth$len[ToothGrowth$dose == 0.5]
  group_1 <- ToothGrowth$len[ToothGrowth$dose == 1]
  canonical_wilcox <- stats::wilcox.test(group_05, group_1, exact = FALSE)
  delta <- mean(sign(outer(group_05, group_1, "-")))
  labels <- vapply(p$fields$columns, function(col) col$label, character(1))

  expect_equal(p$type, "pairwise_comparisons")
  expect_equal(p$fields$method, "Wilcoxon rank-sum test")
  expect_equal(p$fields$source, "rstatix::wilcox_test")
  expect_equal(row$contrast, "0.5 - 1")
  expect_false(isTRUE(row$paired))
  expect_equal(row$n_1, 20L)
  expect_equal(row$n_2, 20L)
  expect_equal(row$median_1, stats::median(group_05))
  expect_equal(row$median_2, stats::median(group_1))
  expect_equal(row$median_difference, stats::median(group_05) - stats::median(group_1))
  expect_equal(row$statistic_label, "W")
  expect_equal(row$statistic, fit$statistic[[1]])
  expect_equal(as.numeric(row$p_raw), as.numeric(canonical_wilcox$p.value), tolerance = 1e-6)
  expect_equal(row$p_value, fit$p.adj[[1]])
  expect_equal(row$effect_size_name, "cliffs_delta")
  expect_equal(row$effect_size, delta)
  expect_true("n 1" %in% labels)
  expect_true("n 2" %in% labels)
  expect_true("Cliff's delta" %in% labels)
  expect_match(p$fields$note, "Cliff's delta")
})

test_that("rstatix pairwise_t_test payload enriches pooled t rows", {
  skip_if_not_installed("rstatix")

  fit <- rstatix::pairwise_t_test(
    ToothGrowth,
    len ~ dose,
    p.adjust.method = "holm"
  )
  p <- mellio_payload(fit)
  row <- p$fields$rows[[1]]
  group_05 <- ToothGrowth$len[ToothGrowth$dose == 0.5]
  group_1 <- ToothGrowth$len[ToothGrowth$dose == 1]
  by_dose <- split(ToothGrowth$len, ToothGrowth$dose)
  group_n <- vapply(by_dose, length, integer(1))
  group_sd <- vapply(by_dose, stats::sd, numeric(1))
  pooled_sd <- sqrt(sum((group_n - 1) * group_sd^2) / sum(group_n - 1))
  labels <- vapply(p$fields$columns, function(col) col$label, character(1))

  expect_equal(p$type, "pairwise_comparisons")
  expect_equal(p$fields$method, "t tests with pooled SD")
  expect_equal(p$fields$source, "rstatix::t_test")
  expect_match(p$fields$grouping_note, "dose is numeric")
  expect_equal(row$contrast, "0.5 - 1")
  expect_false(isTRUE(row$paired))
  expect_equal(row$n_1, 20L)
  expect_equal(row$n_2, 20L)
  expect_equal(row$mean_1, mean(group_05))
  expect_equal(row$mean_2, mean(group_1))
  expect_equal(row$estimate, mean(group_05) - mean(group_1))
  expect_equal(row$statistic_label, "t")
  expect_true(is.numeric(row$statistic))
  expect_true(is.numeric(row$df))
  expect_equal(row$p_raw, fit$p[[1]])
  expect_equal(row$p_value, fit$p.adj[[1]])
  expect_equal(row$effect_size_name, "cohens_d")
  expect_equal(row$effect_size_method, "all_groups_pooled_sd")
  expect_equal(row$effect_size, (mean(group_05) - mean(group_1)) / pooled_sd)
  expect_true("M 1" %in% labels)
  expect_true("M 2" %in% labels)
  expect_true("Cohen's d" %in% labels)
  expect_match(p$fields$note, "Cohen's d")
})

test_that("rstatix t_test payload keeps t statistics and confidence intervals", {
  skip_if_not_installed("rstatix")

  fit <- rstatix::t_test(ToothGrowth, len ~ supp)
  p <- mellio_payload(fit)
  row <- p$fields$rows[[1]]

  expect_equal(p$type, "pairwise_comparisons")
  expect_equal(p$fields$method, "Welch t test")
  expect_equal(row$contrast, "OJ - VC")
  expect_equal(row$statistic_label, "t")
  expect_true(is.numeric(row$statistic))
  expect_true(is.numeric(row$df))
  expect_true(is.numeric(row$ci_lower))
  expect_true(is.numeric(row$ci_upper))
  expect_equal(row$effect_size_name, "cohens_d")
  expect_equal(row$effect_size_method, "pairwise_pooled_sd")
  expect_equal(p$metadata$available_figures[[1]]$type, "pairwise_forest")
})

test_that("rstatix paired pairwise_t_test rows use paired counts and dz", {
  skip_if_not_installed("rstatix")

  paired_df <- data.frame(
    id = rep(1:8, times = 3),
    condition = factor(
      rep(c("before", "mid", "after"), each = 8),
      levels = c("before", "mid", "after")
    ),
    score = c(
      4.1, 5.2, 5.0, 6.4, 7.1, 7.8, 8.3, 9.1,
      5.0, 6.7, 6.1, 8.0, 8.3, 8.9, 9.7, 10.0,
      6.4, 7.1, 7.9, 8.3, 9.5, 10.4, 10.8, 11.6
    )
  )
  fit <- rstatix::pairwise_t_test(
    paired_df,
    score ~ condition,
    paired = TRUE,
    p.adjust.method = "holm"
  )
  p <- mellio_payload(fit)
  row <- p$fields$rows[[1]]
  before <- paired_df$score[paired_df$condition == "before"]
  mid <- paired_df$score[paired_df$condition == "mid"]
  labels <- vapply(p$fields$columns, function(col) col$label, character(1))

  expect_equal(p$fields$method, "paired t test")
  expect_true(isTRUE(row$paired))
  expect_equal(row$n_pairs, 8L)
  expect_null(row$n_1)
  expect_null(row$n_2)
  expect_equal(row$mean_1, mean(before))
  expect_equal(row$mean_2, mean(mid))
  expect_equal(row$effect_size_name, "cohens_dz")
  expect_equal(row$effect_size, mean(before - mid) / stats::sd(before - mid))
  expect_true("n pairs" %in% labels)
  expect_false("n 1" %in% labels)
  expect_match(p$fields$note, "Cohen's dz")
})

test_that("rstatix pairwise_wilcox_test infers paired Friedman follow-ups", {
  skip_if_not_installed("rstatix")

  paired_df <- data.frame(
    id = rep(1:8, times = 3),
    condition = factor(
      rep(c("before", "mid", "after"), each = 8),
      levels = c("before", "mid", "after")
    ),
    score = c(
      4, 5, 5, 6, 7, 7, 8, 9,
      5, 6, 6, 7, 8, 8, 9, 10,
      6, 7, 7, 8, 9, 9, 10, 11
    )
  )
  fit <- rstatix::pairwise_wilcox_test(
    paired_df,
    score ~ condition,
    paired = TRUE,
    p.adjust.method = "holm",
    exact = FALSE
  )
  p <- mellio_payload(fit)
  row <- p$fields$rows[[1]]
  before <- paired_df$score[paired_df$condition == "before"]
  mid <- paired_df$score[paired_df$condition == "mid"]
  labels <- vapply(p$fields$columns, function(col) col$label, character(1))

  expect_equal(p$fields$method, "Wilcoxon signed-rank test")
  expect_equal(row$contrast, "before - mid")
  expect_true(isTRUE(row$paired))
  expect_equal(row$n_pairs, 8L)
  expect_null(row$n_1)
  expect_null(row$n_2)
  expect_equal(row$median_difference, stats::median(before - mid))
  expect_equal(row$statistic_label, "V")
  expect_equal(row$statistic, fit$statistic[[1]])
  expect_equal(row$effect_size_name, "rank_biserial")
  expect_equal(row$effect_size, ms_signed_rank_biserial(before - mid))
  expect_true("n pairs" %in% labels)
  expect_false("n 1" %in% labels)
  expect_false("n 2" %in% labels)
  expect_match(p$fields$note, "median of paired differences")
  expect_match(p$fields$note, "rank-biserial")
})

test_that("mellio_open routes rstatix Wilcoxon tables to Stats", {
  skip_if_not_installed("rstatix")

  paired_df <- data.frame(
    id = rep(1:8, times = 2),
    condition = factor(rep(c("before", "after"), each = 8), levels = c("before", "after")),
    score = c(4, 5, 5, 6, 7, 7, 8, 9, 5, 6, 6, 7, 8, 8, 9, 10)
  )
  url <- withr::with_options(
    list(mellio.editor_url = "https://example.com"),
    mellio_open(
      rstatix::pairwise_wilcox_test(
        paired_df,
        score ~ condition,
        paired = TRUE,
        p.adjust.method = "holm",
        exact = FALSE
      ),
      browse = FALSE
    )
  )
  b64 <- sub("&.*$", "", sub(".*payload=", "", url))
  std <- chartr("-_", "+/", b64)
  pad <- (4 - nchar(std) %% 4) %% 4
  if (pad > 0) std <- paste0(std, strrep("=", pad))
  parsed <- jsonlite::fromJSON(rawToChar(jsonlite::base64_dec(std)), simplifyVector = FALSE)

  expect_equal(parsed$type, "pairwise_comparisons")
  expect_equal(parsed$fields$method, "Wilcoxon signed-rank test")
  expect_equal(parsed$fields$rows[[1]]$statistic_label, "V")
  expect_equal(parsed$fields$rows[[1]]$effect_size_name, "rank_biserial")
})

test_that("rstatix dunn_test payload produces Dunn post hoc rows", {
  skip_if_not_installed("rstatix")

  dunn_df <- data.frame(
    score = c(1.2, 1.5, 1.7, 2.1, 2.4,
              2.6, 2.8, 3.2, 3.5, 3.9,
              4.1, 4.4, 4.7, 5.0, 5.3),
    group = factor(rep(c("A", "B", "C"), each = 5))
  )
  fit <- rstatix::dunn_test(
    dunn_df,
    score ~ group,
    p.adjust.method = "holm",
    detailed = TRUE
  )
  p <- mellio_payload(fit)
  row <- p$fields$rows[[1]]
  group_b <- dunn_df$score[dunn_df$group == "B"]
  group_a <- dunn_df$score[dunn_df$group == "A"]
  labels <- vapply(p$fields$columns, function(col) col$label, character(1))

  expect_equal(p$type, "pairwise_comparisons")
  expect_equal(p$fields$method, "Dunn's Kruskal-Wallis post hoc test")
  expect_equal(p$fields$adjustment_method, "holm")
  expect_equal(p$fields$source, "rstatix::dunn_test")
  expect_equal(row$contrast, "A - B")
  expect_equal(row$n_1, 5L)
  expect_equal(row$n_2, 5L)
  expect_equal(row$median_1, stats::median(group_a))
  expect_equal(row$median_2, stats::median(group_b))
  expect_equal(row$median_difference, stats::median(group_a) - stats::median(group_b))
  expect_equal(row$rank_mean_difference, -fit$estimate[[1]])
  expect_equal(row$statistic_label, "z")
  expect_equal(row$statistic, -fit$statistic[[1]])
  expect_equal(row$p_raw, fit$p[[1]])
  expect_equal(row$p_value, fit$p.adj[[1]])
  expect_equal(row$effect_size_name, "cliffs_delta")
  expect_equal(row$effect_size, -1)
  expect_true("Mean rank difference" %in% labels)
  expect_true("p (raw)" %in% labels)
  expect_true("p (adjusted)" %in% labels)
  expect_true("Cliff's delta" %in% labels)
  expect_match(p$fields$note, "Dunn z statistics")
})

test_that("FSA dunnTest-like objects are normalized as pairwise rows", {
  fit <- structure(
    list(
      method = "Holm adjustment",
      res = data.frame(
        Comparison = c("A - B", "A - C"),
        Z = c(-1.25, -2.50),
        P.unadj = c(.211, .012),
        P.adj = c(.211, .024),
        stringsAsFactors = FALSE
      )
    ),
    class = "dunnTest"
  )
  p <- mellio_payload(fit)
  row <- p$fields$rows[[1]]

  expect_equal(p$type, "pairwise_comparisons")
  expect_equal(p$fields$source, "FSA::dunnTest")
  expect_equal(p$fields$adjustment_method, "holm")
  expect_equal(row$contrast, "A - B")
  expect_equal(row$statistic_label, "z")
  expect_equal(row$statistic, -1.25)
  expect_equal(row$p_raw, .211)
  expect_equal(row$p_value, .211)
})

test_that("dunn.test objects are normalized as pairwise rows", {
  fit <- structure(
    list(
      method = "bonferroni",
      comparisons = c("B - A", "C - A"),
      Z = c(1.25, 2.50),
      P = c(.211, .012),
      P.adjusted = c(.422, .024)
    ),
    class = "dunn.test"
  )
  p <- mellio_payload(fit)
  row <- p$fields$rows[[2]]

  expect_equal(p$type, "pairwise_comparisons")
  expect_equal(p$fields$source, "dunn.test::dunn.test")
  expect_equal(p$fields$adjustment_method, "bonferroni")
  expect_equal(row$contrast, "A - C")
  expect_equal(row$statistic_label, "z")
  expect_equal(row$statistic, -2.50)
  expect_equal(row$p_raw, .012)
  expect_equal(row$p_value, .024)
})

test_that("emmeans contrasts payload produces pairwise comparison rows", {
  skip_if_not_installed("emmeans")

  fit <- aov(mpg ~ factor(cyl), data = mtcars)
  contrasts <- emmeans::emmeans(fit, pairwise ~ cyl)$contrasts
  p <- mellio_payload(contrasts)

  expect_equal(p$type, "pairwise_comparisons")
  expect_equal(p$fields$method, "estimated marginal means contrasts")
  expect_equal(p$fields$adjustment_method, "tukey")
  expect_match(p$fields$note, "estimated marginal means \\(EMMs\\)")
  expect_equal(p$metadata$available_figures[[1]]$type, "pairwise_forest")
  expect_equal(p$figure_data$pairwise_forest$rows[[1]]$contrast, "4 - 6")
  expect_equal(length(p$fields$rows), 3L)
  expect_equal(p$fields$rows[[1]]$contrast, "4 - 6")
  expect_equal(p$fields$rows[[1]]$level_1, "4")
  expect_equal(p$fields$rows[[1]]$level_2, "6")
  expect_equal(p$fields$rows[[1]]$estimate, 6.920779, tolerance = 1e-6)
  expect_equal(p$fields$rows[[1]]$statistic, 4.441099, tolerance = 1e-6)
  expect_true(is.numeric(p$fields$rows[[1]]$ci_lower))
  expect_true(is.numeric(p$fields$rows[[1]]$p_value))

  p_alias <- mellio_payload(contrasts, .call = "em")
  expect_equal(p_alias$call, "emmeans::pairs(...)")
})

test_that("emmeans pairwise lists route to contrast payloads", {
  skip_if_not_installed("emmeans")

  fit <- aov(mpg ~ factor(cyl), data = mtcars)
  em <- emmeans::emmeans(fit, pairwise ~ cyl)
  p <- mellio_payload(em, .call = "emmeans::emmeans(fit, pairwise ~ cyl)")

  expect_equal(p$type, "pairwise_comparisons")
  expect_equal(p$call, "emmeans::emmeans(fit, pairwise ~ cyl)")
  expect_equal(p$fields$source, "emmeans::emm_list")
  expect_equal(p$fields$method, "estimated marginal means contrasts")
  expect_equal(p$fields$adjustment_method, "tukey")
  expect_equal(length(p$fields$rows), 3L)
  expect_equal(p$fields$rows[[1]]$contrast, "4 - 6")
  expect_equal(p$fields$rows[[1]]$level_1, "4")
  expect_equal(p$fields$rows[[1]]$level_2, "6")
  expect_equal(length(p$fields$emmeans_rows), 3L)
  expect_equal(p$fields$emmeans_rows[[1]]$cyl, "4")
  expect_true(is.numeric(p$fields$emmeans_rows[[1]]$emmean))
  expect_match(p$fields$emmeans_note, "Confidence level")
  expect_match(p$fields$note, "Family-wise 95% confidence intervals")
  expect_match(p$raw_output, "\\$emmeans")
  expect_match(p$raw_output, "\\$contrasts")
})

test_that("emmeans pairwise lists keep conditional mean labels", {
  skip_if_not_installed("emmeans")

  tg <- transform(ToothGrowth, dose = factor(dose))
  fit <- aov(len ~ dose * supp, data = tg)

  marginal <- suppressMessages(emmeans::emmeans(fit, pairwise ~ dose))
  p_marginal <- mellio_payload(marginal)
  expect_match(p_marginal$fields$emmeans_note, "misleading due to involvement in interactions")

  by_supp <- emmeans::emmeans(fit, pairwise ~ dose | supp)
  p_by_supp <- mellio_payload(by_supp)
  by_supp_cols <- vapply(p_by_supp$fields$emmeans_columns, function(col) col$key, character(1))
  expect_true(all(c("supp", "dose") %in% by_supp_cols))
  expect_equal(p_by_supp$fields$emmeans_rows[[1]]$supp, "OJ")
  expect_equal(p_by_supp$fields$emmeans_rows[[1]]$dose, "0.5")
  expect_equal(p_by_supp$fields$adjustment_method, "tukey")

  by_dose <- emmeans::emmeans(fit, pairwise ~ supp | dose)
  p_by_dose <- mellio_payload(by_dose)
  by_dose_cols <- vapply(p_by_dose$fields$emmeans_columns, function(col) col$key, character(1))
  expect_true(all(c("dose", "supp") %in% by_dose_cols))
  expect_equal(p_by_dose$fields$emmeans_rows[[1]]$dose, "0.5")
  expect_equal(p_by_dose$fields$emmeans_rows[[1]]$supp, "OJ")
  expect_equal(p_by_dose$fields$adjustment_method, "none")
  expect_match(p_by_dose$fields$note, "p values are unadjusted")
})

test_that("multcomp glht payload produces pairwise comparison rows", {
  skip_if_not_installed("multcomp")

  fit <- aov(mpg ~ factor(cyl), data = mtcars)
  g <- multcomp::glht(fit, linfct = multcomp::mcp(`factor(cyl)` = "Tukey"))
  p <- mellio_payload(g)
  p_alias <- mellio_payload(g, .call = "g")

  expect_equal(p$type, "pairwise_comparisons")
  expect_equal(p$fields$method, "Tukey contrasts")
  expect_equal(p$fields$adjustment_method, "single-step")
  expect_match(p$fields$note, "single-step method for Tukey contrasts")
  expect_equal(p$metadata$available_figures[[1]]$type, "pairwise_forest")
  expect_equal(p$figure_data$pairwise_forest$rows[[1]]$contrast, "4 - 6")
  expect_equal(length(p$fields$rows), 3L)
  expect_equal(p$fields$rows[[1]]$contrast, "4 - 6")
  expect_equal(p$fields$rows[[1]]$level_1, "4")
  expect_equal(p$fields$rows[[1]]$level_2, "6")
  expect_equal(p$fields$rows[[1]]$df, 29)
  expect_equal(p$fields$rows[[1]]$statistic_label, "t")
  expect_true(is.numeric(p$fields$rows[[1]]$estimate))
  expect_true("95% CI lower" %in% vapply(p$fields$columns, function(col) col$label, character(1)))
  expect_equal(p_alias$call, "multcomp::glht(...)")

  ps <- mellio_payload(summary(g))
  expect_equal(ps$type, "pairwise_comparisons")
  expect_equal(ps$fields$method, "Tukey contrasts")
  expect_equal(ps$metadata$available_figures[[1]]$type, "pairwise_forest")
  expect_equal(ps$figure_data$pairwise_forest$rows[[1]]$contrast, "4 - 6")
  expect_equal(length(ps$fields$rows), 3L)
  expect_equal(ps$fields$rows[[1]]$df, 29)
})
