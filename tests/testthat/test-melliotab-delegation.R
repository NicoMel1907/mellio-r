test_that("melliotab delegates TukeyHSD payloads to table projection", {
  fit <- aov(mpg ~ factor(cyl), data = mtcars)
  tbl <- melliotab(TukeyHSD(fit))

  expect_s3_class(tbl, "melliotab")
  expect_equal(tbl$title, "Pairwise comparisons")
  expect_match(tbl$note, "Tukey")
  expect_true("Contrast" %in% names(tbl$raw_data))
  expect_equal(nrow(tbl$raw_data), 3L)
})

test_that("melliotab delegates pairwise.htest payloads to table projection", {
  tbl <- melliotab(pairwise.t.test(mtcars$mpg, mtcars$cyl))

  expect_s3_class(tbl, "melliotab")
  expect_match(tbl$note, "Holm")
  expect_true("Contrast" %in% names(tbl$raw_data))
  expect_true("p (adjusted)" %in% names(tbl$raw_data))
  expect_equal(detect_column_types(names(tbl$raw_data))[match("p (adjusted)", names(tbl$raw_data))], "pvalue")
  expect_equal(nrow(tbl$raw_data), 3L)
})

test_that("melliotab delegates emmeans and glht pairwise payloads when available", {
  skip_if_not_installed("emmeans")
  skip_if_not_installed("multcomp")

  fit <- aov(mpg ~ factor(cyl), data = mtcars)
  em <- emmeans::emmeans(fit, pairwise ~ cyl)$contrasts
  tbl_em <- melliotab(em)
  expect_s3_class(tbl_em, "melliotab")
  expect_true("Contrast" %in% names(tbl_em$raw_data))

  g <- multcomp::glht(fit, linfct = multcomp::mcp(`factor(cyl)` = "Tukey"))
  tbl_g <- melliotab(g)
  expect_s3_class(tbl_g, "melliotab")
  expect_true("Contrast" %in% names(tbl_g$raw_data))
})

test_that("melliotab can project inline payloads", {
  payload <- mellio_payload(t.test(extra ~ group, data = sleep))
  tbl <- melliotab(payload)

  expect_s3_class(tbl, "melliotab")
  expect_equal(tbl$title, payload$type_label)
  expect_true("p" %in% names(tbl$raw_data))
})

test_that("melliotab projects hierarchical comparisons with df pairs and full predictors", {
  skip_if_not_installed("broom")

  m1 <- lm(mpg ~ wt + hp, data = mtcars)
  m2 <- lm(mpg ~ wt + hp + cyl, data = mtcars)
  payload <- mellio_compare(m1, m2)
  tbl <- melliotab(payload)

  expect_s3_class(tbl, "melliotab")
  # Normalize the encoding marker before comparison. In C/POSIX locale
  # sessions, the production-side column names may end up tagged "UTF-8"
  # while the test literals here end up tagged "unknown" (or vice versa),
  # even when the bytes are identical. Encoding<- changes the marker
  # without touching bytes, so this comparison is locale-tolerant.
  actual_names <- names(tbl$raw_data)
  expected_names <- c("Model", "Predictors", "Terms entered", "R\u00B2",
                      "Adjusted R\u00B2", "\u0394R\u00B2", "F change",
                      "df", "p", "n")
  Encoding(actual_names) <- "UTF-8"
  Encoding(expected_names) <- "UTF-8"
  expect_equal(actual_names, expected_names)
  expect_equal(tbl$raw_data$Predictors[[1]], "wt, hp")
  expect_equal(tbl$raw_data$Predictors[[2]], "wt, hp, cyl")
  expect_equal(tbl$raw_data[["Terms entered"]][[1]], "wt, hp")
  expect_equal(tbl$raw_data[["Terms entered"]][[2]], "cyl")
  expect_equal(tbl$raw_data$df[[2]], "(1, 28)")
  expect_true(is.numeric(tbl$raw_data[["Adjusted R\u00B2"]]))
  expect_true(is.numeric(tbl$raw_data[["\u0394R\u00B2"]]))
  expect_match(tbl$note, "Terms entered")
})

test_that("melliotab projects all rows from raw anova model comparisons", {
  m1 <- lm(mpg ~ wt, data = mtcars)
  m2 <- lm(mpg ~ wt + hp, data = mtcars)
  m3 <- lm(mpg ~ wt + hp + factor(am), data = mtcars)

  p <- mellio_payload(anova(m1, m2, m3))
  tbl <- melliotab(p)

  expect_equal(p$fields$table_type, "model_comparison")
  expect_length(p$fields$rows, 2)
  expect_equal(vapply(p$fields$rows, `[[`, character(1), "comparison"),
               c("model 1 vs. model 2", "model 2 vs. model 3"))
  expect_true("Residual SS" %in% names(tbl$raw_data))
  expect_equal(nrow(tbl$raw_data), 2L)
})

test_that("structural payload projection lists and extracts table sections", {
  payload <- structure(
    list(
      card_kind = "structural",
      type = "example_structural",
      type_label = "Example structural model",
      call = "sem(...)",
      fields = list(
        report_zone = list(
          fit_indices = list(
            list(name = "CFI", value = 0.95),
            list(name = "RMSEA", value = 0.04, ci = I(c(0.01, 0.07)), ci_level = 0.90)
          )
        ),
        inspection_zone = list(
          parameters = list(
            list(lhs = "visual", op = "=~", rhs = "x1", estimate = 0.7,
                 std_error = 0.1, statistic = 7, p_value = 0.001,
                 ci_lower = 0.5, ci_upper = 0.9, std_estimate = 0.8),
            list(lhs = "y", op = "~", rhs = "visual", estimate = 0.4,
                 std_error = 0.1, statistic = 4, p_value = 0.002)
          ),
          reliability = list(
            list(factor = "visual", omega = 0.82, ave = 0.61, n_indicators = 3L)
          )
        )
      )
    ),
    class = c("mellio_payload", "list")
  )

  expect_error(melliotab(payload), "several tables")
  expect_error(melliotab(payload), "loadings")

  loadings <- melliotab(payload, section = "loadings")
  expect_equal(loadings$title, "Factor loadings")
  expect_equal(nrow(loadings$raw_data), 1L)
  expect_equal(loadings$raw_data$Parameter, "visual -> x1")

  fit <- melliotab(payload, section = "fit")
  expect_equal(fit$title, "Model fit indices")
  expect_equal(nrow(fit$raw_data), 2L)
  expect_true("Fit index" %in% names(fit$raw_data))

  fit_alias <- melliotab(payload, what = "fit_indices")
  expect_equal(fit_alias$title, "Model fit indices")
  expect_equal(nrow(fit_alias$raw_data), 2L)

  rel <- melliotab(payload, section = "reliability")
  expect_equal(rel$title, "Reliability estimates")
  expect_true("\u03C9" %in% names(rel$raw_data))
})
