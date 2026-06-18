test_that("melliotab.data.frame works", {
  df <- data.frame(
    Variable = c("Age", "Gender"),
    B = c(0.45, -1.23),
    SE = c(0.12, 0.34),
    t = c(3.75, -3.62),
    p = c(0.0003, 0.0004)
  )
  tbl <- melliotab(df, style = "apa7", title = "Test")
  expect_s3_class(tbl, "melliotab")
  expect_equal(tbl$style, "apa7")
  expect_equal(tbl$title, "Test")
  expect_equal(ncol(tbl$data), 5)
})

test_that("melliotab.NULL explains printed-but-not-returned results", {
  expect_error(melliotab(NULL), "received")
  expect_error(melliotab(NULL), "printed output")
  expect_error(melliotab(NULL), "data.frame")
  expect_error(melliotab(NULL), "capture.output")
  err <- tryCatch(melliotab(NULL), error = conditionMessage)
  expect_false(grepl("mediation", err, fixed = TRUE))
})

test_that("melliotab.character parses simple captured key-value output", {
  lines <- c(
    "motive_enjoyment       r = +0.019, p = 0.5573",
    "motive_ltb             r = +0.076, p = 0.0180",
    "motive_appearance      r = +0.078, p = 0.0149"
  )
  tbl <- melliotab(lines, style = "apa7", title = "SES correlations")

  expect_s3_class(tbl, "melliotab")
  expect_equal(tbl$title, "SES correlations")
  expect_equal(tbl$source, "capture.output()")
  expect_equal(names(tbl$raw_data), c("Variable", "r", "p"))
  expect_equal(tbl$raw_data$Variable[[1]], "motive_enjoyment")
  expect_equal(tbl$raw_data$r[[2]], 0.076)
  expect_equal(tbl$raw_data$p[[3]], 0.0149)
})

test_that("melliotab.character rejects prose that is not table-shaped", {
  expect_error(melliotab("this is just a sentence"), "could not parse")
})

test_that("melliotab.summaryDefault turns base summary output into a table", {
  s <- summary(c(1, 2, 3, 4, 5, NA))
  tbl <- melliotab(s, style = "apa7", title = "SES")

  expect_s3_class(tbl, "melliotab")
  expect_equal(tbl$title, "SES")
  expect_equal(tbl$source, "base R summary()")
  expect_true(all(c("Min", "Q1", "Median", "Mean", "Q3", "Max", "Missing") %in% names(tbl$raw_data)))
  expect_equal(tbl$raw_data$Missing[[1]], 1L)
  expect_equal(tbl$data$Missing[[1]], "1")
})

test_that("melliotab.table turns base frequency tables into table cards", {
  x <- table(SES_band = c("Under $50k", "$50k-$99k", "Under $50k"),
             useNA = "ifany")
  tbl <- melliotab(x, style = "apa7")

  expect_s3_class(tbl, "melliotab")
  expect_equal(tbl$title, "Frequency table")
  expect_equal(tbl$source, "base R table()")
  expect_true(all(c("SES_band", "n") %in% names(tbl$raw_data)))
  expect_equal(sum(tbl$raw_data$n), 3L)
  expect_true(all(grepl("^[0-9]+$", tbl$data$n)))
})

test_that("melliotab.table names unnamed one-way tables as Category", {
  tbl <- melliotab(table(c("a", "b", "a")), style = "apa7")

  expect_true("Category" %in% names(tbl$raw_data))
  expect_equal(sum(tbl$raw_data$n), 3L)
})

test_that("melliotab.lm works", {
  model <- lm(mpg ~ wt + hp, data = mtcars)
  tbl <- melliotab(model, style = "apa7", title = "Regression")
  expect_s3_class(tbl, "melliotab")
  expect_true(!is.null(tbl$note))
  expect_true(grepl("R\u00B2", tbl$note))
})

test_that("melliotab compares multiple models side by side", {
  m1 <- lm(mpg ~ wt, data = mtcars)
  m2 <- lm(mpg ~ wt + hp, data = mtcars)

  tbl <- melliotab(m1, m2, labels = c("Step 1", "Step 2"))

  expect_s3_class(tbl, "melliotab")
  expect_equal(tbl$title, "Model comparison (2 models)")
  expect_true(all(c("Step 1", "Step 2") %in% names(tbl$raw_data)))
  expect_true("hp" %in% tbl$raw_data[[" "]])
  expect_true(isTRUE(tbl$options$is_comparison))
})

test_that("model comparison tables carry stargazer-style structure", {
  skip_if_not_installed("broom")
  skip_if_not_installed("gt")

  m1 <- lm(Ozone ~ Temp, data = airquality)
  m2 <- lm(Ozone ~ Temp + Wind, data = airquality)
  m3 <- lm(Ozone ~ Temp + Wind + Solar.R, data = airquality)

  tbl <- melliotab(
    m1, m2, m3,
    labels = c("Step 1", "Step 2", "Step 3"),
    dep.var.labels = "Ozone concentration"
  )
  ed <- mellio:::build_editor_data(tbl)
  stat_row <- sum(ed$rowRoles == "spanheader") +
    tbl$options$comparison_stat_start

  expect_equal(
    vapply(tbl$spanners, `[[`, character(1), "label"),
    c("Dependent variable:", "Ozone concentration")
  )
  expect_equal(ed$rowRoles[seq_len(3)], c("spanheader", "spanheader", "header"))
  expect_true(isTRUE(ed$cellFormats[["0,1"]]$italic))
  expect_false(ed$cellFormats[["0,1"]]$bold)
  expect_true(isTRUE(ed$cellBorders[["0,1"]]$bottom))
  expect_true(isTRUE(ed$cellBorders[[paste0(stat_row, ",0")]]$top))

  gt_tbl <- mt_as_gt(tbl)
  boxhead <- gt_tbl$`_boxhead`
  model_widths <- boxhead$column_width[
    match(c("Step 1", "Step 2", "Step 3"), boxhead$var)
  ]
  model_aligns <- boxhead$column_align[
    match(c("Step 1", "Step 2", "Step 3"), boxhead$var)
  ]
  expect_length(unique(model_widths), 1)
  expect_equal(unique(model_aligns), "center")

  html <- gt::as_raw_html(gt_tbl)
  expect_match(html, "id=\"mellio_dependent_variable\"", fixed = TRUE)
  expect_match(html, "width: 72%;", fixed = TRUE)
  expect_match(html, "id=\"mellio_dependent_variable_label\"", fixed = TRUE)
  expect_match(html, "font-weight: bold", fixed = TRUE)
})

test_that("melliotab.aov works", {
  aov_model <- aov(mpg ~ factor(cyl), data = mtcars)
  tbl <- melliotab(aov_model, style = "apa7")
  expect_s3_class(tbl, "melliotab")
  expect_true("\u03B7\u00B2" %in% names(tbl$data))
})

test_that("melliotab.mediate works", {
  skip_if_not_installed("mediation")
  set.seed(1)
  d <- data.frame(
    x = rnorm(60),
    m = rnorm(60),
    y = rnorm(60)
  )
  d$m <- 0.4 * d$x + d$m
  d$y <- 0.2 * d$x + 0.5 * d$m + d$y

  res <- mediation::mediate(
    lm(m ~ x, data = d),
    lm(y ~ x + m, data = d),
    treat = "x",
    mediator = "m",
    boot = FALSE,
    sims = 20
  )
  tbl <- melliotab(res, style = "apa7", title = "Mediation")
  expect_s3_class(tbl, "melliotab")
  expect_equal(tbl$title, "Mediation")
  expect_true("Effect" %in% names(tbl$data))
  expect_true(grepl("mediation::mediate", tbl$note))
})

test_that("melliotab.matrix works for correlation", {
  cor_mat <- cor(mtcars[, c("mpg", "wt", "hp")])
  tbl <- melliotab(cor_mat, style = "apa7", diagonal = "dash", triangle = "lower")
  expect_s3_class(tbl, "melliotab")
  expect_true(tbl$options$is_correlation)
  # Diagonal should be em-dash
  expect_equal(tbl$data[1, 2], "\u2014")
  # Above diagonal should be blank
  expect_equal(tbl$data[1, 3], "")
})

test_that("mt_set_style changes formatting", {
  df <- data.frame(Variable = "X", p = 0.045)
  tbl <- melliotab(df, style = "apa7")
  # APA removes leading zeros
  expect_equal(tbl$data$p[1], ".045")

  tbl_ieee <- mt_set_style(tbl, "ieee")
  # IEEE keeps leading zeros
  expect_equal(tbl_ieee$data$p[1], "0.045")
})

test_that("mt_sig_stars adds stars", {
  df <- data.frame(
    Variable = c("A", "B", "C"),
    B = c(1.5, 2.3, 0.1),
    p = c(0.001, 0.03, 0.5)
  )
  tbl <- melliotab(df, style = "apa7") |> mt_sig_stars()
  # B column should have stars, p column removed
  expect_false("p" %in% names(tbl$data))
  expect_true(grepl("\\*", tbl$data$B[1]))
})

test_that("format_apa_number handles p-values", {
  expect_equal(format_apa_number("0.0002", "pvalue"), "< .001")
  expect_equal(format_apa_number("0.045", "pvalue"), ".045")
  expect_equal(format_apa_number("0.597", "pvalue"), ".597")
  expect_equal(format_apa_number(NA, "pvalue"), "")
})

test_that("format_apa_number respects custom p_decimals", {
  expect_equal(
    format_apa_number("0.045", "pvalue", p_decimals = 4L),
    ".0450"
  )
  expect_equal(
    format_apa_number("0.045", "pvalue", p_decimals = 2L),
    ".04"
  )
})

test_that("format_apa_number handles estimates", {
  expect_equal(format_apa_number("3.756", "estimate"), "3.76")
  expect_equal(format_apa_number("-0.123", "statistic"), "-0.12")
})

test_that("detect_column_types identifies known headers", {
  types <- detect_column_types(c("Variable", "B", "SE", "t", "p"))
  expect_equal(types, c("stub", "estimate", "estimate", "statistic", "pvalue"))

  r2_types <- detect_column_types(c("R\u00B2", "Adjusted R\u00B2", "\u0394R\u00B2"))
  expect_equal(r2_types, c("estimate", "estimate", "estimate"))
  expect_true(all(detect_leading_zero_cols(c("R\u00B2", "Adjusted R\u00B2", "\u0394R\u00B2"))))

  adjusted_p_headers <- c("p (adjusted)", "p.adj", "p.adjust", "adjusted p")
  expect_equal(detect_column_types(adjusted_p_headers),
               rep("pvalue", length(adjusted_p_headers)))
  expect_true(all(detect_leading_zero_cols(adjusted_p_headers)))
  expect_equal(detect_p_value_col(c("Contrast", "p (adjusted)")), 2L)

  dunn_headers <- c(".y.", "group1", "group2", "n1", "n2", "estimate")
  expect_equal(detect_column_types(dunn_headers),
               c("stub", "stub", "stub", "integer", "integer", "estimate"))
  expect_equal(detect_column_types(c("Contrast", "n pairs", "p (adjusted)")),
               c("stub", "integer", "pvalue"))
  expect_false(any(detect_leading_zero_cols(c("n1", "n_2", "n pairs"))))
})

test_that("is_stat_symbol includes SE for APA header italics", {
  expect_true(is_stat_symbol("SE"))
  expect_true(is_stat_symbol("p (adjusted)"))
  expect_true(is_stat_symbol("n1"))
  expect_true(is_stat_symbol("n pairs"))
})

test_that("exports produce output", {
  model <- lm(mpg ~ wt, data = mtcars)
  tbl <- melliotab(model, style = "apa7", title = "Test") |> mt_number(1)

  md <- mt_as_markdown(tbl)
  expect_true(nchar(md) > 50)
  expect_true(grepl("\\|", md))

  ltx <- mt_as_latex(tbl)
  expect_true(grepl("\\\\begin\\{table\\}", ltx))
  expect_true(grepl("\\\\toprule", ltx))

  html <- mt_as_html(tbl)
  expect_true(nchar(html) > 100)
  expect_true(grepl("<table", html))
})

test_that("IEEE tables use Roman-numbered captions above the table", {
  tbl <- melliotab(
    data.frame(Variable = "X", Value = 1),
    style = "ieee",
    title = "Sample Characteristics"
  ) |> mt_number(1)

  expect_true(tbl$style_config$table_title$separate_line)
  expect_equal(tbl$style_config$table_label$numbering, "roman")

  md <- mt_as_markdown(tbl)
  md_lines <- strsplit(md, "\n", fixed = TRUE)[[1]]
  expect_equal(md_lines[1], "**TABLE I**")
  expect_equal(md_lines[2], "SAMPLE CHARACTERISTICS")

  ltx <- mt_as_latex(tbl)
  expect_true(grepl("TABLE I", ltx, fixed = TRUE))
  expect_true(grepl("SAMPLE CHARACTERISTICS", ltx, fixed = TRUE))

  title_html <- build_title_html(tbl)
  expect_true(grepl("TABLE I", title_html, fixed = TRUE))
  expect_true(grepl("SAMPLE CHARACTERISTICS", title_html, fixed = TRUE))
})

test_that("APA table notes start prose after Note with a capital", {
  tbl <- melliotab(
    data.frame(Variable = "X", Value = 1),
    style = "apa7",
    note = "this note should read like a sentence"
  )
  md <- mt_as_markdown(tbl)
  expect_true(grepl("*Note.* This note should read like a sentence", md, fixed = TRUE))

  stat_tbl <- melliotab(
    data.frame(Variable = "X", Value = 1),
    style = "apa7",
    note = "p < .05"
  )
  expect_true(grepl("*Note.* p < .05", mt_as_markdown(stat_tbl), fixed = TRUE))
})

test_that("mt_save works for supported formats", {
  model <- lm(mpg ~ wt, data = mtcars)
  tbl <- melliotab(model, style = "apa7", title = "Test") |> mt_number(1)

  tmp_html <- tempfile(fileext = ".html")
  tmp_md <- tempfile(fileext = ".md")
  tmp_tex <- tempfile(fileext = ".tex")

  mt_save(tbl, tmp_html)
  mt_save(tbl, tmp_md)
  mt_save(tbl, tmp_tex)

  expect_true(file.exists(tmp_html))
  expect_true(file.exists(tmp_md))
  expect_true(file.exists(tmp_tex))
  expect_error(mt_save(tbl, tempfile(fileext = ".docx")), "Unsupported")
})

test_that("available styles produce valid HTML", {
  model <- lm(mpg ~ wt, data = mtcars)
  tbl <- melliotab(model, style = "apa7", title = "Test") |> mt_number(1)

  for (s in list_styles()) {
    tbl_s <- mt_set_style(tbl, s)
    html <- mt_as_html(tbl_s)
    expect_true(nchar(html) > 100, label = paste("HTML for", s))
  }
})

test_that("table_to_tsv produces tab-delimited text", {
  df <- data.frame(Variable = c("A", "B"), B = c(1.5, 2.3), p = c(0.01, 0.05))
  tsv <- table_to_tsv(df)
  lines <- strsplit(tsv, "\n")[[1]]
  expect_equal(length(lines), 3)  # header + 2 rows
  expect_equal(lines[1], "Variable\tB\tp")
  expect_true(grepl("\t", lines[2]))
})

test_that("viewer copy replaces fragile table edge borders with Word-safe rules", {
  html <- mellio_toolbar_html(mode = "table")
  expect_true(grepl("ruleRow", html, fixed = TRUE))
  expect_true(grepl("clearTopBorders", html, fixed = TRUE))
  expect_true(grepl("clearBottomBorders", html, fixed = TRUE))
  expect_true(grepl("thead.insertBefore", html, fixed = TRUE))
  expect_true(grepl("tbody.appendChild", html, fixed = TRUE))
  expect_false(grepl('borderTopStyle="none";c.style.borderBottomStyle="none"', html, fixed = TRUE))
})

test_that("viewer toolbar stays local and quiet", {
  html <- mellio_toolbar_html(mode = "table", edit_url = "https://mellio.app/edit?id=test")
  expect_true(grepl("mellio-tab-logo", html, fixed = TRUE))
  expect_true(grepl("Open in Mellio", html, fixed = TRUE))
  expect_true(grepl("data-tooltip", html, fixed = TRUE))
  expect_true(grepl("background: #ffffff", html, fixed = TRUE))
  expect_false(grepl("fonts.googleapis.com", html, fixed = TRUE))
  expect_false(grepl("mellio-copy-btn", html, fixed = TRUE))
  expect_false(grepl("mellio-wordmark", html, fixed = TRUE))
  expect_false(grepl("mellio-open-label", html, fixed = TRUE))
  expect_false(grepl(">HTML<", html, fixed = TRUE))
  expect_false(grepl(">LaTeX<", html, fixed = TRUE))
  expect_false(grepl(">Markdown<", html, fixed = TRUE))
})

test_that("viewer resolver respects configured non-RStudio Viewer functions", {
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      isTRUE(rstudioapi::isAvailable())) {
    skip("RStudio viewer takes precedence in the IDE.")
  }
  custom_viewer <- function(url) url
  withr::local_options(list(viewer = custom_viewer))
  expect_identical(mellio_viewer_function(), custom_viewer)
})
