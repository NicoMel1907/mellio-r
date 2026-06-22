test_that("only unified bridge names are exported", {
  exported <- getNamespaceExports("mellio")

  expect_true(all(c(
    "mellio_open",
    "mellio_payload",
    "mellio_to_json",
    "mellio_compare",
    "mellio_addin_send"
  ) %in% exported))

  expect_false(any(c(
    "melliofig",
    "mf_as_html",
    "mf_as_latex",
    "mf_as_markdown",
    "mf_copy",
    "ms_edit",
    "mt_edit",
    "mf_edit",
    "ms_payload",
    "ms_to_json",
    "ms_compare",
    "mt_compare",
    "ms_addin_send"
  ) %in% exported))
})

test_that("unified payload helpers smoke test", {
  payload <- mellio_payload(t.test(extra ~ group, data = sleep), .call = "payload_call")
  json <- mellio_to_json(payload)

  expect_s3_class(payload, "mellio_payload")
  expect_type(json, "character")
  expect_equal(jsonlite::fromJSON(json, simplifyVector = FALSE)$call, "payload_call")
})

test_that("unified comparison helper emits Stats payloads", {
  m1 <- lm(mpg ~ wt, data = mtcars)
  m2 <- lm(mpg ~ wt + hp, data = mtcars)
  payload <- mellio_compare(m1, m2, labels = c("Step 1", "Step 2"), .call = "compare_call")

  expect_s3_class(payload, "mellio_payload")
  expect_equal(payload$type, "hierarchical_regression_comparison")
  expect_equal(payload$call, "compare_call")
})

test_that("RStudio addin binding uses the unified public name", {
  dcf_path <- system.file("rstudio/addins.dcf", package = "mellio")
  if (!nzchar(dcf_path)) {
    dcf_path <- testthat::test_path("../../inst/rstudio/addins.dcf")
  }
  addin <- read.dcf(dcf_path)

  expect_equal(unname(addin[1, "Binding"]), "mellio_addin_send")
})
