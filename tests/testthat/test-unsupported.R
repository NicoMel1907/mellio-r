test_that("mellio_payload.default returns an unsupported card for unknown classes", {
  x <- structure(list(a = 1), class = "definitely_not_a_known_mellio_class")

  p <- mellio_payload(x)

  expect_s3_class(p, "mellio_payload")
  expect_equal(p$card_kind, "unsupported")
  expect_equal(p$type, "unsupported")
  expect_equal(p$type_label, "Unrecognized object")
  expect_true("definitely_not_a_known_mellio_class" %in% p$fields$class)
  expect_true(nchar(p$fields$printed) > 0)
  expect_match(p$fields$message, "could not safely extract")
  expect_null(p$fields$suggestions)
  expect_identical(p$raw_output, p$fields$printed)
})

test_that("unsupported payloads suggest close supported alternatives", {
  x <- structure(list(a = 1), class = c("summary.mystery_model", "list"))

  p <- mellio_payload(x, .call = "summary(fit)")

  expect_equal(p$card_kind, "unsupported")
  expect_true(length(p$fields$suggestions) >= 1L)
  expect_equal(p$fields$suggestions[[1]]$title, "Try the original model object")
  expect_equal(p$fields$suggestions[[1]]$code, "mellio_open(fit)")
})

test_that("unsupported suggestion helper recognizes ANOVA and emmeans near misses", {
  anova_text <- paste(
    "            Df Sum Sq Mean Sq F value Pr(>F)",
    "dose         2 2426.4  1213.2    92.0 <2e-16",
    "Residuals   54  712.1    13.2",
    sep = "\n"
  )
  emmeans_text <- paste(
    "$emmeans",
    " dose emmean SE df lower.CL upper.CL",
    "$contrasts",
    " contrast estimate SE df t.ratio p.value",
    sep = "\n"
  )

  anova_suggestions <- ms_unsupported_suggestions(
    cls = "mystery_anova",
    call_str = NA_character_,
    printed_text = anova_text
  )
  emmeans_suggestions <- ms_unsupported_suggestions(
    cls = "mystery_emmeans",
    call_str = NA_character_,
    printed_text = emmeans_text
  )

  expect_true(any(vapply(anova_suggestions, function(item) {
    identical(item$code, "mellio_open(anova(fit))")
  }, logical(1))))
  expect_true(any(vapply(emmeans_suggestions, function(item) {
    identical(item$code, "mellio_open(result$contrasts)")
  }, logical(1))))
})

test_that("mellio_payload.default partially extracts unknown broom-supported classes", {
  result <- nls(
    rate ~ Vm * conc / (K + conc),
    data = Puromycin,
    start = list(Vm = 200, K = 0.05)
  )

  p <- mellio_payload(result, .call = "nls(...)")

  expect_s3_class(p, "mellio_payload")
  expect_equal(p$card_kind, "table")
  expect_equal(p$type, "partial_broom_summary")
  expect_equal(p$fields$parse_state, "partial")
  expect_equal(p$fields$adapter, "broom::tidy")
  expect_true("nls" %in% p$fields$class)
  expect_true(length(p$fields$rows) >= 2L)
  expect_true(any(vapply(p$fields$columns, function(col) {
    identical(col$key, "term")
  }, logical(1))))
  expect_match(p$fields$message, "dedicated adapter")
  expect_match(p$fields$table_note, "verify")
  expect_identical(p$raw_output, p$fields$printed)
})

test_that("melliotab uses the partial payload fallback for broom-supported classes", {
  result <- nls(
    rate ~ Vm * conc / (K + conc),
    data = Puromycin,
    start = list(Vm = 200, K = 0.05)
  )

  tbl <- melliotab(result)

  expect_s3_class(tbl, "melliotab")
  expect_equal(tbl$title, "Partially parsed R result")
  expect_s3_class(tbl$model, "mellio_payload")
  expect_equal(tbl$model$type, "partial_broom_summary")
  expect_true("term" %in% names(tbl$raw_data))
})

test_that("unsupported payload serialises to JSON", {
  json <- mellio_to_json(mellio_payload(structure(list(), class = "foo")))
  parsed <- jsonlite::fromJSON(as.character(json), simplifyVector = FALSE)

  expect_true(nchar(as.character(json)) > 0)
  expect_equal(parsed$card_kind, "unsupported")
  expect_equal(parsed$type, "unsupported")
  expect_equal(parsed$fields$class[[1]], "foo")
  expect_type(parsed$fields$printed, "character")
})

test_that("unsupported table defaults include paste-flow tips", {
  expect_error(melliotab(structure(list(), class = "not_a_table")), "Tables in Mellio")
})
