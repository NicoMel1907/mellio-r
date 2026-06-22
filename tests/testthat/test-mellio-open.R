decode_stats_payload_url <- function(url) {
  b64 <- sub(".*payload=", "", url)
  b64 <- sub("&.*$", "", b64)
  std <- chartr("-_", "+/", b64)
  pad <- (4 - nchar(std) %% 4) %% 4
  if (pad > 0) std <- paste0(std, strrep("=", pad))
  jsonlite::fromJSON(rawToChar(jsonlite::base64_dec(std)),
                     simplifyVector = FALSE)
}

stable_payload_bits <- function(payload) {
  payload[c("schema_version", "card_kind", "type", "type_label", "call", "fields")]
}

decode_hash_params <- function(url) {
  hash <- sub("^[^#]*#?", "", url)
  parts <- strsplit(hash, "&", fixed = TRUE)[[1]]
  out <- list()
  for (part in parts) {
    bits <- strsplit(part, "=", fixed = TRUE)[[1]]
    key <- bits[[1]]
    value <- if (length(bits) > 1L) {
      paste(bits[-1], collapse = "=")
    } else {
      ""
    }
    out[[key]] <- utils::URLdecode(value)
  }
  out
}

decode_base64_json_param <- function(url, key) {
  hash <- sub("^[^#]*#?", "", url)
  pattern <- paste0("(?:^|&)", key, "=([^&]+)")
  hit <- regmatches(hash, regexec(pattern, hash, perl = TRUE))[[1]]
  value <- if (length(hit) >= 2L) hit[[2]] else NULL
  if (is.null(value) || !nzchar(value)) return(NULL)
  value <- utils::URLdecode(value)
  jsonlite::fromJSON(rawToChar(jsonlite::base64_dec(value)),
                     simplifyVector = FALSE)
}

test_that("mellio_open builds stats payload URLs", {
  withr::with_options(list(
    mellio.editor_url = "https://example.com"
  ), {
    url <- mellio_open(t.test(extra ~ group, data = sleep), browse = FALSE)
    payload <- decode_stats_payload_url(url)

    expect_match(url, "^https://example.com/#stats/payload=")
    expect_equal(payload$type, "welch_t_test")
    expect_equal(payload$call, "t.test(extra ~ group, data = sleep)")
    expect_equal(stable_payload_bits(payload)$card_kind, "inline")
  })
})

test_that("custom Mellio URLs are validated", {
  withr::local_options(list(mellio.editor_url = NULL))
  withr::local_envvar(c(MELLIO_URL = ""))

  withr::with_options(list(mellio.editor_url = "javascript:alert(1)"), {
    expect_error(
      mellio_open(t.test(extra ~ group, data = sleep), browse = FALSE),
      "must start with"
    )
  })

  withr::with_options(list(mellio.editor_url = "https://example.com/#x"), {
    expect_error(
      mellio_open(t.test(extra ~ group, data = sleep), browse = FALSE),
      "must not include"
    )
  })

  withr::with_options(list(mellio.editor_url = "https://example.com/?x=1"), {
    expect_error(
      mellio_open(mtcars[1:2, 1:2], browse = FALSE),
      "must not include"
    )
  })
})

test_that("mellio_open on a bare data.frame routes to Tables, not Stats", {
  url <- mellio_open(mtcars[1:3, 1:3], title = "MTCars rows", browse = FALSE)

  expect_match(url, "#data=", fixed = TRUE)
  expect_no_match(url, "#stats/payload=", fixed = TRUE)
  expect_match(url, "provenance=", fixed = TRUE)
  provenance <- decode_base64_json_param(url, "provenance")
  expect_equal(provenance$r_version, R.version.string)
})

test_that("mellio_open on melliotab built from a data.frame routes to Tables", {
  tbl <- melliotab(mtcars[1:2, 1:2], title = "Prepared table")
  url <- mellio_open(tbl, browse = FALSE)

  expect_match(url, "#data=", fixed = TRUE)
  expect_no_match(url, "#stats/payload=", fixed = TRUE)
  expect_match(url, "provenance=", fixed = TRUE)
})

test_that("mellio_open compares multiple models in Tables", {
  m1 <- lm(mpg ~ wt, data = mtcars)
  m2 <- lm(mpg ~ wt + hp, data = mtcars)

  url <- mellio_open(
    m1, m2,
    labels = c("Step 1", "Step 2"),
    dep.var.labels = "Fuel efficiency",
    browse = FALSE
  )
  formats <- decode_base64_json_param(url, "cellFormats")
  borders <- decode_base64_json_param(url, "cellBorders")

  expect_match(url, "#data=", fixed = TRUE)
  expect_no_match(url, "#stats/payload=", fixed = TRUE)
  expect_match(url, "title=Model%20comparison%20%282%20models%29", fixed = TRUE)
  expect_true(isTRUE(formats[["0,1"]]$italic))
  expect_false(formats[["0,1"]]$bold)
  expect_true(isTRUE(borders[["0,1"]]$bottom))
})

test_that("mellio_open keeps positional browse compatibility", {
  url <- mellio_open(mtcars[1:2, 1:2], FALSE)

  expect_match(url, "#data=", fixed = TRUE)
  expect_no_match(url, "#stats/payload=", fixed = TRUE)
})

test_that("mellio_open on melliotab built from an lm model routes to Stats", {
  fit <- lm(mpg ~ wt, data = mtcars)
  tbl <- melliotab(fit)
  url <- mellio_open(tbl, browse = FALSE)

  expect_match(url, "#stats/payload=", fixed = TRUE)
  expect_no_match(url, "#data=", fixed = TRUE)
})

test_that("direct Tables URLs preserve metadata and provenance", {
  tbl <- melliotab(
    mtcars[1:2, 1:2],
    title = "Prepared table",
    note = "A compact note.",
    source = "R"
  )
  tbl$number <- 4
  tbl$provenance <- list(
    r_version = "R version test",
    working_dir = getwd(),
    sender = list(user = "tester", host = "test-host")
  )
  url <- mellio:::mellio_table_url(tbl)
  params <- decode_hash_params(url)
  provenance <- decode_base64_json_param(url, "provenance")

  expect_equal(params$title, "Prepared table")
  expect_equal(params$note, "A compact note.")
  expect_equal(params$source, "R")
  expect_equal(params$tableNumber, "4")
  expect_equal(params$style, "APA7")
  expect_equal(provenance$r_version, "R version test")
  expect_equal(provenance$sender$user, "tester")
})

test_that("mellio_open on a raw stats::anova table routes to Stats without figures", {
  fit <- lm(mpg ~ factor(cyl) * factor(am), data = mtcars)
  url <- mellio_open(anova(fit), browse = FALSE)
  payload <- decode_stats_payload_url(url)

  expect_match(url, "#stats/payload=", fixed = TRUE)
  expect_no_match(url, "#data=", fixed = TRUE)
  expect_equal(payload$type, "anova_single_model")
  expect_equal(
    vapply(payload$fields$all_terms, function(row) row$term, character(1)),
    c("factor(cyl)", "factor(am)", "factor(cyl):factor(am)")
  )
  expect_null(payload$metadata$available_figures)
  expect_null(payload$figure_data)
})

test_that("mellio_open on car::Anova Type II and III tables routes to Stats", {
  skip_if_not_installed("car")

  fit <- lm(mpg ~ factor(cyl) + wt + hp, data = mtcars)
  type_ii_url <- mellio_open(car::Anova(fit, type = 2), browse = FALSE)
  type_ii_payload <- decode_stats_payload_url(type_ii_url)

  expect_match(type_ii_url, "#stats/payload=", fixed = TRUE)
  expect_no_match(type_ii_url, "#data=", fixed = TRUE)
  expect_equal(type_ii_payload$fields$ss_type, "type_ii")
  expect_equal(type_ii_payload$fields$model_kind, "ancova")
  expect_equal(type_ii_payload$fields$term, "factor(cyl)")

  withr::local_options(list(contrasts = c("contr.sum", "contr.poly")))
  fit_type_iii <- lm(mpg ~ factor(cyl) + wt + hp, data = mtcars)
  type_iii_url <- mellio_open(car::Anova(fit_type_iii, type = 3), browse = FALSE)
  type_iii_payload <- decode_stats_payload_url(type_iii_url)

  expect_match(type_iii_url, "#stats/payload=", fixed = TRUE)
  expect_no_match(type_iii_url, "#data=", fixed = TRUE)
  expect_equal(type_iii_payload$fields$ss_type, "type_iii")
  expect_equal(type_iii_payload$fields$model_kind, "ancova")
  expect_equal(type_iii_payload$fields$term, "factor(cyl)")
})

test_that("mellio_open on a numeric matrix routes to Tables", {
  m <- matrix(1:12, nrow = 3, dimnames = list(letters[1:3], LETTERS[1:4]))
  url <- mellio_open(m, browse = FALSE)

  expect_match(url, "#data=", fixed = TRUE)
  expect_no_match(url, "#stats/payload=", fixed = TRUE)
})

test_that("mellio_open on a contingency table (xtabs/table) routes to Tables", {
  url <- mellio_open(table(mtcars$cyl, mtcars$gear), browse = FALSE)

  expect_match(url, "#data=", fixed = TRUE)
  expect_no_match(url, "#stats/payload=", fixed = TRUE)
})

# ── Regression guards: unchanged Stats paths ────────────────────────────
# These pin the behavior that statistical inputs still route to Stats,
# guarding against accidental misrouting if the dispatch changes again.

test_that("mellio_open on an lm still routes to Stats", {
  url <- mellio_open(lm(mpg ~ wt, data = mtcars), browse = FALSE)
  expect_match(url, "#stats/payload=", fixed = TRUE)
})

test_that("mellio_open on a t.test still routes to Stats", {
  url <- mellio_open(t.test(extra ~ group, data = sleep), browse = FALSE)
  expect_match(url, "#stats/payload=", fixed = TRUE)
})

test_that("mellio_open on an aov still routes to Stats", {
  url <- mellio_open(aov(mpg ~ factor(cyl), data = mtcars), browse = FALSE)
  expect_match(url, "#stats/payload=", fixed = TRUE)
})

test_that("mellio_open sends prebuilt payloads through Stats", {
  payload <- mellio_payload(cor.test(mtcars$mpg, mtcars$wt))
  url <- mellio_open(payload, browse = FALSE)
  decoded <- decode_stats_payload_url(url)

  expect_equal(decoded$type, payload$type)
  expect_equal(decoded$fields$estimate$value, payload$fields$estimate$value)
})

test_that("mellio_open sends image files through Figures", {
  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  grDevices::png(tmp, width = 200, height = 150)
  plot(1:10, main = "test")
  grDevices::dev.off()

  url <- mellio_open(tmp, title = "Test figure", browse = FALSE)

  expect_match(url, "#mode=figure", fixed = TRUE)
  expect_match(url, "imageData=", fixed = TRUE)
  expect_match(url, "figTitle=Test%20figure", fixed = TRUE)
})

test_that("mellio_open sends ggplot objects through Figures", {
  testthat::skip_if_not_installed("ggplot2")

  p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
    ggplot2::geom_point()
  url <- mellio_open(p, title = "Weight and mpg", browse = FALSE)

  expect_match(url, "#mode=figure", fixed = TRUE)
  expect_match(url, "imageData=", fixed = TRUE)
  expect_match(url, "figTitle=Weight%20and%20mpg", fixed = TRUE)
})
