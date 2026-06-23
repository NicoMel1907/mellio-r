close_graphics_devices <- function() {
  while (!is.null(grDevices::dev.list())) {
    grDevices::dev.off()
  }
}

expect_figure_url <- function(url, title = NULL) {
  expect_match(url, "#mode=figure", fixed = TRUE)
  expect_match(url, "imageData=", fixed = TRUE)

  if (!is.null(title)) {
    expect_match(
      url,
      paste0("figTitle=", utils::URLencode(title, reserved = TRUE)),
      fixed = TRUE
    )
  }
}

decode_figure_payload <- function(url) {
  value <- sub(".*[?&#]figurePayload=([^&]+).*", "\\1", url)
  expect_false(identical(value, url))
  value <- utils::URLdecode(value)
  std <- chartr("-_", "+/", value)
  pad <- (4 - nchar(std) %% 4) %% 4
  if (pad) std <- paste0(std, paste(rep("=", pad), collapse = ""))
  jsonlite::fromJSON(rawToChar(jsonlite::base64_dec(std)),
                     simplifyVector = FALSE)
}

skip_if_no_webshot_browser <- function() {
  testthat::skip_if_not_installed("webshot2")

  probe_html <- tempfile(fileext = ".html")
  probe_png <- tempfile(fileext = ".png")
  on.exit(unlink(c(probe_html, probe_png), force = TRUE), add = TRUE)
  writeLines("<!doctype html><html><body></body></html>", probe_html)

  browser_available <- tryCatch({
    suppressMessages(suppressWarnings(
      webshot2::webshot(probe_html, file = probe_png, vwidth = 10, vheight = 10)
    ))
    file.exists(probe_png)
  }, error = function(e) FALSE)

  testthat::skip_if(
    !isTRUE(browser_available),
    "Chrome/Chromium is not available"
  )
}

test_that("mellio_open sends lattice trellis objects through Figures", {
  testthat::skip_if_not_installed("lattice")
  withr::local_options(list(mellio.editor_url = "https://example.com"))

  p <- lattice::xyplot(mpg ~ wt, mtcars)
  url <- mellio_open(p, title = "Lattice plot", browse = FALSE)

  expect_figure_url(url, "Lattice plot")
})

test_that("mellio_open sends simple ggplots with editable raw scatter data", {
  testthat::skip_if_not_installed("ggplot2")
  withr::local_options(list(mellio.editor_url = "https://example.com"))

  p <- ggplot2::ggplot(
    mtcars,
    ggplot2::aes(wt, mpg, colour = factor(cyl))
  ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(
      x = "Weight (1000 lbs)",
      y = "Miles per gallon",
      colour = "Cylinders"
    ) +
    ggplot2::theme_minimal()

  url <- mellio_open(
    p,
    title = "Fuel efficiency by weight",
    number = 1,
    note = "Cylinders shown by color.",
    browse = FALSE
  )

  expect_figure_url(url, "Fuel efficiency by weight")
  expect_match(url, "figurePayload=", fixed = TRUE)
  expect_match(url, "figureType=raw_scatter", fixed = TRUE)
  expect_lt(regexpr("figurePayload=", url, fixed = TRUE)[[1]],
            regexpr("imageData=", url, fixed = TRUE)[[1]])

  payload <- decode_figure_payload(url)
  scatter <- payload$figure_data$raw_scatter
  expect_equal(payload$metadata$source, "r_ggplot")
  expect_equal(payload$metadata$available_figures[[1]]$type, "raw_scatter")
  expect_equal(scatter$x$label, "Weight (1000 lbs)")
  expect_equal(scatter$y$label, "Miles per gallon")
  expect_equal(scatter$group$label, "Cylinders")
  expect_equal(vapply(scatter$group$levels, `[[`, character(1), "value"),
               c("4", "6", "8"))
  expect_length(scatter$points, nrow(mtcars))
})

test_that("mellio_open sends ggplot count bars as editable bar charts", {
  testthat::skip_if_not_installed("ggplot2")
  withr::local_options(list(mellio.editor_url = "https://example.com"))

  p <- ggplot2::ggplot(mtcars, ggplot2::aes(factor(cyl))) +
    ggplot2::geom_bar() +
    ggplot2::labs(x = "Cylinders")

  url <- mellio_open(p, title = "Cars by cylinder count", browse = FALSE)

  expect_figure_url(url, "Cars by cylinder count")
  expect_match(url, "figureType=bar_chart", fixed = TRUE)
  payload <- decode_figure_payload(url)
  bars <- payload$figure_data$bar_chart
  expect_equal(payload$metadata$available_figures[[1]]$type, "bar_chart")
  expect_equal(bars$variable, "Cylinders")
  expect_equal(bars$y_label, "Count")
  expect_equal(vapply(bars$categories, `[[`, character(1), "label"),
               c("4", "6", "8"))
  expect_equal(vapply(bars$categories, `[[`, integer(1), "count"),
               as.integer(c(11, 7, 14)))
})

test_that("mellio_open preserves ggplot coord_flip bars as horizontal", {
  testthat::skip_if_not_installed("ggplot2")
  withr::local_options(list(mellio.editor_url = "https://example.com"))

  p <- ggplot2::ggplot(mtcars, ggplot2::aes(factor(cyl))) +
    ggplot2::geom_bar() +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Cylinders")

  url <- mellio_open(p, title = "Horizontal cylinder counts", browse = FALSE)

  expect_figure_url(url, "Horizontal cylinder counts")
  payload <- decode_figure_payload(url)
  expect_equal(payload$metadata$available_figures[[1]]$type, "bar_chart")
  expect_equal(payload$metadata$figure_settings$barOrientation, "horizontal")
})

test_that("mellio_open sends simple ggplot value bars as editable bar charts", {
  testthat::skip_if_not_installed("ggplot2")
  withr::local_options(list(mellio.editor_url = "https://example.com"))

  df <- data.frame(group = c("Control", "Treatment"), mean = c(3.2, 4.7))
  p <- ggplot2::ggplot(df, ggplot2::aes(group, mean)) +
    ggplot2::geom_col() +
    ggplot2::labs(x = "Group", y = "Mean score")

  url <- mellio_open(p, title = "Mean score by group", browse = FALSE)

  expect_figure_url(url, "Mean score by group")
  expect_match(url, "figureType=bar_chart", fixed = TRUE)
  payload <- decode_figure_payload(url)
  bars <- payload$figure_data$bar_chart
  expect_equal(bars$y_label, "Mean score")
  expect_equal(vapply(bars$categories, `[[`, character(1), "label"),
               c("Control", "Treatment"))
  expect_equal(vapply(bars$categories, `[[`, numeric(1), "count"),
               c(3.2, 4.7))
})

test_that("mellio_open sends ggplot lm smooths as editable regression scatterplots", {
  testthat::skip_if_not_installed("ggplot2")
  withr::local_options(list(mellio.editor_url = "https://example.com"))

  p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "lm") +
    ggplot2::labs(x = "Weight (1000 lbs)", y = "Miles per gallon")

  url <- mellio_open(p, title = "Fuel efficiency regression", browse = FALSE)

  expect_figure_url(url, "Fuel efficiency regression")
  expect_match(url, "figureType=scatter_plot", fixed = TRUE)
  payload <- decode_figure_payload(url)
  scatter <- payload$figure_data$scatter_plot
  expect_equal(payload$metadata$available_figures[[1]]$type, "scatter_plot")
  expect_equal(scatter$x_label, "Weight (1000 lbs)")
  expect_equal(scatter$y_label, "Miles per gallon")
  expect_length(scatter$observations, nrow(mtcars))
  expect_true(length(scatter$fit$line) >= 2)
  expect_true(all(vapply(scatter$fit$line, function(point) {
    !is.null(point$lower) && !is.null(point$upper)
  }, logical(1))))
})

test_that("mellio_open keeps complex ggplots as static figure fallbacks", {
  testthat::skip_if_not_installed("ggplot2")
  withr::local_options(list(mellio.editor_url = "https://example.com"))

  p <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg)) +
    ggplot2::geom_histogram(binwidth = 5)

  url <- mellio_open(p, title = "MPG histogram", browse = FALSE)

  expect_figure_url(url, "MPG histogram")
  expect_no_match(url, "figurePayload=", fixed = TRUE)
})

test_that("mellio_open sends htmlwidgets through Figures", {
  testthat::skip_if_not_installed("plotly")
  skip_if_no_webshot_browser()
  withr::local_options(list(mellio.editor_url = "https://example.com"))

  p <- plotly::plot_ly(mtcars, x = ~wt, y = ~mpg, type = "scatter", mode = "markers")
  url <- mellio_open(p, title = "Interactive plot", browse = FALSE)

  expect_figure_url(url, "Interactive plot")
})

test_that("mellio_open sends ggExtra plots through Figures", {
  testthat::skip_if_not_installed("ggplot2")
  testthat::skip_if_not_installed("ggExtra")
  withr::local_options(list(mellio.editor_url = "https://example.com"))

  p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
    ggplot2::geom_point()
  marginal <- ggExtra::ggMarginal(p)

  url <- mellio_open(marginal, title = "Marginal plot", browse = FALSE)

  expect_figure_url(url, "Marginal plot")
})

test_that("mellio_capture errors when no graphics device is active", {
  close_graphics_devices()
  on.exit(close_graphics_devices(), add = TRUE)

  expect_error(
    mellio_capture(browse = FALSE),
    "No active graphics device to capture"
  )
})

test_that("mellio_open points base R plotting calls to mellio_capture", {
  tmp <- tempfile(fileext = ".png")
  on.exit({
    close_graphics_devices()
    unlink(tmp)
  }, add = TRUE)

  grDevices::png(tmp, width = 600, height = 400)

  expect_error(
    mellio_open(plot(mtcars$wt, mtcars$mpg), browse = FALSE),
    "mellio_capture"
  )
})

test_that("mellio_capture sends the current base plot through Figures", {
  withr::local_options(list(mellio.editor_url = "https://example.com"))
  tmp <- tempfile(fileext = ".png")
  on.exit({
    close_graphics_devices()
    unlink(tmp)
  }, add = TRUE)

  grDevices::png(tmp, width = 600, height = 400)
  grDevices::dev.control("enable")
  plot(mtcars$wt, mtcars$mpg)

  url <- mellio_capture(title = "Base plot", browse = FALSE)

  expect_figure_url(url, "Base plot")
})
