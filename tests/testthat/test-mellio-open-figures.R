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
