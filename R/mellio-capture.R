#' Capture the current graphics device and open it in Mellio
#'
#' Convenience wrapper for base R plots. Call after producing a base plot
#' with `plot()`, `hist()`, `barplot()`, or another graphics function.
#'
#' @param ... Arguments passed to [mellio_open()], such as `title`,
#'   `number`, `note`, and `browse`.
#' @return Invisibly, the Mellio URL string.
#' @export
#' @family R bridge
#' @examples
#' \dontrun{
#' plot(mtcars$wt, mtcars$mpg)
#' mellio_capture(title = "Fuel efficiency by weight", number = 1)
#' }
mellio_capture <- function(...) {
  if (is.null(grDevices::dev.list())) {
    cli::cli_abort(c(
      "No active graphics device to capture.",
      "i" = "Run a plotting function such as {.code plot()} or {.code hist()} first."
    ))
  }

  mellio_open(grDevices::recordPlot(), ...)
}
