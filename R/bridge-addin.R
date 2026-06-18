# R bridge — RStudio addin.
#
# Adds a "Send to Mellio" entry under RStudio's Addins menu. When the
# user highlights an R expression (e.g. `t.test(score ~ group, data = df)`)
# and runs the addin, Mellio receives the evaluated result as a Result
# Card — the original selection text becomes the displayed call so
# users don't lose the expression even though we evaluated it.
#
# Addin declaration: inst/rstudio/addins.dcf
# Requires rstudioapi (Suggests) — fails with a clear message if absent.

#' RStudio addin — Send the selected R expression to Mellio
#'
#' Reads the active selection in the RStudio editor, evaluates it in
#' the global environment, and opens the result as a Result Card in
#' the Mellio web app. The selected text is preserved as the displayed
#' call so the saved card still says e.g. "t.test(score ~ group, data
#' = df)" rather than the evaluated variable.
#'
#' Registered via `inst/rstudio/addins.dcf`. Appears in RStudio's
#' **Addins** menu once the package is loaded.
#'
#' @return Invisibly `NULL`. Side effect: opens the browser at Mellio.
#' @export
#' @family R bridge
#' @examples
#' \dontrun{
#' # Highlight   t.test(extra ~ group, data = sleep)   in RStudio
#' # → Addins menu → Send to Mellio
#' mellio_addin_send()
#' }
mellio_addin_send <- function() {
  rlang::check_installed("rstudioapi",
    reason = "for the Mellio RStudio addin")

  if (!rstudioapi::isAvailable()) {
    cli::cli_abort(c(
      "The Mellio addin requires RStudio.",
      "i" = "Use {.fn mellio_open} from the console outside RStudio."
    ))
  }

  ctx <- rstudioapi::getActiveDocumentContext()
  selected <- NA_character_
  if (!is.null(ctx$selection) && length(ctx$selection) > 0L) {
    selected <- ctx$selection[[1]]$text
  }

  if (is.na(selected) || !nzchar(selected)) {
    rstudioapi::showDialog(
      "Mellio",
      paste0(
        "Highlight an R expression first ",
        "(e.g. <code>t.test(score ~ group, data = df)</code>), ",
        "then run the addin again."
      )
    )
    return(invisible(NULL))
  }

  # Parse and evaluate the selection. The eval uses the global env so
  # references to user-defined data frames / variables resolve.
  expr <- tryCatch(parse(text = selected),
                   error = function(e) NULL)
  if (is.null(expr)) {
    rstudioapi::showDialog(
      "Mellio",
      paste0("The selected text isn't valid R:\n", selected)
    )
    return(invisible(NULL))
  }

  result <- tryCatch(eval(expr, envir = globalenv()),
                     error = function(e) {
                       rstudioapi::showDialog(
                         "Mellio",
                         paste0("Evaluating the selection failed:\n",
                                conditionMessage(e))
                       )
                       NULL
                     })
  if (is.null(result)) return(invisible(NULL))

  call_str <- trimws(gsub("\\s+", " ", selected))

  url <- tryCatch(
    mellio_open(result, .call = call_str),
    error = function(e) {
      rstudioapi::showDialog(
        "Mellio",
        paste0(
          "Could not send to Mellio:\n",
          conditionMessage(e), "\n\n",
          "Tip: mellio_payload() supports htest, lm, anova, lmer, glmer, ",
          "lavaan, mediation::mediate, psych::alpha(), psych::corr.test(), ",
          "numeric vectors, summary() results, and result-shaped data frames."
        )
      )
      NULL
    }
  )
  invisible(url)
}
