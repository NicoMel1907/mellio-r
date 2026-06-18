# Save to file

#' Save a melliotab table to file
#'
#' Auto-detects the format from the file extension.
#'
#' @param x A melliotab object
#' @param filename Output file path. Supported extensions: .html, .tex, .md
#' @param ... Additional arguments passed to format-specific functions
#' @return Invisible file path
#' @export
#'
#' @examples
#' \dontrun{
#' model <- lm(mpg ~ wt + hp, data = mtcars)
#' tbl <- melliotab(model, style = "apa7", title = "Results")
#'
#' mt_save(tbl, "table.html")
#' mt_save(tbl, "table.tex")
#' mt_save(tbl, "table.md")
#' }
mt_save <- function(x, filename, ...) {
  if (!inherits(x, "melliotab")) {
    cli::cli_abort("{.arg x} must be a melliotab object.")
  }

  ext <- tolower(tools::file_ext(filename))

  switch(ext,
    "html" = {
      html_str <- mt_as_html(x)
      writeLines(html_str, filename)
    },
    "tex" = {
      latex_str <- mt_as_latex(x)
      writeLines(latex_str, filename)
    },
    "md" = {
      md_str <- mt_as_markdown(x)
      writeLines(md_str, filename)
    },
    cli::cli_abort("Unsupported file extension {.val {ext}}. Use .html, .tex, or .md.")
  )

  cli::cli_inform("Saved to {.file {filename}}")
  invisible(filename)
}
