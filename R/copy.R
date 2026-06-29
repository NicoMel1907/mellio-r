#' Copy a melliotab table to the system clipboard
#'
#' Copies the formatted table to the clipboard as rich HTML so it can be
#' pasted into Word with formatting preserved.
#'
#' @param x A melliotab object
#' @return Invisible melliotab object (for piping)
#' @export
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   model <- lm(mpg ~ wt + hp, data = mtcars)
#'   tab <- melliotab(model, title = "Regression Results")
#'   mt_copy(tab)
#'   # Now paste into Word
#' }
#' }
mt_copy <- function(x) {
  if (!inherits(x, "melliotab")) {
    cli::cli_abort("{.arg x} must be a melliotab object.")
  }

  os <- Sys.info()[["sysname"]]

  full_html <- mt_as_html(x)

  if (os == "Darwin") {
    tmp <- tempfile(fileext = ".html")
    writeLines(full_html, tmp)
    jxa <- paste0(
      'ObjC.import("AppKit");',
      'var pb=$.NSPasteboard.generalPasteboard;',
      'pb.clearContents;',
      'var html=$.NSString.alloc',
      '.initWithContentsOfFileEncodingError("', tmp, '",',
      '$.NSUTF8StringEncoding,null);',
      'pb.setStringForType(html,$.NSPasteboardTypeHTML);'
    )
    result <- system2("osascript",
      c("-l", "JavaScript", "-e", shQuote(jxa)),
      stdout = TRUE, stderr = TRUE
    )
    if (!is.null(attr(result, "status")) &&
        attr(result, "status") != 0) {
      cli::cli_abort("Failed to copy table to clipboard.")
    }
    cli::cli_inform("Table copied to clipboard.")
  } else if (os == "Windows") {
    tmp <- tempfile(fileext = ".html")
    writeLines(full_html, tmp)
    ps_cmd <- paste0(
      'Add-Type -AssemblyName System.Windows.Forms; ',
      '$html = [System.IO.File]::ReadAllText("',
      gsub("/", "\\\\", tmp), '"); ',
      '$data = New-Object System.Windows.Forms.DataObject; ',
      '$data.SetData(',
      '[System.Windows.Forms.DataFormats]::Html, $html); ',
      '[System.Windows.Forms.Clipboard]::',
      'SetDataObject($data, $true)'
    )
    result <- system2("powershell", c("-command", ps_cmd),
                      stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(result, "status")) &&
        attr(result, "status") != 0) {
      cli::cli_abort("Failed to copy table to clipboard.")
    }
    cli::cli_inform("Table copied to clipboard.")
  } else {
    cli::cli_inform(c(
      "!" = "Clipboard copy is not supported on {os}.",
      "i" = "Use {.fun mt_save} to save the table to a file instead."
    ))
  }

  invisible(x)
}
