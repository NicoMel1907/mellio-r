#' Apply significance stars to a melliotab object
#'
#' Adds significance stars to the target column based on p-values,
#' and optionally removes the p-value column.
#'
#' @param mt A melliotab object
#' @param target Column index or name to add stars to ("auto" to detect)
#' @param remove_p Whether to remove the p-value column
#' @param levels Named character vector of significance levels
#' @return Modified melliotab object
#' @keywords internal
apply_sig_stars <- function(mt, target = "auto", remove_p = TRUE,
                             levels = c("*" = 0.05, "**" = 0.01, "***" = 0.001)) {
  headers <- names(mt$data)

  p_col <- detect_p_value_col(headers)
  if (is.null(p_col)) {
    cli::cli_warn("No p-value column detected. Significance stars not applied.")
    return(mt)
  }

  if (identical(target, "auto")) {
    target_col <- detect_estimate_col(headers, p_col)
  } else if (is.character(target)) {
    target_col <- match(target, headers)
    if (is.na(target_col)) {
      cli::cli_warn("Target column {.val {target}} not found.")
      return(mt)
    }
  } else {
    target_col <- as.integer(target)
  }

  p_vals <- mt$raw_data[[p_col]]
  target_vals <- mt$data[[target_col]]

  for (i in seq_along(p_vals)) {
    stars <- get_sig_stars(p_vals[i])
    if (nzchar(stars)) {
      target_vals[i] <- paste0(target_vals[i], stars)
    }
  }
  mt$data[[target_col]] <- target_vals

  if (remove_p) {
    mt$data <- mt$data[, -p_col, drop = FALSE]
    mt$raw_data <- mt$raw_data[, -p_col, drop = FALSE]
    mt$column_types <- mt$column_types[-p_col]
  }

  mt$options$sig_stars <- TRUE
  mt
}
