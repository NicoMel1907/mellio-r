# Correlation matrix detection and formatting
# Ported from apa-generator.js lines 1407-1469

#' Check if a matrix looks like a correlation matrix
#'
#' @param x A matrix
#' @return Logical
#' @keywords internal
is_correlation_matrix <- function(x) {
  if (!is.matrix(x)) return(FALSE)
  if (nrow(x) != ncol(x)) return(FALSE)

  # Check diagonal is all 1s (or very close)
  diag_vals <- diag(x)
  if (all(abs(diag_vals - 1) < 0.001, na.rm = TRUE)) return(TRUE)

  FALSE
}

#' Apply correlation-specific formatting to a melliotab object
#'
#' @param mt A melliotab object
#' @return Modified melliotab object
#' @keywords internal
apply_correlation_formatting <- function(mt) {
  data <- mt$data
  diagonal <- mt$options$diagonal_mode
  triangle <- mt$options$triangle

  # Find where numeric columns start (skip stub column)
  start_col <- 1L
  for (j in seq_along(mt$column_types)) {
    if (mt$column_types[j] != "stub") {
      start_col <- j
      break
    }
  }

  n_vars <- ncol(data) - start_col + 1
  n_rows <- nrow(data)

  # Format diagonal
  if (diagonal == "dash") {
    for (i in seq_len(min(n_rows, n_vars))) {
      col_idx <- start_col + i - 1
      if (col_idx <= ncol(data)) {
        data[i, col_idx] <- "\u2014"  # em-dash
      }
    }
  } else if (diagonal == "blank") {
    for (i in seq_len(min(n_rows, n_vars))) {
      col_idx <- start_col + i - 1
      if (col_idx <= ncol(data)) {
        data[i, col_idx] <- ""
      }
    }
  }
  # "one" mode: keep the formatted value (already formatted as 1.00 etc.)

  # Blank above or below diagonal
  if (triangle == "lower") {
    # Blank above diagonal
    for (i in seq_len(n_rows)) {
      for (j in seq_len(n_vars)) {
        col_idx <- start_col + j - 1
        if (j > i && col_idx <= ncol(data)) {
          data[i, col_idx] <- ""
        }
      }
    }
  } else if (triangle == "upper") {
    # Blank below diagonal
    for (i in seq_len(n_rows)) {
      for (j in seq_len(n_vars)) {
        col_idx <- start_col + j - 1
        if (j < i && col_idx <= ncol(data)) {
          data[i, col_idx] <- ""
        }
      }
    }
  }

  mt$data <- data
  mt
}
