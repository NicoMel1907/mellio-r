#' Create a formatted statistical table
#'
#' The main entry point for melliotab. Accepts data frames, model objects,
#' or correlation matrices and formats them according to the specified
#' citation style. Passing two or more model objects creates a side-by-side
#' model comparison table.
#'
#' @param x A data.frame, model object (lm, glm, aov, htest), or correlation matrix
#' @param section Optional section selector for multi-section tables. Examples
#'   include `"fit"`, `"loadings"`, `"paths"`, `"covariances"`, `"defined"`,
#'   `"reliability"`, or `"modification_indices"` for structural model
#'   payloads; `"loadings"`, `"variance"`, or `"fit"` for EFA; and
#'   `"comparison"`, `"fit"`, or `"diff"` for FitDiff objects.
#' @param what Backward-compatible alias for `section` in methods that
#'   previously used `what`. Prefer `section`.
#' @param ... Additional arguments passed to methods
#' @return A melliotab object
#' @export
#' @section Supported inputs:
#' `melliotab()` is intended for table output inside R. It supports plain data
#' frames, matrices, base `table` objects, correlation matrices, common model
#' objects, hypothesis-test objects, and payloads created by `mellio_payload()`.
#' For objects that can produce several tables, call `melliotab(x)` once to
#' see the available `section` choices, then request one explicitly, for
#' example `melliotab(fit, section = "loadings")`.
#'
#' @section Common modifiers:
#' The `mt_*()` helpers are optional table modifiers. They follow the same
#' pattern as many R table packages: create a table once, then add formatting
#' only where needed. Common helpers include `mt_title()`, `mt_note()`,
#' `mt_decimals()`, `mt_format_ci()`, `mt_remove_leading_zeros()`,
#' `mt_sig_stars()`, `mt_spanner()`, `mt_indent()`, and
#' `mt_section_title()`.
#'
#' Significance stars are never added by default. Use `mt_sig_stars()` only
#' when that convention is appropriate for your manuscript, course, or journal.
#'
#' @examples
#' # From a data frame
#' df <- data.frame(
#'   Variable = c("Age", "Gender"),
#'   B = c(0.45, -1.23),
#'   SE = c(0.12, 0.34),
#'   t = c(3.75, -3.62),
#'   p = c(0.0003, 0.0004)
#' )
#' melliotab(df, style = "apa7", title = "Regression Results")
#'
#' # From a linear model
#' model <- lm(Ozone ~ Temp + Wind, data = airquality)
#' melliotab(model, style = "apa7", title = "Predictors of ozone concentration")
#'
#' # Compare multiple models side by side
#' m1 <- lm(Ozone ~ Temp, data = airquality)
#' m2 <- lm(Ozone ~ Temp + Wind, data = airquality)
#' m3 <- lm(Ozone ~ Temp + Wind + Solar.R, data = airquality)
#' melliotab(
#'   m1, m2, m3,
#'   labels = c("Step 1", "Step 2", "Step 3"),
#'   dep.var.labels = "Ozone concentration"
#' )
melliotab <- function(x, ...) {
  dots <- list(...)
  if (mellio_should_build_model_comparison_table(x, dots)) {
    return(mellio_model_comparison_table(x, dots))
  }

  UseMethod("melliotab")
}

mellio_should_build_model_comparison_table <- function(x, dots) {
  if (!mellio_is_model_comparison_object(x) || !length(dots)) return(FALSE)

  dot_names <- names(dots) %||% rep("", length(dots))
  unnamed <- !nzchar(dot_names)
  any(unnamed & vapply(dots, mellio_is_model_comparison_object, logical(1)))
}

mellio_model_comparison_table <- function(x, dots) {
  model_count <- 1L + sum(vapply(dots, mellio_is_model_comparison_object, logical(1)))
  if (!is.null(dots$labels)) {
    if (is.null(dots$column.labels)) dots$column.labels <- dots$labels
    dots$labels <- NULL
  }
  if (is.null(dots$title)) {
    dots$title <- paste0("Model comparison (", model_count, " models)")
  }
  do.call(mt_compare, c(list(x), dots))
}

mellio_is_model_comparison_object <- function(x) {
  if (inherits(x, c(
    "lm", "glm", "aov", "coxph", "survreg", "lme", "gls",
    "merMod", "lmerMod", "glmerMod", "polr", "multinom",
    "gam", "gee", "geeglm", "negbin", "zeroinfl"
  ))) {
    return(TRUE)
  }

  if (inherits(x, c(
    "formula", "data.frame", "matrix", "table", "htest", "character",
    "melliotab", "mellio_payload"
  ))) {
    return(FALSE)
  }

  f <- tryCatch(stats::formula(x), error = function(e) NULL)
  !is.null(f) && length(f) >= 2L
}

#' @rdname melliotab
#' @export
melliotab.default <- function(x, ..., section = NULL) {
  if (is.data.frame(x)) {
    return(melliotab.data.frame(x, ...))
  }

  cls <- class(x)
  payload_method <- NULL
  for (class_name in cls) {
    payload_method <- utils::getS3method("mellio_payload", class_name, optional = TRUE)
    if (!is.null(payload_method)) break
  }

  if (!is.null(payload_method)) {
    payload <- mellio_payload(x, ...)
    return(melliotab_from_payload(payload, section = section, ...))
  }

  fallback_payload <- tryCatch(
    mellio_payload(x, ...),
    error = function(e) NULL
  )
  if (inherits(fallback_payload, "mellio_payload") &&
      !fallback_payload$card_kind %in% c("unsupported", "raw_text")) {
    return(melliotab_from_payload(fallback_payload, section = section, ...))
  }

  cli::cli_abort(c(
    "melliotab does not know how to handle objects of class {.cls {class(x)}}.",
    "i" = "Tip: paste the printed output into Tables in Mellio."
  ))
}

#' @rdname melliotab
#' @export
melliotab.NULL <- function(x, ...) {
  cli::cli_abort(c(
    "{.fun melliotab} received {.code NULL}.",
    "i" = "This usually means the R code printed output but did not return a table or model object.",
    "i" = "Return a {.cls data.frame}, {.cls table}, {.fun summary}, or supported model object instead of only calling {.fun print} or {.fun cat}.",
    "i" = "To keep simple printed loop output, wrap the loop in {.fun capture.output} and pass those lines to {.fun mellio_open}."
  ))
}

#' @rdname melliotab
#' @export
melliotab.character <- function(x, style = "apa7", title = NULL,
                                number = NULL, note = NULL,
                                source = NULL, decimals = 2L,
                                p_decimals = 3L, ...) {
  df <- console_key_value_table(x)
  if (is.null(df)) {
    cli::cli_abort(c(
      "{.fun melliotab} could not parse character input as a table.",
      "i" = "Pass a {.cls data.frame} when possible.",
      "i" = "For printed loop output, use lines like {.code variable r = .12, p = .034}."
    ))
  }
  if (is.null(title)) title <- "Captured console output"
  if (is.null(source)) source <- "capture.output()"

  result <- melliotab.data.frame(
    df, style = style, title = title,
    number = number, note = note,
    source = source, decimals = decimals,
    p_decimals = p_decimals, ...
  )
  result$model <- x
  result
}

#' @rdname melliotab
#' @export
melliotab.summaryDefault <- function(x, style = "apa7", title = NULL,
                                     number = NULL, note = NULL,
                                     source = NULL, decimals = 2L,
                                     p_decimals = 3L, ...) {
  df <- summary_default_table(x)
  if (is.null(title)) title <- "Summary statistics"
  if (is.null(source)) source <- "base R summary()"

  result <- melliotab.data.frame(
    df, style = style, title = title,
    number = number, note = note,
    source = source, decimals = decimals,
    p_decimals = p_decimals, ...
  )
  result$model <- x
  result
}

#' @rdname melliotab
#' @export
melliotab.table <- function(x, style = "apa7", title = NULL,
                            number = NULL, note = NULL,
                            source = NULL, decimals = 2L,
                            p_decimals = 3L, ...) {
  df <- frequency_table_data(x)
  if (is.null(title)) title <- "Frequency table"
  if (is.null(source)) source <- "base R table()"

  result <- melliotab.data.frame(
    df, style = style, title = title,
    number = number, note = note,
    source = source, decimals = decimals,
    p_decimals = p_decimals, ...
  )
  result$model <- x
  result
}

#' @rdname melliotab
#' @param style Citation style: "apa7" or "ieee"
#' @param title Table title
#' @param number Table number (integer or character)
#' @param note Table note text
#' @param source Table source text
#' @param decimals Decimal places for estimates/statistics (1-4)
#' @param p_decimals Decimal places for p-values (2-4)
#' @export
melliotab.data.frame <- function(x, style = "apa7", title = NULL,
                                  number = NULL, note = NULL,
                                  source = NULL, decimals = 2L,
                                  p_decimals = 3L, ...) {
  style <- match.arg(style, list_styles())
  style_config <- get_style(style)

  # Store original data for re-styling

  raw_data <- x

  # Detect column types
  col_types <- detect_column_types(names(x))

  # Determine if leading zero removal applies based on style
  remove_lz <- isTRUE(style_config$remove_leading_zeros)
  lz_cols <- detect_leading_zero_cols(names(x))

  # Format the data
  formatted <- format_table_data(x, col_types, lz_cols,
                                  decimals = decimals,
                                  p_decimals = p_decimals,
                                  remove_lz = remove_lz,
                                  style_config = style_config)

  # Build the melliotab object
  structure(
    list(
      data = formatted,
      raw_data = raw_data,
      column_types = col_types,
      style = style,
      style_config = style_config,
      title = title,
      number = number,
      note = note,
      source = source,
      spanners = list(),
      section_titles = list(),
      merged_regions = list(),
      indent_levels = integer(0),
      decimals = as.integer(decimals),
      p_decimals = as.integer(p_decimals),
      options = list(
        remove_leading_zeros = remove_lz,
        italic_stat_headers = isTRUE(style_config$italic_stat_headers),
        bold_section_titles = isTRUE(style_config$bold_section_titles),
        sig_stars = FALSE,
        format_ci = FALSE,
        diagonal_mode = "all",
        triangle = "all",
        is_correlation = FALSE
      ),
      model = NULL,
      model_summary = NULL
    ),
    class = "melliotab"
  )
}

summary_default_table <- function(x) {
  values <- as.numeric(unname(x))
  labels <- names(x)
  if (is.null(labels) || length(labels) != length(values)) {
    labels <- paste0("Value ", seq_along(values))
  }

  clean <- vapply(labels, function(label) {
    if (ms_summary_missing_label(label)) return("Missing")
    switch(label,
      "Min." = "Min",
      "1st Qu." = "Q1",
      "Median" = "Median",
      "Mean" = "Mean",
      "3rd Qu." = "Q3",
      "Max." = "Max",
      label
    )
  }, character(1), USE.NAMES = FALSE)

  df <- as.data.frame(as.list(values), check.names = FALSE)
  names(df) <- clean
  if ("Missing" %in% names(df)) df[["Missing"]] <- as.integer(df[["Missing"]])
  df
}

frequency_table_data <- function(x) {
  df <- as.data.frame(x, responseName = "n", stringsAsFactors = FALSE)
  nms <- names(df)
  dim_names <- names(dimnames(x))
  value_cols <- setdiff(seq_along(nms), length(nms))

  for (i in value_cols) {
    dim_label <- if (length(dim_names) >= i && nzchar(dim_names[[i]])) {
      dim_names[[i]]
    } else if (length(value_cols) == 1L) {
      "Category"
    } else {
      paste0("Variable ", i)
    }
    names(df)[[i]] <- dim_label
  }

  df[["n"]] <- as.integer(df[["n"]])
  df
}

console_key_value_table <- function(x) {
  lines <- unlist(strsplit(paste(x, collapse = "\n"), "\n", fixed = TRUE),
                  use.names = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0L) return(NULL)

  parsed <- lapply(lines, console_key_value_row)
  if (any(vapply(parsed, is.null, logical(1)))) return(NULL)

  keys <- unique(unlist(lapply(parsed, function(row) names(row$values)),
                        use.names = FALSE))
  if (length(keys) == 0L) return(NULL)

  out <- data.frame(Variable = vapply(parsed, `[[`, character(1), "label"),
                    stringsAsFactors = FALSE, check.names = FALSE)
  for (key in keys) {
    out[[key]] <- vapply(parsed, function(row) {
      row$values[[key]] %||% NA_real_
    }, numeric(1))
  }
  out
}

console_key_value_row <- function(line) {
  number <- "[+-]?(?:\\d+\\.?\\d*|\\.\\d+)(?:[eE][+-]?\\d+)?"
  pattern <- paste0("([[:alpha:].][[:alnum:]_.]*)\\s*=\\s*(", number, ")")
  hits <- gregexpr(pattern, line, perl = TRUE)[[1]]
  if (identical(hits[[1]], -1L)) return(NULL)

  starts <- as.integer(hits)
  lens <- attr(hits, "match.length")
  first <- starts[[1]]
  label <- trimws(substr(line, 1L, first - 1L))
  label <- trimws(gsub("[,;:]+$", "", label))
  if (!nzchar(label)) return(NULL)

  values <- list()
  for (i in seq_along(starts)) {
    txt <- substr(line, starts[[i]], starts[[i]] + lens[[i]] - 1L)
    parts <- regmatches(txt, regexec(pattern, txt, perl = TRUE))[[1]]
    key <- parts[[2]]
    value <- suppressWarnings(as.numeric(parts[[3]]))
    if (is.na(value)) return(NULL)
    values[[key]] <- value
  }

  list(label = label, values = values)
}

#' Format table data according to APA/style rules
#'
#' @keywords internal
format_table_data <- function(data, col_types, lz_cols, decimals, p_decimals,
                               remove_lz, style_config) {
  formatted <- data

  for (j in seq_along(col_types)) {
    ct <- col_types[j]
    if (ct == "stub" || ct == "default") {
      if (is.numeric(formatted[[j]])) {
        formatted[[j]] <- formatC(as.numeric(formatted[[j]]),
                                   format = "f", digits = decimals)
      }
      next
    }

    # Convert column to character for formatting
    vals <- as.character(formatted[[j]])

    # Apply APA number formatting
    vals <- vapply(vals, function(v) {
      format_apa_number(v, ct, decimals = decimals, p_decimals = p_decimals,
                         remove_lz = remove_lz)
    }, character(1), USE.NAMES = FALSE)

    formatted[[j]] <- vals
  }

  # Apply leading zero removal to bounded stat columns
  if (remove_lz) {
    for (j in seq_along(lz_cols)) {
      if (lz_cols[j] && col_types[j] != "pvalue") {
        formatted[[j]] <- vapply(as.character(formatted[[j]]),
                                  remove_leading_zero,
                                  character(1), USE.NAMES = FALSE)
      }
    }
  }

  # Optional thousands/parenthetical number formatting
  if (isTRUE(style_config$thousands_separator)) {
    use_parens <- isTRUE(style_config$parenthetical_negatives)
    for (j in seq_along(col_types)) {
      if (col_types[j] == "stub") next
      formatted[[j]] <- vapply(as.character(formatted[[j]]),
                                function(v) format_business_number(v, use_parens),
                                character(1), USE.NAMES = FALSE)
    }
  }

  # Ensure all columns are character
  for (j in seq_along(formatted)) {
    formatted[[j]] <- as.character(formatted[[j]])
  }

  formatted
}
