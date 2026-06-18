# Number formatting functions
# Ported from apa-generator.js lines 560-662

#' Format a number according to APA rules
#'
#' @param value Character string containing the value
#' @param col_type Column type: "pvalue", "estimate", "statistic", "integer", "default"
#' @param decimals Number of decimal places for estimates/statistics
#' @param p_decimals Number of decimal places for p-values
#' @param remove_lz Whether to remove leading zeros (style-dependent)
#' @return Formatted character string
#' @keywords internal
format_apa_number <- function(value, col_type, decimals = 2L, p_decimals = 3L,
                               remove_lz = TRUE) {
  if (is.na(value)) return("")
  trimmed <- trimws(as.character(value))
  if (!nzchar(trimmed)) return(trimmed)

  # Already formatted "< .001"
  if (grepl("^<\\s", trimmed) || trimmed == "< .001") {
    ann <- extract_annotation(trimmed)
    return(paste0("< .001", ann$annotation))
  }

  ann <- extract_annotation(trimmed)
  core <- ann$core
  suffix <- ann$annotation

  # Remove thousand separators for parsing
  no_commas <- gsub(",(?=\\d{3})", "", core, perl = TRUE)
  num <- suppressWarnings(as.numeric(no_commas))
  if (is.na(num)) return(trimmed)

  # Pure integer in a default column — don't touch
  is_integer <- grepl("^-?\\d+$", no_commas)
  if (is_integer && col_type == "default") return(trimmed)
  if (col_type == "integer") return(paste0(as.character(round(num)), suffix))

  if (col_type == "pvalue") {
    if (num < 0.001 && num > 0) return(paste0("< .001", suffix))
    if (num == 0) return(paste0("< .001", suffix))
    formatted <- formatC(num, digits = p_decimals, format = "f")
    if (remove_lz) {
      formatted <- sub("^0\\.", ".", formatted)
      formatted <- sub("^-0\\.", "-.", formatted)
    }
    return(paste0(formatted, suffix))
  }

  if (col_type == "estimate" || col_type == "statistic") {
    formatted <- formatC(num, digits = decimals, format = "f")
    return(paste0(formatted, suffix))
  }

  trimmed
}

#' Remove leading zero from a bounded statistic value
#'
#' @param value Character string
#' @return Character string with leading zero removed if applicable
#' @keywords internal
remove_leading_zero <- function(value) {
  if (is.na(value)) return("")
  trimmed <- trimws(as.character(value))
  if (!nzchar(trimmed)) return(trimmed)

  ann <- extract_annotation(trimmed)
  core <- ann$core
  no_commas <- gsub(",(?=\\d{3})", "", core, perl = TRUE)
  num <- suppressWarnings(as.numeric(no_commas))
  if (is.na(num) || abs(num) >= 1) return(trimmed)

  result <- trimmed
  result <- sub("^0\\.", ".", result)
  result <- sub("^-0\\.", "-.", result)
  result
}

#' Get significance stars for a p-value
#'
#' @param p_value Character string or numeric p-value
#' @return Character string: "", "*", "**", or "***"
#' @keywords internal
get_sig_stars <- function(p_value) {
  if (is.character(p_value)) {
    trimmed <- trimws(p_value)
    has_lt <- grepl("^<", trimmed)
    ann <- extract_annotation(trimmed)
    num <- suppressWarnings(as.numeric(gsub("[<\\s]", "", ann$core)))
  } else {
    num <- p_value
    has_lt <- FALSE
  }

  if (is.na(num)) return("")

  if (has_lt) {
    if (num <= .001) return("***")
    if (num <= .01)  return("**")
    if (num <= .05)  return("*")
    return("")
  }

  if (num < .001) return("***")
  if (num < .01)  return("**")
  if (num < .05)  return("*")
  ""
}

#' Format a confidence interval value
#'
#' @param value Character string like "\[0.12, 0.45\]" or "0.12, 0.45"
#' @param decimals Number of decimal places
#' @return Formatted CI string "\[ll, ul\]"
#' @keywords internal
format_ci_value <- function(value, decimals = 2L) {
  trimmed <- trimws(value)
  m <- regmatches(trimmed, regexec("^\\[?\\s*(-?[\\d.]+)\\s*,\\s*(-?[\\d.]+)\\s*\\]?$", trimmed))[[1]]
  if (length(m) < 3) return(trimmed)

  ll <- as.numeric(m[2])
  ul <- as.numeric(m[3])
  if (is.na(ll) || is.na(ul)) return(trimmed)

  paste0("[", formatC(ll, digits = decimals, format = "f"), ", ",
         formatC(ul, digits = decimals, format = "f"), "]")
}

#' Format a number with thousands separators and parenthetical negatives
#'
#' @param value Character string containing the number
#' @param use_parenthetical Whether to use parentheses for negatives
#' @return Formatted character string
#' @keywords internal
format_business_number <- function(value, use_parenthetical = TRUE) {
  if (is.na(value)) return("")
  trimmed <- trimws(as.character(value))
  if (!nzchar(trimmed)) return(trimmed)

  ann <- extract_annotation(trimmed)
  core <- ann$core
  suffix <- ann$annotation

  no_commas <- gsub(",(?=\\d{3})", "", core, perl = TRUE)
  num <- suppressWarnings(as.numeric(no_commas))
  if (is.na(num)) return(trimmed)

  is_neg <- num < 0
  abs_val <- abs(num)

  # Detect original decimal places (minimum 2)
  dec_match <- regmatches(no_commas, regexpr("\\.\\d+$", no_commas))
  dec_places <- if (length(dec_match) > 0) max(nchar(dec_match) - 1, 2) else 2

  formatted <- formatC(abs_val, digits = dec_places, format = "f", big.mark = ",")

  if (is_neg && use_parenthetical) {
    formatted <- paste0("(", formatted, ")")
  } else if (is_neg) {
    formatted <- paste0("-", formatted)
  }

  paste0(formatted, suffix)
}
