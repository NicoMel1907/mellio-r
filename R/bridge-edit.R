ms_warn_missing_script <- function(payload) {
  if (!identical(ms_provenance_mode(), "full")) return(invisible())
  if (is.null(payload$provenance)) return(invisible())

  script <- payload$provenance$script

  if (!is.null(script) && is.na(script$line %||% NA)) {
    cli::cli_alert_info(
      "Script captured ({.file {basename(script$file)}}), but no line number.",
      " " = "For exact {.code Open in R} later, run {.fn mellio_open} from a \\
             saved RStudio script line; sourcing the file also works."
    )
    return(invisible())
  }

  if (!is.null(script)) return(invisible())

  in_rstudio <- requireNamespace("rstudioapi", quietly = TRUE) &&
                isTRUE(tryCatch(rstudioapi::isAvailable(),
                                error = function(e) FALSE))

  if (!in_rstudio) {
    cli::cli_alert_info(c(
      "Script source not captured (running outside RStudio).",
      " " = "Run {.fn mellio_open} from a saved RStudio script to capture \\
             file/line provenance for {.code Open in R} later."
    ))
    return(invisible())
  }

  cli::cli_alert_info(c(
    "Script source not captured (no saved file in the active editor).",
    " " = "Save the script, then run {.fn mellio_open} from the editor to \\
           capture the file/line for {.code Open in R} later."
  ))
  invisible()
}

send_payload_to_stats <- function(payload, browse = TRUE) {
  url <- ms_stats_url(payload)

  if (isTRUE(browse) && interactive()) {
    utils::browseURL(url)
  }

  if (interactive()) {
    cli::cli_alert_success("{payload$type_label}")
    cli::cli_inform(if (isTRUE(browse)) "  Opened in Mellio." else "  Mellio URL prepared.")
    ms_warn_missing_script(payload)
  }

  invisible(url)
}

ms_stats_url <- function(payload) {
  rlang::check_installed("jsonlite", reason = "to encode R bridge payloads")

  base <- ms_mellio_base_url()
  json <- jsonlite::toJSON(
    payload,
    auto_unbox = TRUE,
    na = "null",
    null = "null",
    digits = NA
  )
  encoded <- ms_base64url_encode(as.character(json))

  url <- paste0(base, "/#stats/payload=", encoded,
                mellio_launch_hash_param())

  if (nchar(url) > ms_url_size_limit) {
    cli::cli_warn(c(
      "Payload too large for URL transfer ({nchar(url)} chars).",
      "i" = "Use {.fn ms_send} (signed-in cloud save) when it lands in v0.2,",
      "i" = "or save the payload offline with {.fn mellio_to_json}."
    ))
  }

  url
}

ms_mellio_base_url <- function() {
  opt_url <- getOption("mellio.editor_url", default = NULL)
  if (!is.null(opt_url) && nzchar(opt_url)) {
    return(ms_validate_mellio_base_url(opt_url, "options('mellio.editor_url')"))
  }

  "https://www.mellioapp.com"
}

ms_validate_mellio_base_url <- function(url, source = "Mellio URL") {
  if (length(url) != 1L || is.na(url)) {
    cli::cli_abort("{source} must be a single URL.")
  }

  url <- trimws(as.character(url))
  if (!nzchar(url)) {
    cli::cli_abort("{source} must not be empty.")
  }
  if (grepl("[[:space:]]", url) || grepl("[[:cntrl:]]", url)) {
    cli::cli_abort("{source} must not contain whitespace or control characters.")
  }
  if (!grepl("^https?://", url, ignore.case = TRUE)) {
    cli::cli_abort("{source} must start with {.val https://} or {.val http://}.")
  }
  if (grepl("[?#]", url)) {
    cli::cli_abort("{source} must not include a query string or URL fragment.")
  }

  authority <- sub("^https?://", "", url, ignore.case = TRUE)
  authority <- sub("/.*$", "", authority)
  if (!nzchar(authority) || grepl("@", authority, fixed = TRUE)) {
    cli::cli_abort("{source} must include a plain host name.")
  }

  sub("/+$", "", url)
}

ms_url_size_limit <- 1900000L

ms_base64url_encode <- function(x) {
  raw_bytes <- charToRaw(enc2utf8(x))
  b64 <- ms_base64_encode_raw(raw_bytes)
  b64 <- chartr("+/", "-_", b64)
  sub("=+$", "", b64)
}

ms_base64_encode_raw <- function(bytes) {
  alpha <- c(LETTERS, letters, as.character(0:9), "+", "/")
  n <- length(bytes)
  if (n == 0) return("")

  pad <- (3 - n %% 3) %% 3
  if (pad > 0) bytes <- c(bytes, as.raw(rep(0, pad)))

  ints <- as.integer(bytes)
  ng <- length(ints) / 3
  out <- character(ng * 4)
  for (i in seq_len(ng)) {
    j <- (i - 1L) * 3L
    b1 <- ints[j + 1L]; b2 <- ints[j + 2L]; b3 <- ints[j + 3L]
    out[(i - 1L) * 4L + 1L] <- alpha[bitwShiftR(b1, 2) + 1L]
    out[(i - 1L) * 4L + 2L] <- alpha[bitwOr(bitwShiftL(bitwAnd(b1, 0x03), 4), bitwShiftR(b2, 4)) + 1L]
    out[(i - 1L) * 4L + 3L] <- alpha[bitwOr(bitwShiftL(bitwAnd(b2, 0x0F), 2), bitwShiftR(b3, 6)) + 1L]
    out[(i - 1L) * 4L + 4L] <- alpha[bitwAnd(b3, 0x3F) + 1L]
  }

  s <- paste(out, collapse = "")
  if (pad > 0) {
    s <- paste0(substr(s, 1, nchar(s) - pad), strrep("=", pad))
  }
  s
}
