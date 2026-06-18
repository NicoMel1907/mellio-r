#' Convert a data frame to tab-delimited text
#'
#' @param data A data frame
#' @return A single character string with tab-separated values
#' @keywords internal
table_to_tsv <- function(data) {
  header <- paste(names(data), collapse = "\t")
  rows <- apply(data, 1, function(row) {
    paste(as.character(row), collapse = "\t")
  })
  paste(c(header, rows), collapse = "\n")
}

mellio_launch_hash_param <- function() {
  nonce <- paste0(format(Sys.time(), "%Y%m%d%H%M%OS6", tz = "UTC"),
                  "-", Sys.getpid())
  paste0("&mellio_launch=", utils::URLencode(nonce, reserved = TRUE))
}

#' Build extended editor data with metadata
#'
#' Constructs a TSV grid that includes spanner rows above the header,
#' plus rowRoles and mergedRegions arrays for PaperEditor.
#'
#' @param x A melliotab object
#' @return List with tsv (string), rowRoles (character vector),
#'   mergedRegions (list of lists with row/col/rowSpan/colSpan),
#'   cellFormats, and cellBorders.
#' @keywords internal
build_editor_data <- function(x) {
  tbl_data <- x$data
  n_cols <- ncol(tbl_data)
  opts <- x$options %||% list()

  # Build grid as list of character vectors: header + data rows
  header_row <- as.character(names(tbl_data))
  data_rows <- lapply(seq_len(nrow(tbl_data)), function(i) {
    as.character(tbl_data[i, ])
  })

  all_rows <- c(list(header_row), data_rows)
  roles <- c("header", rep("data", nrow(tbl_data)))
  merged <- list()
  cell_formats <- list()
  cell_borders <- list()

  # Insert spanner rows above header (process in reverse so indices stay correct)
  if (length(x$spanners) > 0) {
    for (sp in rev(x$spanners)) {
      span_row <- rep("", n_cols)
      start_col <- min(sp$columns)
      span_row[start_col] <- sp$label

      # Prepend spanner row
      all_rows <- c(list(span_row), all_rows)
      roles <- c("spanheader", roles)

      # Shift existing merged regions down by 1
      merged <- lapply(merged, function(m) {
        m$row <- m$row + 1L
        m
      })

      # Add merged region for this spanner (0-indexed)
      merged <- c(merged, list(list(
        row = 0L,
        col = as.integer(start_col - 1),
        rowSpan = 1L,
        colSpan = as.integer(length(sp$columns))
      )))
    }
  }

  # Apply indent levels: prefix stub cell with em-spaces
  if (length(x$indent_levels) > 0) {
    n_span_rows <- sum(roles == "spanheader")
    header_offset <- n_span_rows + 1L  # spanner rows + header row

    for (row_str in names(x$indent_levels)) {
      row_idx <- as.integer(row_str)
      level <- x$indent_levels[row_str]
      actual_idx <- row_idx + header_offset  # 1-based data row → grid index
      if (actual_idx >= 1 && actual_idx <= length(all_rows)) {
        all_rows[[actual_idx]][1] <- paste0(
          paste(rep("\u2003", level), collapse = ""),
          all_rows[[actual_idx]][1]
        )
      }
    }
  }

  set_cell_format <- function(row, col, values) {
    key <- paste(row, col, sep = ",")
    cell_formats[[key]] <<- utils::modifyList(cell_formats[[key]] %||% list(),
                                              values)
  }

  set_cell_border <- function(row, col, side) {
    key <- paste(row, col, sep = ",")
    cell_borders[[key]] <<- cell_borders[[key]] %||% list()
    cell_borders[[key]][[side]] <<- TRUE
  }

  if (isTRUE(opts$is_comparison)) {
    dep_rows <- which(vapply(all_rows, function(row) {
      any(row == "Dependent variable:")
    }, logical(1)))

    for (row_idx in dep_rows) {
      col_idx <- which(all_rows[[row_idx]] == "Dependent variable:")[[1]]
      set_cell_format(row_idx - 1L, col_idx - 1L,
                      list(italic = TRUE, bold = FALSE, align = "center"))
      set_cell_border(row_idx - 1L, col_idx - 1L, "bottom")
    }

    stat_start <- opts$comparison_stat_start
    if (!is.null(stat_start) && length(stat_start) == 1L &&
        !is.na(stat_start) && stat_start > 0L) {
      n_span_rows <- sum(roles == "spanheader")
      stat_row <- as.integer(n_span_rows + stat_start)
      if (stat_row >= 0L && stat_row < length(all_rows)) {
        for (col in seq_len(n_cols) - 1L) {
          set_cell_border(stat_row, col, "top")
        }
      }
    }
  }

  # Convert grid to TSV
  tsv_lines <- vapply(all_rows, function(row) {
    paste(row, collapse = "\t")
  }, character(1))
  tsv <- paste(tsv_lines, collapse = "\n")

  list(
    tsv = tsv,
    rowRoles = roles,
    mergedRegions = merged,
    cellFormats = cell_formats,
    cellBorders = cell_borders
  )
}

send_figure_to_mellio <- function(x, browse = TRUE) {
  if (is.null(x$provenance)) {
    x$provenance <- ms_provenance_basic()
  }
  url <- mellio_figure_edit_url(x)

  if (isTRUE(browse) && interactive()) {
    cli::cli_inform(c(
      "i" = "Opening figure in Mellio..."
    ))

    utils::browseURL(url)
  }

  invisible(url)
}

mellio_figure_edit_url <- function(x) {
  # Build URL with figure mode
  base <- ms_mellio_base_url()
  url <- paste0(sub("/+$", "", base), "/#mode=figure",
                mellio_launch_hash_param())

  # Try to include compressed image data
  img_b64 <- compress_image_for_url(x)
  if (!is.null(img_b64)) {
    mime <- attr(img_b64, "mime") %||% "image/jpeg"
    img_data_uri <- paste0("data:", mime, ";base64,", img_b64)
    encoded_img <- utils::URLencode(img_data_uri, reserved = TRUE)
    # Check URL length — browsers handle ~2MB URLs
    if (nchar(encoded_img) < 2000000) {
      url <- paste0(url, "&imageData=", encoded_img)
    } else {
      cli::cli_inform(c(
        "i" = "Image too large for URL transfer. Please upload it in Mellio."
      ))
    }
  }

  # Add metadata
  if (!is.null(x$title)) {
    url <- paste0(url, "&figTitle=",
                  utils::URLencode(x$title, reserved = TRUE))
  }
  if (!is.null(x$number)) {
    url <- paste0(url, "&figNumber=",
                  utils::URLencode(as.character(x$number), reserved = TRUE))
  }
  if (!is.null(x$note)) {
    url <- paste0(url, "&figNote=",
                  utils::URLencode(x$note, reserved = TRUE))
  }

  # Pass citation style
  if (!is.null(x$style)) {
    style_map <- c(apa7 = "APA7", ieee = "IEEE")
    web_style <- style_map[x$style]
    if (!is.na(web_style)) {
      url <- paste0(url, "&style=", web_style)
    }
  }

  # Pass R source provenance as base64-encoded JSON
  if (!is.null(x$provenance)) {
    prov_json <- jsonlite::toJSON(x$provenance, auto_unbox = TRUE,
                                  null = "null", na = "null")
    url <- paste0(url, "&provenance=",
      utils::URLencode(
        mellio_base64_compact(charToRaw(as.character(prov_json))),
        reserved = TRUE
      ))
  }

  url
}

send_table_to_mellio <- function(x, browse = TRUE) {
  if (inherits(x, "melliotab") && is.null(x$provenance)) {
    x$provenance <- ms_provenance_basic()
  }
  url <- mellio_table_url(x)

  if (isTRUE(browse) && interactive()) {
    cli::cli_inform(c(
      "i" = "Opening table in Mellio..."
    ))

    utils::browseURL(url)
  }

  invisible(url)
}

mellio_table_url <- function(x) {
  base <- ms_mellio_base_url()
  opts <- x$options %||% list()
  has_comparison_rule <- !is.null(opts$comparison_stat_start) &&
    length(opts$comparison_stat_start) == 1L &&
    !is.na(opts$comparison_stat_start)

  has_structure <- length(x$spanners %||% list()) > 0L ||
                   length(x$indent_levels %||% integer(0)) > 0L ||
                   length(x$merged_regions %||% list()) > 0L ||
                   has_comparison_rule

  if (has_structure) {
    ed <- build_editor_data(x)
    tsv <- ed$tsv
    row_roles_json <- jsonlite::toJSON(ed$rowRoles, auto_unbox = FALSE)
    merged_regions_json <- if (length(ed$mergedRegions) > 0L) {
      jsonlite::toJSON(ed$mergedRegions, auto_unbox = TRUE)
    } else {
      NULL
    }
    cell_formats_json <- if (length(ed$cellFormats) > 0L) {
      jsonlite::toJSON(ed$cellFormats, auto_unbox = TRUE)
    } else {
      NULL
    }
    cell_borders_json <- if (length(ed$cellBorders) > 0L) {
      jsonlite::toJSON(ed$cellBorders, auto_unbox = TRUE)
    } else {
      NULL
    }
  } else {
    tsv <- table_to_tsv(x$data)
    row_roles_json <- NULL
    merged_regions_json <- NULL
    cell_formats_json <- NULL
    cell_borders_json <- NULL
  }

  tsv_b64 <- ms_base64_encode_raw(charToRaw(enc2utf8(tsv)))
  base <- sub("/+$", "", base)
  url <- paste0(base, "/#data=", utils::URLencode(tsv_b64, reserved = TRUE),
                mellio_launch_hash_param())

  if (!is.null(row_roles_json)) {
    rr_b64 <- ms_base64_encode_raw(
      charToRaw(enc2utf8(as.character(row_roles_json)))
    )
    url <- paste0(url, "&rowRoles=",
                  utils::URLencode(rr_b64, reserved = TRUE))
  }

  if (!is.null(merged_regions_json)) {
    merged_b64 <- ms_base64_encode_raw(
      charToRaw(enc2utf8(as.character(merged_regions_json)))
    )
    url <- paste0(url, "&mergedRegions=",
                  utils::URLencode(merged_b64, reserved = TRUE))
  }

  if (!is.null(cell_formats_json)) {
    formats_b64 <- ms_base64_encode_raw(
      charToRaw(enc2utf8(as.character(cell_formats_json)))
    )
    url <- paste0(url, "&cellFormats=",
                  utils::URLencode(formats_b64, reserved = TRUE))
  }

  if (!is.null(cell_borders_json)) {
    borders_b64 <- ms_base64_encode_raw(
      charToRaw(enc2utf8(as.character(cell_borders_json)))
    )
    url <- paste0(url, "&cellBorders=",
                  utils::URLencode(borders_b64, reserved = TRUE))
  }

  if (!is.null(x$title)) {
    url <- paste0(url, "&title=",
                  utils::URLencode(x$title, reserved = TRUE))
  }
  if (!is.null(x$note)) {
    url <- paste0(url, "&note=",
                  utils::URLencode(x$note, reserved = TRUE))
  }
  if (!is.null(x$number)) {
    url <- paste0(url, "&tableNumber=",
                  utils::URLencode(as.character(x$number), reserved = TRUE))
  }
  if (!is.null(x$style)) {
    style_map <- c(apa7 = "APA7", ieee = "IEEE")
    web_style <- style_map[x$style]
    if (!is.na(web_style)) {
      url <- paste0(url, "&style=", web_style)
    }
  }
  if (!is.null(x$source)) {
    url <- paste0(url, "&source=",
                  utils::URLencode(x$source, reserved = TRUE))
  }
  if (!is.null(x$provenance)) {
    prov_json <- jsonlite::toJSON(x$provenance, auto_unbox = TRUE,
                                  null = "null", na = "null")
    prov_b64 <- ms_base64_encode_raw(
      charToRaw(enc2utf8(as.character(prov_json)))
    )
    url <- paste0(url, "&provenance=",
                  utils::URLencode(prov_b64, reserved = TRUE))
  }

  url
}

#' Compress image for URL transfer
#'
#' @param fig Internal Mellio plot handoff object.
#' @return Base64-encoded string or NULL if compression fails
#' @keywords internal
compress_image_for_url <- function(fig) {
  if (!is.null(fig$plot) && inherits(fig$plot, "gg")) {
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      w <- fig$width_in %||% 6
      h <- fig$height_in %||% 4

      tmp_png <- tempfile(fileext = ".png")
      on.exit(unlink(tmp_png), add = TRUE)
      tryCatch({
        ggplot2::ggsave(tmp_png, plot = fig$plot,
                        width = w, height = h,
                        dpi = 150, device = "png", bg = "white")
        raw_data <- readBin(tmp_png, "raw", file.info(tmp_png)$size)
        encoded <- mellio_base64_compact(raw_data)
        if (nchar(encoded) < 1400000) {
          attr(encoded, "mime") <- "image/png"
          return(encoded)
        }
      }, error = function(e) NULL)

      tmp_jpg <- tempfile(fileext = ".jpg")
      on.exit(unlink(tmp_jpg), add = TRUE)
      tryCatch({
        ggplot2::ggsave(tmp_jpg, plot = fig$plot,
                        width = w, height = h,
                        dpi = 150, device = "jpeg", bg = "white")
        raw_data <- readBin(tmp_jpg, "raw", file.info(tmp_jpg)$size)
        encoded <- mellio_base64_compact(raw_data)
        attr(encoded, "mime") <- "image/jpeg"
        return(encoded)
      }, error = function(e) NULL)
    }
  }

  if (!is.null(fig$image_raw)) {
    raw_size <- length(fig$image_raw)
    if (raw_size < 500000) {
      encoded <- mellio_base64_compact(fig$image_raw)
      mime <- switch(fig$image_format %||% "png",
        "png" = "image/png", "jpeg" = "image/jpeg",
        "image/png")
      attr(encoded, "mime") <- mime
      return(encoded)
    }
  }

  NULL
}
