# HTML export via gt
# Style-aware conversion from melliotab -> gt table object

#' Convert a melliotab object to a gt table
#'
#' @param x A melliotab object
#' @return A gt table object
#' @export
mt_as_gt <- function(x) {
  if (!inherits(x, "melliotab")) {
    cli::cli_abort("{.arg x} must be a melliotab object.")
  }
  rlang::check_installed("gt", reason = "for HTML table rendering")

  data <- x$data
  sc <- x$style_config
  col_types <- x$column_types

  # Build gt table
  gt_tbl <- gt::gt(data)

  # ── Title and table number ──
  title_text <- build_table_label(x)
  subtitle_text <- NULL

  if (!is.null(x$title)) {
    if (!is.null(title_text) && isTRUE(sc$table_title$separate_line)) {
      # Table label exists -> title goes on separate line as subtitle
      subtitle_text <- x$title
    } else if (!is.null(title_text)) {
      # Table label exists, same line
      inline_title <- apply_case(x$title, sc$table_title$case)
      if (isTRUE(sc$table_title$italic)) {
        inline_title <- paste0("*", inline_title, "*")
      }
      title_text <- paste(title_text, inline_title)
    } else {
      # No table label -> title is the main title
      title_text <- x$title
    }
  }

  if (!is.null(title_text) || !is.null(subtitle_text)) {
    # Apply case/formatting per style
    if (!is.null(subtitle_text)) {
      subtitle_text <- apply_case(subtitle_text, sc$table_title$case)
    }

    # When title is standalone (no table label), apply case/italic
    if (!is.null(title_text) && is.null(x$number)) {
      title_text <- apply_case(title_text, sc$table_title$case)
      if (isTRUE(sc$table_title$italic)) {
        title_text <- paste0("*", title_text, "*")
      }
    }

    gt_tbl <- gt::tab_header(
      gt_tbl,
      title = if (!is.null(title_text)) gt::md(title_text),
      subtitle = if (!is.null(subtitle_text)) {
        if (isTRUE(sc$table_title$italic)) {
          gt::md(paste0("*", subtitle_text, "*"))
        } else {
          subtitle_text
        }
      }
    )
  }

  # ── Font and sizing ──
  gt_tbl <- gt::tab_options(
    gt_tbl,
    table.width = gt::pct(80),
    table.align = "center",
    table.font.names = sc$font,
    table.font.size = gt::px(sc$font_size),
    data_row.padding = gt::px(sc$cell_padding[1]),
    column_labels.padding = gt::px(sc$cell_padding[1])
  )

  # ── Line spacing ──
  if (sc$spacing == "double") {
    gt_tbl <- gt::tab_options(gt_tbl, data_row.padding = gt::px(8))
  }

  # ── Borders ──
  # Build all border options in one list to avoid override conflicts
  border_opts <- list(
    table_body.hlines.style = "hidden",
    heading.border.bottom.style = "hidden"
  )

  has_heading <- !is.null(title_text) || !is.null(subtitle_text)

  if (isTRUE(sc$borders$top)) {
    border_opts$table.border.top.width <- gt::px(2)
    border_opts$table.border.top.style <- "solid"
    border_opts$table.border.top.color <- "black"
    if (!has_heading) {
      # No heading: column labels are the top, they need the border
      border_opts$column_labels.border.top.width <- gt::px(2)
      border_opts$column_labels.border.top.style <- "solid"
      border_opts$column_labels.border.top.color <- "black"
    } else {
      # Heading present: table.border.top is the top rule,
      # no extra border between heading and column labels
      border_opts$column_labels.border.top.style <- "hidden"
    }
  } else {
    border_opts$table.border.top.style <- "hidden"
    border_opts$column_labels.border.top.style <- "hidden"
  }

  if (isTRUE(sc$borders$header_bottom)) {
    border_opts$column_labels.border.bottom.width <- gt::px(1)
    border_opts$column_labels.border.bottom.style <- "solid"
    border_opts$column_labels.border.bottom.color <- "black"
    border_opts$table_body.border.top.width <- gt::px(1)
    border_opts$table_body.border.top.style <- "solid"
    border_opts$table_body.border.top.color <- "black"
  } else {
    border_opts$column_labels.border.bottom.style <- "hidden"
    border_opts$table_body.border.top.style <- "hidden"
  }

  if (isTRUE(sc$borders$bottom)) {
    border_opts$table_body.border.bottom.width <- gt::px(2)
    border_opts$table_body.border.bottom.style <- "solid"
    border_opts$table_body.border.bottom.color <- "black"
    border_opts$table.border.bottom.style <- "none"
  } else {
    border_opts$table_body.border.bottom.style <- "hidden"
    border_opts$table.border.bottom.style <- "none"
  }

  gt_tbl <- do.call(gt::tab_options, c(list(data = gt_tbl), border_opts))

  # Midrule: add border-top on first body row (<td> cells render in
  # RStudio Viewer, unlike <th> borders which get suppressed)
  if (isTRUE(sc$borders$header_bottom)) {
    gt_tbl <- gt::tab_style(
      gt_tbl,
      style = gt::cell_borders(
        sides = "top",
        color = "#000000",
        weight = gt::px(1)
      ),
      locations = gt::cells_body(rows = 1)
    )
  }

  # Bottom rule: same Viewer workaround — border-bottom on last body row
  if (isTRUE(sc$borders$bottom)) {
    gt_tbl <- gt::tab_style(
      gt_tbl,
      style = gt::cell_borders(
        sides = "bottom",
        color = "#000000",
        weight = gt::px(2)
      ),
      locations = gt::cells_body(rows = nrow(data))
    )
  }

  if (isTRUE(sc$borders$internal_rows)) {
    gt_tbl <- gt::tab_options(
      gt_tbl,
      table_body.hlines.width = gt::px(1),
      table_body.hlines.style = "solid",
      table_body.hlines.color = "black"
    )
  }

  if (isTRUE(sc$borders$vertical)) {
    gt_tbl <- gt::tab_options(
      gt_tbl,
      table_body.vlines.width = gt::px(1),
      table_body.vlines.style = "solid",
      table_body.vlines.color = "black"
    )
    # Also add vertical borders to header row
    gt_tbl <- gt::tab_style(
      gt_tbl,
      style = gt::cell_borders(
        sides = c("left", "right"),
        color = "black",
        weight = gt::px(1)
      ),
      locations = gt::cells_column_labels()
    )
  }

  # ── Alignment ──
  # First column (stub) left-aligned, rest right-aligned
  if (ncol(data) > 1) {
    gt_tbl <- gt::cols_align(gt_tbl, align = "left", columns = 1)
    gt_tbl <- gt::cols_align(gt_tbl, align = "right",
                              columns = seq(2, ncol(data)))
  }

  # ── Italic stat headers (APA style) ──
  if (isTRUE(sc$italic_stat_headers)) {
    italic_cols <- which(vapply(names(data), is_stat_symbol, logical(1)))
    if (length(italic_cols) > 0) {
      gt_tbl <- gt::tab_style(
        gt_tbl,
        style = gt::cell_text(style = "italic"),
        locations = gt::cells_column_labels(columns = italic_cols)
      )
    }
  }

  # ── Optional header background ──
  if (!is.null(sc$header_bg)) {
    gt_tbl <- gt::tab_style(
      gt_tbl,
      style = list(
        gt::cell_fill(color = sc$header_bg),
        gt::cell_text(color = sc$header_color %||% "#ffffff", weight = "bold")
      ),
      locations = gt::cells_column_labels()
    )
  }

  # ── Optional zebra striping ──
  if (isTRUE(sc$zebra_striping)) {
    gt_tbl <- gt::opt_row_striping(gt_tbl)
  }

  # ── Spanners ──
  opts <- x$options %||% list()
  is_comparison <- isTRUE(opts$is_comparison)
  dependent_variable_label <- opts$dependent_variable_label
  has_dependent_variable_spanner <- FALSE
  has_dependent_variable_label_spanner <- FALSE
  for (sp in x$spanners) {
    spanner_label <- sp$label
    spanner_id <- sp$label
    if (is_comparison && identical(sp$label, "Dependent variable:")) {
      spanner_label <- gt::md("*Dependent variable:*")
      spanner_id <- "mellio_dependent_variable"
      has_dependent_variable_spanner <- TRUE
    }
    if (is_comparison && !is.null(dependent_variable_label) &&
        identical(sp$label, dependent_variable_label)) {
      spanner_id <- "mellio_dependent_variable_label"
      has_dependent_variable_label_spanner <- TRUE
    }
    gt_tbl <- gt::tab_spanner(
      gt_tbl,
      label = spanner_label,
      columns = sp$columns,
      level = sp$level %||% 1L,
      id = spanner_id
    )
  }

  if (is_comparison) {
    model_cols <- names(data)[-1]
    if (length(model_cols) > 0L) {
      stub_width <- 170
      model_width <- 160
      table_width <- stub_width + model_width * length(model_cols)
      gt_tbl <- gt::tab_options(
        gt_tbl,
        table.width = gt::px(table_width)
      )
      gt_tbl <- gt::cols_width(
        gt_tbl,
        .list = list(
          rlang::new_formula(names(data)[1], gt::px(stub_width)),
          rlang::new_formula(model_cols, gt::px(model_width))
        )
      )
      gt_tbl <- gt::cols_align(gt_tbl, align = "center", columns = model_cols)
    }

    if (has_dependent_variable_spanner) {
      gt_tbl <- gt::tab_style(
        gt_tbl,
        style = gt::cell_text(style = "italic", weight = "normal"),
        locations = gt::cells_column_spanners(
          spanners = "mellio_dependent_variable"
        )
      )
    }
    if (has_dependent_variable_label_spanner) {
      gt_tbl <- gt::opt_css(
        gt_tbl,
        css = paste(
          "#mellio_dependent_variable .gt_column_spanner {",
          "width: 72% !important;",
          "}"
        )
      )
      gt_tbl <- gt::tab_style(
        gt_tbl,
        style = gt::cell_text(weight = "bold"),
        locations = gt::cells_column_spanners(
          spanners = "mellio_dependent_variable_label"
        )
      )
      gt_tbl <- gt::opt_css(
        gt_tbl,
        css = paste(
          "#mellio_dependent_variable_label,",
          "#mellio_dependent_variable_label .gt_column_spanner {",
          "border-bottom-style: none !important;",
          "border-bottom-width: 0 !important;",
          "}",
          "#mellio_dependent_variable_label .gt_column_spanner {",
          "font-weight: bold !important;",
          "}"
        )
      )
    }

    stat_start <- opts$comparison_stat_start
    if (!is.null(stat_start) && length(stat_start) == 1L &&
        !is.na(stat_start) && stat_start > 1L) {
      gt_tbl <- gt::tab_style(
        gt_tbl,
        style = gt::cell_borders(
          sides = "top",
          color = "#000000",
          weight = gt::px(1)
        ),
        locations = gt::cells_body(rows = as.integer(stat_start))
      )
    }
  }

  # ── Notes ──
  if (!is.null(x$note)) {
    note_text <- capitalize_note_after_label(x$note, sc$notes$general_label)
    if (isTRUE(sc$notes$general_label_italic) && !is.null(sc$notes$general_label)) {
      note_text <- paste0("*", sc$notes$general_label, "* ", note_text)
    } else if (!is.null(sc$notes$general_label)) {
      note_text <- paste0(sc$notes$general_label, " ", note_text)
    }
    gt_tbl <- gt::tab_source_note(gt_tbl, source_note = gt::md(note_text))
  }

  if (!is.null(x$source) && !is.null(sc$notes$source_label)) {
    source_text <- paste(sc$notes$source_label, x$source)
    gt_tbl <- gt::tab_source_note(gt_tbl, source_note = source_text)
  }

  # ── Label alignment ──
  if (!is.null(sc$label_align) && sc$label_align == "center") {
    gt_tbl <- gt::tab_options(gt_tbl, heading.align = "center")
  } else {
    gt_tbl <- gt::tab_options(gt_tbl, heading.align = "left")
  }

  gt_tbl
}

#' Build table label text (e.g., "Table 1" or "TABLE I")
#' @keywords internal
build_table_label <- function(mt) {
  sc <- mt$style_config
  if (is.null(mt$number)) return(NULL)

  num <- mt$number
  if (isTRUE(sc$table_label$numbering == "roman")) {
    num <- to_roman(as.integer(num))
  }

  separator <- sc$table_label$separator %||% ""
  label <- paste0(sc$table_label$prefix, " ", num, separator)

  # Apply case
  if (isTRUE(sc$table_label$case == "upper")) {
    label <- toupper(label)
  }

  # Bold
  if (isTRUE(sc$table_label$bold)) {
    label <- paste0("**", label, "**")
  }

  label
}

#' Apply case transformation
#' @keywords internal
apply_case <- function(text, case_type) {
  switch(case_type %||% "title",
    "upper" = toupper(text),
    "sentence" = to_sentence_case(text),
    "title" = to_title_case(text),
    text
  )
}

#' Get HTML string from a melliotab object
#'
#' Produces a standalone HTML fragment with the title rendered above the table
#' (outside the gt borders), matching the viewer layout.
#'
#' @param x A melliotab object
#' @return Character string of HTML
#' @export
mt_as_html <- function(x) {
  # Build title HTML separately (above the table, no borders)
  title_html <- build_title_html(x)

  # Build gt table WITHOUT title/number
  x_notitle <- x
  x_notitle$title <- NULL
  x_notitle$number <- NULL
  gt_tbl <- mt_as_gt(x_notitle)

  # Hide table-level top border when title is present
  if (nzchar(title_html)) {
    gt_tbl <- gt::tab_options(gt_tbl,
      table.border.top.width = gt::px(0),
      table.border.top.style = "none"
    )
  }

  table_html <- gt::as_raw_html(gt_tbl)
  paste0(title_html, as.character(table_html))
}
