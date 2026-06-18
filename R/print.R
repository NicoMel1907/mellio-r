# Print and knit_print methods

#' Print a melliotab object
#'
#' Renders the table in the RStudio Viewer pane (or browser).
#'
#' @param x A melliotab object
#' @param ... Ignored
#' @return Invisible melliotab object
#' @export
print.melliotab <- function(x, ...) {
  # Non-interactive (R CMD check, Rscript, knit): print a text summary
  # instead of launching the viewer/browser. CRAN policy forbids spawning
  # external software during examples; the viewer fallback would be
  # utils::browseURL outside RStudio.
  if (!interactive()) {
    cat("<melliotab>",
        if (!is.null(x$title)) paste0(" - ", x$title) else "",
        "\n", sep = "")
    print(x$data, ...)
    return(invisible(x))
  }

  # Build title HTML separately (above the table, no borders)
  title_html <- build_title_html(x)

  # Build gt table WITHOUT title/number (title is rendered above)
  x_notitle <- x
  x_notitle$title <- NULL
  x_notitle$number <- NULL
  gt_tbl <- mt_as_gt(x_notitle)

  # Hide table-level top border; column_labels.border.top is the top rule
  if (nzchar(title_html)) {
    gt_tbl <- gt::tab_options(gt_tbl,
      table.border.top.width = gt::px(0),
      table.border.top.style = "none"
    )
  }
  html <- as.character(gt::as_raw_html(gt_tbl))

  # Pre-generate export data for toolbar
  html_export <- mt_as_html(x)
  latex_export <- mt_as_latex(x)
  md_export <- mt_as_markdown(x)

  # Build filename label
  label <- if (!is.null(x$number)) paste0("table", x$number) else "table"

  # Build Mellio URL for print/open helpers.
  edit_url <- build_edit_url(x, mode = "table")

  www_dir <- tempfile("mellio_tab")
  dir.create(www_dir)

  html_file <- file.path(www_dir, "index.html")
  writeLines(c(
    "<!DOCTYPE html>",
    '<html lang="en"><head><meta charset="utf-8"/>',
    '<meta name="viewport" content="width=device-width, initial-scale=1.0"/>',
    "<style>body{background-color:white;margin:0;padding:0;}</style>",
    "</head><body>",
    mellio_toolbar_html(
      mode = "table",
      downloads = list(html = html_export, latex = latex_export,
                       markdown = md_export),
      label = label,
      edit_url = edit_url
    ),
    '<div id="mellio-content">', title_html, html, "</div>",
    mellio_copy_script(),
    "</body></html>"
  ), html_file)

  viewer <- mellio_viewer_function()
  viewer(html_file)
  mellio_copy_tip()
  invisible(x)
}

#' Build title HTML string for export (used by mt_copy)
#' @keywords internal
build_title_html <- function(x) {
  sc <- x$style_config
  font <- sc$font %||% "Times New Roman"
  align <- sc$table_title$align %||% "left"
  lines <- character(0)
  label <- NULL

  # Table number line (e.g., "Table 1")
  if (!is.null(x$number)) {
    label <- build_table_label(x)
    if (!is.null(label)) {
      # Convert markdown bold **text** to HTML <b>text</b>
      label <- gsub("\\*\\*(.+?)\\*\\*", "<b>\\1</b>", label)
    }
  }

  title <- NULL
  # Title line (with top margin for gap from table number)
  if (!is.null(x$title)) {
    title <- apply_case(x$title, sc$table_title$case)
    if (isTRUE(sc$table_title$italic)) {
      title <- paste0("<em>", title, "</em>")
    }
  }

  if (!is.null(label) && !is.null(title) &&
      isFALSE(sc$table_title$separate_line)) {
    joiner <- if (nzchar(sc$table_label$separator %||% "")) " " else ": "
    lines <- c(lines, paste0("<div>", label, joiner, title, "</div>"))
  } else {
    if (!is.null(label)) {
      lines <- c(lines, paste0("<div>", label, "</div>"))
    }
    if (!is.null(title)) {
      gap <- if (length(lines) > 0) ' style="margin-top: 0.5em;"' else ""
      lines <- c(lines, paste0("<div", gap, ">", title, "</div>"))
    }
  }

  if (length(lines) == 0) return("")

  paste0(
    '<div style="font-family: \'', font, '\'; ',
    'font-size: 12pt; ',
    'text-align: ', align, '; ',
    'padding-bottom: 2px; ',
    'width: 80%; margin: 0 auto;">',
    paste(lines, collapse = ""),
    '</div>'
  )
}

#' Compact base64-encode an object for URL transport
#'
#' Strips line breaks from `jsonlite::base64_enc()` output so the result
#' can be appended to a URL hash without breaking the URL.
#'
#' @param x A character string to encode.
#' @return Base64-encoded string with no embedded newlines.
#' @keywords internal
mellio_base64_compact <- function(x) {
  gsub("[\r\n]", "", jsonlite::base64_enc(x))
}

mellio_viewer_function <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      isTRUE(rstudioapi::isAvailable())) {
    return(function(url) {
      reset_file <- tempfile("mellio_viewer_reset", fileext = ".html")
      writeLines(c(
        "<!DOCTYPE html>",
        '<html lang="en"><head><meta charset="utf-8"/>',
        "</head><body></body></html>"
      ), reset_file)
      try(rstudioapi::viewer(reset_file), silent = TRUE)
      Sys.sleep(0.05)
      rstudioapi::viewer(url)
    })
  }

  opt_viewer <- getOption("viewer", NULL)
  if (is.function(opt_viewer)) {
    return(opt_viewer)
  }

  utils::browseURL
}

build_edit_url <- function(x, mode = "table") {
  base <- ms_mellio_base_url()
  style_map <- c(
    apa7 = "APA7", ieee = "IEEE"
  )

  if (mode == "table" && inherits(x, "melliotab")) {
    # Build extended grid with spanner rows and metadata
    editor_data <- build_editor_data(x)
    encoded <- mellio_base64_compact(charToRaw(editor_data$tsv))
    url <- paste0(
      sub("/+$", "", base),
      "/#data=", utils::URLencode(encoded, reserved = TRUE),
      mellio_launch_hash_param()
    )
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
        utils::URLencode(as.character(x$number),
                         reserved = TRUE))
    }
    if (!is.null(x$style)) {
      web_style <- style_map[x$style]
      if (!is.na(web_style)) {
        url <- paste0(url, "&style=", web_style)
      }
    }
    # Pass rowRoles as base64-encoded JSON
    if (length(editor_data$rowRoles) > 0) {
      roles_json <- jsonlite::toJSON(editor_data$rowRoles, auto_unbox = TRUE)
      url <- paste0(url, "&rowRoles=",
        utils::URLencode(
          mellio_base64_compact(charToRaw(as.character(roles_json))),
          reserved = TRUE
        ))
    }
    # Pass mergedRegions as base64-encoded JSON
    if (length(editor_data$mergedRegions) > 0) {
      merged_json <- jsonlite::toJSON(editor_data$mergedRegions,
                                       auto_unbox = TRUE)
      url <- paste0(url, "&mergedRegions=",
        utils::URLencode(
          mellio_base64_compact(charToRaw(as.character(merged_json))),
          reserved = TRUE
        ))
    }
    # Pass cell-level formatting/borders as base64-encoded JSON
    if (length(editor_data$cellFormats) > 0) {
      formats_json <- jsonlite::toJSON(editor_data$cellFormats,
                                       auto_unbox = TRUE)
      url <- paste0(url, "&cellFormats=",
        utils::URLencode(
          mellio_base64_compact(charToRaw(as.character(formats_json))),
          reserved = TRUE
        ))
    }
    if (length(editor_data$cellBorders) > 0) {
      borders_json <- jsonlite::toJSON(editor_data$cellBorders,
                                       auto_unbox = TRUE)
      url <- paste0(url, "&cellBorders=",
        utils::URLencode(
          mellio_base64_compact(charToRaw(as.character(borders_json))),
          reserved = TRUE
        ))
    }
    # Pass source text
    if (!is.null(x$source)) {
      url <- paste0(url, "&source=",
        utils::URLencode(x$source, reserved = TRUE))
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
    return(url)
  }

  base
}

# ── Viewer copy helpers ──

#' Generate JavaScript for Cmd+C / Ctrl+C clipboard copy in Viewer
#' @param figure If TRUE, swaps img src with base64 before copying
#' @return HTML script tag string
#' @keywords internal
mellio_copy_script <- function(figure = FALSE) {
  js <- "
document.addEventListener('keydown', function(e) {
  if ((e.ctrlKey || e.metaKey) && e.key === 'c') {
    e.preventDefault();
    if (typeof mellioCopy === 'function') mellioCopy();
  }
});"
  paste0("<script>", js, "</script>")
}

#' Show once-per-session copy tip in console
#' @keywords internal
mellio_copy_tip <- function() {
  if (!isTRUE(getOption("mellio.copy_tip_shown"))) {
    cli::cli_inform(c(
      "i" = "Use the Copy button in the Viewer for Word-ready copy, or the Mo logo to open in Mellio."
    ))
    options(mellio.copy_tip_shown = TRUE)
  }
}

# ── Viewer toolbar ──

#' Generate toolbar HTML for viewer
#'
#' Creates a quiet local toolbar with an optional Mellio link.
#'
#' @param mode "table" or "figure"
#' @param copy_img_b64 Base64 data URI for image copy (figures only)
#' @param downloads Named list: html, latex, markdown (text), png_b64 (raw base64)
#' @param label Filename prefix for downloads (e.g., "table1", "figure1")
#' @return Character string of HTML (style + toolbar + scripts)
#' @keywords internal
mellio_toolbar_html <- function(mode, copy_img_b64 = NULL,
                                downloads = list(), label = "mellio",
                                edit_url = NULL) {

  # --- CSS ---
  css <- '
<style>
.mellio-toolbar {
  display: flex; align-items: center; padding: 3px 18px;
  border-bottom: 1px solid #eee9df; background: #ffffff;
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", system-ui, sans-serif;
  box-sizing: border-box; font-size: 12px; height: 38px;
  position: sticky; top: 0; z-index: 100;
}
.mellio-toolbar .mellio-tool-button {
  width: 30px; height: 30px; display: inline-flex; align-items: center;
  justify-content: center; border: 0; border-radius: 7px; background: transparent;
  color: #17140f; cursor: pointer; position: relative; padding: 0;
  transition: background-color .12s, box-shadow .12s, color .12s;
}
.mellio-toolbar .mellio-tool-button:hover {
  background: #fbfaf7; box-shadow: 0 1px 2px rgba(31, 25, 18, .06);
}
.mellio-toolbar .mellio-tool-button.mellio-copied {
  color: #1E4D3A; background: #edf6f1;
}
.mellio-toolbar .mellio-tool-button::after {
  content: attr(data-tooltip); position: absolute; left: 0; top: calc(100% + 7px);
  padding: 5px 8px; border-radius: 6px; background: #17140f; color: #ffffff;
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", system-ui, sans-serif;
  font-size: 11px; font-weight: 500; line-height: 1; white-space: nowrap;
  opacity: 0; transform: translateY(-2px); pointer-events: none;
  transition: opacity .12s, transform .12s;
}
.mellio-toolbar .mellio-tool-button:hover::after,
.mellio-toolbar .mellio-tool-button:focus-visible::after {
  opacity: 1; transform: translateY(0);
}
.mellio-toolbar .mellio-tool-button svg {
  width: 17px; height: 17px; display: block;
}
.mellio-toolbar .mellio-logo-link {
  width: 32px; height: 32px; margin-left: auto; display: inline-flex;
  align-items: center; justify-content: center; border-radius: 8px;
  color: #17140f; text-decoration: none; cursor: pointer; position: relative;
  transition: background-color .12s, box-shadow .12s;
}
.mellio-toolbar .mellio-logo-link:hover {
  background: #fbfaf7; box-shadow: 0 1px 2px rgba(31, 25, 18, .06);
}
.mellio-toolbar .mellio-logo-link::after {
  content: attr(data-tooltip); position: absolute; right: 0; top: calc(100% + 7px);
  padding: 5px 8px; border-radius: 6px; background: #17140f; color: #ffffff;
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", system-ui, sans-serif;
  font-size: 11px; font-weight: 500; line-height: 1; white-space: nowrap;
  opacity: 0; transform: translateY(-2px); pointer-events: none;
  transition: opacity .12s, transform .12s;
}
.mellio-toolbar .mellio-logo-link:hover::after,
.mellio-toolbar .mellio-logo-link:focus-visible::after {
  opacity: 1; transform: translateY(0);
}
.mellio-toolbar .mellio-tab-logo {
  width: 30px; height: 30px; display: block;
}
.mellio-toolbar .mellio-tab-logo text {
  font-family: "Newsreader", "Source Serif 4", Georgia, serif;
}
</style>'

  # --- Buttons ---
  logo_svg <- paste0(
    '<svg class="mellio-tab-logo" viewBox="0 0 220 220" fill="none" ',
    'xmlns="http://www.w3.org/2000/svg" aria-hidden="true">',
    '<text x="110" y="130" text-anchor="middle" font-size="144" ',
    'font-weight="600" letter-spacing="-6" fill="#1A1814">M',
    '<tspan dx="2" font-size="100" font-style="italic" ',
    'font-weight="500" fill="#1E4D3A">o</tspan></text>',
    '<line x1="48" y1="152" x2="172" y2="152" stroke="#1A1814" ',
    'stroke-width="3" stroke-linecap="round" stroke-opacity="0.42"/>',
    '</svg>'
  )

  copy_svg <- paste0(
    '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" ',
    'stroke-width="2" stroke-linecap="round" stroke-linejoin="round" ',
    'aria-hidden="true">',
    '<rect x="9" y="9" width="13" height="13" rx="2" ry="2"></rect>',
    '<path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9',
    'a2 2 0 0 1 2 2v1"></path>',
    '</svg>'
  )
  copy_button <- paste0(
    '<button class="mellio-tool-button mellio-copy-button" type="button" ',
    'onclick="mellioCopy(this)" title="Copy for Word" ',
    'data-tooltip="Copy for Word" aria-label="Copy for Word">',
    copy_svg,
    '</button>'
  )

  if (!is.null(edit_url)) {
    brand <- paste0(
      '<a class="mellio-logo-link" href="', edit_url,
      '" target="_blank" title="Open in Mellio" ',
      'data-tooltip="Open in Mellio" ',
      'aria-label="Open in Mellio">',
      logo_svg, '</a>'
    )
  } else {
    brand <- paste0(
      '<span class="mellio-logo-link" ',
      'title="Mellio" data-tooltip="Mellio">',
      logo_svg, '</span>'
    )
  }
  toolbar_div <- paste0(
    '<div class="mellio-toolbar">',
    copy_button,
    brand,
    '</div>'
  )

  # --- JS data variables ---
  js_escape <- function(s) {
    s <- gsub("\\\\", "\\\\\\\\", s)
    s <- gsub("`", "\\\\`", s)
    s <- gsub("\\$\\{", "\\\\${", s)
    s
  }

  js_vars <- paste0(
    "<script>\n",
    "window.MELLIO_MODE='", mode, "';\n",
    "window.MELLIO_LABEL='", label, "';\n"
  )

  if (!is.null(downloads$html)) {
    js_vars <- paste0(js_vars, "window.MELLIO_HTML=`", js_escape(downloads$html), "`;\n")
  }
  if (!is.null(downloads$latex)) {
    js_vars <- paste0(js_vars, "window.MELLIO_LATEX=`", js_escape(downloads$latex), "`;\n")
  }
  if (!is.null(downloads$markdown)) {
    js_vars <- paste0(js_vars, "window.MELLIO_MARKDOWN=`", js_escape(downloads$markdown), "`;\n")
  }
  if (!is.null(copy_img_b64)) {
    js_vars <- paste0(js_vars, "window.MELLIO_IMG_B64='", copy_img_b64, "';\n")
  }
  if (!is.null(downloads$png_b64)) {
    js_vars <- paste0(js_vars, "window.MELLIO_PNG_B64='", downloads$png_b64, "';\n")
  }

  js_vars <- paste0(js_vars, "</script>")

  # --- JS functions ---
  js_functions <- '<script>
function mellioCopy(b){
  b=b||null;
  if(window.MELLIO_MODE==="table"){
    var e=document.getElementById("mellio-content");
    var cl=e.cloneNode(true);
    var countCols=function(row){var n=0;
      Array.prototype.forEach.call(row.children,function(c){n+=parseInt(c.getAttribute("colspan")||"1",10)||1});
      return n;};
    var ncols=1;
    Array.prototype.forEach.call(cl.querySelectorAll("thead tr,tbody tr"),function(row){
      ncols=Math.max(ncols,countCols(row));
    });
    var clearTopBorders=function(nodes){
      Array.prototype.forEach.call(nodes,function(c){
        c.style.borderTopStyle="none";c.style.borderTopWidth="0";
      });
    };
    var clearBottomBorders=function(nodes){
      Array.prototype.forEach.call(nodes,function(c){
        c.style.borderBottomStyle="none";c.style.borderBottomWidth="0";
      });
    };
    var ruleRow=function(){var tr=document.createElement("tr");
      var td=document.createElement("td");
      td.setAttribute("colspan",ncols);
      td.style.cssText="border-bottom:2px solid black;height:0;padding:0;font-size:0;line-height:0";
      tr.appendChild(td);return tr;};
    var thead=cl.querySelector("thead");
    if(thead){
      clearTopBorders(thead.querySelectorAll("tr:first-child th,tr:first-child td"));
      thead.insertBefore(ruleRow(),thead.firstChild);
    }
    var tbody=cl.querySelector("tbody");
    if(tbody){
      clearBottomBorders(tbody.querySelectorAll("tr:last-child th,tr:last-child td"));
      tbody.appendChild(ruleRow());
    }
    var t=document.createElement("div");
    t.style.cssText="position:fixed;left:-9999px";
    t.appendChild(cl);
    document.body.appendChild(t);
    var r=document.createRange();
    r.selectNodeContents(t);
    var s=window.getSelection();
    s.removeAllRanges();
    s.addRange(r);
    document.execCommand("copy");
    s.removeAllRanges();
    document.body.removeChild(t);
    if(b) showCopied(b);
  } else {
    var e=document.getElementById("mellio-content");
    var cl=e.cloneNode(true);
    var im=cl.querySelector("img");
    if(im&&window.MELLIO_IMG_B64) im.src=window.MELLIO_IMG_B64;
    var t=document.createElement("div");
    t.style.cssText="position:fixed;left:-9999px";
    t.appendChild(cl);
    document.body.appendChild(t);
    var r=document.createRange();
    r.selectNodeContents(t);
    var s=window.getSelection();
    s.removeAllRanges();
    s.addRange(r);
    document.execCommand("copy");
    s.removeAllRanges();
    document.body.removeChild(t);
    if(b) showCopied(b);
  }
}
function showCopied(b){
  if(!b) return;
  var tooltip=b.getAttribute("data-tooltip")||"Copy for Word";
  var title=b.getAttribute("title")||tooltip;
  var aria=b.getAttribute("aria-label")||title;
  b.setAttribute("data-tooltip","Copied");
  b.setAttribute("title","Copied");
  b.setAttribute("aria-label","Copied");
  b.classList.add("mellio-copied");
  setTimeout(function(){
    b.setAttribute("data-tooltip",tooltip);
    b.setAttribute("title",title);
    b.setAttribute("aria-label",aria);
    b.classList.remove("mellio-copied");
  },1400);
}
function mellioDownload(f){
  var l=window.MELLIO_LABEL||"mellio";
  var c,n,m;
  switch(f){
    case "html":c=window.MELLIO_HTML;n=l+".html";m="text/html";break;
    case "tex":c=window.MELLIO_LATEX;n=l+".tex";m="text/x-latex";break;
    case "md":c=window.MELLIO_MARKDOWN;n=l+".md";m="text/markdown";break;
    case "png":
      var r=atob(window.MELLIO_PNG_B64);
      var a=new Uint8Array(r.length);
      for(var i=0;i<r.length;i++) a[i]=r.charCodeAt(i);
      triggerDownload(new Blob([a],{type:"image/png"}),l+".png");
      return;
  }
  triggerDownload(new Blob([c],{type:m+";charset=utf-8"}),n);
}
function triggerDownload(b,f){
  var u=URL.createObjectURL(b);
  var a=document.createElement("a");
  a.href=u;a.download=f;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(u);
}
</script>'

  paste0(css, "\n", toolbar_div, "\n", js_vars, "\n", js_functions)
}

#' Knit print method for R Markdown / Quarto
#'
#' Automatically detects the output format and uses the appropriate backend.
#'
#' @param x A melliotab object
#' @param ... Additional arguments
#' @return knitr output
#' @export
knit_print.melliotab <- function(x, ...) {
  output_format <- knitr::opts_knit$get("rmarkdown.pandoc.to")

  if (identical(output_format, "latex")) {
    latex_str <- mt_as_latex(x)
    knitr::asis_output(latex_str)
  } else if (identical(output_format, "docx")) {
    knitr::asis_output(mt_as_markdown(x))
  } else {
    gt_tbl <- mt_as_gt(x)
    knitr::knit_print(gt_tbl, ...)
  }
}
