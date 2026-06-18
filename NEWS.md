# mellio 0.9.0

Public beta release.

## Unified R bridge

* Added `mellio_open()` as the unified web-editor entry point for statistical
  results, tables, and figures.
* Renamed the stats payload API: `ms_payload()` is now `mellio_payload()`,
  and `ms_to_json()` is now `mellio_to_json()`.
* Renamed hierarchical model comparison for Stats cards from `ms_compare()` to
  `mellio_compare()`. `mt_compare()` is unchanged because it produces
  manuscript tables.
* Renamed the RStudio addin binding from `ms_addin_send()` to
  `mellio_addin_send()`.
* Removed the old editor entry points (`ms_edit()`, `mt_edit()`, `mf_edit()`)
  so R-to-Mellio handoff is canalized through `mellio_open()`.
* `mellio_open()` now routes by intent. Statistical inputs (`lm`, `glm`,
  `aov`, `htest`, `lavaan`, model-shaped `melliotab` carrying p-value
  columns, etc.) go to the Stats workspace so they can be narrated.
  Pure tabular inputs (bare `data.frame`, `matrix`, `table`, and
  `melliotab` built from a data.frame without p-values) now route
  directly to the Tables workspace via the existing `#data=` handoff
  — no more round-trip through Stats. The "Open in Tables" button on
  Stats cards remains as the escape valve for promoting statistical
  results to a manuscript table.
* Added `mellio_capture()` for sending the current base R plotting device to
  Mellio as a figure.
