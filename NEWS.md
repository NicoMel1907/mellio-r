# mellio 1.1.0

* Added Johnson-Neyman and simple-slopes support: `mellio_open()` and
  `melliotab()` now accept `johnson_neyman()` and `sim_slopes()` results from
  the 'interactions' package. Results open as a Stats card with the interval
  narration and a key-point table, plus a native, fully editable
  Johnson-Neyman figure in the web app (significance-colored slope and
  confidence band, interval bounds, observed moderator range, legend).
* Long Mellio URLs now open through a small local redirect page. Some URL
  launchers (notably inside RStudio) silently truncate very long URLs, which
  dropped the result payload; the redirect file sidesteps the limit.

# mellio 1.0.2

* Fixed lavaan bootstrap metadata extraction so Mellio supports both the
  existing numeric `bootstrap` option and the list-shaped `bootstrap` option
  introduced in lavaan 0.7-1.

# mellio 1.0.1

CRAN resubmission and figure handoff improvements.

* Added editable figure handoff for simple `ggplot2` point-range and
  point-plus-errorbar plots.
* Added editable figure handoff for simple and grouped `ggplot2` box plots.
* Increased the editable raw scatter handoff limit for `ggplot2` plots from
  1,000 to 1,500 points.
* Updated package/software/API name formatting in `DESCRIPTION` for CRAN.
* Revised examples in R help pages to follow CRAN guidance.

# mellio 1.0.0

Initial public release.

## R-to-Mellio handoff

* `mellio_open()` opens supported R objects in the Mellio web app, including
  hypothesis tests, model objects, model comparisons, tabular data, plots, and
  image files.
* Statistical results route to the Stats workspace, tabular data routes to the
  Tables workspace, and supported plots/images route to Figures.
* R handoff data includes package-version metadata so Mellio can notify users
  when an installed R package is behind the current release.

## Tables in R

* `melliotab()` creates polished, editable tables directly in R for data frames,
  model summaries, correlation matrices, and side-by-side model comparisons.
* Table helpers support APA 7 and IEEE styling, p-value formatting,
  significance stars, column spanners, table notes, and HTML/LaTeX/Markdown
  output.
* SEM/CFA results can be converted to specific table sections such as
  loadings, fit indices, paths, reliability, and observed-variable summaries.
* EFA results default to factor loadings, with optional variance and fit tables.
