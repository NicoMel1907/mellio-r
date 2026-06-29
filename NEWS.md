# mellio 1.0.1

Documentation-only CRAN resubmission of the initial public release.

* Updated package/software/API name formatting in `DESCRIPTION` for CRAN.
* Revised examples in R help pages to follow CRAN guidance.
* No user-facing behavior changes.

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
