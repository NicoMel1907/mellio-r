## Submission

This is a resubmission. In this version I have:

* wrapped package, software, and API names in single quotes in `DESCRIPTION`;
* replaced `\dontrun{}` examples with `\donttest{}` examples or executable
  temporary-file examples, as appropriate; and
* bumped the package version from 1.0.0 to 1.0.1 for this resubmission.

## Test environments

* local macOS (aarch64-apple-darwin20), R 4.4.0
  (`R CMD check --no-manual --ignore-vignettes --no-build-vignettes`,
  with `_R_CHECK_FORCE_SUGGESTS_=false`)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking package dependencies ... NOTE
  Packages suggested but not available for checking:
    'webshot2', 'ggExtra', 'plotly', 'randomForest', 'brms', 'rstanarm'

This note is local-only. The packages are optional integrations guarded by
`requireNamespace()`, and the check was run with `_R_CHECK_FORCE_SUGGESTS_=false`
because these optional packages are not installed in this local environment.

## Notes for the reviewer

* `mellio_open()` hands an R object to the Mellio web app by encoding it in a
  URL fragment and opening the browser. The fragment is not sent in the HTTP
  request. Examples that would open a browser are wrapped in `\donttest{}` with
  an `interactive()` guard, so no example launches a browser during checks.
* Several packages in Suggests are used conditionally for optional model and
  plot inputs, each behind `requireNamespace()`.
