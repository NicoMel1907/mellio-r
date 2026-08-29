## Submission

This is an urgent maintenance release requested by the `lavaan` maintainers.

In this version I have:

* fixed `lavaan` bootstrap metadata extraction so `mellio` supports both the
  existing numeric `bootstrap` option and the list-shaped `bootstrap` option
  introduced in `lavaan` 0.7-1; and
* bumped the package version from 1.0.1 to 1.0.2.

The short interval since the previous CRAN release is intentional because this
release unblocks `lavaan` reverse-dependency checks.

## Test environments

* local macOS (aarch64-apple-darwin20), R 4.4.0
  (`R CMD check --as-cran`, with `_R_CHECK_FORCE_SUGGESTS_=false`)
* win-builder, R-oldrelease 4.5.3 (2026-03-11 ucrt)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Days since last update: 2

This note is expected for this urgent maintenance release.

## Notes for the reviewer

* `mellio_open()` hands an R object to the Mellio web app by encoding it in a
  URL fragment and opening the browser. The fragment is not sent in the HTTP
  request. Examples that would open a browser are wrapped in `\donttest{}` with
  an `interactive()` guard, so no example launches a browser during checks.
* Several packages in Suggests are used conditionally for optional model and
  plot inputs, each behind `requireNamespace()`.
