## Submission

This is a new submission.

## Test environments

* local macOS (aarch64-apple-darwin20), R 4.4.0
* win-builder, R-devel and R-release

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Melih Sahin <nicomelpro@pm.me>'
  New submission

This note is expected for a first submission.

## Notes for the reviewer

* `mellio_open()` hands an R object to the Mellio web app by encoding it in a
  URL fragment and opening the browser. The fragment is not sent in the HTTP
  request. Examples that would open a browser are wrapped in `\donttest{}` with
  an `interactive()` guard, so no example launches a browser during checks.
* Several packages in Suggests are used conditionally for optional model and
  plot inputs, each behind `requireNamespace()`.
