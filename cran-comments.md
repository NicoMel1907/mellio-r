## Submission

This is a resubmission. In this version I have:

* wrapped package, software, and API names in single quotes in `DESCRIPTION`;
* replaced `\dontrun{}` examples with `\donttest{}` examples or executable
  temporary-file examples, as appropriate; and
* bumped the package version from 1.0.0 to 1.0.1 for this resubmission.

The 1.0.1 release also includes small `ggplot2` figure handoff improvements
documented in `NEWS.md`.

## Test environments

* local macOS (aarch64-apple-darwin20), R 4.4.0
  (`R CMD check --no-manual --ignore-vignettes --no-build-vignettes`,
  with `_R_CHECK_FORCE_SUGGESTS_=false`)
* win-builder, R-release 4.6.1 (2026-06-24 ucrt)
* win-builder, R-devel (2026-06-26 r90195 ucrt)

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Melih Sahin <nicomelpro@pm.me>'
  New submission
  Possibly misspelled words in DESCRIPTION:
    APA (11:55)

This note is expected for a first CRAN submission. `APA` is a standard acronym,
and the package/software/API names in `DESCRIPTION` have been wrapped in single
quotes as requested in the previous review.

## Notes for the reviewer

* `mellio_open()` hands an R object to the Mellio web app by encoding it in a
  URL fragment and opening the browser. The fragment is not sent in the HTTP
  request. Examples that would open a browser are wrapped in `\donttest{}` with
  an `interactive()` guard, so no example launches a browser during checks.
* Several packages in Suggests are used conditionally for optional model and
  plot inputs, each behind `requireNamespace()`.
