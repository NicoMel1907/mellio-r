## Submission

This is a feature release.

In this version I have:

* added support for `johnson_neyman()` and `sim_slopes()` results from the
  'interactions' package, opening as editable result cards and figures in the
  'Mellio' web app ('interactions' added to Suggests);
* routed long 'Mellio' URLs through a small local redirect page, because some
  IDE URL launchers silently truncate very long URLs; and
* bumped the package version from 1.0.2 to 1.1.0.

## Test environments

* local macOS (aarch64-apple-darwin20), R 4.4.0
  (`R CMD check --as-cran`, with `_R_CHECK_FORCE_SUGGESTS_=false`)
* win-builder, R-release 4.6.1 (2026-06-24 ucrt) — Status: OK
* win-builder, R-devel (2026-08-27 r90452 ucrt) — Status: OK

## R CMD check results

0 errors | 0 warnings | 2 notes (local only)

* Both local NOTEs are environment artifacts: several optional Suggests
  packages are not installed on the check machine, and "unable to verify
  current time" reflects the local clock service. Neither reproduces on
  win-builder.
