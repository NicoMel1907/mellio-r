.onLoad <- function(libname, pkgname) {
  # Register knit_print methods
  if (requireNamespace("knitr", quietly = TRUE)) {
    registerS3method("knit_print", "melliotab", knit_print.melliotab,
                     envir = asNamespace("knitr"))
  }
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("mellio: use mellio_open() for Mellio, melliotab() for R tables")
}
