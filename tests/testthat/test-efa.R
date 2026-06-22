# Tests for the psych EFA / PCA bridge extractor (bridge-extract-efa.R).

test_that("mellio_payload.fa builds an EFA loadings table card", {
  skip_if_not_installed("psych")
  vars <- c("mpg", "hp", "wt", "disp", "drat", "qsec")
  fit <- suppressWarnings(
    psych::fa(mtcars[, vars], nfactors = 2, rotate = "oblimin")
  )
  p <- mellio_payload(fit)

  expect_equal(p$card_kind, "table")
  expect_equal(p$type, "efa_loadings")
  expect_equal(p$type_label, "Exploratory factor analysis")
  expect_equal(p$fields$table_type, "factor_loadings")
  expect_equal(p$fields$n_factors, 2L)

  # One row per measured variable
  expect_equal(length(p$fields$rows), length(vars))

  # Columns: Variable + 2 factor columns + h2 + u2
  col_keys <- vapply(p$fields$columns, function(col) col$key, character(1))
  expect_equal(col_keys[[1]], "variable")
  expect_true(all(c("factor_1", "factor_2", "h2", "u2") %in% col_keys))

  col_labels <- vapply(p$fields$columns, function(col) col$label, character(1))
  expect_true(all(c("Factor 1", "Factor 2") %in% col_labels))

  # Rows carry the variable name and numeric loadings
  row1 <- p$fields$rows[[1]]
  expect_true(is.character(row1$variable))
  expect_true(is.numeric(row1$factor_1))

  # Communality + uniqueness sum to 1 for every variable
  for (row in p$fields$rows) {
    expect_equal(row$h2 + row$u2, 1, tolerance = 1e-6)
  }

  # EFA metadata rides along in fields
  expect_equal(p$fields$rotation, "oblimin")
  expect_true(is.character(p$fields$method) && nzchar(p$fields$method))
  expect_true(is.numeric(p$fields$variance_explained))
})

test_that("mellio_payload.principal builds a PCA loadings table card", {
  skip_if_not_installed("psych")
  vars <- c("mpg", "hp", "wt", "disp", "drat", "qsec")
  fit <- suppressWarnings(
    psych::principal(mtcars[, vars], nfactors = 2, rotate = "varimax")
  )
  p <- mellio_payload(fit)

  expect_equal(p$card_kind, "table")
  expect_equal(p$type, "pca_loadings")
  expect_equal(p$type_label, "Principal component analysis")
  expect_equal(p$fields$table_type, "component_loadings")
  expect_equal(p$fields$n_factors, 2L)
  expect_equal(length(p$fields$rows), length(vars))

  # PCA labels its dimensions "Component", not "Factor"
  col_keys <- vapply(p$fields$columns, function(col) col$key, character(1))
  expect_true(all(c("component_1", "component_2") %in% col_keys))

  col_labels <- vapply(p$fields$columns, function(col) col$label, character(1))
  expect_true(all(c("Component 1", "Component 2") %in% col_labels))
})

test_that("an EFA payload serialises to JSON", {
  skip_if_not_installed("psych")
  fit <- suppressWarnings(
    psych::fa(mtcars[, c("mpg", "hp", "wt", "disp")], nfactors = 1)
  )
  json <- mellio_to_json(mellio_payload(fit))
  expect_true(nchar(as.character(json)) > 0)
})

test_that("melliotab.fa uses section with what as a compatibility alias", {
  skip_if_not_installed("psych")
  fit <- suppressWarnings(
    psych::fa(mtcars[, c("mpg", "hp", "wt", "disp")], nfactors = 1)
  )

  default_tab <- melliotab(fit)
  variance_tab <- melliotab(fit, section = "variance")
  legacy_variance_tab <- melliotab(fit, what = "variance")

  expect_s3_class(default_tab, "melliotab")
  expect_s3_class(variance_tab, "melliotab")
  expect_equal(variance_tab$data, legacy_variance_tab$data)
  expect_error(
    melliotab(fit, section = "fit", what = "variance"),
    "different table sections"
  )
})

test_that("mellio_payload.factanal builds an EFA loadings table card", {
  vars <- c("mpg", "hp", "wt", "disp", "drat", "qsec")
  fit <- suppressWarnings(factanal(mtcars[, vars], factors = 2))
  p <- mellio_payload(fit)

  expect_equal(p$card_kind, "table")
  expect_equal(p$type, "efa_loadings")
  expect_equal(p$type_label, "Exploratory factor analysis")
  expect_equal(p$fields$table_type, "factor_loadings")
  expect_equal(p$fields$n_factors, 2L)
  expect_equal(length(p$fields$rows), length(vars))

  col_keys <- vapply(p$fields$columns, function(col) col$key, character(1))
  expect_true(all(c("variable", "factor_1", "factor_2", "h2", "u2") %in% col_keys))

  # h2 = 1 - u2 for factanal
  for (row in p$fields$rows) {
    expect_equal(row$h2 + row$u2, 1, tolerance = 1e-6)
  }

  expect_equal(p$fields$method, "Maximum Likelihood")
  expect_equal(p$fields$rotation, "varimax")  # factanal's documented default
  expect_true(is.numeric(p$fields$variance_explained))
})

test_that("mellio_payload.prcomp builds a PCA loadings card without h2/u2", {
  fit <- prcomp(mtcars[, 1:6], scale. = TRUE)
  p <- mellio_payload(fit)

  expect_equal(p$card_kind, "table")
  expect_equal(p$type, "pca_loadings")
  expect_equal(p$type_label, "Principal component analysis")
  expect_equal(p$fields$table_type, "component_loadings")
  expect_equal(length(p$fields$rows), 6L)

  col_keys <- vapply(p$fields$columns, function(col) col$key, character(1))
  expect_true(all(c("variable", "component_1") %in% col_keys))
  # prcomp produces no communality/uniqueness columns
  expect_false("h2" %in% col_keys)
  expect_false("u2" %in% col_keys)

  expect_equal(p$fields$method, "Principal components")
  # All PCs retained by default -> variance_explained == 1
  expect_equal(p$fields$variance_explained, 1, tolerance = 1e-6)
})
