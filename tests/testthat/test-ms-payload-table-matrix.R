test_that("mellio_payload.table emits rows and columns for one- and two-way tables", {
  one_way <- mellio_payload(table(mtcars$cyl))
  two_way <- mellio_payload(table(mtcars$cyl, mtcars$am))

  expect_equal(one_way$type, "frequency_table")
  expect_equal(one_way$card_kind, "table")
  expect_true(length(one_way$fields$rows) > 0L)
  expect_true(length(one_way$fields$columns) > 0L)
  expect_equal(one_way$fields$columns[[1]]$label, "Category")
  expect_equal(one_way$fields$total, nrow(mtcars))

  expect_equal(two_way$type, "frequency_table")
  expect_equal(two_way$card_kind, "table")
  expect_true(length(two_way$fields$rows) > 0L)
  expect_true(length(two_way$fields$columns) > 0L)
  expect_equal(two_way$fields$n_dimensions, 2L)
})

test_that("mellio_payload.matrix detects correlation matrices", {
  p <- mellio_payload(cor(mtcars[, c("mpg", "wt", "hp")]))

  expect_equal(p$type, "correlation_matrix")
  expect_equal(p$card_kind, "table")
  expect_equal(p$fields$table_type, "correlation_matrix")
  expect_equal(length(p$fields$rows), 3L)
  expect_equal(length(p$fields$columns), 4L)
  expect_equal(p$fields$columns[[1]]$key, "Variable")
  expect_equal(p$metadata$available_figures[[1]]$type, "correlation_heatmap")
  expect_true(isTRUE(p$metadata$available_figures[[1]]$default))
  expect_equal(p$metadata$available_figures[[2]]$type, "correlation_forest")
  expect_false(isTRUE(p$metadata$available_figures[[2]]$default))

  fig <- p$figure_data$correlation_heatmap
  expect_equal(fig$variables, c("mpg", "wt", "hp"))
  expect_equal(length(fig$r), 3L)
  expect_equal(fig$r[[1]][[1]], 1)
  expect_equal(fig$default_triangle, "lower")

  forest <- p$figure_data$correlation_forest
  expect_equal(forest$variables, c("mpg", "wt", "hp"))
  expect_equal(length(forest$pairs), 3L)
  expect_equal(forest$pairs[[1]]$x, "mpg")
  expect_equal(forest$pairs[[1]]$y, "wt")
})

test_that("mellio_payload.matrix delegates generic matrices to data-frame payloads", {
  p <- mellio_payload(matrix(1:9, nrow = 3))

  expect_s3_class(p, "mellio_payload")
  expect_false(identical(p$type, "correlation_matrix"))
})

test_that("new table payloads use supported card kinds and include rows and columns", {
  supported_kinds <- c("inline", "structural", "table", "raw_text", "multi", "unsupported")
  payloads <- list(
    mellio_payload(cor(mtcars[, c("mpg", "wt", "hp")])),
    mellio_payload(table(mtcars$cyl)),
    mellio_payload(table(mtcars$cyl, mtcars$am)),
    mellio_payload(melliotab(mtcars[1:3, 1:3]))
  )

  for (p in payloads) {
    expect_true(p$card_kind %in% supported_kinds, info = p$type)
    if (identical(p$card_kind, "table")) {
      expect_true(length(p$fields$rows) > 0L, info = paste("rows missing for", p$type))
      expect_true(length(p$fields$columns) > 0L, info = paste("columns missing for", p$type))
    }
  }
})
