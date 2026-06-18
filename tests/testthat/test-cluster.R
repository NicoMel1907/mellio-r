test_that("kmeans payload produces cluster-center rows", {
  set.seed(42)
  fit <- kmeans(iris[, 1:4], centers = 3)
  p <- mellio_payload(fit)

  expect_equal(p$card_kind, "table")
  expect_equal(p$type, "cluster_summary")
  expect_equal(p$fields$table_type, "cluster_centers")
  expect_equal(p$fields$k, 3L)
  expect_equal(length(p$fields$rows), 3L)
  expect_true(is.numeric(p$fields$rows[[1]]$size))
})

test_that("hclust payload produces hierarchy rows", {
  fit <- hclust(dist(USArrests))
  p <- mellio_payload(fit)

  expect_equal(p$type, "cluster_summary")
  expect_equal(p$fields$table_type, "cluster_hierarchy")
  expect_equal(p$fields$n_observations, nrow(USArrests))
  expect_equal(length(p$fields$rows), nrow(USArrests) - 1L)
  expect_true(is.numeric(p$fields$rows[[1]]$height))
})

test_that("randomForest-like payload extracts confusion matrix", {
  confusion <- matrix(
    c(10, 2, 0.1667,
      1, 12, 0.0769),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      c("setosa", "versicolor"),
      c("setosa", "versicolor", "class.error")
    )
  )
  err <- matrix(0.12, nrow = 1, dimnames = list(NULL, "OOB"))
  fit <- structure(
    list(
      confusion = confusion,
      type = "classification",
      ntree = 100L,
      mtry = 2L,
      err.rate = err
    ),
    class = "randomForest"
  )

  p <- mellio_payload(fit)

  expect_equal(p$type, "classification_model")
  expect_equal(p$fields$table_type, "classification_confusion_matrix")
  expect_equal(p$fields$ntree, 100L)
  expect_equal(p$fields$oob_error, 0.12)
  expect_equal(length(p$fields$rows), 2L)
  expect_true("class_error" %in% names(p$fields$rows[[1]]))
})
