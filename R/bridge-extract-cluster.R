# R bridge -- clustering and classification extractors.

#' @rdname mellio_payload
#' @export
mellio_payload.kmeans <- function(x, ..., .call = NULL) {
  centers <- as.matrix(x$centers %||% matrix(nrow = 0L, ncol = 0L))
  if (!nrow(centers)) {
    stop("kmeans object does not contain cluster centers.", call. = FALSE)
  }

  variable_names <- colnames(centers) %||% paste0("V", seq_len(ncol(centers)))
  center_keys <- ms_safe_keys(variable_names, "center")
  sizes <- x$size %||% rep(NA_integer_, nrow(centers))

  rows <- lapply(seq_len(nrow(centers)), function(i) {
    row <- list(
      cluster = ms_matrix_row_name(centers, i, "cluster"),
      size = as.integer(sizes[[i]] %||% NA_integer_)
    )
    for (j in seq_len(ncol(centers))) {
      row[[center_keys[[j]]]] <- ms_safe_numeric(centers[i, j])
    }
    row
  })

  columns <- list(
    list(key = "cluster", label = "Cluster", format = "text"),
    list(key = "size", label = "n", format = "integer")
  )
  for (j in seq_along(variable_names)) {
    columns <- c(columns, list(list(
      key = center_keys[[j]],
      label = variable_names[[j]],
      format = "number"
    )))
  }

  fields <- list(
    table_type = "cluster_centers",
    method = "k-means",
    k = nrow(centers),
    columns = columns,
    rows = rows,
    total_withinss = ms_safe_numeric(x$tot.withinss %||% NA_real_),
    betweenss = ms_safe_numeric(x$betweenss %||% NA_real_),
    totss = ms_safe_numeric(x$totss %||% NA_real_),
    iterations = as.integer(x$iter %||% NA_integer_),
    ifault = as.integer(x$ifault %||% NA_integer_),
    source = "stats::kmeans"
  )
  n <- length(x$cluster %||% integer(0))
  if (n > 0L) fields$n <- as.integer(n)

  ms_build_envelope(
    type = "cluster_summary",
    type_label = "Cluster summary",
    call = ms_cluster_call(x, .call),
    fields = fields,
    raw_output = ms_capture_output(x),
    packages = ms_packages_basic(),
    card_kind = "table"
  )
}

#' @rdname mellio_payload
#' @export
mellio_payload.hclust <- function(x, ..., .call = NULL) {
  merge <- as.matrix(x$merge %||% matrix(nrow = 0L, ncol = 2L))
  if (!nrow(merge)) {
    stop("hclust object does not contain merge steps.", call. = FALSE)
  }

  max_rows <- 500L
  keep <- seq_len(min(nrow(merge), max_rows))
  heights <- x$height %||% rep(NA_real_, nrow(merge))
  labels <- x$labels %||% character(0)

  rows <- lapply(keep, function(i) {
    list(
      step = as.integer(i),
      left = ms_hclust_node_label(merge[i, 1], labels),
      right = ms_hclust_node_label(merge[i, 2], labels),
      height = ms_safe_numeric(heights[[i]] %||% NA_real_)
    )
  })

  fields <- list(
    table_type = "cluster_hierarchy",
    method = as.character(x$method %||% "hierarchical clustering"),
    dist_method = as.character(x$dist.method %||% ""),
    n_observations = if (length(labels)) length(labels) else nrow(merge) + 1L,
    n_merges = nrow(merge),
    truncated = nrow(merge) > max_rows,
    columns = list(
      list(key = "step", label = "Step", format = "integer"),
      list(key = "left", label = "Left branch", format = "text"),
      list(key = "right", label = "Right branch", format = "text"),
      list(key = "height", label = "Height", format = "number")
    ),
    rows = rows,
    source = "stats::hclust"
  )

  ms_build_envelope(
    type = "cluster_summary",
    type_label = "Cluster summary",
    call = ms_cluster_call(x, .call),
    fields = fields,
    raw_output = ms_capture_output(x),
    packages = ms_packages_basic(),
    card_kind = "table"
  )
}

#' @rdname mellio_payload
#' @export
mellio_payload.randomForest <- function(x, ..., .call = NULL) {
  if (!is.null(x$confusion)) {
    return(ms_random_forest_confusion_payload(x, .call = .call))
  }
  importance <- x$importance %||% NULL
  if (!is.null(importance)) {
    return(ms_random_forest_importance_payload(x, .call = .call))
  }
  stop("randomForest object does not contain a confusion matrix or importance table.",
       call. = FALSE)
}

ms_random_forest_confusion_payload <- function(x, .call = NULL) {
  confusion <- as.matrix(x$confusion)
  actual <- rownames(confusion) %||% paste0("class_", seq_len(nrow(confusion)))
  class_error_col <- grep("^class\\.error$|^class_error$", colnames(confusion), ignore.case = TRUE)
  predicted_cols <- setdiff(seq_len(ncol(confusion)), class_error_col)
  predicted_names <- colnames(confusion)[predicted_cols] %||% paste0("class_", predicted_cols)
  predicted_keys <- ms_safe_keys(predicted_names, "predicted")

  rows <- lapply(seq_len(nrow(confusion)), function(i) {
    row <- list(actual = actual[[i]])
    for (j in seq_along(predicted_cols)) {
      row[[predicted_keys[[j]]]] <- ms_safe_numeric(confusion[i, predicted_cols[[j]]])
    }
    if (length(class_error_col)) {
      row$class_error <- ms_safe_numeric(confusion[i, class_error_col[[1]]])
    }
    row
  })

  columns <- list(list(key = "actual", label = "Actual class", format = "text"))
  for (j in seq_along(predicted_names)) {
    columns <- c(columns, list(list(
      key = predicted_keys[[j]],
      label = paste("Predicted", predicted_names[[j]]),
      format = "integer"
    )))
  }
  if (length(class_error_col)) {
    columns <- c(columns, list(list(key = "class_error", label = "Class error", format = "number")))
  }

  ms_build_envelope(
    type = "classification_model",
    type_label = "Classification model",
    call = ms_cluster_call(x, .call),
    fields = c(
      ms_random_forest_common_fields(x),
      list(
        table_type = "classification_confusion_matrix",
        columns = columns,
        rows = rows,
        classes = I(as.character(actual)),
        source = "randomForest::randomForest"
      )
    ),
    raw_output = ms_capture_output(x),
    packages = ms_packages_basic(extras = "randomForest"),
    card_kind = "table"
  )
}

ms_random_forest_importance_payload <- function(x, .call = NULL) {
  importance <- as.matrix(x$importance)
  vars <- rownames(importance) %||% paste0("variable_", seq_len(nrow(importance)))
  metric_names <- colnames(importance) %||% paste0("metric_", seq_len(ncol(importance)))
  metric_keys <- ms_safe_keys(metric_names, "metric")

  rows <- lapply(seq_len(nrow(importance)), function(i) {
    row <- list(variable = vars[[i]])
    for (j in seq_len(ncol(importance))) {
      row[[metric_keys[[j]]]] <- ms_safe_numeric(importance[i, j])
    }
    row
  })

  columns <- list(list(key = "variable", label = "Variable", format = "text"))
  for (j in seq_along(metric_names)) {
    columns <- c(columns, list(list(
      key = metric_keys[[j]],
      label = metric_names[[j]],
      format = "number"
    )))
  }

  model_type <- tolower(as.character(x$type %||% ""))
  payload_type <- if (identical(model_type, "regression")) "random_forest_model" else "classification_model"
  payload_label <- if (identical(model_type, "regression")) "Random forest model" else "Classification model"

  ms_build_envelope(
    type = payload_type,
    type_label = payload_label,
    call = ms_cluster_call(x, .call),
    fields = c(
      ms_random_forest_common_fields(x),
      list(
        table_type = "variable_importance",
        columns = columns,
        rows = rows,
        source = "randomForest::randomForest"
      )
    ),
    raw_output = ms_capture_output(x),
    packages = ms_packages_basic(extras = "randomForest"),
    card_kind = "table"
  )
}

ms_random_forest_common_fields <- function(x) {
  out <- list(
    method = "random forest",
    problem_type = as.character(x$type %||% ""),
    ntree = as.integer(x$ntree %||% NA_integer_),
    mtry = as.integer(x$mtry %||% NA_integer_)
  )
  err <- x$err.rate %||% NULL
  if (!is.null(err)) {
    err <- as.matrix(err)
    if (nrow(err) && ncol(err)) {
      col <- which(tolower(colnames(err) %||% "") == "oob")
      if (!length(col)) col <- 1L
      out$oob_error <- unname(ms_safe_numeric(err[nrow(err), col[[1]]]))
    }
  }
  out
}

ms_hclust_node_label <- function(value, labels) {
  value <- as.integer(value)
  if (is.na(value)) return(NA_character_)
  if (value < 0L) {
    idx <- abs(value)
    if (length(labels) >= idx && nzchar(labels[[idx]])) return(labels[[idx]])
    return(paste0("Observation ", idx))
  }
  paste0("Merge ", value)
}

ms_cluster_call <- function(x, .call = NULL) {
  if (!is.null(.call)) return(trimws(gsub("\\s+", " ", .call)))
  call_obj <- tryCatch(stats::getCall(x), error = function(e) NULL)
  if (!is.null(call_obj)) return(ms_deparse_call(call_obj))
  call_obj <- x$call %||% NULL
  if (!is.null(call_obj)) return(ms_deparse_call(call_obj))
  NA_character_
}

ms_safe_keys <- function(labels, prefix) {
  labels <- as.character(labels %||% character(0))
  keys <- tolower(gsub("[^A-Za-z0-9]+", "_", labels))
  keys <- gsub("^_+|_+$", "", keys)
  keys[!nzchar(keys)] <- paste0(prefix, "_", which(!nzchar(keys)))
  duplicated_any <- duplicated(keys) | duplicated(keys, fromLast = TRUE)
  if (any(duplicated_any)) {
    keys[duplicated_any] <- paste0(keys[duplicated_any], "_", seq_along(keys)[duplicated_any])
  }
  keys
}

ms_matrix_row_name <- function(x, i, prefix) {
  rn <- rownames(x)
  if (!is.null(rn) && length(rn) >= i && nzchar(rn[[i]])) return(as.character(rn[[i]]))
  paste0(prefix, "_", i)
}
