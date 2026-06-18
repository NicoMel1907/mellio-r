# R bridge -- psych::fa() / psych::principal() loadings extractor.
#
# Both produce a rotated loadings matrix, so each maps to a
# `card_kind: "table"` payload: one row per measured variable, one
# column per factor/component, plus communality (h2) and uniqueness
# (u2). Rotation, factoring method (fa only), N, the factor count,
# and total variance explained ride along in `fields`.
#
# Requires the psych package (Suggests).

#' @rdname mellio_payload
#' @export
mellio_payload.fa <- function(x, ..., .call = NULL) {
  ms_efa_payload(x, .call = .call, kind = "fa")
}

#' @rdname mellio_payload
#' @export
mellio_payload.principal <- function(x, ..., .call = NULL) {
  ms_efa_payload(x, .call = .call, kind = "principal")
}

#' @rdname mellio_payload
#' @export
mellio_payload.factanal <- function(x, ..., .call = NULL) {
  loadings <- tryCatch(as.matrix(unclass(x$loadings)),
                       error = function(e) NULL)
  if (is.null(loadings) || nrow(loadings) == 0L || ncol(loadings) == 0L) {
    stop("This factanal object has no loadings matrix to report.",
         call. = FALSE)
  }
  uniqueness  <- x$uniquenesses
  # factanal exposes uniquenesses; communality = 1 - uniquenesses.
  communality <- if (!is.null(uniqueness)) 1 - uniqueness else NULL

  ms_efa_build_table_card(
    loadings    = loadings,
    communality = communality,
    uniqueness  = uniqueness,
    unit_label  = "Factor",
    key_prefix  = "factor_",
    type        = "efa_loadings",
    type_label  = "Exploratory factor analysis",
    table_type  = "factor_loadings",
    n_obs       = ms_safe_numeric(x$n.obs),
    rotation    = ms_factanal_rotation(x),
    method      = "Maximum Likelihood",
    variance    = ms_efa_variance_from_loadings(loadings),
    call_str    = ms_psych_call(x, .call, "factanal(...)"),
    raw_output  = ms_capture_output(x),
    packages    = ms_packages_basic()
  )
}

#' @rdname mellio_payload
#' @export
mellio_payload.prcomp <- function(x, ..., .call = NULL) {
  # prcomp stores loadings in $rotation (the eigenvector matrix);
  # the name overloads the EFA "rotation" concept (prcomp does not
  # rotate post-hoc).
  loadings <- tryCatch(as.matrix(x$rotation), error = function(e) NULL)
  if (is.null(loadings) || nrow(loadings) == 0L || ncol(loadings) == 0L) {
    stop("This prcomp object has no rotation/loadings matrix to report.",
         call. = FALSE)
  }
  n_obs    <- if (!is.null(x$x)) ms_safe_numeric(nrow(x$x)) else NA_real_
  variance <- ms_prcomp_variance_explained(x, n_retained = ncol(loadings))

  ms_efa_build_table_card(
    loadings    = loadings,
    communality = NULL,
    uniqueness  = NULL,
    unit_label  = "Component",
    key_prefix  = "component_",
    type        = "pca_loadings",
    type_label  = "Principal component analysis",
    table_type  = "component_loadings",
    n_obs       = n_obs,
    rotation    = NA_character_,
    method      = "Principal components",
    variance    = variance,
    call_str    = ms_psych_call(x, .call, "prcomp(...)"),
    raw_output  = ms_capture_output(x),
    packages    = ms_packages_basic()
  )
}

# Shared worker. psych::fa and psych::principal share a loadings
# structure; only the payload type, the user-facing labels, and the
# loadings column header ("Factor" vs "Component") differ.
ms_efa_payload <- function(x, .call = NULL, kind = c("fa", "principal")) {
  kind <- match.arg(kind)
  rlang::check_installed("psych", reason = "to extract psych factor analyses")

  cfg <- if (kind == "principal") {
    list(unit = "Component", key = "component_",
         table_type = "component_loadings", type = "pca_loadings",
         label = "Principal component analysis",
         fallback = "psych::principal(...)")
  } else {
    list(unit = "Factor", key = "factor_",
         table_type = "factor_loadings", type = "efa_loadings",
         label = "Exploratory factor analysis",
         fallback = "psych::fa(...)")
  }

  loadings <- tryCatch(as.matrix(unclass(x$loadings)),
                       error = function(e) NULL)
  if (is.null(loadings) || nrow(loadings) == 0L || ncol(loadings) == 0L) {
    stop("This psych object has no loadings matrix to report.",
         call. = FALSE)
  }

  n_factors  <- ncol(loadings)
  item_names <- rownames(loadings)
  if (is.null(item_names)) {
    item_names <- paste0("V", seq_len(nrow(loadings)))
  }
  factor_keys   <- paste0(cfg$key, seq_len(n_factors))
  factor_labels <- paste(cfg$unit, seq_len(n_factors))

  comm <- if (!is.null(x$communality)) x$communality[item_names] else NULL
  uniq <- if (!is.null(x$uniquenesses)) x$uniquenesses[item_names] else NULL
  has_h2 <- !is.null(comm) && any(is.finite(comm))
  has_u2 <- !is.null(uniq) && any(is.finite(uniq))

  # ── Columns: Variable | <unit> 1..k | h2 | u2 ──────────────────────
  columns <- list(list(key = "variable", label = "Variable", format = "text"))
  for (j in seq_len(n_factors)) {
    columns[[length(columns) + 1L]] <- list(
      key = factor_keys[[j]], label = factor_labels[[j]], format = "bounded"
    )
  }
  if (has_h2) {
    columns[[length(columns) + 1L]] <-
      list(key = "h2", label = "h\u00b2", format = "bounded")
  }
  if (has_u2) {
    columns[[length(columns) + 1L]] <-
      list(key = "u2", label = "u\u00b2", format = "bounded")
  }

  # ── Rows: one per measured variable ────────────────────────────────
  rows <- lapply(seq_len(nrow(loadings)), function(i) {
    row <- list(variable = item_names[[i]])
    for (j in seq_len(n_factors)) {
      row[[factor_keys[[j]]]] <- ms_safe_numeric(loadings[i, j])
    }
    if (has_h2) row$h2 <- ms_safe_numeric(comm[[i]])
    if (has_u2) row$u2 <- ms_safe_numeric(uniq[[i]])
    row
  })

  fields <- list(
    table_type = cfg$table_type,
    columns    = columns,
    rows       = rows,
    n_factors  = as.integer(n_factors)
  )

  rotation <- x$rotation %||% NA_character_
  if (length(rotation) == 1L && !is.na(rotation) && nzchar(rotation)) {
    fields$rotation <- rotation
  }

  fm <- x$fm %||% NA_character_
  if (length(fm) == 1L && !is.na(fm) && nzchar(fm)) {
    fields$method <- ms_efa_method_label(fm)
  }

  n_obs <- ms_safe_numeric(x$n.obs)
  if (!is.na(n_obs)) fields$n <- as.integer(n_obs)

  variance <- ms_efa_variance_explained(x)
  if (!is.na(variance)) fields$variance_explained <- variance

  ms_build_envelope(
    type       = cfg$type,
    type_label = cfg$label,
    call       = trimws(gsub("\\s+", " ",
                             ms_psych_call(x, .call, cfg$fallback))),
    fields     = fields,
    raw_output = ms_capture_output(x),
    packages   = ms_packages_basic("psych"),
    card_kind  = "table"
  )
}

# ── Helpers ───────────────────────────────────────────────────────────

# Map psych's terse factoring-method codes to readable labels. Mirrors
# the labels melliotab.fa() uses so the Tables and Stats paths agree.
ms_efa_method_label <- function(fm) {
  switch(
    fm,
    pa     = "Principal Axis",
    minres = "Minimum Residual",
    ml     = "Maximum Likelihood",
    wls    = "Weighted Least Squares",
    gls    = "Generalized Least Squares",
    fm
  )
}

# Total variance explained by the retained factors, read from the
# psych object's Vaccounted matrix. Orthogonal solutions expose a
# "Cumulative Var" row; oblique solutions may only carry "Proportion
# Var", so sum that as a fallback. Returns NA when neither is present.
ms_efa_variance_explained <- function(x) {
  vacc <- x$Vaccounted
  if (is.null(vacc) || !is.matrix(vacc) || ncol(vacc) == 0L) {
    return(NA_real_)
  }
  rn <- rownames(vacc) %||% character(0)
  if ("Cumulative Var" %in% rn) {
    return(ms_safe_numeric(vacc["Cumulative Var", ncol(vacc)]))
  }
  if ("Proportion Var" %in% rn) {
    return(ms_safe_numeric(sum(vacc["Proportion Var", ], na.rm = TRUE)))
  }
  NA_real_
}

# Shared table-card builder for base R loadings producers (factanal,
# prcomp). Builds the columns / rows / envelope from already-extracted
# loadings, communality, uniqueness, and metadata.
ms_efa_build_table_card <- function(loadings, communality, uniqueness,
                                    unit_label, key_prefix,
                                    type, type_label, table_type,
                                    n_obs, rotation, method, variance,
                                    call_str, raw_output, packages) {
  n_factors  <- ncol(loadings)
  item_names <- rownames(loadings)
  if (is.null(item_names)) {
    item_names <- paste0("V", seq_len(nrow(loadings)))
  }
  factor_keys   <- paste0(key_prefix, seq_len(n_factors))
  factor_labels <- paste(unit_label, seq_len(n_factors))

  align <- function(v) {
    if (is.null(v)) return(NULL)
    if (!is.null(names(v)) && all(item_names %in% names(v))) v[item_names] else v
  }
  comm <- align(communality)
  uniq <- align(uniqueness)
  has_h2 <- !is.null(comm) && any(is.finite(comm))
  has_u2 <- !is.null(uniq) && any(is.finite(uniq))

  columns <- list(list(key = "variable", label = "Variable", format = "text"))
  for (j in seq_len(n_factors)) {
    columns[[length(columns) + 1L]] <- list(
      key = factor_keys[[j]], label = factor_labels[[j]], format = "bounded"
    )
  }
  if (has_h2) {
    columns[[length(columns) + 1L]] <-
      list(key = "h2", label = "h\u00b2", format = "bounded")
  }
  if (has_u2) {
    columns[[length(columns) + 1L]] <-
      list(key = "u2", label = "u\u00b2", format = "bounded")
  }

  rows <- lapply(seq_len(nrow(loadings)), function(i) {
    row <- list(variable = item_names[[i]])
    for (j in seq_len(n_factors)) {
      row[[factor_keys[[j]]]] <- ms_safe_numeric(loadings[i, j])
    }
    if (has_h2) row$h2 <- ms_safe_numeric(comm[[i]])
    if (has_u2) row$u2 <- ms_safe_numeric(uniq[[i]])
    row
  })

  fields <- list(
    table_type = table_type,
    columns    = columns,
    rows       = rows,
    n_factors  = as.integer(n_factors)
  )
  if (length(rotation) == 1L && !is.na(rotation) && nzchar(rotation)) {
    fields$rotation <- rotation
  }
  if (!is.null(method) && length(method) == 1L && nzchar(method)) {
    fields$method <- method
  }
  if (!is.na(n_obs)) fields$n <- as.integer(n_obs)
  if (!is.null(variance) && !is.na(variance)) {
    fields$variance_explained <- variance
  }

  ms_build_envelope(
    type       = type,
    type_label = type_label,
    call       = trimws(gsub("\\s+", " ", call_str)),
    fields     = fields,
    raw_output = raw_output,
    packages   = packages,
    card_kind  = "table"
  )
}

# Extract the rotation method from a factanal object. Reads the user's
# call; factanal's rotation defaults to "varimax" when omitted.
ms_factanal_rotation <- function(x) {
  rotation <- tryCatch({
    call_obj <- x$call
    if (!is.null(call_obj) && "rotation" %in% names(call_obj)) {
      as.character(call_obj$rotation)
    } else {
      "varimax"
    }
  }, error = function(e) NA_character_)
  if (length(rotation) != 1L) NA_character_ else rotation
}

# Variance explained when no Vaccounted matrix exists: sum of squared
# loadings across retained factors divided by the number of items.
# Equals psych's "Cumulative Var" for the last factor in orthogonal
# solutions.
ms_efa_variance_from_loadings <- function(loadings) {
  if (is.null(loadings) || nrow(loadings) == 0L || ncol(loadings) == 0L) {
    return(NA_real_)
  }
  ss <- colSums(loadings^2, na.rm = TRUE)
  ms_safe_numeric(sum(ss) / nrow(loadings))
}

# Variance explained by the retained PCs in a prcomp object.
# var of each PC = sdev^2; retained = the first n_retained; total =
# the full sdev^2 sum.
ms_prcomp_variance_explained <- function(x, n_retained) {
  if (is.null(x$sdev) || length(x$sdev) == 0L || n_retained < 1L) {
    return(NA_real_)
  }
  vars <- x$sdev^2
  n_retained <- min(n_retained, length(vars))
  total <- sum(vars, na.rm = TRUE)
  if (!is.finite(total) || total <= 0) return(NA_real_)
  ms_safe_numeric(sum(vars[seq_len(n_retained)], na.rm = TRUE) / total)
}
