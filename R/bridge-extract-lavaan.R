# R bridge — lavaan extractor (CFA / SEM / mediation).
#
# Produces a `card_kind: "structural"` payload with two zones:
#
#   report_zone:
#     fit_indices = [chi2(df,p), CFI, TLI, RMSEA(+CI), SRMR, AIC, BIC]
#
#   inspection_zone:
#     parameters = full parameterEstimates(fit, standardized = TRUE, ci = TRUE)
#                  as ParameterRow[]   (lhs/op/rhs + est/se/z/p/ci/std)
#
# Type auto-detection from the model's operators:
#   supported mediation signature → lavaan_mediation
#   "~"  present → lavaan_sem         (structural paths)
#   "=~" only    → lavaan_cfa         (measurement model only)
#   anything else → lavaan_sem        (safe fallback)
#
# Requires the `lavaan` package (Suggests). For standardised estimates
# we read `std.all` (default in lavaan) so users don't have to refit.
# Schema reference: docs/STATS-R-BRIDGE-SCHEMA.md

#' @rdname mellio_payload
#' @param focal Character vector of parameter labels to mark as focal
#'   in the rendered card. Format: `"lhs op rhs"`, whitespace
#'   insensitive (e.g. `"dem60 ~ ind60"` or `"dem60~ind60"`). Marked
#'   rows surface in the report zone above the fit indices; the rest
#'   stay in the collapsible inspection zone. Only used for
#'   `card_kind: "structural"` payloads (lavaan).
#' @param diagram_omit Optional character vector of parameter keys or
#'   variable names to omit from the generated SEM path diagram while
#'   retaining them in diagram metadata. Use exact path keys like
#'   `"y ~ ideology"` or a variable name like `"ideology"` for hidden
#'   covariate-style paths.
#' @param standardized Include standardized estimates and standardized-solution
#'   metadata when available.
#' @export
mellio_payload.lavaan <- function(x, focal = NULL, diagram_omit = NULL,
                                  standardized = TRUE, ..., .call = NULL) {
  rlang::check_installed("lavaan", reason = "to extract lavaan fits")
  standardized <- isTRUE(standardized)

  # Normalise focal-path keys for comparison (strip whitespace).
  focal_keys <- if (length(focal) > 0L) {
    vapply(focal, function(f) gsub("\\s+", "", as.character(f)), character(1))
  } else character(0)

  call_str <- if (!is.null(.call)) {
    .call
  } else {
    user_call <- match.call()$x
    if (!is.null(user_call) && !identical(user_call, as.name("x"))) {
      paste(deparse(user_call, width.cutoff = 500L), collapse = " ")
    } else NA_character_
  }

  # ── Parameter estimates (inspection zone) ──────────────────────────
  pe_err <- NULL
  pe <- tryCatch(
    lavaan::parameterEstimates(x, standardized = standardized, ci = TRUE),
    error = function(e) { pe_err <<- conditionMessage(e); NULL }
  )
  if (is.null(pe)) {
    stop("Could not extract parameter estimates from lavaan fit: ",
         pe_err %||% "unknown",
         call. = FALSE)
  }

  std_solution <- if (standardized) {
    ms_lavaan_standardized_solution_map(x)
  } else {
    list()
  }
  params <- lapply(seq_len(nrow(pe)), function(i) {
    row <- list(
      lhs       = as.character(pe$lhs[i]),
      op        = as.character(pe$op[i]),
      rhs       = as.character(pe$rhs[i]),
      estimate  = ms_safe_numeric(pe$est[i])
    )
    if (!is.null(pe$se))       row$std_error    <- ms_safe_numeric(pe$se[i])
    if (!is.null(pe$z))        row$statistic    <- ms_safe_numeric(pe$z[i])
    if (!is.null(pe$pvalue))   row$p_value      <- ms_safe_numeric(pe$pvalue[i])
    if (!is.null(pe$ci.lower)) row$ci_lower     <- ms_safe_numeric(pe$ci.lower[i])
    if (!is.null(pe$ci.upper)) row$ci_upper     <- ms_safe_numeric(pe$ci.upper[i])
    if (!is.null(pe$std.all))  row$std_estimate <- ms_safe_numeric(pe$std.all[i])
    if (!is.null(pe$free))     row$free         <- ms_safe_numeric(pe$free[i])
    if (!is.null(pe$label)) {
      label <- as.character(pe$label[i])
      if (!is.na(label) && nzchar(label)) row$label <- label
    }
    if (!is.null(pe$group))    row$group        <- as.character(pe$group[i])

    std_stats <- ms_lavaan_standardized_solution_lookup(std_solution, row)
    if (length(std_stats) > 0L) {
      for (nm in names(std_stats)) {
        if (identical(nm, "std_estimate") && !is.null(row$std_estimate)) next
        row[[nm]] <- std_stats[[nm]]
      }
    }

    # Focal mark — match on whitespace-stripped "lhs op rhs" key.
    if (length(focal_keys) > 0L) {
      row_key <- gsub("\\s+", "",
                      paste0(row$lhs, row$op, row$rhs))
      if (row_key %in% focal_keys) row$is_focal <- TRUE
    }
    row
  })
  params <- ms_lavaan_attach_loading_rsquare(params, x)
  params <- ms_lavaan_mark_diagram_omissions(params, diagram_omit)

  model_type <- tryCatch(x@Options$model.type, error = function(e) NULL)
  type <- ms_lavaan_payload_type(pe, model_type)
  model_class <- ms_lavaan_model_class(pe, model_type)

  # ── Fit indices (report zone) ──────────────────────────────────────
  # Robust estimators (WLSMV, MLR, MLM, etc.) carry separate "*.scaled"
  # fit indices that adjust for the non-normality / categorical-data
  # correction; reporting the unscaled values for these estimators is
  # methodologically wrong. ms_lavaan_fit_indices prefers the scaled
  # version when fitMeasures emits one.
  estimator   <- ms_lavaan_estimator(x)
  fit_indices <- ms_lavaan_fit_indices(x, estimator = estimator)
  nobs_info <- ms_lavaan_nobs_info(x)
  n_total <- nobs_info$total
  if (is.na(n_total)) {
    n_total <- tryCatch(
      ms_safe_numeric(unname(lavaan::fitMeasures(x, "ntotal"))),
      error = function(e) NA_real_
    )
  }
  inference_meta <- ms_lavaan_inference_meta(x, pe)
  structural_r_squared <- ms_lavaan_structural_rsquare(params, x)
  observed_variables <- ms_lavaan_observed_variable_summaries(x)
  model_info <- ms_lavaan_model_info(
    x,
    pe = pe,
    params = params,
    payload_type = type,
    model_class = model_class,
    estimator = estimator,
    nobs_info = nobs_info,
    n_total = n_total
  )
  diagram <- ms_lavaan_diagram_schema(params, model_class = model_class)
  diagram_figure_type <- ms_lavaan_diagram_figure_type(model_class)
  diagram_figure_data <- ms_lavaan_diagram_figure_data(
    diagram,
    diagram_figure_type
  )
  diagram_available_figures <- ms_lavaan_diagram_available_figures(
    diagram_figure_type
  )

  # Focal paths copied into report_zone for prominent display.
  focal_paths <- Filter(function(p) isTRUE(p$is_focal), params)

  # Auto-promote defined parameters (:=) when the user hasn't marked any
  # paths as focal explicitly. The whole point of writing `indirect := a*b`
  # in a lavaan model is to surface the derived effect; burying it among
  # variances and intercepts in the inspection zone defeats that. Most
  # commonly this fires for mediation models (indirect, total,
  # prop_mediated), but it also benefits any structural model whose
  # author bothered to define derived quantities.
  if (length(focal_paths) == 0L) {
    defined_paths <- Filter(function(p) identical(p$op, ":="), params)
    if (length(defined_paths) > 0L) {
      defined_paths <- lapply(defined_paths, function(p) {
        p$is_focal <- TRUE
        p$auto_focal <- TRUE
        p
      })
      focal_paths <- defined_paths
    }
  }

  rz <- list(fit_indices = fit_indices)
  if (!is.null(estimator) && nzchar(estimator)) rz$estimator <- estimator
  if (length(focal_paths) > 0L) rz$focal_paths <- focal_paths

  # P3b: standardized loading range. lavaan emits each measurement
  # relation as a "=~" row in parameterEstimates(); std.all gives the
  # standardized loading. Summarising the range (min, max) lets the
  # paragraph close with "Standardized loadings ranged from .61 to .83"
  # without dumping the full table. Skip when there are no =~ rows
  # (path-only models) or no std.all column (rare).
  loadings <- vapply(params, function(p) {
    if (identical(p$op, "=~") &&
        !is.null(p$std_estimate) &&
        !is.na(p$std_estimate) &&
        is.finite(p$std_estimate)) {
      p$std_estimate
    } else {
      NA_real_
    }
  }, numeric(1))
  loadings <- loadings[is.finite(loadings)]
  if (length(loadings) > 0L) {
    rz$loadings_range <- list(
      min = ms_safe_numeric(min(loadings)),
      max = ms_safe_numeric(max(loadings)),
      n   = as.integer(length(loadings))
    )
  }

  reliability <- ms_lavaan_reliability(params)
  modification_indices <- ms_lavaan_modification_indices(x, top = 10L)

  inspection_zone <- list(parameters = params)
  if (length(reliability) > 0L) {
    inspection_zone$reliability <- reliability
  }
  if (length(modification_indices) > 0L) {
    inspection_zone$modification_indices <- modification_indices
  }

  fields <- list(
    report_zone     = rz,
    inspection_zone = inspection_zone,
    n               = n_total,
    model_class     = model_class,
    model_info      = model_info
  )
  if (length(diagram) > 0L) {
    fields$diagram <- diagram
  }
  if (length(structural_r_squared) > 0L) {
    fields$structural_r_squared <- structural_r_squared
  }
  if (length(observed_variables) > 0L) {
    fields$observed_variables <- observed_variables
  }
  for (nm in names(inference_meta)) {
    fields[[nm]] <- inference_meta[[nm]]
  }

  # ── Type detection ─────────────────────────────────────────────────
  # Use lavaan's own model.type when it's growth() — it preserves the
  # constructor identity. Otherwise infer from operators.
  model_type <- tryCatch(x@Options$model.type, error = function(e) NULL)
  type <- ms_lavaan_payload_type(pe, model_type = model_type)

  type_label <- switch(
    type,
    lavaan_cfa       = "Confirmatory factor analysis",
    lavaan_sem       = "Structural equation model",
    lavaan_mediation = "Mediation model",
    lavaan_growth    = "Latent growth model",
    "Structural equation model"
  )

  # Multi-group hint in the label
  ng <- tryCatch(lavaan::lavInspect(x, "ngroups"), error = function(e) 1L)
  if (!is.null(ng) && ng > 1L) {
    group_labels <- tryCatch(
      lavaan::lavInspect(x, "group.label"),
      error = function(e) NULL
    )
    type_label <- if (!is.null(group_labels) && length(group_labels) > 0) {
      paste0(type_label, " \u2014 ", ng, " groups (",
             paste(group_labels, collapse = ", "), ")")
    } else {
      paste0(type_label, " \u2014 ", ng, " groups")
    }
  }

  figure_data <- NULL
  available_figures <- NULL
  if (!identical(type, "lavaan_mediation")) {
    structural_diagram <- ms_lavaan_structural_diagram(
      params,
      type = type,
      type_label = type_label,
      n = n_total
    )
    if (!is.null(structural_diagram)) {
      figure_data <- list(structural_path_diagram = structural_diagram)
      available_figures <- list(list(
        type = "structural_path_diagram",
        label = if (identical(type, "lavaan_cfa")) "CFA path diagram" else "Path diagram",
        default = TRUE
      ))
    }
  }

  payload <- ms_build_envelope(
    type       = type,
    type_label = type_label,
    call       = trimws(gsub("\\s+", " ", call_str)),
    fields     = fields,
    raw_output = ms_lavaan_raw_output(x, standardized = standardized),
    packages   = ms_packages_basic(extras = "lavaan"),
    card_kind  = "structural",
    figure_data = figure_data,
    available_figures = available_figures
  )
  payload$metadata <- payload$metadata %||% list()
  payload$metadata$available_tables <- payload_structural_table_options(payload)
  payload
}

# ── Helpers ───────────────────────────────────────────────────────────

ms_lavaan_payload_type <- function(pe, model_type = NULL) {
  if (identical(model_type, "growth")) return("lavaan_growth")
  if (is.null(pe) || !is.data.frame(pe) || !nrow(pe) || !"op" %in% names(pe)) {
    return("lavaan_sem")
  }
  ops <- unique(as.character(pe$op))
  if (ms_lavaan_has_mediation_signature(pe)) return("lavaan_mediation")
  if ("~" %in% ops) return("lavaan_sem")
  if ("=~" %in% ops) return("lavaan_cfa")
  "lavaan_sem"
}

ms_lavaan_model_class <- function(pe, model_type = NULL) {
  type <- ms_lavaan_payload_type(pe, model_type = model_type)
  if (identical(type, "lavaan_growth")) return("growth")
  if (identical(type, "lavaan_mediation")) return("mediation")
  if (is.null(pe) || !is.data.frame(pe) || !nrow(pe) || !"op" %in% names(pe)) {
    return("observed_sem")
  }
  ops <- unique(as.character(pe$op))
  if ("=~" %in% ops && "~" %in% ops) return("latent_sem")
  if ("=~" %in% ops) return("cfa")
  if ("~" %in% ops) return("observed_sem")
  "sem"
}

ms_lavaan_has_mediation_signature <- function(pe) {
  if (is.null(pe) || !is.data.frame(pe) || !nrow(pe) || !"op" %in% names(pe)) {
    return(FALSE)
  }
  ops <- as.character(pe$op %||% character(0))
  if (!(":=" %in% ops) || !("~" %in% ops)) return(FALSE)
  if (!ms_lavaan_has_regression_chain(pe)) return(FALSE)

  defined <- pe[ops == ":=", , drop = FALSE]
  if (!nrow(defined)) return(FALSE)
  lhs <- if ("lhs" %in% names(defined)) tolower(as.character(defined$lhs)) else character(nrow(defined))
  rhs <- if ("rhs" %in% names(defined)) tolower(as.character(defined$rhs)) else character(nrow(defined))
  name_hits <- grepl(paste(c(
    "^(ind|indirect)([._-]|$)",
    "([._-])(ind|indirect)([._-]|$)",
    "serial",
    "mediat",
    "prop[._-]?med",
    "index[._-]?modmed",
    "modmed"
  ), collapse = "|"), lhs)
  formula_hits <- grepl("\\*", rhs) & grepl("[[:alpha:]]", rhs)
  any(name_hits | formula_hits, na.rm = TRUE)
}

ms_lavaan_has_regression_chain <- function(pe) {
  if (is.null(pe) || !is.data.frame(pe) || !nrow(pe) ||
      !all(c("lhs", "op", "rhs") %in% names(pe))) {
    return(FALSE)
  }
  regs <- pe[as.character(pe$op) == "~", c("lhs", "rhs"), drop = FALSE]
  if (nrow(regs) < 2L) return(FALSE)
  from <- as.character(regs$rhs)
  to <- as.character(regs$lhs)
  valid <- nzchar(from) & nzchar(to)
  from <- from[valid]
  to <- to[valid]
  if (length(from) < 2L) return(FALSE)
  for (i in seq_along(from)) {
    mediator <- to[[i]]
    next_steps <- which(from == mediator)
    if (!length(next_steps)) next
    for (j in next_steps) {
      if (!identical(from[[i]], to[[j]]) &&
          !identical(from[[i]], mediator) &&
          !identical(mediator, to[[j]])) {
        return(TRUE)
      }
    }
  }
  FALSE
}

ms_lavaan_mark_diagram_omissions <- function(params, diagram_omit = NULL) {
  if (!is.list(params) || !length(params) || !length(diagram_omit)) return(params)
  omit <- trimws(as.character(diagram_omit))
  omit <- omit[nzchar(omit) & !is.na(omit)]
  if (!length(omit)) return(params)
  omit_keys <- gsub("\\s+", "", omit)
  omit_vars <- omit[!grepl("=~|~~|~", omit)]

  lapply(params, function(p) {
    lhs <- as.character(p$lhs %||% "")
    op <- as.character(p$op %||% "")
    rhs <- as.character(p$rhs %||% "")
    key <- gsub("\\s+", "", paste0(lhs, op, rhs))
    var_hit <- length(omit_vars) > 0L &&
      (lhs %in% omit_vars || rhs %in% omit_vars)
    if (key %in% omit_keys || var_hit) {
      p$diagram_hidden <- TRUE
      p$diagram_hidden_reason <- if (var_hit && identical(op, "~") &&
                                     rhs %in% omit_vars) {
        "hidden_covariate"
      } else {
        "omitted_path"
      }
    }
    p
  })
}

ms_lavaan_structural_diagram <- function(params, type = "lavaan_sem",
                                         type_label = NULL, n = NULL) {
  if (!is.list(params) || !length(params)) return(NULL)

  group_values <- unique(vapply(params, function(p) {
    as.character(p$group %||% "")
  }, character(1)))
  group_values <- group_values[nzchar(group_values)]
  if (length(group_values) > 1L) return(NULL)

  measurement <- Filter(function(p) identical(p$op, "=~") &&
                          nzchar(as.character(p$lhs %||% "")) &&
                          nzchar(as.character(p$rhs %||% "")), params)
  structural <- Filter(function(p) identical(p$op, "~") &&
                         nzchar(as.character(p$lhs %||% "")) &&
                         nzchar(as.character(p$rhs %||% "")), params)
  covariances <- Filter(function(p) identical(p$op, "~~") &&
                          nzchar(as.character(p$lhs %||% "")) &&
                          nzchar(as.character(p$rhs %||% "")) &&
                          !identical(as.character(p$lhs), as.character(p$rhs)), params)
  hidden_params <- Filter(ms_lavaan_param_hidden_in_diagram,
                          c(measurement, structural, covariances))
  measurement <- Filter(function(p) !ms_lavaan_param_hidden_in_diagram(p),
                        measurement)
  structural <- Filter(function(p) !ms_lavaan_param_hidden_in_diagram(p),
                       structural)
  covariances <- Filter(function(p) !ms_lavaan_param_hidden_in_diagram(p),
                        covariances)
  if (!length(measurement) && !length(structural)) return(NULL)

  latent <- unique(vapply(measurement, function(p) as.character(p$lhs), character(1)))
  indicators <- unique(vapply(measurement, function(p) as.character(p$rhs), character(1)))
  endogenous <- unique(vapply(structural, function(p) as.character(p$lhs), character(1)))
  exogenous <- unique(vapply(structural, function(p) as.character(p$rhs), character(1)))
  all_vars <- unique(c(
    latent,
    indicators,
    endogenous,
    exogenous,
    vapply(covariances, function(p) as.character(p$lhs), character(1)),
    vapply(covariances, function(p) as.character(p$rhs), character(1))
  ))
  all_vars <- all_vars[nzchar(all_vars)]
  if (!length(all_vars)) return(NULL)

  nodes <- lapply(all_vars, function(name) {
    role <- if (name %in% latent) {
      "latent"
    } else if (name %in% indicators) {
      "indicator"
    } else if (name %in% endogenous) {
      "endogenous"
    } else if (name %in% exogenous) {
      "exogenous"
    } else {
      "observed"
    }
    list(
      id = ms_lavaan_diagram_id(name),
      label = name,
      variable = name,
      role = role,
      observed = !(name %in% latent)
    )
  })

  edges <- list()
  add_edges <- function(rows, edge_type) {
    for (p in rows) {
      source <- if (identical(edge_type, "structural")) p$rhs else p$lhs
      target <- if (identical(edge_type, "structural")) p$lhs else p$rhs
      edge <- list(
        id = paste0("e", length(edges) + 1L),
        source = ms_lavaan_diagram_id(source),
        target = ms_lavaan_diagram_id(target),
        source_label = as.character(source),
        target_label = as.character(target),
        type = edge_type,
        op = as.character(p$op %||% ""),
        estimate = ms_safe_numeric(p$estimate %||% NA_real_),
        std_error = ms_safe_numeric(p$std_error %||% NA_real_),
        statistic = ms_safe_numeric(p$statistic %||% NA_real_),
        p_value = ms_safe_numeric(p$p_value %||% NA_real_),
        ci_lower = ms_safe_numeric(p$ci_lower %||% NA_real_),
        ci_upper = ms_safe_numeric(p$ci_upper %||% NA_real_),
        std_estimate = ms_safe_numeric(p$std_estimate %||% NA_real_)
      )
      if (!is.null(p$free)) edge$free <- ms_safe_numeric(p$free)
      label <- as.character(p$label %||% "")
      if (nzchar(label)) edge$label <- label
      edges[[length(edges) + 1L]] <<- edge
    }
  }
  add_edges(measurement, "measurement")
  add_edges(structural, "structural")
  add_edges(covariances, "covariance")
  if (!length(edges)) return(NULL)

  omitted_paths <- ms_lavaan_omitted_diagram_paths(hidden_params)
  hidden_covariates <- ms_lavaan_hidden_covariates(hidden_params)
  out <- list(
    source = "lavaan",
    type = "structural_path_diagram",
    model_type = type,
    model_class = ms_lavaan_model_class(ms_lavaan_params_to_data_frame(params)),
    title = type_label %||% "Structural equation model",
    nodes = nodes,
    edges = edges,
    coefficient_scale = if (identical(type, "lavaan_cfa")) "standardized" else "unstandardized",
    layout = "auto"
  )
  n_value <- ms_safe_numeric(n)
  if (!is.na(n_value) && is.finite(n_value)) out$n <- n_value
  if (length(omitted_paths)) {
    out$omitted_paths <- omitted_paths
    out$omitted_path_count <- length(omitted_paths)
  }
  if (length(hidden_covariates)) {
    out$hidden_covariates <- hidden_covariates
  }
  out
}

ms_lavaan_param_hidden_in_diagram <- function(p) {
  isTRUE(p$diagram_hidden) || isTRUE(p$hidden) || isTRUE(p$omitted)
}

ms_lavaan_omitted_diagram_paths <- function(params) {
  if (!is.list(params) || !length(params)) return(list())
  out <- list()
  for (p in params) {
    op <- as.character(p$op %||% "")
    lhs <- as.character(p$lhs %||% "")
    rhs <- as.character(p$rhs %||% "")
    if (!nzchar(lhs) || !nzchar(rhs) || !(op %in% c("=~", "~", "~~"))) next
    source <- if (identical(op, "~")) rhs else lhs
    target <- if (identical(op, "~")) lhs else rhs
    edge_type <- if (identical(op, "=~")) {
      "measurement"
    } else if (identical(op, "~~")) {
      "covariance"
    } else {
      "structural"
    }
    row <- list(
      id = paste0("omitted", length(out) + 1L),
      source = ms_lavaan_diagram_id(source),
      target = ms_lavaan_diagram_id(target),
      source_label = source,
      target_label = target,
      type = edge_type,
      op = op,
      hidden = TRUE,
      reason = as.character(p$diagram_hidden_reason %||% "omitted_path"),
      estimate = ms_safe_numeric(p$estimate %||% NA_real_),
      std_error = ms_safe_numeric(p$std_error %||% NA_real_),
      statistic = ms_safe_numeric(p$statistic %||% NA_real_),
      p_value = ms_safe_numeric(p$p_value %||% NA_real_),
      ci_lower = ms_safe_numeric(p$ci_lower %||% NA_real_),
      ci_upper = ms_safe_numeric(p$ci_upper %||% NA_real_),
      std_estimate = ms_safe_numeric(p$std_estimate %||% NA_real_)
    )
    label <- as.character(p$label %||% "")
    if (nzchar(label)) row$label <- label
    out[[length(out) + 1L]] <- row
  }
  out
}

ms_lavaan_hidden_covariates <- function(params) {
  structural <- Filter(function(p) {
    identical(as.character(p$op %||% ""), "~") &&
      identical(as.character(p$diagram_hidden_reason %||% ""), "hidden_covariate")
  }, params)
  if (!length(structural)) return(list())
  vars <- unique(vapply(structural, function(p) as.character(p$rhs %||% ""),
                        character(1)))
  vars <- vars[nzchar(vars)]
  lapply(vars, function(name) {
    targets <- unique(vapply(structural, function(p) {
      if (identical(as.character(p$rhs %||% ""), name)) {
        as.character(p$lhs %||% "")
      } else {
        ""
      }
    }, character(1)))
    targets <- targets[nzchar(targets)]
    list(
      id = ms_lavaan_diagram_id(name),
      label = name,
      variable = name,
      path_count = length(targets),
      targets = as.list(targets)
    )
  })
}

ms_lavaan_params_to_data_frame <- function(params) {
  if (!is.list(params) || !length(params)) return(data.frame())
  data.frame(
    lhs = vapply(params, function(p) as.character(p$lhs %||% ""), character(1)),
    op = vapply(params, function(p) as.character(p$op %||% ""), character(1)),
    rhs = vapply(params, function(p) as.character(p$rhs %||% ""), character(1)),
    stringsAsFactors = FALSE
  )
}

ms_lavaan_diagram_id <- function(name) {
  name <- gsub("[^A-Za-z0-9_]+", "_", as.character(name %||% ""))
  name <- gsub("^_+|_+$", "", name)
  if (!nzchar(name)) name <- "node"
  name
}

ms_lavaan_raw_output <- function(x, standardized = TRUE) {
  out <- tryCatch(
    utils::capture.output(lavaan::summary(
      x,
      standardized = isTRUE(standardized),
      fit.measures = TRUE
    )),
    error = function(e) NULL
  )
  if (!is.null(out) && length(out) > 0L) {
    txt <- paste(out, collapse = "\n")
    if (!grepl("^Length\\s+Class\\s+Mode$", trimws(out[[1]]))) {
      return(txt)
    }
  }
  ms_capture_output(x)
}

ms_lavaan_standardized_solution_map <- function(fit) {
  ss <- tryCatch(
    lavaan::standardizedSolution(fit, ci = TRUE, output = "data.frame"),
    error = function(e) NULL
  )
  if (is.null(ss) || !is.data.frame(ss) || !nrow(ss)) return(list())
  if (!all(c("lhs", "op", "rhs") %in% names(ss))) return(list())

  out <- list()
  groups <- if ("group" %in% names(ss)) {
    as.character(ss$group)
  } else {
    rep("", nrow(ss))
  }
  non_empty_groups <- unique(groups[!is.na(groups) & nzchar(groups)])
  single_group <- length(non_empty_groups) <= 1L
  for (i in seq_len(nrow(ss))) {
    row <- list()
    if ("est.std" %in% names(ss)) {
      value <- ms_safe_numeric(ss$est.std[i])
      if (!is.na(value)) row$std_estimate <- value
    }
    if ("se" %in% names(ss)) {
      value <- ms_safe_numeric(ss$se[i])
      if (!is.na(value)) row$std_std_error <- value
    }
    if ("ci.lower" %in% names(ss)) {
      value <- ms_safe_numeric(ss$ci.lower[i])
      if (!is.na(value)) row$std_ci_lower <- value
    }
    if ("ci.upper" %in% names(ss)) {
      value <- ms_safe_numeric(ss$ci.upper[i])
      if (!is.na(value)) row$std_ci_upper <- value
    }
    if (!length(row)) next

    group <- groups[[i]]
    key <- ms_lavaan_parameter_key(ss$lhs[i], ss$op[i], ss$rhs[i], group)
    out[[key]] <- row
    if (isTRUE(single_group) && nzchar(as.character(group %||% ""))) {
      fallback_key <- ms_lavaan_parameter_key(ss$lhs[i], ss$op[i], ss$rhs[i], "")
      out[[fallback_key]] <- row
    }
  }
  out
}

ms_lavaan_standardized_solution_lookup <- function(map, row) {
  if (!length(map) || !is.list(row)) return(list())
  key <- ms_lavaan_parameter_key(row$lhs, row$op, row$rhs, row$group %||% "")
  if (key %in% names(map)) return(map[[key]])
  fallback <- ms_lavaan_parameter_key(row$lhs, row$op, row$rhs, "")
  if (fallback %in% names(map)) return(map[[fallback]])
  list()
}

ms_lavaan_parameter_key <- function(lhs, op, rhs, group = "") {
  paste(
    as.character(group %||% ""),
    as.character(lhs %||% ""),
    as.character(op %||% ""),
    as.character(rhs %||% ""),
    sep = "\r"
  )
}

ms_lavaan_payload_type <- function(pe, model_type = NULL) {
  ops <- unique(as.character(pe$op %||% character(0)))
  if (identical(model_type, "growth")) return("lavaan_growth")
  if (ms_lavaan_has_mediation_signature(pe)) return("lavaan_mediation")
  if ("~" %in% ops) return("lavaan_sem")
  if ("=~" %in% ops) return("lavaan_cfa")
  "lavaan_sem"
}

ms_lavaan_model_class <- function(pe, model_type = NULL) {
  ops <- unique(as.character(pe$op %||% character(0)))
  has_measurement <- "=~" %in% ops
  has_structural <- "~" %in% ops
  has_mediation <- ms_lavaan_has_mediation_signature(pe)

  if (identical(model_type, "growth")) return("growth")
  if (has_mediation) return("mediation")
  if (has_measurement && has_structural) return("latent_sem")
  if (has_structural) return("observed_sem")
  if (has_measurement) return("cfa")
  "structural_model"
}

ms_lavaan_has_mediation_signature <- function(pe) {
  if (is.null(pe) || !is.data.frame(pe) || !nrow(pe) ||
      !"op" %in% names(pe)) {
    return(FALSE)
  }
  ops <- as.character(pe$op %||% character(0))
  if (!(":=" %in% ops) || !("~" %in% ops)) return(FALSE)
  if (!ms_lavaan_has_regression_chain(pe)) return(FALSE)

  defined <- pe[ops == ":=", , drop = FALSE]
  if (!nrow(defined)) return(FALSE)

  lhs <- if ("lhs" %in% names(defined)) {
    tolower(as.character(defined$lhs))
  } else {
    character(nrow(defined))
  }
  rhs <- if ("rhs" %in% names(defined)) {
    tolower(as.character(defined$rhs))
  } else {
    character(nrow(defined))
  }

  name_hits <- grepl(
    paste(
      c(
        "^(ind|indirect)([._-]|$)",
        "([._-])(ind|indirect)([._-]|$)",
        "serial",
        "mediat",
        "prop[._-]?med",
        "index[._-]?modmed",
        "modmed"
      ),
      collapse = "|"
    ),
    lhs
  )
  formula_hits <- grepl("\\*", rhs) &
    grepl("[[:alpha:]]", rhs)
  any(name_hits | formula_hits, na.rm = TRUE)
}

ms_lavaan_has_regression_chain <- function(pe) {
  if (is.null(pe) || !is.data.frame(pe) || !nrow(pe) ||
      !all(c("lhs", "op", "rhs") %in% names(pe))) {
    return(FALSE)
  }
  regs <- pe[as.character(pe$op) == "~", c("lhs", "rhs"), drop = FALSE]
  if (nrow(regs) < 2L) return(FALSE)
  from <- as.character(regs$rhs)
  to <- as.character(regs$lhs)
  valid <- nzchar(from) & nzchar(to)
  from <- from[valid]
  to <- to[valid]
  if (length(from) < 2L) return(FALSE)
  for (i in seq_along(from)) {
    mediator <- to[[i]]
    next_steps <- which(from == mediator)
    if (!length(next_steps)) next
    for (j in next_steps) {
      if (!identical(from[[i]], to[[j]]) &&
          !identical(from[[i]], mediator) &&
          !identical(mediator, to[[j]])) {
        return(TRUE)
      }
    }
  }
  FALSE
}

ms_lavaan_model_info <- function(fit, pe, params, payload_type,
                                 model_class, estimator, nobs_info,
                                 n_total) {
  opts <- tryCatch(fit@Options, error = function(e) NULL)
  ops <- unique(as.character(pe$op %||% character(0)))
  latent <- unique(vapply(
    Filter(function(p) identical(p$op, "=~") && !is.null(p$lhs), params),
    function(p) as.character(p$lhs),
    character(1)
  ))

  out <- list(
    payload_type = payload_type,
    model_class = model_class,
    has_measurement = "=~" %in% ops,
    has_structural_paths = "~" %in% ops,
    has_defined_parameters = ":=" %in% ops
  )
  if (length(latent) > 0L) out$latent_variables <- as.list(latent)
  if (!is.null(estimator) && nzchar(estimator)) out$estimator <- estimator
  if (!is.na(n_total) && is.finite(n_total)) {
    out$n <- as.integer(round(n_total))
  }
  if (length(nobs_info$groups %||% list()) > 0L) {
    out$nobs <- nobs_info$groups
  }

  if (!is.null(opts)) {
    se <- ms_lavaan_option_label(opts$se)
    test <- ms_lavaan_option_label(opts$test)
    missing <- ms_lavaan_option_label(opts$missing)
    if (!is.na(se) && nzchar(se)) out$se <- se
    if (!is.na(test) && nzchar(test)) out$test <- test
    if (!is.na(missing) && nzchar(missing)) out$missing <- missing
    if (!is.null(opts$fixed.x)) out$fixed_x <- isTRUE(opts$fixed.x)
    if (!is.null(opts$std.lv)) out$std_lv <- isTRUE(opts$std.lv)
  }

  group_info <- ms_lavaan_group_info(fit, nobs_info)
  if (length(group_info) > 0L) out$groups <- group_info
  out
}

ms_lavaan_option_label <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NA_character_)
  vals <- unique(trimws(as.character(x)))
  vals <- vals[!is.na(vals) & nzchar(vals)]
  if (!length(vals)) return(NA_character_)
  paste(vals, collapse = ", ")
}

ms_lavaan_nobs_info <- function(fit) {
  nobs <- tryCatch(lavaan::lavInspect(fit, "nobs"),
                   error = function(e) NULL)
  if (is.null(nobs)) return(list(total = NA_real_, groups = list()))

  values <- suppressWarnings(as.numeric(nobs))
  values <- values[is.finite(values)]
  if (!length(values)) return(list(total = NA_real_, groups = list()))

  group_labels <- tryCatch(lavaan::lavInspect(fit, "group.label"),
                           error = function(e) NULL)
  names <- names(nobs)
  groups <- lapply(seq_along(values), function(i) {
    label <- NULL
    if (!is.null(group_labels) && length(group_labels) >= i &&
        nzchar(as.character(group_labels[[i]]))) {
      label <- as.character(group_labels[[i]])
    } else if (!is.null(names) && length(names) >= i &&
               nzchar(as.character(names[[i]]))) {
      label <- as.character(names[[i]])
    } else if (length(values) > 1L) {
      label <- as.character(i)
    }
    row <- list(n = as.integer(round(values[[i]])))
    if (!is.null(label) && nzchar(label)) row$group <- label
    row
  })
  list(total = sum(values), groups = groups)
}

ms_lavaan_group_info <- function(fit, nobs_info = NULL) {
  ng <- tryCatch(lavaan::lavInspect(fit, "ngroups"),
                 error = function(e) 1L)
  if (is.null(ng) || ng <= 1L) return(list())
  labels <- tryCatch(lavaan::lavInspect(fit, "group.label"),
                     error = function(e) NULL)
  nobs_groups <- (nobs_info %||% list())$groups %||% list()
  lapply(seq_len(ng), function(i) {
    label <- if (!is.null(labels) && length(labels) >= i &&
                 nzchar(as.character(labels[[i]]))) {
      as.character(labels[[i]])
    } else {
      as.character(i)
    }
    row <- list(group = label)
    if (length(nobs_groups) >= i && !is.null(nobs_groups[[i]]$n)) {
      row$n <- nobs_groups[[i]]$n
    }
    row
  })
}

ms_lavaan_diagram_schema <- function(params, model_class = NULL) {
  params <- Filter(function(p) is.list(p) && !is.null(p$op), params)
  if (!length(params)) return(list())

  latent <- unique(vapply(
    Filter(function(p) identical(p$op, "=~") && !is.null(p$lhs), params),
    function(p) as.character(p$lhs),
    character(1)
  ))
  reg_lhs <- unique(vapply(
    Filter(function(p) identical(p$op, "~") && !is.null(p$lhs), params),
    function(p) as.character(p$lhs),
    character(1)
  ))
  reg_rhs <- unique(vapply(
    Filter(function(p) identical(p$op, "~") && !is.null(p$rhs), params),
    function(p) as.character(p$rhs),
    character(1)
  ))
  indicators <- unique(vapply(
    Filter(function(p) identical(p$op, "=~") && !is.null(p$rhs), params),
    function(p) as.character(p$rhs),
    character(1)
  ))

  nodes <- list()
  node_keys <- character(0)
  edges <- list()
  defined <- list()

  node_type_for <- function(name) {
    name <- as.character(name %||% "")
    if (name %in% latent) return("latent")
    if (name %in% reg_rhs && !(name %in% reg_lhs) && !(name %in% indicators)) {
      return("covariate")
    }
    "observed"
  }

  add_node <- function(name, type = NULL, group = "") {
    name <- as.character(name %||% "")
    if (!nzchar(name)) return(NULL)
    group <- as.character(group %||% "")
    type <- type %||% node_type_for(name)
    key <- paste(type, group, name, sep = "\r")
    if (key %in% node_keys) {
      idx <- match(key, node_keys)
      return(nodes[[idx]])
    }
    node <- list(
      id = ms_lavaan_diagram_node_id(type, name, group),
      label = name,
      type = type
    )
    if (nzchar(group)) node$group <- group
    if (identical(type, "residual")) node$hidden <- TRUE
    nodes[[length(nodes) + 1L]] <<- node
    node_keys <<- c(node_keys, key)
    node
  }

  add_edge <- function(row, type, from, to, bidirectional = FALSE,
                       hidden = FALSE) {
    group <- as.character(row$group %||% "")
    from_node <- add_node(from, NULL, group)
    to_node <- add_node(to, NULL, group)
    if (is.null(from_node) || is.null(to_node)) return()
    edge <- ms_lavaan_diagram_edge(row, type, from_node$id, to_node$id,
                                   bidirectional = bidirectional,
                                   hidden = hidden,
                                   from_label = from,
                                   to_label = to)
    edges[[length(edges) + 1L]] <<- edge
  }

  for (p in params) {
    op <- as.character(p$op %||% "")
    group <- as.character(p$group %||% "")
    if (identical(op, "=~")) {
      add_node(p$lhs, "latent", group)
      add_node(p$rhs, "observed", group)
      add_edge(p, "loading", p$lhs, p$rhs)
    } else if (identical(op, "~")) {
      add_node(p$rhs, node_type_for(p$rhs), group)
      add_node(p$lhs, node_type_for(p$lhs), group)
      add_edge(p, "regression", p$rhs, p$lhs)
    } else if (identical(op, "~~")) {
      if (identical(as.character(p$lhs), as.character(p$rhs))) {
        target <- add_node(p$lhs, node_type_for(p$lhs), group)
        if (!is.null(target)) {
          residual_name <- paste0("e_", as.character(p$lhs))
          residual <- add_node(residual_name, "residual", group)
          edge <- ms_lavaan_diagram_edge(
            p, "residual", residual$id, target$id,
            bidirectional = FALSE,
            hidden = TRUE,
            from_label = residual_name,
            to_label = p$lhs
          )
          edge$subtype <- "variance"
          edges[[length(edges) + 1L]] <- edge
        }
      } else {
        add_node(p$lhs, node_type_for(p$lhs), group)
        add_node(p$rhs, node_type_for(p$rhs), group)
        add_edge(p, "covariance", p$lhs, p$rhs, bidirectional = TRUE)
      }
    } else if (identical(op, ":=")) {
      item <- ms_lavaan_diagram_defined_effect(p)
      if (length(item) > 0L) defined[[length(defined) + 1L]] <- item
    }
  }

  if (!length(nodes) || !length(edges)) return(list())
  nodes <- ms_lavaan_diagram_layout(nodes, edges, params, model_class)

  out <- list(
    schema_version = "mellio_sem_diagram_v1",
    model_class = model_class %||% "structural_model",
    nodes = nodes,
    edges = edges,
    layout = list(engine = "mellio_simple_v1", editable = TRUE),
    editor = list(
      node_overrides = list(),
      edge_overrides = list(),
      hidden_nodes = list(),
      hidden_edges = list(),
      label_overrides = list()
    )
  )
  if (length(defined) > 0L) out$defined_effects <- defined
  out
}

ms_lavaan_diagram_figure_type <- function(model_class) {
  model_class <- as.character(model_class %||% "")
  if (identical(model_class, "cfa")) return("cfa_path_diagram")
  if (model_class %in% c("observed_sem", "latent_sem")) return("sem_path_diagram")
  NULL
}

ms_lavaan_diagram_figure_data <- function(diagram, figure_type) {
  if (!length(diagram) || is.null(figure_type) || !nzchar(figure_type)) {
    return(NULL)
  }
  out <- list()
  out[[figure_type]] <- diagram
  out
}

ms_lavaan_diagram_available_figures <- function(figure_type) {
  if (is.null(figure_type) || !nzchar(figure_type)) return(NULL)
  label <- if (identical(figure_type, "cfa_path_diagram")) {
    "CFA path diagram"
  } else {
    "SEM path diagram"
  }
  list(list(type = figure_type, label = label, default = TRUE))
}

ms_lavaan_diagram_node_id <- function(type, name, group = "") {
  parts <- c("node", type, ms_lavaan_id_component(name))
  group <- as.character(group %||% "")
  if (nzchar(group)) parts <- c(parts, paste0("g", ms_lavaan_id_component(group)))
  paste(parts, collapse = ":")
}

ms_lavaan_diagram_edge_id <- function(type, row, from, to) {
  group <- as.character(row$group %||% "")
  parts <- c(
    "edge",
    type,
    ms_lavaan_id_component(from),
    ms_lavaan_id_component(to)
  )
  if (nzchar(group)) parts <- c(parts, paste0("g", ms_lavaan_id_component(group)))
  paste(parts, collapse = ":")
}

ms_lavaan_id_component <- function(value) {
  value <- as.character(value %||% "")
  if (!nzchar(value)) return("empty")
  utils::URLencode(value, reserved = TRUE)
}

ms_lavaan_diagram_edge <- function(row, type, from_id, to_id,
                                   bidirectional = FALSE, hidden = FALSE,
                                   from_label = NULL, to_label = NULL) {
  from_label <- from_label %||% row$rhs %||% row$lhs
  to_label <- to_label %||% row$lhs %||% row$rhs
  edge <- list(
    id = ms_lavaan_diagram_edge_id(type, row, from_label, to_label),
    type = type,
    from = from_id,
    to = to_id,
    op = as.character(row$op %||% "")
  )
  if (isTRUE(bidirectional)) edge$bidirectional <- TRUE
  if (isTRUE(hidden)) edge$hidden <- TRUE
  if (!is.null(row$lhs)) edge$lhs <- as.character(row$lhs)
  if (!is.null(row$rhs)) edge$rhs <- as.character(row$rhs)
  if (!is.null(row$group) && nzchar(as.character(row$group))) {
    edge$group <- as.character(row$group)
  }
  if (!is.null(row$estimate)) edge$estimate <- ms_safe_numeric(row$estimate)
  if (!is.null(row$free)) edge$free <- ms_safe_numeric(row$free)
  if (!is.null(row$std_estimate)) {
    edge$std_estimate <- ms_safe_numeric(row$std_estimate)
  }
  if (!is.null(row$std_std_error)) {
    edge$std_std_error <- ms_safe_numeric(row$std_std_error)
  }
  if (!is.null(row$std_error)) edge$std_error <- ms_safe_numeric(row$std_error)
  if (!is.null(row$p_value)) {
    edge$p_value <- ms_safe_numeric(row$p_value)
    stars <- ms_lavaan_p_stars(edge$p_value)
    if (nzchar(stars)) edge$stars <- stars
  }
  if (!is.null(row$ci_lower) && !is.null(row$ci_upper)) {
    lo <- ms_safe_numeric(row$ci_lower)
    hi <- ms_safe_numeric(row$ci_upper)
    if (!is.na(lo) && !is.na(hi)) edge$ci <- I(c(lo, hi))
  }
  if (!is.null(row$std_ci_lower) && !is.null(row$std_ci_upper)) {
    lo <- ms_safe_numeric(row$std_ci_lower)
    hi <- ms_safe_numeric(row$std_ci_upper)
    if (!is.na(lo) && !is.na(hi)) edge$std_ci <- I(c(lo, hi))
  }
  label <- as.character(row$label %||% "")
  if (nzchar(label)) edge$parameter_label <- label
  edge
}

ms_lavaan_diagram_defined_effect <- function(row) {
  if (is.null(row$lhs)) return(list())
  item <- list(
    id = paste("defined", ms_lavaan_id_component(row$lhs), sep = ":"),
    type = "defined",
    label = as.character(row$lhs),
    op = ":="
  )
  if (!is.null(row$rhs) && nzchar(as.character(row$rhs))) {
    item$formula <- as.character(row$rhs)
  }
  if (!is.null(row$estimate)) item$estimate <- ms_safe_numeric(row$estimate)
  if (!is.null(row$std_error)) item$std_error <- ms_safe_numeric(row$std_error)
  if (!is.null(row$p_value)) {
    item$p_value <- ms_safe_numeric(row$p_value)
    stars <- ms_lavaan_p_stars(item$p_value)
    if (nzchar(stars)) item$stars <- stars
  }
  if (!is.null(row$ci_lower) && !is.null(row$ci_upper)) {
    lo <- ms_safe_numeric(row$ci_lower)
    hi <- ms_safe_numeric(row$ci_upper)
    if (!is.na(lo) && !is.na(hi)) item$ci <- I(c(lo, hi))
  }
  item
}

ms_lavaan_p_stars <- function(p) {
  p <- ms_safe_numeric(p)
  if (is.na(p)) return("")
  if (p < .001) return("***")
  if (p < .01) return("**")
  if (p < .05) return("*")
  ""
}

ms_lavaan_diagram_layout <- function(nodes, edges, params, model_class = NULL) {
  if (!length(nodes)) return(nodes)
  if (identical(model_class, "cfa")) {
    return(ms_lavaan_cfa_diagram_layout(nodes, params))
  }

  ids <- vapply(nodes, function(n) n$id, character(1))
  level <- setNames(rep(0L, length(ids)), ids)
  for (iter in seq_along(ids)) {
    changed <- FALSE
    for (edge in edges) {
      if (!identical(edge$type, "regression") &&
          !identical(edge$type, "loading")) next
      from <- as.character(edge$from %||% "")
      to <- as.character(edge$to %||% "")
      if (!from %in% names(level) || !to %in% names(level)) next
      next_level <- level[[from]] + 1L
      if (next_level > level[[to]]) {
        level[[to]] <- next_level
        changed <- TRUE
      }
    }
    if (!changed) break
  }

  split_ids <- split(ids, level[ids])
  y_step <- 92
  x_step <- 220
  for (lvl in names(split_ids)) {
    group_ids <- split_ids[[lvl]]
    group_ids <- group_ids[order(vapply(group_ids, function(id) {
      node <- nodes[[match(id, ids)]]
      paste(node$type, node$label, sep = "\r")
    }, character(1)))]
    for (i in seq_along(group_ids)) {
      idx <- match(group_ids[[i]], ids)
      y0 <- 80 + (i - 1L) * y_step
      if (identical(nodes[[idx]]$type, "residual")) y0 <- y0 + 34
      nodes[[idx]]$position <- list(
        x = 120 + as.integer(lvl) * x_step,
        y = y0
      )
    }
  }
  nodes
}

ms_lavaan_cfa_diagram_layout <- function(nodes, params) {
  matching_node_indices <- function(type, label) {
    which(vapply(nodes, function(n) {
      identical(n$type, type) && identical(as.character(n$label), as.character(label))
    }, logical(1)))
  }
  latent_names <- unique(vapply(
    Filter(function(p) identical(p$op, "=~") && !is.null(p$lhs), params),
    function(p) as.character(p$lhs),
    character(1)
  ))
  for (f_idx in seq_along(latent_names)) {
    factor <- latent_names[[f_idx]]
    factor_y <- 90 + (f_idx - 1L) * 180
    for (node_idx in matching_node_indices("latent", factor)) {
      nodes[[node_idx]]$position <- list(x = 120, y = factor_y)
    }
    indicators <- vapply(
      Filter(function(p) identical(p$op, "=~") &&
               identical(as.character(p$lhs), factor) &&
               !is.null(p$rhs), params),
      function(p) as.character(p$rhs),
      character(1)
    )
    if (!length(indicators)) next
    offsets <- (seq_along(indicators) - (length(indicators) + 1) / 2) * 54
    for (i in seq_along(indicators)) {
      for (ind_idx in matching_node_indices("observed", indicators[[i]])) {
        nodes[[ind_idx]]$position <- list(x = 390, y = factor_y + offsets[[i]])
      }
    }
  }
  for (i in seq_along(nodes)) {
    if (!is.null(nodes[[i]]$position)) next
    nodes[[i]]$position <- list(x = 650, y = 80 + (i - 1L) * 48)
  }
  nodes
}

ms_lavaan_attach_loading_rsquare <- function(params, fit) {
  if (length(params) == 0L) return(params)
  rsq <- ms_lavaan_rsquare_map(fit)

  loading_keys <- vapply(params, function(p) {
    if (!identical(p$op, "=~") || is.null(p$rhs)) return(NA_character_)
    ms_lavaan_rsquare_key(p$rhs, p$group %||% "")
  }, character(1))
  loading_counts <- table(loading_keys[!is.na(loading_keys)])

  lapply(params, function(p) {
    if (!identical(p$op, "=~") || is.null(p$rhs)) return(p)
    value <- ms_lavaan_rsquare_lookup(rsq, p$rhs, p$group %||% "")

    # Fallback for simple one-factor indicators and older lavaan objects
    # where lavInspect("rsquare") is unavailable.
    key <- ms_lavaan_rsquare_key(p$rhs, p$group %||% "")
    if (is.na(value) &&
        key %in% names(loading_counts) &&
        identical(as.integer(loading_counts[[key]]), 1L) &&
        !is.null(p$std_estimate) &&
        is.finite(p$std_estimate)) {
      value <- as.numeric(p$std_estimate)^2
    }

    if (!is.na(value) && is.finite(value)) {
      p$r_squared <- ms_safe_numeric(value)
    }
    p
  })
}

ms_lavaan_rsquare_map <- function(fit) {
  rsq <- tryCatch(
    lavaan::lavInspect(fit, "rsquare"),
    error = function(e) NULL
  )
  out <- list()
  add_values <- function(values, group = "") {
    if (is.null(values)) return()
    if (is.matrix(values) || is.data.frame(values)) values <- as.vector(values)
    if (!is.numeric(values)) return()
    names <- names(values)
    if (is.null(names)) return()
    for (i in seq_along(values)) {
      name <- names[[i]]
      value <- ms_safe_numeric(values[[i]])
      if (!nzchar(name) || is.na(value) || !is.finite(value)) next
      out[[ms_lavaan_rsquare_key(name, group)]] <<- value
    }
  }

  if (is.list(rsq) && !is.data.frame(rsq) && !is.numeric(rsq)) {
    for (i in seq_along(rsq)) {
      add_values(rsq[[i]], as.character(i))
    }
  } else {
    add_values(rsq, "")
  }
  out
}

ms_lavaan_rsquare_key <- function(name, group = "") {
  paste0(as.character(group %||% ""), "\r", as.character(name %||% ""))
}

ms_lavaan_rsquare_lookup <- function(map, name, group = "") {
  if (!length(map) || is.null(name)) return(NA_real_)
  key <- ms_lavaan_rsquare_key(name, group)
  if (key %in% names(map)) return(map[[key]])
  fallback <- ms_lavaan_rsquare_key(name, "")
  if (fallback %in% names(map)) return(map[[fallback]])
  NA_real_
}

ms_lavaan_structural_rsquare <- function(params, fit) {
  if (length(params) == 0L) return(list())
  rsq <- ms_lavaan_rsquare_map(fit)
  if (!length(rsq)) return(list())

  seen <- character(0)
  out <- list()
  for (p in params) {
    if (!identical(p$op, "~") || is.null(p$lhs)) next
    group <- p$group %||% ""
    key <- ms_lavaan_rsquare_key(p$lhs, group)
    if (key %in% seen) next
    value <- ms_lavaan_rsquare_lookup(rsq, p$lhs, group)
    if (is.na(value) || !is.finite(value)) next
    seen <- c(seen, key)
    row <- list(
      variable = as.character(p$lhs),
      r_squared = ms_safe_numeric(value)
    )
    if (!is.null(p$group) && nzchar(as.character(p$group))) {
      row$group <- as.character(p$group)
    }
    out[[length(out) + 1L]] <- row
  }
  out
}

ms_lavaan_observed_variable_summaries <- function(fit) {
  out <- list()
  seen <- character(0)

  add_row <- function(variable, values = NULL, mean_value = NA_real_,
                      group = "") {
    variable <- as.character(variable %||% "")
    if (!nzchar(variable)) return()
    group <- as.character(group %||% "")
    key <- paste0(group, "\r", variable)
    if (key %in% seen) return()

    if (!is.null(values)) {
      values <- suppressWarnings(as.numeric(values))
      values <- values[is.finite(values)]
      if (!length(values)) return()
      mean_value <- mean(values)
      row <- list(
        variable = variable,
        mean = ms_safe_numeric(mean_value),
        n = as.integer(length(values))
      )
      if (length(values) > 1L) {
        row$sd <- ms_safe_numeric(stats::sd(values))
      }
    } else {
      mean_value <- ms_safe_numeric(mean_value)
      if (is.na(mean_value) || !is.finite(mean_value)) return()
      row <- list(
        variable = variable,
        mean = mean_value
      )
    }

    if (nzchar(group)) row$group <- group
    seen <<- c(seen, key)
    out[[length(out) + 1L]] <<- row
  }

  add_data_frame <- function(dat, group = "") {
    if (is.null(dat)) return()
    dat <- tryCatch(as.data.frame(dat, optional = TRUE),
                    error = function(e) NULL)
    if (is.null(dat) || !length(dat)) return()
    for (nm in names(dat)) {
      col <- dat[[nm]]
      if (!is.numeric(col) && !is.integer(col) && !is.logical(col)) next
      add_row(nm, values = col, group = group)
    }
  }

  group_labels <- tryCatch(lavaan::lavInspect(fit, "group.label"),
                           error = function(e) NULL)
  group_label <- function(i) {
    if (!is.null(group_labels) && length(group_labels) >= i &&
        nzchar(as.character(group_labels[[i]]))) {
      return(as.character(group_labels[[i]]))
    }
    if (length(group_labels) > 0L) return(as.character(i))
    ""
  }

  raw_data <- tryCatch(lavaan::lavInspect(fit, "data"),
                       error = function(e) NULL)
  if (is.list(raw_data) && !is.data.frame(raw_data)) {
    for (i in seq_along(raw_data)) add_data_frame(raw_data[[i]], group_label(i))
  } else {
    add_data_frame(raw_data, "")
  }
  if (length(out) > 0L) return(out)

  add_sample_means <- function(stat, group = "") {
    if (is.null(stat) || !is.list(stat)) return()
    means <- stat$mean %||% stat$mean.ov %||% stat$mean.x
    if (is.null(means) || !is.numeric(means)) return()
    nms <- names(means)
    if (is.null(nms)) return()
    for (i in seq_along(means)) {
      add_row(nms[[i]], mean_value = means[[i]], group = group)
    }
  }

  sampstat <- tryCatch(lavaan::lavInspect(fit, "sampstat"),
                       error = function(e) NULL)
  if (is.list(sampstat) && length(sampstat) > 0L &&
      !is.null(sampstat[[1L]]) && is.list(sampstat[[1L]]) &&
      is.null(sampstat$mean)) {
    for (i in seq_along(sampstat)) add_sample_means(sampstat[[i]], group_label(i))
  } else {
    add_sample_means(sampstat, "")
  }
  out
}

ms_lavaan_reliability <- function(params) {
  value_or_empty <- function(p, key) {
    value <- p[[key]]
    if (is.null(value) || length(value) == 0L) return("")
    value <- value[[1L]]
    if (is.null(value) || is.na(value)) "" else as.character(value)
  }
  param_group <- function(p) value_or_empty(p, "group")
  param_key <- function(group, name) paste0(group, "\r", name)
  is_free_parameter <- function(p) {
    free <- ms_safe_numeric(p$free)
    !is.finite(free) || free != 0
  }

  loadings <- Filter(function(p) {
    identical(p$op, "=~") &&
      !is.null(p$lhs) &&
      !is.null(p$std_estimate) &&
      is.finite(p$std_estimate)
  }, params)
  if (length(loadings) == 0L) return(list())

  residual_variances <- list()
  for (p in params) {
    if (!identical(p$op, "~~")) next
    lhs <- value_or_empty(p, "lhs")
    rhs <- value_or_empty(p, "rhs")
    if (!nzchar(lhs) || !identical(lhs, rhs)) next
    theta <- ms_safe_numeric(p$std_estimate)
    if (!is.finite(theta)) next
    residual_variances[[param_key(param_group(p), lhs)]] <- max(0, as.numeric(theta))
  }

  factor_index <- list()
  for (p in loadings) {
    factor <- value_or_empty(p, "lhs")
    group <- param_group(p)
    key <- param_key(group, factor)
    if (is.null(factor_index[[key]])) {
      factor_index[[key]] <- list(factor = factor, group = group)
    }
  }

  out <- lapply(factor_index, function(entry) {
    factor <- entry$factor
    group <- entry$group
    rows <- Filter(function(p) {
      identical(value_or_empty(p, "lhs"), factor) &&
        identical(param_group(p), group)
    }, loadings)
    rows <- Filter(function(p) {
      lambda <- ms_safe_numeric(p$std_estimate)
      is.finite(lambda) && nzchar(value_or_empty(p, "rhs"))
    }, rows)
    lambdas <- vapply(rows, function(p) as.numeric(p$std_estimate), numeric(1))
    if (length(lambdas) < 2L) return(NULL)
    indicators <- vapply(rows, function(p) value_or_empty(p, "rhs"), character(1))
    squared <- lambdas^2
    uniqueness <- vapply(seq_along(lambdas), function(i) {
      indicator_key <- param_key(group, indicators[[i]])
      theta <- residual_variances[[indicator_key]]
      if (is.null(theta) || !is.finite(theta)) theta <- 1 - squared[[i]]
      max(0, as.numeric(theta))
    }, numeric(1))
    names(uniqueness) <- indicators

    residual_covariance_sum <- 0
    residual_covariance_count <- 0L
    for (p in params) {
      if (!identical(p$op, "~~") || !is_free_parameter(p)) next
      if (!identical(param_group(p), group)) next
      lhs <- value_or_empty(p, "lhs")
      rhs <- value_or_empty(p, "rhs")
      if (!nzchar(lhs) || !nzchar(rhs) || identical(lhs, rhs)) next
      if (!(lhs %in% indicators) || !(rhs %in% indicators)) next
      r_theta <- ms_safe_numeric(p$std_estimate)
      if (!is.finite(r_theta)) next
      theta_lhs <- uniqueness[[lhs]]
      theta_rhs <- uniqueness[[rhs]]
      residual_covariance_sum <- residual_covariance_sum +
        as.numeric(r_theta) * sqrt(max(0, theta_lhs) * max(0, theta_rhs))
      residual_covariance_count <- residual_covariance_count + 1L
    }

    numerator <- sum(abs(lambdas))^2
    denominator <- numerator + sum(uniqueness) + 2 * residual_covariance_sum
    omega <- if (is.finite(denominator) && denominator > 0) numerator / denominator else NA_real_
    ave <- mean(squared)
    row <- list(
      factor = factor,
      omega = ms_safe_numeric(omega),
      ave = ms_safe_numeric(ave),
      n_indicators = as.integer(length(lambdas))
    )
    if (nzchar(group)) row$group <- group
    if (residual_covariance_count > 0L) {
      row$residual_covariances <- as.integer(residual_covariance_count)
      row$omega_residual_covariance_adjusted <- TRUE
    }
    row
  })
  Filter(Negate(is.null), out)
}

ms_lavaan_modification_indices <- function(fit, top = 10L) {
  mi <- tryCatch(
    lavaan::modindices(fit, sort. = TRUE),
    error = function(e) NULL
  )
  if (is.null(mi) || !nrow(mi) || !"mi" %in% names(mi)) return(list())
  mi <- mi[is.finite(mi$mi), , drop = FALSE]
  if (!nrow(mi)) return(list())
  mi <- utils::head(mi, max(1L, as.integer(top)))
  lapply(seq_len(nrow(mi)), function(i) {
    row <- list(
      lhs = as.character(mi$lhs[i]),
      op = as.character(mi$op[i]),
      rhs = as.character(mi$rhs[i]),
      mi = ms_safe_numeric(mi$mi[i])
    )
    if ("epc" %in% names(mi)) {
      row$epc <- ms_safe_numeric(mi$epc[i])
    }
    if ("sepc.all" %in% names(mi)) {
      row$std_epc <- ms_safe_numeric(mi$sepc.all[i])
    }
    row
  })
}

ms_lavaan_inference_meta <- function(fit, pe = NULL) {
  opts <- tryCatch(fit@Options, error = function(e) NULL)
  if (is.null(opts)) return(list())

  se <- tolower(as.character(opts$se %||% ""))
  boot <- any(grepl("boot", se, fixed = TRUE))
  out <- list(boot = isTRUE(boot))

  sims <- ms_safe_numeric(opts$bootstrap %||% NA_real_)
  if (isTRUE(boot) && !is.na(sims) && is.finite(sims) && sims > 0) {
    out$sims <- as.integer(round(sims))
  }

  ci_type <- ms_chr_or_na(opts$bootstrap.ci.type %||% opts$boot.ci.type)
  if (isTRUE(boot) && !is.na(ci_type) && nzchar(ci_type)) {
    out$ci_type <- ci_type
  }

  level <- ms_safe_numeric(opts$level %||% opts$conf.level %||% opts$conf_level)
  if (is.na(level) && !is.null(pe) &&
      all(c("ci.lower", "ci.upper") %in% names(pe))) {
    level <- 0.95
  }
  if (!is.na(level) && is.finite(level) && level > 0) {
    out$conf_level <- level
  }

  out
}

# Identify the user-facing estimator label (ML / MLM / MLMV / MLR /
# WLSMV / etc.). lavaan stores the user's "MLR" as a combination —
# `@Options$estimator = "ML"` plus `@Options$test = c("standard",
# "yuan.bentler.mplus")` — so we have to look at the test correction
# to recover the original label.
ms_lavaan_estimator <- function(fit) {
  opts <- tryCatch(fit@Options, error = function(e) NULL)
  if (is.null(opts)) return(NULL)
  base <- toupper(trimws(as.character(opts$estimator %||% "")))
  tests <- as.character(opts$test %||% character(0))
  correction <- tests[tolower(tests) != "standard"]
  if (length(correction) == 0L) {
    return(if (nzchar(base)) base else NULL)
  }
  variant <- tolower(correction[[1L]])
  suffix <- switch(
    variant,
    "yuan.bentler"       = "R",
    "yuan.bentler.mplus" = "R",
    "satorra.bentler"    = "M",
    "scaled.shifted"     = "MV",
    "mean.var.adjusted"  = "MVS",
    NULL
  )
  # Yuan-Bentler is always reported as "MLR" regardless of base — it
  # only applies to ML. Other suffixes glue onto whatever base lavaan
  # selected (ML/DWLS/ULS).
  label <- if (identical(suffix, "R")) {
    "MLR"
  } else if (!is.null(suffix) && nzchar(base)) {
    paste0(base, suffix)
  } else if (nzchar(base)) {
    base
  } else {
    NULL
  }
  label
}

# Pull the canonical APA-style fit-indices set from lavaan::fitMeasures.
# Returns a list of FitIndex objects matching the schema. When lavaan
# emits "*.scaled" keys (Satorra-Bentler / Yuan-Bentler / WLSMV-style
# corrections), those values take precedence — reporting the uncorrected
# chi-square / RMSEA for a robust fit would misrepresent the fit.
ms_lavaan_fit_indices <- function(fit, estimator = NULL) {
  fm <- tryCatch(lavaan::fitMeasures(fit), error = function(e) NULL)
  if (is.null(fm)) return(list())

  get_num <- function(key) {
    if (key %in% names(fm)) ms_safe_numeric(fm[[key]]) else NA_real_
  }

  # Presence of `chisq.scaled` in fitMeasures is the source of truth
  # for whether lavaan applied a scaling correction (covers MLM, MLMV,
  # MLR, WLSM, WLSMV, etc.). The estimator string argument is kept for
  # callers that want it but no longer drives the gate.
  scaled <- "chisq.scaled" %in% names(fm)
  pick <- function(unscaled_key, scaled_key) {
    if (scaled && scaled_key %in% names(fm)) {
      v <- ms_safe_numeric(fm[[scaled_key]])
      if (!is.na(v)) return(list(value = v, scaled = TRUE))
    }
    list(value = get_num(unscaled_key), scaled = FALSE)
  }

  out <- list()

  # chi-square test of model fit
  chi <- pick("chisq", "chisq.scaled")
  if (!is.na(chi$value)) {
    df_choice <- if (chi$scaled) pick("df", "df.scaled") else list(value = get_num("df"))
    p_choice  <- if (chi$scaled) pick("pvalue", "pvalue.scaled") else list(value = get_num("pvalue"))
    entry <- list(
      name    = "chi2",
      value   = chi$value,
      df      = df_choice$value,
      p_value = p_choice$value
    )
    if (chi$scaled) entry$scaled <- TRUE
    out[[length(out) + 1L]] <- entry
  }

  # CFI, TLI — comparative fit
  for (k in c("cfi", "tli")) {
    ch <- pick(k, paste0(k, ".scaled"))
    if (!is.na(ch$value)) {
      entry <- list(
        name      = toupper(k),
        value     = ch$value,
        threshold = ms_threshold_higher(ch$value, good = 0.95, acceptable = 0.90)
      )
      if (ch$scaled) entry$scaled <- TRUE
      out[[length(out) + 1L]] <- entry
    }
  }

  # RMSEA with 90% CI
  rmsea_ch <- pick("rmsea", "rmsea.scaled")
  if (!is.na(rmsea_ch$value)) {
    ci_lo <- if (rmsea_ch$scaled) get_num("rmsea.ci.lower.scaled") else get_num("rmsea.ci.lower")
    ci_hi <- if (rmsea_ch$scaled) get_num("rmsea.ci.upper.scaled") else get_num("rmsea.ci.upper")
    # Fall back to unscaled CI bounds when the .scaled keys are absent
    # but the .scaled point estimate is present (rare but possible).
    if (rmsea_ch$scaled && (is.na(ci_lo) || is.na(ci_hi))) {
      ci_lo <- get_num("rmsea.ci.lower")
      ci_hi <- get_num("rmsea.ci.upper")
    }
    entry <- list(
      name      = "RMSEA",
      value     = rmsea_ch$value,
      ci_level  = 0.90,
      threshold = ms_lavaan_rmsea_threshold(rmsea_ch$value, ci_hi)
    )
    if (!is.na(ci_lo) && !is.na(ci_hi)) {
      entry$ci <- I(c(ci_lo, ci_hi))
    }
    if (rmsea_ch$scaled) entry$scaled <- TRUE
    out[[length(out) + 1L]] <- entry
  }

  # SRMR — lavaan uses "srmr_bentler" for the WLSMV/robust variants and
  # "srmr" for the standard ML path. Prefer the Bentler form when scaled.
  srmr_v <- NA_real_
  srmr_scaled <- FALSE
  if (scaled && "srmr_bentler" %in% names(fm)) {
    srmr_v <- ms_safe_numeric(fm[["srmr_bentler"]])
    srmr_scaled <- !is.na(srmr_v)
  }
  if (is.na(srmr_v)) srmr_v <- get_num("srmr")
  if (!is.na(srmr_v)) {
    entry <- list(
      name      = "SRMR",
      value     = srmr_v,
      threshold = ms_threshold_lower(srmr_v, good = 0.05, acceptable = 0.08)
    )
    if (srmr_scaled) entry$scaled <- TRUE
    out[[length(out) + 1L]] <- entry
  }

  # AIC, BIC — information criteria (no threshold). No scaled variants.
  for (k in c("aic", "bic")) {
    v <- get_num(k)
    if (!is.na(v)) {
      out[[length(out) + 1L]] <- list(
        name  = toupper(k),
        value = v
      )
    }
  }

  out
}

ms_lavaan_rmsea_threshold <- function(value, ci_upper = NA_real_) {
  threshold <- ms_threshold_lower(value, good = 0.05, acceptable = 0.08)
  if (identical(threshold, "good") &&
      !is.na(ci_upper) && is.finite(ci_upper) && ci_upper > 0.08) {
    return("acceptable")
  }
  threshold
}

# Threshold helpers — bucket a fit index value into good/acceptable/poor.
# For "higher is better" indices (CFI, TLI): >= good → good, >= acceptable → acceptable, else poor
ms_threshold_higher <- function(v, good, acceptable) {
  if (is.na(v)) return(NULL)
  if (v >= good)       "good"
  else if (v >= acceptable) "acceptable"
  else                      "poor"
}

# For "lower is better" indices (RMSEA, SRMR): <= good → good, <= acceptable → acceptable, else poor
ms_threshold_lower <- function(v, good, acceptable) {
  if (is.na(v)) return(NULL)
  if (v <= good)       "good"
  else if (v <= acceptable) "acceptable"
  else                      "poor"
}
