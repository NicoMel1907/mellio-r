# R bridge + table extractors for bruceR::PROCESS objects.
#
# bruceR::PROCESS returns a plain list rather than a distinct S3 class.
# Detection therefore has to be conservative: only lists with the PROCESS
# metadata and the expected result blocks are handled here.

ms_is_bruce_process <- function(x) {
  is.list(x) &&
    !is.null(x$process.id) &&
    length(x$process.id) > 0L &&
    !is.null(x$process.type) &&
    length(x$process.type) > 0L &&
    is.character(x$process.type) &&
    is.list(x$results)
}

ms_bruce_process_payload <- function(x, .call = NULL) {
  process_id <- ms_safe_numeric(ms_bruce_first(x$process.id, NA_real_))
  if (identical(process_id, 4) || isTRUE(process_id == 4)) {
    payload <- ms_bruce_process_model4_payload(x, .call = .call)
    if (!is.null(payload)) return(payload)
  }
  if (identical(process_id, 6) || isTRUE(process_id == 6)) {
    payload <- ms_bruce_process_model6_payload(x, .call = .call)
    if (!is.null(payload)) return(payload)
  }
  if (identical(process_id, 7) || isTRUE(process_id == 7)) {
    payload <- ms_bruce_process_model7_payload(x, .call = .call)
    if (!is.null(payload)) return(payload)
  }
  if (identical(process_id, 8) || isTRUE(process_id == 8)) {
    payload <- ms_bruce_process_model8_payload(x, .call = .call)
    if (!is.null(payload)) return(payload)
  }
  if (identical(process_id, 14) || isTRUE(process_id == 14)) {
    payload <- ms_bruce_process_model14_payload(x, .call = .call)
    if (!is.null(payload)) return(payload)
  }
  if (identical(process_id, 15) || isTRUE(process_id == 15)) {
    payload <- ms_bruce_process_model15_payload(x, .call = .call)
    if (!is.null(payload)) return(payload)
  }
  ms_bruce_process_unsupported_payload(x, .call = .call)
}

ms_bruce_process_model4_payload <- function(x, .call = NULL) {
  model_m_list <- ms_bruce_process_model_m_list(x)
  if (length(model_m_list) > 1L) {
    return(ms_bruce_process_parallel_payload(x, .call = .call))
  }

  med_df <- ms_bruce_process_mediation_df(x)
  model_m <- ms_bruce_process_model_m(x)
  model_y <- x$model.y
  if (!is.data.frame(med_df) ||
      !inherits(model_m, "lm") ||
      !inherits(model_y, "lm")) {
    return(NULL)
  }

  vars <- ms_bruce_process_model4_vars(x, med_df, model_m, model_y)
  rows <- ms_bruce_process_mediation_rows(med_df)
  if (!length(rows)) return(NULL)

  ci_meta <- ms_bruce_process_ci_meta(x, med_df)
  paths <- ms_bruce_process_model4_paths(model_m, model_y, vars)
  path_rows <- ms_bruce_process_model4_path_rows(model_m, model_y, vars)
  n <- ms_safe_numeric(stats::nobs(model_y))

  fields <- list(
    table_type = "mediation",
    source = "bruceR::PROCESS",
    process_id = as.integer(ms_safe_numeric(ms_bruce_first(x$process.id, NA_real_))),
    process_type = as.character(ms_bruce_first(x$process.type, "")),
    rows = rows,
    columns = ms_bruce_process_effect_columns(),
    note = ms_bruce_process_summary_note(x, n),
    effect_note = ms_bruce_process_effect_note(x, n),
    path_note = ms_bruce_process_path_note(model_m, model_y, n),
    treatment = vars$predictor,
    mediator = vars$mediator,
    outcome = vars$outcome,
    controls = I(vars$covariates),
    n = as.integer(round(n)),
    boot = ci_meta$boot,
    ci_type = ci_meta$ci_type,
    conf_level = ci_meta$conf_level,
    sims = ci_meta$sims
  )
  if (!is.null(paths) && length(paths)) fields$paths <- paths
  if (!is.null(path_rows) && length(path_rows)) {
    fields$path_rows <- path_rows
    fields$path_columns <- ms_bruce_process_path_columns()
  }

  ms_build_envelope(
    type = "mediation_mediate",
    type_label = "Mediation analysis (bruceR PROCESS Model 4)",
    call = ms_bruce_process_call(.call, "bruceR::PROCESS(..., med.type = \"boot\")"),
    fields = fields,
    raw_output = ms_capture_output(x),
    packages = ms_packages_basic(extras = c("bruceR", "mediation")),
    card_kind = "table"
  )
}

ms_bruce_process_parallel_payload <- function(x, .call = NULL) {
  med_dfs <- ms_bruce_process_mediation_dfs(x)
  model_m <- ms_bruce_process_model_m_list(x)
  model_y <- x$model.y
  if (!length(med_dfs) ||
      length(model_m) < 2L ||
      !inherits(model_y, "lm") ||
      !all(vapply(model_m, inherits, logical(1), what = "lm"))) {
    return(NULL)
  }

  vars <- ms_bruce_process_parallel_vars(x, med_dfs, model_m, model_y)
  if (is.na(vars$predictor) || !nzchar(vars$predictor) ||
      length(vars$mediators) < 2L ||
      any(is.na(vars$mediators)) ||
      any(!nzchar(vars$mediators)) ||
      is.na(vars$outcome) || !nzchar(vars$outcome)) {
    return(NULL)
  }

  rows <- ms_bruce_process_parallel_rows(med_dfs, vars$mediators)
  if (!length(rows)) return(NULL)

  ci_meta <- ms_bruce_process_ci_meta(x, med_dfs[[1L]])
  paths <- ms_bruce_process_parallel_paths(model_m, model_y, vars)
  path_rows <- ms_bruce_process_parallel_path_rows(model_m, model_y, vars)
  n <- ms_safe_numeric(stats::nobs(model_y))

  fields <- list(
    table_type = "mediation",
    source = "bruceR::PROCESS",
    process_id = as.integer(ms_safe_numeric(ms_bruce_first(x$process.id, NA_real_))),
    process_type = as.character(ms_bruce_first(x$process.type, "")),
    topology = "parallel",
    hayes_model = "4",
    rows = rows,
    columns = ms_bruce_process_effect_columns(),
    note = ms_bruce_process_parallel_summary_note(x, n),
    effect_note = ms_bruce_process_parallel_effect_note(x, n),
    path_note = ms_bruce_process_model6_path_note(model_m, model_y, n),
    treatment = vars$predictor,
    mediator = vars$mediators[[1L]],
    mediators = I(vars$mediators),
    outcome = vars$outcome,
    controls = I(vars$covariates),
    n = as.integer(round(n)),
    boot = ci_meta$boot,
    ci_type = ci_meta$ci_type,
    conf_level = ci_meta$conf_level,
    sims = ci_meta$sims
  )
  if (!is.null(paths) && length(paths)) fields$paths <- paths
  if (!is.null(path_rows) && length(path_rows)) {
    fields$path_rows <- path_rows
    fields$path_columns <- ms_bruce_process_path_columns()
  }

  ms_build_envelope(
    type = "mediation_mediate",
    type_label = "Parallel mediation analysis (bruceR PROCESS Model 4)",
    call = ms_bruce_process_call(.call, "bruceR::PROCESS(..., med.type = \"parallel\")"),
    fields = fields,
    raw_output = ms_capture_output(x),
    packages = ms_packages_basic(extras = c("bruceR", "mediation")),
    card_kind = "table"
  )
}

ms_bruce_process_model6_payload <- function(x, .call = NULL) {
  med_df <- ms_bruce_process_lavaan_mediation_df(x)
  model_m <- ms_bruce_process_model_m_list(x)
  model_y <- x$model.y
  if (!is.data.frame(med_df) ||
      length(model_m) != 2L ||
      !inherits(model_y, "lm") ||
      !all(vapply(model_m, inherits, logical(1), what = "lm"))) {
    return(NULL)
  }

  vars <- ms_bruce_process_model6_vars(x, med_df, model_m, model_y)
  if (is.na(vars$predictor) || !nzchar(vars$predictor) ||
      length(vars$mediators) != 2L ||
      any(is.na(vars$mediators)) ||
      any(!nzchar(vars$mediators)) ||
      is.na(vars$outcome) || !nzchar(vars$outcome)) {
    return(NULL)
  }

  rows <- ms_bruce_process_model6_rows(med_df, vars$mediators)
  if (!length(rows)) return(NULL)

  ci_meta <- ms_bruce_process_ci_meta(x, med_df)
  paths <- ms_bruce_process_model6_paths(model_m, model_y, vars)
  path_rows <- ms_bruce_process_model6_path_rows(model_m, model_y, vars)
  n <- ms_safe_numeric(stats::nobs(model_y))

  fields <- list(
    table_type = "mediation",
    source = "bruceR::PROCESS",
    process_id = as.integer(ms_safe_numeric(ms_bruce_first(x$process.id, NA_real_))),
    process_type = as.character(ms_bruce_first(x$process.type, "")),
    topology = "serial",
    hayes_model = "6",
    rows = rows,
    columns = ms_bruce_process_effect_columns(),
    note = ms_bruce_process_model6_summary_note(x, n),
    effect_note = ms_bruce_process_model6_effect_note(x, n),
    path_note = ms_bruce_process_model6_path_note(model_m, model_y, n),
    treatment = vars$predictor,
    mediator = vars$mediators[[1L]],
    mediators = I(vars$mediators),
    outcome = vars$outcome,
    controls = I(vars$covariates),
    serial_edges = I(list(list(from = vars$mediators[[1L]],
                               to = vars$mediators[[2L]]))),
    n = as.integer(round(n)),
    boot = ci_meta$boot,
    ci_type = ci_meta$ci_type,
    conf_level = ci_meta$conf_level,
    sims = ci_meta$sims
  )
  if (!is.null(paths) && length(paths)) fields$paths <- paths
  if (!is.null(path_rows) && length(path_rows)) {
    fields$path_rows <- path_rows
    fields$path_columns <- ms_bruce_process_path_columns()
  }

  ms_build_envelope(
    type = "mediation_mediate",
    type_label = "Serial mediation analysis (bruceR PROCESS Model 6)",
    call = ms_bruce_process_call(.call, "bruceR::PROCESS(..., med.type = \"serial\")"),
    fields = fields,
    raw_output = ms_capture_output(x),
    packages = ms_packages_basic(extras = c("bruceR", "lavaan")),
    card_kind = "table"
  )
}

ms_bruce_process_model7_payload <- function(x, .call = NULL) {
  med_df <- ms_bruce_process_mediation_df(x)
  model_m <- ms_bruce_process_model_m(x)
  model_y <- x$model.y
  if (!is.data.frame(med_df) ||
      !inherits(model_m, "lm") ||
      !inherits(model_y, "lm")) {
    return(NULL)
  }

  vars <- ms_bruce_process_model7_vars(x, med_df, model_m, model_y)
  if (is.na(vars$predictor) || !nzchar(vars$predictor) ||
      is.na(vars$moderator) || !nzchar(vars$moderator) ||
      is.na(vars$interaction) || !nzchar(vars$interaction)) {
    return(NULL)
  }

  rows <- ms_bruce_process_conditional_indirect_rows(med_df, vars$moderator)
  if (!length(rows)) return(NULL)

  ci_meta <- ms_bruce_process_ci_meta(x, med_df)
  paths <- ms_bruce_process_model7_paths(model_m, model_y, vars)
  path_rows <- ms_bruce_process_model7_path_rows(model_m, model_y, vars)
  index_row <- ms_bruce_process_model7_index_row(model_m, model_y, vars)
  if (!is.null(index_row)) rows[[length(rows) + 1L]] <- index_row
  n <- ms_safe_numeric(stats::nobs(model_y))

  fields <- list(
    table_type = "mediation",
    source = "bruceR::PROCESS",
    process_id = as.integer(ms_safe_numeric(ms_bruce_first(x$process.id, NA_real_))),
    process_type = as.character(ms_bruce_first(x$process.type, "")),
    topology = "moderated_mediation",
    hayes_model = "7",
    moderation = list(
      hayes_model = "7",
      moderator = vars$moderator,
      moderated_path = "a",
      path = "X_to_M",
      interaction_term = vars$interaction
    ),
    rows = rows,
    columns = ms_bruce_process_bootstrap_indirect_columns(),
    note = ms_bruce_process_summary_note(x, n),
    effect_note = ms_bruce_process_model7_effect_note(
      x, n, has_index = !is.null(index_row)),
    path_note = ms_bruce_process_path_note(
      model_m, model_y, n, has_interaction = TRUE),
    treatment = vars$predictor,
    mediator = vars$mediator,
    mediators = I(c(vars$mediator)),
    moderator = vars$moderator,
    outcome = vars$outcome,
    controls = I(vars$covariates),
    n = as.integer(round(n)),
    boot = ci_meta$boot,
    ci_type = ci_meta$ci_type,
    conf_level = ci_meta$conf_level,
    sims = ci_meta$sims
  )
  if (!is.null(paths) && length(paths)) fields$paths <- paths
  if (!is.null(path_rows) && length(path_rows)) {
    fields$path_rows <- path_rows
    fields$path_columns <- ms_bruce_process_path_columns()
  }

  ms_build_envelope(
    type = "mediation_mediate",
    type_label = "Moderated mediation analysis (bruceR PROCESS Model 7)",
    call = ms_bruce_process_call(.call, "bruceR::PROCESS(..., mod.path = \"x-m\")"),
    fields = fields,
    raw_output = ms_capture_output(x),
    packages = ms_packages_basic(extras = c("bruceR", "mediation", "interactions")),
    card_kind = "table"
  )
}

ms_bruce_process_model8_payload <- function(x, .call = NULL) {
  med_df <- ms_bruce_process_mediation_df(x)
  model_m <- ms_bruce_process_model_m(x)
  model_y <- x$model.y
  if (!is.data.frame(med_df) ||
      !inherits(model_m, "lm") ||
      !inherits(model_y, "lm")) {
    return(NULL)
  }

  vars <- ms_bruce_process_model8_vars(x, med_df, model_m, model_y)
  if (is.na(vars$predictor) || !nzchar(vars$predictor) ||
      is.na(vars$moderator) || !nzchar(vars$moderator) ||
      is.na(vars$interaction) || !nzchar(vars$interaction) ||
      is.na(vars$direct_interaction) || !nzchar(vars$direct_interaction)) {
    return(NULL)
  }

  indirect_rows <- ms_bruce_process_conditional_indirect_rows(
    med_df, vars$moderator)
  if (!length(indirect_rows)) return(NULL)

  direct_rows <- ms_bruce_process_conditional_direct_rows(x, vars, model_y)
  rows <- indirect_rows
  index_row <- ms_bruce_process_model7_index_row(model_m, model_y, vars)
  if (!is.null(index_row)) rows[[length(rows) + 1L]] <- index_row
  if (length(direct_rows)) rows <- c(rows, direct_rows)

  ci_meta <- ms_bruce_process_ci_meta(x, med_df)
  paths <- ms_bruce_process_model8_paths(model_m, model_y, vars)
  path_rows <- ms_bruce_process_model8_path_rows(model_m, model_y, vars)
  n <- ms_safe_numeric(stats::nobs(model_y))

  fields <- list(
    table_type = "mediation",
    source = "bruceR::PROCESS",
    process_id = as.integer(ms_safe_numeric(ms_bruce_first(x$process.id, NA_real_))),
    process_type = as.character(ms_bruce_first(x$process.type, "")),
    topology = "moderated_mediation",
    hayes_model = "8",
    moderation = list(
      hayes_model = "8",
      moderator = vars$moderator,
      moderated_path = "a",
      path = "X_to_M",
      interaction_term = vars$interaction,
      direct_moderated_path = "X_to_Y",
      direct_interaction_term = vars$direct_interaction
    ),
    rows = rows,
    columns = ms_bruce_process_model8_effect_columns(),
    note = ms_bruce_process_summary_note(x, n),
    effect_note = ms_bruce_process_model8_effect_note(
      x, n,
      has_index = !is.null(index_row),
      has_directs = length(direct_rows) > 0L),
    path_note = ms_bruce_process_path_note(
      model_m, model_y, n, has_interaction = TRUE),
    treatment = vars$predictor,
    mediator = vars$mediator,
    mediators = I(c(vars$mediator)),
    moderator = vars$moderator,
    outcome = vars$outcome,
    controls = I(vars$covariates),
    n = as.integer(round(n)),
    boot = ci_meta$boot,
    ci_type = ci_meta$ci_type,
    conf_level = ci_meta$conf_level,
    sims = ci_meta$sims
  )
  if (!is.null(paths) && length(paths)) fields$paths <- paths
  if (!is.null(path_rows) && length(path_rows)) {
    fields$path_rows <- path_rows
    fields$path_columns <- ms_bruce_process_path_columns()
  }

  ms_build_envelope(
    type = "mediation_mediate",
    type_label = "Moderated mediation analysis (bruceR PROCESS Model 8)",
    call = ms_bruce_process_call(.call, "bruceR::PROCESS(..., mod.path = c(\"x-m\", \"x-y\"))"),
    fields = fields,
    raw_output = ms_capture_output(x),
    packages = ms_packages_basic(extras = c("bruceR", "mediation", "interactions")),
    card_kind = "table"
  )
}

ms_bruce_process_model14_payload <- function(x, .call = NULL) {
  med_df <- ms_bruce_process_mediation_df(x)
  model_m <- ms_bruce_process_model_m(x)
  model_y <- x$model.y
  if (!is.data.frame(med_df) ||
      !inherits(model_m, "lm") ||
      !inherits(model_y, "lm")) {
    return(NULL)
  }

  vars <- ms_bruce_process_model14_vars(x, med_df, model_m, model_y)
  if (is.na(vars$predictor) || !nzchar(vars$predictor) ||
      is.na(vars$moderator) || !nzchar(vars$moderator) ||
      is.na(vars$interaction) || !nzchar(vars$interaction)) {
    return(NULL)
  }

  rows <- ms_bruce_process_conditional_indirect_rows(med_df, vars$moderator)
  if (!length(rows)) return(NULL)

  ci_meta <- ms_bruce_process_ci_meta(x, med_df)
  paths <- ms_bruce_process_model14_paths(model_m, model_y, vars)
  path_rows <- ms_bruce_process_model14_path_rows(model_m, model_y, vars)
  index_row <- ms_bruce_process_model14_index_row(model_m, model_y, vars)
  if (!is.null(index_row)) rows[[length(rows) + 1L]] <- index_row
  n <- ms_safe_numeric(stats::nobs(model_y))

  fields <- list(
    table_type = "mediation",
    source = "bruceR::PROCESS",
    process_id = as.integer(ms_safe_numeric(ms_bruce_first(x$process.id, NA_real_))),
    process_type = as.character(ms_bruce_first(x$process.type, "")),
    topology = "moderated_mediation",
    hayes_model = "14",
    moderation = list(
      hayes_model = "14",
      moderator = vars$moderator,
      moderated_path = "b",
      path = "M_to_Y",
      interaction_term = vars$interaction
    ),
    rows = rows,
    columns = ms_bruce_process_bootstrap_indirect_columns(),
    note = ms_bruce_process_summary_note(x, n),
    effect_note = ms_bruce_process_model14_effect_note(
      x, n, has_index = !is.null(index_row)),
    path_note = ms_bruce_process_path_note(
      model_m, model_y, n, has_interaction = TRUE),
    treatment = vars$predictor,
    mediator = vars$mediator,
    mediators = I(c(vars$mediator)),
    moderator = vars$moderator,
    outcome = vars$outcome,
    controls = I(vars$covariates),
    n = as.integer(round(n)),
    boot = ci_meta$boot,
    ci_type = ci_meta$ci_type,
    conf_level = ci_meta$conf_level,
    sims = ci_meta$sims
  )
  if (!is.null(paths) && length(paths)) fields$paths <- paths
  if (!is.null(path_rows) && length(path_rows)) {
    fields$path_rows <- path_rows
    fields$path_columns <- ms_bruce_process_path_columns()
  }

  ms_build_envelope(
    type = "mediation_mediate",
    type_label = "Moderated mediation analysis (bruceR PROCESS Model 14)",
    call = ms_bruce_process_call(.call, "bruceR::PROCESS(..., mod.path = \"m-y\")"),
    fields = fields,
    raw_output = ms_capture_output(x),
    packages = ms_packages_basic(extras = c("bruceR", "mediation", "interactions")),
    card_kind = "table"
  )
}

ms_bruce_process_model15_payload <- function(x, .call = NULL) {
  med_df <- ms_bruce_process_mediation_df(x)
  model_m <- ms_bruce_process_model_m(x)
  model_y <- x$model.y
  if (!is.data.frame(med_df) ||
      !inherits(model_m, "lm") ||
      !inherits(model_y, "lm")) {
    return(NULL)
  }

  vars <- ms_bruce_process_model15_vars(x, med_df, model_m, model_y)
  if (is.na(vars$predictor) || !nzchar(vars$predictor) ||
      is.na(vars$moderator) || !nzchar(vars$moderator) ||
      is.na(vars$interaction) || !nzchar(vars$interaction) ||
      is.na(vars$direct_interaction) || !nzchar(vars$direct_interaction)) {
    return(NULL)
  }

  indirect_rows <- ms_bruce_process_conditional_indirect_rows(
    med_df, vars$moderator)
  if (!length(indirect_rows)) return(NULL)

  direct_rows <- ms_bruce_process_conditional_direct_rows(x, vars, model_y)
  rows <- indirect_rows
  index_row <- ms_bruce_process_model14_index_row(model_m, model_y, vars)
  if (!is.null(index_row)) rows[[length(rows) + 1L]] <- index_row
  if (length(direct_rows)) rows <- c(rows, direct_rows)

  ci_meta <- ms_bruce_process_ci_meta(x, med_df)
  paths <- ms_bruce_process_model15_paths(model_m, model_y, vars)
  path_rows <- ms_bruce_process_model15_path_rows(model_m, model_y, vars)
  n <- ms_safe_numeric(stats::nobs(model_y))

  fields <- list(
    table_type = "mediation",
    source = "bruceR::PROCESS",
    process_id = as.integer(ms_safe_numeric(ms_bruce_first(x$process.id, NA_real_))),
    process_type = as.character(ms_bruce_first(x$process.type, "")),
    topology = "moderated_mediation",
    hayes_model = "15",
    moderation = list(
      hayes_model = "15",
      moderator = vars$moderator,
      moderated_path = "b",
      path = "M_to_Y",
      interaction_term = vars$interaction,
      direct_moderated_path = "X_to_Y",
      direct_interaction_term = vars$direct_interaction
    ),
    rows = rows,
    columns = ms_bruce_process_model8_effect_columns(),
    note = ms_bruce_process_summary_note(x, n),
    effect_note = ms_bruce_process_model15_effect_note(
      x, n,
      has_index = !is.null(index_row),
      has_directs = length(direct_rows) > 0L),
    path_note = ms_bruce_process_path_note(
      model_m, model_y, n, has_interaction = TRUE),
    treatment = vars$predictor,
    mediator = vars$mediator,
    mediators = I(c(vars$mediator)),
    moderator = vars$moderator,
    outcome = vars$outcome,
    controls = I(vars$covariates),
    n = as.integer(round(n)),
    boot = ci_meta$boot,
    ci_type = ci_meta$ci_type,
    conf_level = ci_meta$conf_level,
    sims = ci_meta$sims
  )
  if (!is.null(paths) && length(paths)) fields$paths <- paths
  if (!is.null(path_rows) && length(path_rows)) {
    fields$path_rows <- path_rows
    fields$path_columns <- ms_bruce_process_path_columns()
  }

  ms_build_envelope(
    type = "mediation_mediate",
    type_label = "Moderated mediation analysis (bruceR PROCESS Model 15)",
    call = ms_bruce_process_call(.call, "bruceR::PROCESS(..., mod.path = c(\"m-y\", \"x-y\"))"),
    fields = fields,
    raw_output = ms_capture_output(x),
    packages = ms_packages_basic(extras = c("bruceR", "mediation", "interactions")),
    card_kind = "table"
  )
}

ms_bruce_process_unsupported_payload <- function(x, .call = NULL) {
  printed <- ms_capture_output_safe(x)
  process_id <- ms_safe_numeric(ms_bruce_first(x$process.id, NA_real_))
  process_type <- as.character(ms_bruce_first(x$process.type, ""))
  message <- paste(
    "Mellio recognizes this as a bruceR::PROCESS object, but currently",
    "fully formats bruceR PROCESS Models 4, 6, 7, 8, 14, and 15. Supported full",
    "lavaan mediation output includes simple, parallel, serial, and",
    "Hayes/PROCESS Models 7, 8, 14, and 15."
  )
  ms_build_envelope(
    type = "unsupported",
    type_label = "Unsupported bruceR PROCESS model",
    call = ms_bruce_process_call(.call, "bruceR::PROCESS(...)"),
    fields = list(
      class = I(c("bruceR_PROCESS", as.character(class(x)))),
      printed = printed$text,
      message = message,
      truncated = printed$truncated,
      source = "bruceR::PROCESS",
      process_id = process_id,
      process_type = process_type
    ),
    raw_output = printed$text,
    packages = ms_packages_basic(extras = "bruceR"),
    card_kind = "unsupported"
  )
}

ms_bruce_process_mediation_df <- function(x) {
  results <- x$results
  if (!is.list(results) || !length(results)) return(NULL)
  for (block in results) {
    if (is.list(block) && is.data.frame(block$mediation)) {
      return(block$mediation)
    }
  }
  NULL
}

ms_bruce_process_mediation_dfs <- function(x) {
  results <- x$results
  if (!is.list(results) || !length(results)) return(list())
  out <- list()
  for (block in results) {
    if (is.list(block) && is.data.frame(block$mediation)) {
      out[[length(out) + 1L]] <- block$mediation
    }
  }
  out
}

ms_bruce_process_lavaan_mediation_df <- function(x) {
  results <- x$results
  if (!is.list(results) || !length(results)) return(NULL)
  for (block in results) {
    if (is.list(block) && is.data.frame(block$lavaan.mediation)) {
      return(block$lavaan.mediation)
    }
  }
  NULL
}

ms_bruce_process_lavaan_syntax <- function(x) {
  results <- x$results
  if (!is.list(results) || !length(results)) return("")
  for (block in results) {
    syntax <- tryCatch(block$lavaan.syntax, error = function(e) NULL)
    if (length(syntax) && is.character(syntax)) {
      syntax <- paste(syntax, collapse = "\n")
      if (nzchar(trimws(syntax))) return(syntax)
    }
  }
  ""
}

ms_bruce_process_model_m <- function(x) {
  model_m <- x$model.m
  if (inherits(model_m, "lm")) return(model_m)
  if (!is.list(model_m) || !length(model_m)) return(NULL)
  for (model in model_m) {
    if (inherits(model, "lm")) return(model)
  }
  NULL
}

ms_bruce_process_model_m_list <- function(x) {
  model_m <- x$model.m
  if (inherits(model_m, "lm")) return(list(model_m))
  if (!is.list(model_m) || !length(model_m)) return(list())
  Filter(function(model) inherits(model, "lm"), model_m)
}

ms_bruce_process_model4_vars <- function(x, med_df, model_m, model_y) {
  mediator <- ms_bruce_lm_response(model_m)
  outcome <- ms_bruce_lm_response(model_y)
  m_terms <- ms_bruce_lm_terms(model_m)
  y_terms <- ms_bruce_lm_terms(model_y)
  common_terms <- intersect(m_terms, setdiff(y_terms, mediator))

  b <- ms_bruce_lm_coef(model_y, mediator)
  indirect <- ms_bruce_process_row_value(med_df, "Indirect", "Effect")
  predictor <- NA_character_
  if (length(common_terms) && !is.na(b) && !is.na(indirect)) {
    products <- vapply(common_terms, function(term) {
      a <- ms_bruce_lm_coef(model_m, term)
      if (is.na(a)) return(Inf)
      abs((a * b) - indirect)
    }, numeric(1))
    predictor <- common_terms[[which.min(products)]]
  }
  if (is.na(predictor) || !nzchar(predictor)) {
    fallback_terms <- setdiff(common_terms, mediator)
    predictor <- if (length(fallback_terms)) fallback_terms[[1L]] else NA_character_
  }

  covariates <- union(m_terms, y_terms)
  covariates <- setdiff(covariates, c(predictor, mediator))
  covariates <- covariates[nzchar(covariates)]

  list(
    predictor = predictor,
    mediator = mediator,
    outcome = outcome,
    covariates = covariates
  )
}

ms_bruce_process_parallel_vars <- function(x, med_dfs, model_m, model_y) {
  mediators <- vapply(model_m, ms_bruce_lm_response, character(1))
  outcome <- ms_bruce_lm_response(model_y)
  m_terms <- lapply(model_m, ms_bruce_lm_terms)
  y_terms <- ms_bruce_lm_terms(model_y)

  candidates <- Reduce(intersect, c(m_terms, list(y_terms)))
  candidates <- setdiff(candidates, mediators)
  candidates <- candidates[!grepl(":", candidates, fixed = TRUE)]

  predictor <- NA_character_
  usable_n <- min(length(model_m), length(med_dfs), length(mediators))
  if (length(candidates) && usable_n > 0L) {
    scores <- vapply(candidates, function(term) {
      diffs <- vapply(seq_len(usable_n), function(i) {
        a <- ms_bruce_lm_coef(model_m[[i]], term)
        b <- ms_bruce_lm_coef(model_y, mediators[[i]])
        indirect <- ms_bruce_process_row_value(med_dfs[[i]], "Indirect", "Effect")
        if (any(is.na(c(a, b, indirect)))) return(NA_real_)
        abs((a * b) - indirect)
      }, numeric(1))
      diffs <- diffs[is.finite(diffs)]
      if (!length(diffs)) return(Inf)
      mean(diffs)
    }, numeric(1))
    if (any(is.finite(scores))) {
      predictor <- candidates[[which.min(scores)]]
    }
  }

  if (is.na(predictor) || !nzchar(predictor)) {
    predictor <- if (length(candidates)) candidates[[1L]] else NA_character_
  }

  covariates <- unique(unlist(c(m_terms, list(y_terms)), use.names = FALSE))
  covariates <- setdiff(covariates, c(predictor, mediators))
  covariates <- covariates[!grepl(":", covariates, fixed = TRUE)]
  covariates <- covariates[nzchar(covariates)]

  list(
    predictor = predictor,
    mediators = mediators,
    outcome = outcome,
    covariates = covariates
  )
}

ms_bruce_process_model6_vars <- function(x, med_df, model_m, model_y) {
  mediators <- vapply(model_m, ms_bruce_lm_response, character(1))
  outcome <- ms_bruce_lm_response(model_y)
  m1_terms <- ms_bruce_lm_terms(model_m[[1L]])
  m2_terms <- ms_bruce_lm_terms(model_m[[2L]])
  y_terms <- ms_bruce_lm_terms(model_y)

  candidates <- Reduce(intersect, list(m1_terms, m2_terms, y_terms))
  candidates <- setdiff(candidates, mediators)
  candidates <- candidates[!grepl(":", candidates, fixed = TRUE)]

  predictor <- NA_character_
  target <- ms_bruce_process_row_value(med_df, "Indirect_All", "Estimate")
  if (length(candidates) && !is.na(target)) {
    d21 <- ms_bruce_lm_coef(model_m[[2L]], mediators[[1L]])
    b1 <- ms_bruce_lm_coef(model_y, mediators[[1L]])
    b2 <- ms_bruce_lm_coef(model_y, mediators[[2L]])
    scores <- vapply(candidates, function(term) {
      a1 <- ms_bruce_lm_coef(model_m[[1L]], term)
      a2 <- ms_bruce_lm_coef(model_m[[2L]], term)
      if (any(is.na(c(a1, a2, d21, b1, b2)))) return(Inf)
      indirect <- a1 * b1 + a2 * b2 + a1 * d21 * b2
      abs(indirect - target)
    }, numeric(1))
    if (any(is.finite(scores))) {
      predictor <- candidates[[which.min(scores)]]
    }
  }

  if (is.na(predictor) || !nzchar(predictor)) {
    syntax_predictor <- ms_bruce_process_model6_predictor_from_syntax(x, mediators)
    if (!is.na(syntax_predictor) && nzchar(syntax_predictor) &&
        syntax_predictor %in% union(m1_terms, union(m2_terms, y_terms))) {
      predictor <- syntax_predictor
    }
  }

  if (is.na(predictor) || !nzchar(predictor)) {
    predictor <- if (length(candidates)) candidates[[1L]] else NA_character_
  }

  covariates <- union(union(m1_terms, m2_terms), y_terms)
  covariates <- setdiff(covariates, c(predictor, mediators))
  covariates <- covariates[!grepl(":", covariates, fixed = TRUE)]
  covariates <- covariates[nzchar(covariates)]

  list(
    predictor = predictor,
    mediators = mediators,
    outcome = outcome,
    covariates = covariates
  )
}

ms_bruce_process_model6_predictor_from_syntax <- function(x, mediators) {
  syntax <- ms_bruce_process_lavaan_syntax(x)
  if (!nzchar(syntax)) return(NA_character_)
  lines <- unlist(strsplit(syntax, "\n", fixed = TRUE), use.names = FALSE)
  m1 <- as.character(mediators[[1L]] %||% "")
  if (!nzchar(m1)) return(NA_character_)
  for (line in lines) {
    lhs <- trimws(strsplit(line, "~", fixed = TRUE)[[1L]][[1L]] %||% "")
    if (!identical(lhs, m1)) next
    match <- regexec("a1\\s*\\*\\s*`?([^`+\\s]+)`?", line, perl = TRUE)
    hit <- regmatches(line, match)[[1L]]
    if (length(hit) >= 2L && nzchar(hit[[2L]])) return(hit[[2L]])
  }
  NA_character_
}

ms_bruce_process_model7_vars <- function(x, med_df, model_m, model_y) {
  mediator <- ms_bruce_lm_response(model_m)
  outcome <- ms_bruce_lm_response(model_y)
  m_terms <- ms_bruce_lm_terms(model_m)
  y_terms <- ms_bruce_lm_terms(model_y)

  moderator <- ms_bruce_process_condition_column(med_df)
  interactions <- m_terms[grepl(":", m_terms, fixed = TRUE)]
  interaction <- NA_character_
  parts <- character(0)
  if (length(interactions)) {
    if (!is.na(moderator) && nzchar(moderator)) {
      for (term in interactions) {
        term_parts <- ms_bruce_interaction_parts(term)
        if (moderator %in% term_parts) {
          interaction <- term
          parts <- term_parts
          break
        }
      }
    }
    if (is.na(interaction) || !nzchar(interaction)) {
      interaction <- interactions[[1L]]
      parts <- ms_bruce_interaction_parts(interaction)
    }
  }

  if ((is.na(moderator) || !nzchar(moderator)) && length(parts) >= 2L) {
    condition_names <- ms_bruce_process_result_condition_names(x)
    moderator_hit <- intersect(parts, condition_names)
    if (length(moderator_hit)) moderator <- moderator_hit[[1L]]
  }
  if ((is.na(moderator) || !nzchar(moderator)) && length(parts) >= 2L) {
    moderator <- parts[[2L]]
  }

  predictor <- NA_character_
  if (length(parts) >= 2L && !is.na(moderator) && nzchar(moderator)) {
    other_part <- setdiff(parts, moderator)
    if (length(other_part)) predictor <- other_part[[1L]]
  }
  if (is.na(predictor) || !nzchar(predictor)) {
    candidates <- setdiff(intersect(m_terms, y_terms),
                          c(mediator, moderator, interaction))
    predictor <- candidates[[1L]] %||% NA_character_
  }

  covariates <- union(m_terms, y_terms)
  covariates <- setdiff(covariates,
                        c(predictor, mediator, moderator, interaction))
  covariates <- covariates[nzchar(covariates)]

  list(
    predictor = predictor,
    mediator = mediator,
    moderator = moderator,
    interaction = interaction,
    outcome = outcome,
    covariates = covariates
  )
}

ms_bruce_process_model8_vars <- function(x, med_df, model_m, model_y) {
  vars <- ms_bruce_process_model7_vars(x, med_df, model_m, model_y)
  y_terms <- ms_bruce_lm_terms(model_y)
  parts <- ms_bruce_interaction_parts(vars$interaction)
  direct_interaction <- if (length(parts) >= 2L) {
    ms_bruce_find_interaction_term(y_terms, parts)
  } else {
    NA_character_
  }
  if ((is.na(direct_interaction) || !nzchar(direct_interaction)) &&
      !is.na(vars$interaction) && vars$interaction %in% y_terms) {
    direct_interaction <- vars$interaction
  }
  vars$direct_interaction <- direct_interaction
  vars$covariates <- setdiff(vars$covariates, direct_interaction)
  vars
}

ms_bruce_process_model14_vars <- function(x, med_df, model_m, model_y) {
  mediator <- ms_bruce_lm_response(model_m)
  outcome <- ms_bruce_lm_response(model_y)
  m_terms <- ms_bruce_lm_terms(model_m)
  y_terms <- ms_bruce_lm_terms(model_y)

  moderator <- ms_bruce_process_condition_column(med_df)
  interactions <- y_terms[grepl(":", y_terms, fixed = TRUE)]
  interaction <- NA_character_
  parts <- character(0)
  if (length(interactions)) {
    for (term in interactions) {
      term_parts <- ms_bruce_interaction_parts(term)
      has_mediator <- mediator %in% term_parts
      has_moderator <- !is.na(moderator) && nzchar(moderator) &&
        moderator %in% term_parts
      if (has_mediator && (has_moderator || is.na(moderator) || !nzchar(moderator))) {
        interaction <- term
        parts <- term_parts
        break
      }
    }
    if (is.na(interaction) || !nzchar(interaction)) {
      interaction <- interactions[[1L]]
      parts <- ms_bruce_interaction_parts(interaction)
    }
  }

  if ((is.na(moderator) || !nzchar(moderator)) && length(parts) >= 2L) {
    condition_names <- ms_bruce_process_result_condition_names(x)
    moderator_hit <- intersect(setdiff(parts, mediator), condition_names)
    if (length(moderator_hit)) moderator <- moderator_hit[[1L]]
  }
  if ((is.na(moderator) || !nzchar(moderator)) && length(parts) >= 2L) {
    other_part <- setdiff(parts, mediator)
    if (length(other_part)) moderator <- other_part[[1L]]
  }

  b_mean <- ms_bruce_process_conditional_mean_slope(
    x, model_y, mediator, moderator, interaction)
  if (is.na(b_mean)) b_mean <- ms_bruce_lm_coef(model_y, mediator)
  indirect_mean <- ms_bruce_process_condition_estimate(med_df, moderator, "mean")

  candidates <- setdiff(intersect(m_terms, y_terms),
                        c(mediator, moderator, interaction))
  candidates <- candidates[!grepl(":", candidates, fixed = TRUE)]
  predictor <- NA_character_
  if (length(candidates) && !is.na(b_mean) && !is.na(indirect_mean)) {
    scores <- vapply(candidates, function(term) {
      a <- ms_bruce_lm_coef(model_m, term)
      if (is.na(a)) return(Inf)
      abs((a * b_mean) - indirect_mean)
    }, numeric(1))
    predictor <- candidates[[which.min(scores)]]
  }
  if (is.na(predictor) || !nzchar(predictor)) {
    fallback <- setdiff(candidates, c(moderator, mediator))
    predictor <- fallback[[1L]] %||% NA_character_
  }

  covariates <- union(m_terms, y_terms)
  covariates <- setdiff(covariates,
                        c(predictor, mediator, moderator, interaction))
  covariates <- covariates[!grepl(":", covariates, fixed = TRUE)]
  covariates <- covariates[nzchar(covariates)]

  list(
    predictor = predictor,
    mediator = mediator,
    moderator = moderator,
    interaction = interaction,
    outcome = outcome,
    covariates = covariates
  )
}

ms_bruce_process_model15_vars <- function(x, med_df, model_m, model_y) {
  vars <- ms_bruce_process_model14_vars(x, med_df, model_m, model_y)
  y_terms <- ms_bruce_lm_terms(model_y)
  m_terms <- ms_bruce_lm_terms(model_m)
  direct_interaction <- NA_character_
  direct_parts <- character(0)
  interactions <- y_terms[grepl(":", y_terms, fixed = TRUE)]
  if (!is.na(vars$moderator) && nzchar(vars$moderator) &&
      length(interactions)) {
    for (term in interactions) {
      term_parts <- ms_bruce_interaction_parts(term)
      if (vars$moderator %in% term_parts &&
          !(vars$mediator %in% term_parts)) {
        direct_interaction <- term
        direct_parts <- term_parts
        break
      }
    }
  }
  if (!is.na(direct_interaction) && nzchar(direct_interaction)) {
    predictor_hit <- setdiff(direct_parts, vars$moderator)
    predictor_hit <- intersect(predictor_hit, union(m_terms, y_terms))
    if (length(predictor_hit)) vars$predictor <- predictor_hit[[1L]]
  } else if (!is.na(vars$predictor) && nzchar(vars$predictor) &&
             !is.na(vars$moderator) && nzchar(vars$moderator)) {
    direct_interaction <- ms_bruce_find_interaction_term(
      y_terms, c(vars$predictor, vars$moderator))
  }
  vars$direct_interaction <- direct_interaction
  covariates <- union(m_terms, y_terms)
  covariates <- setdiff(covariates,
                        c(vars$predictor, vars$mediator, vars$moderator,
                          vars$interaction, vars$direct_interaction))
  covariates <- covariates[!grepl(":", covariates, fixed = TRUE)]
  vars$covariates <- covariates[nzchar(covariates)]
  vars
}

ms_bruce_process_model4_paths <- function(model_m, model_y, vars) {
  treat <- vars$predictor
  mediator <- vars$mediator
  outcome <- vars$outcome
  if (is.na(treat) || !nzchar(treat) ||
      is.na(mediator) || !nzchar(mediator) ||
      is.na(outcome) || !nzchar(outcome)) {
    return(NULL)
  }

  treat_is_categorical <- ms_mediate_term_is_categorical(model_m, model_y, treat)
  sd_treat <- ms_mediate_term_sd(model_m, treat)
  sd_mediator_in_y <- ms_mediate_term_sd(model_y, mediator)
  sd_outcome_m <- ms_mediate_response_sd(model_m)
  sd_outcome_y <- ms_mediate_response_sd(model_y)

  paths <- list()
  a <- ms_mediate_coef_row(
    model_m, treat, "a",
    paste0(treat, " \u2192 ", mediator),
    sd_predictor = sd_treat,
    sd_outcome = sd_outcome_m,
    skip_std = treat_is_categorical
  )
  if (!is.null(a)) paths[[length(paths) + 1L]] <- a

  b <- ms_mediate_coef_row(
    model_y, mediator, "b",
    paste0(mediator, " \u2192 ", outcome),
    sd_predictor = sd_mediator_in_y,
    sd_outcome = sd_outcome_y,
    skip_std = FALSE
  )
  if (!is.null(b)) paths[[length(paths) + 1L]] <- b

  cprime <- ms_mediate_coef_row(
    model_y, treat, "c_prime",
    paste0(treat, " \u2192 ", outcome, " (direct)"),
    sd_predictor = sd_treat,
    sd_outcome = sd_outcome_y,
    skip_std = treat_is_categorical
  )
  if (!is.null(cprime)) paths[[length(paths) + 1L]] <- cprime

  if (length(paths) == 0L) NULL else paths
}

ms_bruce_process_model6_paths <- function(model_m, model_y, vars) {
  path_rows <- ms_bruce_process_model6_path_rows(model_m, model_y, vars)
  if (!length(path_rows)) return(NULL)
  keep <- vapply(path_rows, function(row) {
    row$path %in% c("a1", "a2", "d21", "b1", "b2", "c_prime")
  }, logical(1))
  out <- path_rows[keep]
  if (length(out) == 0L) NULL else out
}

ms_bruce_process_parallel_paths <- function(model_m, model_y, vars) {
  path_rows <- ms_bruce_process_parallel_path_rows(model_m, model_y, vars)
  if (!length(path_rows)) return(NULL)
  keep <- vapply(path_rows, function(row) {
    grepl("^(?:a|b)[0-9]+$", row$path) || identical(row$path, "c_prime")
  }, logical(1))
  out <- path_rows[keep]
  if (length(out) == 0L) NULL else out
}

ms_bruce_process_model7_paths <- function(model_m, model_y, vars) {
  path_rows <- ms_bruce_process_model7_path_rows(model_m, model_y, vars)
  if (!length(path_rows)) return(NULL)
  keep <- vapply(path_rows, function(row) {
    row$path %in% c("a1", "a2", "a3", "b", "c_prime", "c2")
  }, logical(1))
  out <- path_rows[keep]
  if (length(out) == 0L) NULL else out
}

ms_bruce_process_model8_paths <- function(model_m, model_y, vars) {
  path_rows <- ms_bruce_process_model8_path_rows(model_m, model_y, vars)
  if (!length(path_rows)) return(NULL)
  keep <- vapply(path_rows, function(row) {
    row$path %in% c("a1", "a2", "a3", "b", "c1", "c_prime", "c2", "c3")
  }, logical(1))
  out <- path_rows[keep]
  if (length(out) == 0L) NULL else out
}

ms_bruce_process_model14_paths <- function(model_m, model_y, vars) {
  path_rows <- ms_bruce_process_model14_path_rows(model_m, model_y, vars)
  if (!length(path_rows)) return(NULL)
  keep <- vapply(path_rows, function(row) {
    row$path %in% c("a", "b1", "b2", "b3", "c_prime")
  }, logical(1))
  out <- path_rows[keep]
  if (length(out) == 0L) NULL else out
}

ms_bruce_process_model15_paths <- function(model_m, model_y, vars) {
  path_rows <- ms_bruce_process_model15_path_rows(model_m, model_y, vars)
  if (!length(path_rows)) return(NULL)
  keep <- vapply(path_rows, function(row) {
    row$path %in% c("a", "b1", "b2", "b3", "c1", "c3")
  }, logical(1))
  out <- path_rows[keep]
  if (length(out) == 0L) NULL else out
}

ms_bruce_process_model4_path_rows <- function(model_m, model_y, vars) {
  treat <- vars$predictor
  mediator <- vars$mediator
  outcome <- vars$outcome
  if (is.na(treat) || !nzchar(treat) ||
      is.na(mediator) || !nzchar(mediator) ||
      is.na(outcome) || !nzchar(outcome)) {
    return(NULL)
  }

  rows <- list()
  append_row <- function(model, term, path_id, label, role, target) {
    sd_predictor <- ms_mediate_term_sd(model, term)
    sd_outcome <- ms_mediate_response_sd(model)
    skip_std <- ms_mediate_term_is_categorical(model, NULL, term)
    row <- ms_mediate_coef_row(
      model, term, path_id, label,
      sd_predictor = sd_predictor,
      sd_outcome = sd_outcome,
      skip_std = skip_std
    )
    if (is.null(row)) return()
    row$parameter <- label
    row$role <- role
    row$outcome <- target
    rows[[length(rows) + 1L]] <<- row
  }

  append_row(model_m, treat, "a",
             paste0("a: ", treat, " \u2192 ", mediator),
             "direct", mediator)
  append_row(model_y, mediator, "b",
             paste0("b: ", mediator, " \u2192 ", outcome),
             "direct", outcome)
  append_row(model_y, treat, "c_prime",
             paste0("c\u2032: ", treat, " \u2192 ", outcome),
             "direct", outcome)

  m_controls <- setdiff(ms_bruce_lm_terms(model_m), treat)
  y_controls <- setdiff(ms_bruce_lm_terms(model_y), c(treat, mediator))
  for (term in m_controls) {
    append_row(model_m, term, paste0("control_m_", term),
               paste0(term, " \u2192 ", mediator),
               "additional", mediator)
  }
  for (term in y_controls) {
    append_row(model_y, term, paste0("control_y_", term),
               paste0(term, " \u2192 ", outcome),
               "additional", outcome)
  }

  if (length(rows) == 0L) NULL else rows
}

ms_bruce_process_model6_path_rows <- function(model_m, model_y, vars) {
  treat <- vars$predictor
  mediators <- as.character(vars$mediators %||% character(0))
  outcome <- vars$outcome
  if (is.na(treat) || !nzchar(treat) ||
      length(mediators) != 2L ||
      any(is.na(mediators)) ||
      any(!nzchar(mediators)) ||
      is.na(outcome) || !nzchar(outcome)) {
    return(NULL)
  }
  m1 <- mediators[[1L]]
  m2 <- mediators[[2L]]

  rows <- list()
  append_row <- function(model, term, path_id, label, role, target) {
    sd_predictor <- ms_mediate_term_sd(model, term)
    sd_outcome <- ms_mediate_response_sd(model)
    skip_std <- ms_mediate_term_is_categorical(model, NULL, term)
    row <- ms_mediate_coef_row(
      model, term, path_id, label,
      sd_predictor = sd_predictor,
      sd_outcome = sd_outcome,
      skip_std = skip_std
    )
    if (is.null(row)) return()
    row$parameter <- label
    row$role <- role
    row$outcome <- target
    rows[[length(rows) + 1L]] <<- row
  }

  append_row(model_m[[1L]], treat, "a1",
             paste0("a1: ", treat, " \u2192 ", m1),
             "direct", m1)
  append_row(model_m[[2L]], treat, "a2",
             paste0("a2: ", treat, " \u2192 ", m2),
             "direct", m2)
  append_row(model_m[[2L]], m1, "d21",
             paste0("d21: ", m1, " \u2192 ", m2),
             "direct", m2)
  append_row(model_y, m1, "b1",
             paste0("b1: ", m1, " \u2192 ", outcome),
             "direct", outcome)
  append_row(model_y, m2, "b2",
             paste0("b2: ", m2, " \u2192 ", outcome),
             "direct", outcome)
  append_row(model_y, treat, "c_prime",
             paste0("c\u2032: ", treat, " \u2192 ", outcome),
             "direct", outcome)

  m1_controls <- setdiff(ms_bruce_lm_terms(model_m[[1L]]),
                         c(treat, mediators))
  m2_controls <- setdiff(ms_bruce_lm_terms(model_m[[2L]]),
                         c(treat, mediators))
  y_controls <- setdiff(ms_bruce_lm_terms(model_y),
                        c(treat, mediators))
  for (term in m1_controls) {
    append_row(model_m[[1L]], term, paste0("control_m1_", term),
               paste0(term, " \u2192 ", m1),
               "additional", m1)
  }
  for (term in m2_controls) {
    append_row(model_m[[2L]], term, paste0("control_m2_", term),
               paste0(term, " \u2192 ", m2),
               "additional", m2)
  }
  for (term in y_controls) {
    append_row(model_y, term, paste0("control_y_", term),
               paste0(term, " \u2192 ", outcome),
               "additional", outcome)
  }

  if (length(rows) == 0L) NULL else rows
}

ms_bruce_process_parallel_path_rows <- function(model_m, model_y, vars) {
  treat <- vars$predictor
  mediators <- as.character(vars$mediators %||% character(0))
  outcome <- vars$outcome
  if (is.na(treat) || !nzchar(treat) ||
      length(mediators) < 2L ||
      any(is.na(mediators)) ||
      any(!nzchar(mediators)) ||
      is.na(outcome) || !nzchar(outcome)) {
    return(NULL)
  }

  rows <- list()
  append_row <- function(model, term, path_id, label, role, target) {
    sd_predictor <- ms_mediate_term_sd(model, term)
    sd_outcome <- ms_mediate_response_sd(model)
    skip_std <- ms_mediate_term_is_categorical(model, NULL, term)
    row <- ms_mediate_coef_row(
      model, term, path_id, label,
      sd_predictor = sd_predictor,
      sd_outcome = sd_outcome,
      skip_std = skip_std
    )
    if (is.null(row)) return()
    row$parameter <- label
    row$role <- role
    row$outcome <- target
    rows[[length(rows) + 1L]] <<- row
  }

  for (i in seq_along(mediators)) {
    if (i > length(model_m)) next
    mediator <- mediators[[i]]
    append_row(model_m[[i]], treat, paste0("a", i),
               paste0("a", i, ": ", treat, " \u2192 ", mediator),
               "direct", mediator)
  }
  for (i in seq_along(mediators)) {
    mediator <- mediators[[i]]
    append_row(model_y, mediator, paste0("b", i),
               paste0("b", i, ": ", mediator, " \u2192 ", outcome),
               "direct", outcome)
  }
  append_row(model_y, treat, "c_prime",
             paste0("c\u2032: ", treat, " \u2192 ", outcome),
             "direct", outcome)

  for (i in seq_along(mediators)) {
    if (i > length(model_m)) next
    mediator <- mediators[[i]]
    m_controls <- setdiff(ms_bruce_lm_terms(model_m[[i]]), treat)
    for (term in m_controls) {
      append_row(model_m[[i]], term, paste0("control_m", i, "_", term),
                 paste0(term, " \u2192 ", mediator),
                 "additional", mediator)
    }
  }
  y_controls <- setdiff(ms_bruce_lm_terms(model_y), c(treat, mediators))
  for (term in y_controls) {
    append_row(model_y, term, paste0("control_y_", term),
               paste0(term, " \u2192 ", outcome),
               "additional", outcome)
  }

  if (length(rows) == 0L) NULL else rows
}

ms_bruce_process_model7_path_rows <- function(model_m, model_y, vars) {
  treat <- vars$predictor
  mediator <- vars$mediator
  moderator <- vars$moderator
  interaction <- vars$interaction
  outcome <- vars$outcome
  if (is.na(treat) || !nzchar(treat) ||
      is.na(mediator) || !nzchar(mediator) ||
      is.na(moderator) || !nzchar(moderator) ||
      is.na(interaction) || !nzchar(interaction) ||
      is.na(outcome) || !nzchar(outcome)) {
    return(NULL)
  }

  rows <- list()
  append_row <- function(model, term, path_id, label, role, target) {
    sd_predictor <- ms_mediate_term_sd(model, term)
    sd_outcome <- ms_mediate_response_sd(model)
    skip_std <- ms_mediate_term_is_categorical(model, NULL, term)
    row <- ms_mediate_coef_row(
      model, term, path_id, label,
      sd_predictor = sd_predictor,
      sd_outcome = sd_outcome,
      skip_std = skip_std
    )
    if (is.null(row)) return()
    row$parameter <- label
    row$role <- role
    row$outcome <- target
    rows[[length(rows) + 1L]] <<- row
  }

  append_row(model_m, treat, "a1",
             paste0("a1: ", treat, " \u2192 ", mediator),
             "direct", mediator)
  append_row(model_m, moderator, "a2",
             paste0("a2: ", moderator, " \u2192 ", mediator),
             "moderator", mediator)
  append_row(model_m, interaction, "a3",
             paste0("a3: ", ms_bruce_pretty_ordered_interaction(
                    interaction, treat, moderator),
                    " \u2192 ", mediator),
             "interaction", mediator)
  append_row(model_y, mediator, "b",
             paste0("b: ", mediator, " \u2192 ", outcome),
             "direct", outcome)
  append_row(model_y, treat, "c_prime",
             paste0("c\u2032: ", treat, " \u2192 ", outcome),
             "direct", outcome)
  append_row(model_y, moderator, "c2",
             paste0(moderator, " \u2192 ", outcome),
             "moderator", outcome)

  m_controls <- setdiff(ms_bruce_lm_terms(model_m),
                        c(treat, moderator, interaction))
  y_controls <- setdiff(ms_bruce_lm_terms(model_y),
                        c(treat, mediator, moderator, interaction))
  for (term in m_controls) {
    append_row(model_m, term, paste0("control_m_", term),
               paste0(term, " \u2192 ", mediator),
               "additional", mediator)
  }
  for (term in y_controls) {
    append_row(model_y, term, paste0("control_y_", term),
               paste0(term, " \u2192 ", outcome),
               "additional", outcome)
  }

  if (length(rows) == 0L) NULL else rows
}

ms_bruce_process_model8_path_rows <- function(model_m, model_y, vars) {
  treat <- vars$predictor
  mediator <- vars$mediator
  moderator <- vars$moderator
  interaction <- vars$interaction
  direct_interaction <- vars$direct_interaction
  outcome <- vars$outcome
  if (is.na(treat) || !nzchar(treat) ||
      is.na(mediator) || !nzchar(mediator) ||
      is.na(moderator) || !nzchar(moderator) ||
      is.na(interaction) || !nzchar(interaction) ||
      is.na(direct_interaction) || !nzchar(direct_interaction) ||
      is.na(outcome) || !nzchar(outcome)) {
    return(NULL)
  }

  rows <- list()
  append_row <- function(model, term, path_id, label, role, target) {
    sd_predictor <- ms_mediate_term_sd(model, term)
    sd_outcome <- ms_mediate_response_sd(model)
    skip_std <- identical(role, "interaction") ||
      ms_mediate_term_is_categorical(model, NULL, term)
    row <- ms_mediate_coef_row(
      model, term, path_id, label,
      sd_predictor = sd_predictor,
      sd_outcome = sd_outcome,
      skip_std = skip_std
    )
    if (is.null(row)) return()
    row$parameter <- label
    row$role <- role
    row$outcome <- target
    rows[[length(rows) + 1L]] <<- row
  }

  append_row(model_m, treat, "a1",
             paste0("a1: ", treat, " \u2192 ", mediator),
             "direct", mediator)
  append_row(model_m, moderator, "a2",
             paste0("a2: ", moderator, " \u2192 ", mediator),
             "moderator", mediator)
  append_row(model_m, interaction, "a3",
             paste0("a3: ", ms_bruce_pretty_ordered_interaction(
                    interaction, treat, moderator),
                    " \u2192 ", mediator),
             "interaction", mediator)
  append_row(model_y, mediator, "b",
             paste0("b: ", mediator, " \u2192 ", outcome),
             "direct", outcome)
  append_row(model_y, treat, "c1",
             paste0("c1: ", treat, " \u2192 ", outcome),
             "direct", outcome)
  append_row(model_y, moderator, "c2",
             paste0("c2: ", moderator, " \u2192 ", outcome),
             "moderator", outcome)
  append_row(model_y, direct_interaction, "c3",
             paste0("c3: ", ms_bruce_pretty_ordered_interaction(
                    direct_interaction, treat, moderator),
                    " \u2192 ", outcome),
             "interaction", outcome)

  m_controls <- setdiff(ms_bruce_lm_terms(model_m),
                        c(treat, moderator, interaction))
  y_controls <- setdiff(ms_bruce_lm_terms(model_y),
                        c(treat, mediator, moderator, interaction,
                          direct_interaction))
  for (term in m_controls) {
    append_row(model_m, term, paste0("control_m_", term),
               paste0(term, " \u2192 ", mediator),
               "additional", mediator)
  }
  for (term in y_controls) {
    append_row(model_y, term, paste0("control_y_", term),
               paste0(term, " \u2192 ", outcome),
               "additional", outcome)
  }

  if (length(rows) == 0L) NULL else rows
}

ms_bruce_process_model14_path_rows <- function(model_m, model_y, vars) {
  treat <- vars$predictor
  mediator <- vars$mediator
  moderator <- vars$moderator
  interaction <- vars$interaction
  outcome <- vars$outcome
  if (is.na(treat) || !nzchar(treat) ||
      is.na(mediator) || !nzchar(mediator) ||
      is.na(moderator) || !nzchar(moderator) ||
      is.na(interaction) || !nzchar(interaction) ||
      is.na(outcome) || !nzchar(outcome)) {
    return(NULL)
  }

  rows <- list()
  append_row <- function(model, term, path_id, label, role, target) {
    sd_predictor <- ms_mediate_term_sd(model, term)
    sd_outcome <- ms_mediate_response_sd(model)
    skip_std <- identical(role, "interaction") ||
      ms_mediate_term_is_categorical(model, NULL, term)
    row <- ms_mediate_coef_row(
      model, term, path_id, label,
      sd_predictor = sd_predictor,
      sd_outcome = sd_outcome,
      skip_std = skip_std
    )
    if (is.null(row)) return()
    row$parameter <- label
    row$role <- role
    row$outcome <- target
    rows[[length(rows) + 1L]] <<- row
  }

  append_row(model_m, treat, "a",
             paste0("a: ", treat, " \u2192 ", mediator),
             "direct", mediator)
  append_row(model_y, mediator, "b1",
             paste0("b1: ", mediator, " \u2192 ", outcome),
             "direct", outcome)
  append_row(model_y, moderator, "b2",
             paste0("b2: ", moderator, " \u2192 ", outcome),
             "moderator", outcome)
  append_row(model_y, interaction, "b3",
             paste0("b3: ", ms_bruce_pretty_ordered_interaction(
                    interaction, mediator, moderator),
                    " \u2192 ", outcome),
             "interaction", outcome)
  append_row(model_y, treat, "c_prime",
             paste0("c\u2032: ", treat, " \u2192 ", outcome),
             "direct", outcome)

  m_controls <- setdiff(ms_bruce_lm_terms(model_m), treat)
  y_controls <- setdiff(ms_bruce_lm_terms(model_y),
                        c(treat, mediator, moderator, interaction))
  for (term in m_controls) {
    append_row(model_m, term, paste0("control_m_", term),
               paste0(term, " \u2192 ", mediator),
               "additional", mediator)
  }
  for (term in y_controls) {
    append_row(model_y, term, paste0("control_y_", term),
               paste0(term, " \u2192 ", outcome),
               "additional", outcome)
  }

  if (length(rows) == 0L) NULL else rows
}

ms_bruce_process_model15_path_rows <- function(model_m, model_y, vars) {
  treat <- vars$predictor
  mediator <- vars$mediator
  moderator <- vars$moderator
  interaction <- vars$interaction
  direct_interaction <- vars$direct_interaction
  outcome <- vars$outcome
  if (is.na(treat) || !nzchar(treat) ||
      is.na(mediator) || !nzchar(mediator) ||
      is.na(moderator) || !nzchar(moderator) ||
      is.na(interaction) || !nzchar(interaction) ||
      is.na(direct_interaction) || !nzchar(direct_interaction) ||
      is.na(outcome) || !nzchar(outcome)) {
    return(NULL)
  }

  rows <- list()
  append_row <- function(model, term, path_id, label, role, target) {
    sd_predictor <- ms_mediate_term_sd(model, term)
    sd_outcome <- ms_mediate_response_sd(model)
    skip_std <- identical(role, "interaction") ||
      ms_mediate_term_is_categorical(model, NULL, term)
    row <- ms_mediate_coef_row(
      model, term, path_id, label,
      sd_predictor = sd_predictor,
      sd_outcome = sd_outcome,
      skip_std = skip_std
    )
    if (is.null(row)) return()
    row$parameter <- label
    row$role <- role
    row$outcome <- target
    rows[[length(rows) + 1L]] <<- row
  }

  append_row(model_m, treat, "a",
             paste0("a: ", treat, " \u2192 ", mediator),
             "direct", mediator)
  append_row(model_y, mediator, "b1",
             paste0("b1: ", mediator, " \u2192 ", outcome),
             "direct", outcome)
  append_row(model_y, moderator, "b2",
             paste0("b2: ", moderator, " \u2192 ", outcome),
             "moderator", outcome)
  append_row(model_y, interaction, "b3",
             paste0("b3: ", ms_bruce_pretty_ordered_interaction(
                    interaction, mediator, moderator),
                    " \u2192 ", outcome),
             "interaction", outcome)
  append_row(model_y, treat, "c1",
             paste0("c1: ", treat, " \u2192 ", outcome),
             "direct", outcome)
  append_row(model_y, direct_interaction, "c3",
             paste0("c3: ", ms_bruce_pretty_ordered_interaction(
                    direct_interaction, treat, moderator),
                    " \u2192 ", outcome),
             "interaction", outcome)

  m_controls <- setdiff(ms_bruce_lm_terms(model_m), treat)
  y_controls <- setdiff(ms_bruce_lm_terms(model_y),
                        c(treat, mediator, moderator, interaction,
                          direct_interaction))
  for (term in m_controls) {
    append_row(model_m, term, paste0("control_m_", term),
               paste0(term, " \u2192 ", mediator),
               "additional", mediator)
  }
  for (term in y_controls) {
    append_row(model_y, term, paste0("control_y_", term),
               paste0(term, " \u2192 ", outcome),
               "additional", outcome)
  }

  if (length(rows) == 0L) NULL else rows
}

ms_bruce_process_model7_index_row <- function(model_m, model_y, vars) {
  a3 <- ms_bruce_lm_coef(model_m, vars$interaction)
  b <- ms_bruce_lm_coef(model_y, vars$mediator)
  if (is.na(a3) || is.na(b)) return(NULL)
  list(
    effect = "Index of moderated mediation",
    raw_effect = "index_modmed",
    estimate = a3 * b,
    std_error = NA_real_,
    statistic = NA_real_,
    p_value = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_
  )
}

ms_bruce_process_model14_index_row <- function(model_m, model_y, vars) {
  a <- ms_bruce_lm_coef(model_m, vars$predictor)
  b3 <- ms_bruce_lm_coef(model_y, vars$interaction)
  if (is.na(a) || is.na(b3)) return(NULL)
  list(
    effect = "Index of moderated mediation",
    raw_effect = "index_modmed",
    estimate = a * b3,
    std_error = NA_real_,
    statistic = NA_real_,
    p_value = NA_real_,
    ci_lower = NA_real_,
    ci_upper = NA_real_
  )
}

ms_bruce_process_model6_rows <- function(df, mediators) {
  labels <- rownames(df)
  if (is.null(labels) || !length(labels)) {
    labels <- as.character(df[["Effect"]] %||% seq_len(nrow(df)))
  }
  rows <- lapply(seq_len(nrow(df)), function(i) {
    raw_label <- trimws(as.character(labels[[i]]))
    ms_bruce_process_bootstrap_indirect_row(list(
      effect = ms_bruce_process_model6_effect_label(raw_label, mediators),
      raw_effect = raw_label,
      statistic_type = "z",
      estimate = ms_bruce_df_numeric(df, i, c("Estimate", "Effect")),
      std_error = ms_bruce_df_numeric(df, i, c("S.E.", "SE", "Std. Error")),
      statistic = ms_bruce_df_numeric(df, i, c("z", "Z")),
      p_value = ms_bruce_df_numeric(df, i, c("pval", "p.value", "p")),
      ci_lower = ms_bruce_df_numeric(df, i, c("BootLLCI", "LLCI", "lower", "conf.low")),
      ci_upper = ms_bruce_df_numeric(df, i, c("BootULCI", "ULCI", "upper", "conf.high"))
    ))
  })
  rows <- Filter(function(row) !is.na(row$estimate) && nzchar(row$effect), rows)
  ranks <- vapply(rows, function(row) {
    switch(row$raw_effect,
      Ind_X_M1_Y = 1L,
      Ind_X_M2_Y = 2L,
      Ind_X_M1_M2_Y = 3L,
      Indirect_All = 4L,
      Direct = 5L,
      Total = 6L,
      99L
    )
  }, integer(1))
  rows[order(ranks)]
}

ms_bruce_process_parallel_rows <- function(dfs, mediators) {
  dfs <- dfs[vapply(dfs, is.data.frame, logical(1))]
  mediators <- as.character(mediators %||% character(0))
  if (!length(dfs) || !length(mediators)) return(list())

  rows <- list()
  usable_n <- min(length(dfs), length(mediators))
  for (i in seq_len(usable_n)) {
    df <- dfs[[i]]
    mediator <- mediators[[i]]
    row <- ms_bruce_process_parallel_row(
      df, "^Indirect\\b", paste0("Indirect effect via ", mediator),
      paste0("indirect_", mediator)
    )
    if (!is.null(row)) {
      row$mediator <- mediator
      row$mediator_index <- i - 1L
      rows[[length(rows) + 1L]] <- row
    }
  }

  total_indirect <- NA_real_
  indirect_estimates <- vapply(rows, function(row) {
    ms_safe_numeric(row$estimate)
  }, numeric(1))
  if (length(indirect_estimates) && all(is.finite(indirect_estimates))) {
    total_indirect <- sum(indirect_estimates)
    rows[[length(rows) + 1L]] <- list(
      effect = "Total indirect effect",
      raw_effect = "total_indirect",
      statistic_type = "z",
      estimate = total_indirect,
      std_error = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      ci_note = "not returned",
      synthesised = TRUE
    )
  }

  direct <- ms_bruce_process_parallel_row(
    dfs[[1L]], "^Direct\\b", "Direct effect", "direct"
  )
  if (!is.null(direct)) rows[[length(rows) + 1L]] <- direct

  if (!is.null(direct) && is.finite(total_indirect)) {
    rows[[length(rows) + 1L]] <- list(
      effect = "Total effect",
      raw_effect = "total",
      statistic_type = "z",
      estimate = ms_safe_numeric(direct$estimate) + total_indirect,
      std_error = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      ci_note = "not returned",
      synthesised = TRUE
    )
  }

  Filter(function(row) !is.na(row$estimate) && nzchar(row$effect), rows)
}

ms_bruce_process_parallel_row <- function(df, pattern, effect, raw_effect) {
  if (!is.data.frame(df)) return(NULL)
  labels <- rownames(df) %||% character(0)
  hit <- grep(pattern, labels, ignore.case = TRUE, perl = TRUE)
  if (!length(hit)) return(NULL)
  i <- hit[[1L]]
  estimate <- ms_bruce_df_numeric(df, i, "Effect")
  if (is.na(estimate)) return(NULL)
  ms_bruce_process_bootstrap_indirect_row(list(
    effect = effect,
    raw_effect = raw_effect,
    statistic_type = "z",
    estimate = estimate,
    std_error = ms_bruce_df_numeric(df, i, c("S.E.", "SE", "Std. Error")),
    statistic = ms_bruce_df_numeric(df, i, c("z", "Z")),
    p_value = ms_bruce_df_numeric(df, i, c("pval", "p.value", "p")),
    ci_lower = ms_bruce_df_numeric(df, i, c("LLCI", "lower", "conf.low")),
    ci_upper = ms_bruce_df_numeric(df, i, c("ULCI", "upper", "conf.high"))
  ))
}

ms_bruce_process_model6_effect_label <- function(label, mediators) {
  value <- trimws(as.character(label %||% ""))
  mediators <- as.character(mediators %||% character(0))
  m1 <- if (length(mediators) >= 1L && nzchar(mediators[[1L]])) {
    mediators[[1L]]
  } else {
    "M1"
  }
  m2 <- if (length(mediators) >= 2L && nzchar(mediators[[2L]])) {
    mediators[[2L]]
  } else {
    "M2"
  }
  switch(value,
    Indirect_All = "Total indirect effect",
    Ind_X_M1_Y = paste0("Indirect effect via ", m1),
    Ind_X_M2_Y = paste0("Indirect effect via ", m2),
    Ind_X_M1_M2_Y = "Serial indirect effect",
    Direct = "Direct effect",
    Total = "Total effect",
    value
  )
}

ms_bruce_process_mediation_rows <- function(df) {
  labels <- rownames(df)
  if (is.null(labels) || !length(labels)) {
    labels <- as.character(df[["Effect"]] %||% seq_len(nrow(df)))
  }
  rows <- lapply(seq_len(nrow(df)), function(i) {
    raw_label <- as.character(labels[[i]])
    ms_bruce_process_bootstrap_indirect_row(list(
      effect = ms_bruce_process_effect_label(raw_label),
      raw_effect = raw_label,
      estimate = ms_bruce_df_numeric(df, i, "Effect"),
      std_error = ms_bruce_df_numeric(df, i, c("S.E.", "SE", "Std. Error")),
      statistic = ms_bruce_df_numeric(df, i, c("z", "Z")),
      p_value = ms_bruce_df_numeric(df, i, c("pval", "p.value", "p")),
      ci_lower = ms_bruce_df_numeric(df, i, c("LLCI", "lower", "conf.low")),
      ci_upper = ms_bruce_df_numeric(df, i, c("ULCI", "upper", "conf.high"))
    ))
  })
  Filter(function(row) !is.na(row$estimate), rows)
}

ms_bruce_process_bootstrap_indirect_row <- function(row) {
  if (is.null(row)) return(row)
  label <- tolower(trimws(paste(row$effect %||% "", row$raw_effect %||% "")))
  has_ci <- !is.na(ms_safe_numeric(row$ci_lower)) &&
    !is.na(ms_safe_numeric(row$ci_upper))
  if (grepl("indirect", label, fixed = TRUE) && isTRUE(has_ci)) {
    row$statistic_type <- "bootstrap_ci"
    row$statistic <- NA_real_
    row$p_value <- NA_real_
  }
  row
}

ms_bruce_process_conditional_indirect_rows <- function(df, moderator) {
  condition_col <- ms_bruce_process_condition_column_raw(df, moderator)
  labels <- if (!is.na(condition_col) && condition_col %in% names(df)) {
    as.character(df[[condition_col]])
  } else {
    rownames(df) %||% as.character(seq_len(nrow(df)))
  }
  rows <- lapply(seq_len(nrow(df)), function(i) {
    condition <- ms_bruce_process_condition_value(labels[[i]], i)
    label <- ms_bruce_process_condition_label(labels[[i]], i)
    list(
      effect = paste0("Conditional indirect effect (", label, " ",
                      moderator, ")"),
      raw_effect = paste0("conditional_indirect_", condition),
      condition = condition,
      condition_label = paste0(label, " ", moderator),
      condition_value = ms_bruce_process_condition_numeric(labels[[i]]),
      moderator = moderator,
      statistic_type = "bootstrap_ci",
      estimate = ms_bruce_df_numeric(df, i, "Effect"),
      std_error = ms_bruce_df_numeric(df, i, c("S.E.", "SE", "Std. Error")),
      statistic = NA_real_,
      p_value = NA_real_,
      ci_lower = ms_bruce_df_numeric(df, i, c("LLCI", "lower", "conf.low")),
      ci_upper = ms_bruce_df_numeric(df, i, c("ULCI", "upper", "conf.high"))
    )
  })
  Filter(function(row) !is.na(row$estimate), rows)
}

ms_bruce_process_conditional_direct_rows <- function(x, vars, model_y) {
  slopes <- ms_bruce_process_simple_slopes_df(
    x, model_y, vars$predictor, vars$moderator, vars$direct_interaction)
  if (!is.data.frame(slopes)) return(list())
  condition_col <- ms_bruce_process_condition_column_raw(
    slopes, vars$moderator)
  labels <- if (!is.na(condition_col) && condition_col %in% names(slopes)) {
    as.character(slopes[[condition_col]])
  } else {
    rownames(slopes) %||% as.character(seq_len(nrow(slopes)))
  }
  rows <- lapply(seq_len(nrow(slopes)), function(i) {
    condition <- ms_bruce_process_condition_value(labels[[i]], i)
    label <- ms_bruce_process_condition_label(labels[[i]], i)
    list(
      effect = paste0("Conditional direct effect (", label, " ",
                      vars$moderator, ")"),
      raw_effect = paste0("conditional_direct_", condition),
      condition = condition,
      condition_label = paste0(label, " ", vars$moderator),
      condition_value = ms_bruce_process_condition_numeric(labels[[i]]),
      moderator = vars$moderator,
      statistic_type = "t",
      estimate = ms_bruce_df_numeric(slopes, i, "Effect"),
      std_error = ms_bruce_df_numeric(slopes, i, c("S.E.", "SE", "Std. Error")),
      statistic = ms_bruce_df_numeric(slopes, i, c("t", "T", "z", "Z")),
      p_value = ms_bruce_df_numeric(slopes, i, c("pval", "p.value", "p")),
      ci_lower = ms_bruce_df_numeric(slopes, i, c("LLCI", "lower", "conf.low")),
      ci_upper = ms_bruce_df_numeric(slopes, i, c("ULCI", "upper", "conf.high"))
    )
  })
  Filter(function(row) !is.na(row$estimate), rows)
}

ms_bruce_process_effect_columns <- function() {
  list(
    list(key = "effect", label = "Effect", format = "text"),
    list(key = "estimate", label = "Estimate", format = "number"),
    list(key = "std_error", label = "SE", format = "number"),
    list(key = "statistic", label = "z", format = "statistic"),
    list(key = "p_value", label = "p", format = "pvalue"),
    list(key = "ci", label = "Bootstrap 95% CI", format = "ci")
  )
}

ms_bruce_process_bootstrap_indirect_columns <- function() {
  list(
    list(key = "effect", label = "Effect", format = "text"),
    list(key = "estimate", label = "Estimate", format = "number"),
    list(key = "std_error", label = "SE", format = "number"),
    list(key = "ci", label = "Bootstrap 95% CI", format = "ci")
  )
}

ms_bruce_process_model8_effect_columns <- function() {
  list(
    list(key = "effect", label = "Effect", format = "text"),
    list(key = "estimate", label = "Estimate", format = "number"),
    list(key = "std_error", label = "SE", format = "number"),
    list(key = "statistic", label = "test", format = "statistic"),
    list(key = "p_value", label = "p", format = "pvalue"),
    list(key = "ci", label = "95% CI", format = "ci")
  )
}

ms_bruce_process_path_columns <- function() {
  list(
    list(key = "parameter", label = "Parameter", format = "text"),
    list(key = "estimate", label = "B", format = "number"),
    list(key = "std_error", label = "SE", format = "number"),
    list(key = "statistic", label = "t", format = "statistic"),
    list(key = "p_value", label = "p", format = "pvalue"),
    list(key = "ci", label = "95% CI", format = "ci"),
    list(key = "std_estimate", label = "\u03b2", format = "bounded")
  )
}

ms_bruce_process_summary_note <- function(x, n = NA_real_) {
  bits <- c(
    paste0("Effects are from bruceR::PROCESS Model ",
           ms_bruce_first(x$process.id, "")),
    ms_bruce_process_n_note(n),
    "mediation-effect intervals are reported as provided by bruceR",
    "bootstrap CIs are used for inference for indirect effects; z and p are omitted for indirect-effect rows",
    "path coefficients are extracted from the returned lm models"
  )
  bits <- bits[nzchar(bits)]
  paste0(paste(bits, collapse = "; "), ".")
}

ms_bruce_process_effect_note <- function(x, n = NA_real_) {
  bits <- c(
    paste0("Effects are from bruceR::PROCESS Model ",
           ms_bruce_first(x$process.id, "")),
    ms_bruce_process_n_note(n),
    "mediation-effect intervals are reported as provided by bruceR",
    "bootstrap CIs are used for inference for indirect effects; z and p are omitted for indirect-effect rows"
  )
  bits <- bits[nzchar(bits)]
  paste0(paste(bits, collapse = "; "), ".")
}

ms_bruce_process_parallel_summary_note <- function(x, n = NA_real_) {
  bits <- c(
    paste0("Effects are from bruceR::PROCESS Model ",
           ms_bruce_first(x$process.id, "")),
    ms_bruce_process_n_note(n),
    "specific indirect-effect intervals are reported as provided by bruceR",
    "bootstrap CIs are used for inference for specific indirect effects; z and p are omitted for specific indirect-effect rows",
    "total indirect and total effects are computed as point estimates because bruceR does not return intervals for these totals",
    "path coefficients are extracted from the returned lm models"
  )
  bits <- bits[nzchar(bits)]
  paste0(paste(bits, collapse = "; "), ".")
}

ms_bruce_process_parallel_effect_note <- function(x, n = NA_real_) {
  bits <- c(
    paste0("Parallel mediation effects are from bruceR::PROCESS Model ",
           ms_bruce_first(x$process.id, "")),
    ms_bruce_process_n_note(n),
    "specific indirect-effect intervals are reported as provided by bruceR",
    "bootstrap CIs are used for inference for specific indirect effects; z and p are omitted for specific indirect-effect rows",
    "total indirect and total effects are computed as point estimates because bruceR does not return intervals for these totals"
  )
  bits <- bits[nzchar(bits)]
  paste0(paste(bits, collapse = "; "), ".")
}

ms_bruce_process_model6_summary_note <- function(x, n = NA_real_) {
  bits <- c(
    paste0("Effects are from bruceR::PROCESS Model ",
           ms_bruce_first(x$process.id, "")),
    ms_bruce_process_n_note(n),
    "serial mediation-effect intervals are reported as provided by bruceR/lavaan",
    "bootstrap CIs are used for inference for indirect effects; z and p are omitted for indirect-effect rows",
    "path coefficients are extracted from the returned lm models"
  )
  bits <- bits[nzchar(bits)]
  paste0(paste(bits, collapse = "; "), ".")
}

ms_bruce_process_model6_effect_note <- function(x, n = NA_real_) {
  bits <- c(
    paste0("Serial mediation effects are from bruceR::PROCESS Model ",
           ms_bruce_first(x$process.id, "")),
    ms_bruce_process_n_note(n),
    "bootstrap intervals are reported as provided by bruceR/lavaan",
    "bootstrap CIs are used for inference for indirect effects; z and p are omitted for indirect-effect rows"
  )
  bits <- bits[nzchar(bits)]
  paste0(paste(bits, collapse = "; "), ".")
}

ms_bruce_process_model7_effect_note <- function(x, n = NA_real_,
                                                has_index = FALSE) {
  bits <- c(
    paste0("Conditional indirect effects are from bruceR::PROCESS Model ",
           ms_bruce_first(x$process.id, "")),
    ms_bruce_process_n_note(n),
    "mediation-effect intervals are reported as provided by bruceR"
  )
  bits <- c(bits,
            "bootstrap CIs are used for inference for conditional indirect effects")
  if (isTRUE(has_index)) {
    bits <- c(bits,
              "the index of moderated mediation is computed as a3 \u00d7 b and reported as a point estimate because bruceR does not return an interval for the index")
  }
  bits <- bits[nzchar(bits)]
  paste0(paste(bits, collapse = "; "), ".")
}

ms_bruce_process_model8_effect_note <- function(x, n = NA_real_,
                                                has_index = FALSE,
                                                has_directs = FALSE) {
  bits <- c(
    paste0("Conditional effects are from bruceR::PROCESS Model ",
           ms_bruce_first(x$process.id, "")),
    ms_bruce_process_n_note(n),
    "conditional indirect intervals are bootstrap intervals as provided by bruceR",
    "bootstrap CIs are used for inference for conditional indirect effects"
  )
  if (isTRUE(has_directs)) {
    bits <- c(bits,
              "conditional direct intervals are analytic simple-slope CIs as provided by bruceR",
              "the test and p columns are shown for conditional direct effects")
  }
  if (isTRUE(has_index)) {
    bits <- c(bits,
              "the index of moderated mediation is computed as a3 \u00d7 b and reported as a point estimate because bruceR does not return an interval for the index")
  }
  bits <- bits[nzchar(bits)]
  paste0(paste(bits, collapse = "; "), ".")
}

ms_bruce_process_model14_effect_note <- function(x, n = NA_real_,
                                                 has_index = FALSE) {
  bits <- c(
    paste0("Conditional indirect effects are from bruceR::PROCESS Model ",
           ms_bruce_first(x$process.id, "")),
    ms_bruce_process_n_note(n),
    "mediation-effect intervals are reported as provided by bruceR"
  )
  bits <- c(bits,
            "bootstrap CIs are used for inference for conditional indirect effects")
  if (isTRUE(has_index)) {
    bits <- c(bits,
              "the index of moderated mediation is computed as a \u00d7 b3 and reported as a point estimate because bruceR does not return an interval for the index")
  }
  bits <- bits[nzchar(bits)]
  paste0(paste(bits, collapse = "; "), ".")
}

ms_bruce_process_model15_effect_note <- function(x, n = NA_real_,
                                                 has_index = FALSE,
                                                 has_directs = FALSE) {
  bits <- c(
    paste0("Conditional effects are from bruceR::PROCESS Model ",
           ms_bruce_first(x$process.id, "")),
    ms_bruce_process_n_note(n),
    "conditional indirect intervals are bootstrap intervals as provided by bruceR",
    "bootstrap CIs are used for inference for conditional indirect effects"
  )
  if (isTRUE(has_directs)) {
    bits <- c(bits,
              "conditional direct intervals are analytic simple-slope CIs as provided by bruceR",
              "the test and p columns are shown for conditional direct effects")
  }
  if (isTRUE(has_index)) {
    bits <- c(bits,
              "the index of moderated mediation is computed as a \u00d7 b3 and reported as a point estimate because bruceR does not return an interval for the index")
  }
  bits <- bits[nzchar(bits)]
  paste0(paste(bits, collapse = "; "), ".")
}

ms_bruce_process_model6_path_note <- function(model_m, model_y,
                                              n = NA_real_) {
  bits <- c(
    ms_bruce_process_n_note(n),
    "Path coefficients are extracted from the lm models returned by bruceR::PROCESS",
    "SE, t, p, and path CIs are analytic lm summaries"
  )
  r2 <- ms_bruce_process_r_squared_note_models(c(model_m, list(model_y)))
  if (nzchar(r2)) bits <- c(bits, r2)
  bits <- bits[nzchar(bits)]
  paste0(paste(bits, collapse = "; "), ".")
}

ms_bruce_process_path_note <- function(model_m, model_y, n = NA_real_,
                                       has_interaction = FALSE) {
  bits <- c(
    ms_bruce_process_n_note(n),
    "Path coefficients are extracted from the lm models returned by bruceR::PROCESS",
    "SE, t, p, and path CIs are analytic lm summaries"
  )
  if (isTRUE(has_interaction)) {
    bits <- c(bits,
              "\u03b2 is not reported for interaction paths because standardized product-term coefficients are scale-dependent")
  }
  r2 <- ms_bruce_process_r_squared_note(model_m, model_y)
  if (nzchar(r2)) bits <- c(bits, r2)
  bits <- bits[nzchar(bits)]
  paste0(paste(bits, collapse = "; "), ".")
}

ms_bruce_process_n_note <- function(n) {
  n <- ms_safe_numeric(n)
  if (is.na(n) || n <= 0) return("")
  paste0("N = ", as.integer(round(n)))
}

ms_bruce_process_r_squared_note <- function(model_m, model_y) {
  ms_bruce_process_r_squared_note_models(list(model_m, model_y))
}

ms_bruce_process_r_squared_note_models <- function(models) {
  models <- if (is.list(models)) models else list(models)
  bits <- vapply(models, function(model) {
    if (!inherits(model, "lm")) return("")
    outcome <- ms_bruce_lm_response(model)
    sm <- tryCatch(summary(model), error = function(e) NULL)
    r2 <- if (!is.null(sm)) ms_safe_numeric(sm$r.squared) else NA_real_
    if (is.na(r2) || is.na(outcome) || !nzchar(outcome)) return("")
    paste0(outcome, " = ", ms_bruce_process_r2_value(r2))
  }, character(1))
  bits <- bits[nzchar(bits)]
  if (!length(bits)) return("")
  paste0("R\u00b2 by regression: ", paste(bits, collapse = "; "))
}

ms_bruce_process_r2_value <- function(value) {
  text <- formatC(ms_safe_numeric(value), digits = 2L, format = "f")
  sub("^0\\.", ".", text)
}

ms_bruce_process_ci_meta <- function(x, med_df) {
  cols <- names(med_df) %||% character(0)
  lowered <- tolower(cols)
  has_boot <- any(grepl("boot.*ci|ci.*boot", cols, ignore.case = TRUE)) ||
    all(c("bootllci", "bootulci") %in% lowered) ||
    any(grepl("^boot", lowered)) ||
    all(c("llci", "ulci") %in% lowered)
  list(
    boot = isTRUE(has_boot),
    ci_type = if (isTRUE(has_boot)) "perc" else NA_character_,
    conf_level = ms_bruce_process_conf_level(cols),
    sims = ms_bruce_process_sims(x)
  )
}

ms_bruce_process_conf_level <- function(cols) {
  lowered <- tolower(as.character(cols %||% character(0)))
  match <- regexpr("[0-9]+(?=\\s*%\\s*(?:boot\\s*)?ci)", lowered, perl = TRUE)
  hit <- regmatches(lowered, match)
  hit <- hit[nzchar(hit)]
  if (length(hit)) {
    pct <- ms_safe_numeric(hit[[1L]])
    if (!is.na(pct) && pct > 0 && pct <= 100) return(pct / 100)
  }
  0.95
}

ms_bruce_process_sims <- function(x) {
  keys <- c("sims", "simulations", "nsim", "boot", "R")
  for (key in keys) {
    value <- tryCatch(x[[key]], error = function(e) NULL)
    value <- ms_safe_numeric(ms_bruce_first(value, value))
    if (!is.na(value) && value > 0) return(as.integer(round(value)))
  }
  NA_integer_
}

ms_bruce_process_effect_label <- function(label) {
  value <- trimws(as.character(label %||% ""))
  lower <- tolower(value)
  if (grepl("^indirect\\b", lower)) return("Indirect effect")
  if (grepl("^direct\\b", lower)) return("Direct effect")
  if (grepl("^total\\b", lower)) return("Total effect")
  value
}

ms_bruce_process_row_value <- function(df, row_pattern, col) {
  if (!is.data.frame(df) || !(col %in% names(df))) return(NA_real_)
  rn <- rownames(df) %||% character(0)
  hit <- grep(row_pattern, rn, ignore.case = TRUE)
  if (!length(hit)) return(NA_real_)
  ms_safe_numeric(df[[col]][[hit[[1L]]]])
}

ms_bruce_process_condition_column <- function(df, preferred = NA_character_) {
  raw <- ms_bruce_process_condition_column_raw(df, preferred)
  ms_bruce_clean_process_name(raw)
}

ms_bruce_process_condition_column_raw <- function(df, preferred = NA_character_) {
  if (!is.data.frame(df)) return(NA_character_)
  effect_cols <- tolower(c(
    "effect", "s.e.", "se", "std. error", "llci", "ulci",
    "lower", "upper", "conf.low", "conf.high",
    "p", "pval", "p.value", "z", "t", "f", "df1", "df2",
    "[boot 95% ci]", "[95% ci]"
  ))
  names_df <- names(df) %||% character(0)
  if (!is.na(preferred) && nzchar(preferred)) {
    clean_names <- vapply(names_df, ms_bruce_clean_process_name, character(1))
    hit <- which(clean_names == preferred)
    if (length(hit)) return(names_df[[hit[[1L]]]])
  }
  candidates <- names_df[!(tolower(names_df) %in% effect_cols)]
  if (!length(candidates)) return(NA_character_)
  candidates[[1L]]
}

ms_bruce_process_result_condition_names <- function(x) {
  results <- x$results
  if (!is.list(results)) return(character(0))
  out <- character(0)
  for (block in results) {
    if (!is.list(block)) next
    conditional <- block$conditional
    if (is.data.frame(conditional)) {
      out <- union(out, vapply(names(conditional) %||% character(0),
                               ms_bruce_clean_process_name,
                               character(1)))
    }
  }
  out
}

ms_bruce_process_simple_slopes_df <- function(x, model, predictor,
                                              moderator, interaction) {
  results <- x$results
  if (!is.list(results)) return(NULL)
  candidates <- list()
  for (block in results) {
    if (is.list(block) && is.data.frame(block$simple.slopes)) {
      candidates[[length(candidates) + 1L]] <- block$simple.slopes
    }
  }
  if (!length(candidates)) return(NULL)
  if (length(candidates) == 1L) return(candidates[[1L]])
  scores <- vapply(candidates, ms_bruce_process_simple_slope_score,
                   numeric(1),
                   model = model,
                   predictor = predictor,
                   moderator = moderator,
                   interaction = interaction)
  finite <- is.finite(scores)
  if (!any(finite)) return(candidates[[1L]])
  candidates[[which.min(scores)]]
}

ms_bruce_process_simple_slope_score <- function(df, model, predictor,
                                                moderator, interaction) {
  if (!is.data.frame(df) || !inherits(model, "lm")) return(Inf)
  base <- ms_bruce_lm_coef(model, predictor)
  mod <- ms_bruce_lm_coef(model, interaction)
  if (is.na(base) || is.na(mod)) return(Inf)
  condition_col <- ms_bruce_process_condition_column_raw(df, moderator)
  if (is.na(condition_col) || !(condition_col %in% names(df))) return(Inf)
  labels <- as.character(df[[condition_col]])
  w <- vapply(labels, ms_bruce_process_condition_numeric, numeric(1))
  effects <- vapply(seq_len(nrow(df)), function(i) {
    ms_bruce_df_numeric(df, i, "Effect")
  }, numeric(1))
  conditions <- vapply(seq_len(nrow(df)), function(i) {
    ms_bruce_process_condition_value(labels[[i]], i)
  }, character(1))
  mean_idx <- which(conditions == "mean" & is.finite(w))
  center <- if (length(mean_idx)) w[[mean_idx[[1L]]]] else 0
  expected <- base + mod * (w - center)
  keep <- is.finite(expected) & is.finite(effects)
  if (!any(keep)) return(Inf)
  mean(abs(effects[keep] - expected[keep]))
}

ms_bruce_process_conditional_mean_slope <- function(x, model, predictor,
                                                    moderator, interaction) {
  slopes <- ms_bruce_process_simple_slopes_df(
    x, model, predictor, moderator, interaction)
  if (!is.data.frame(slopes)) return(NA_real_)
  ms_bruce_process_condition_estimate(slopes, moderator, "mean")
}

ms_bruce_process_condition_estimate <- function(df, moderator,
                                                condition = "mean") {
  if (!is.data.frame(df) || !nrow(df)) return(NA_real_)
  condition_col <- ms_bruce_process_condition_column_raw(df, moderator)
  labels <- if (!is.na(condition_col) && condition_col %in% names(df)) {
    as.character(df[[condition_col]])
  } else {
    rownames(df) %||% as.character(seq_len(nrow(df)))
  }
  conditions <- vapply(seq_along(labels), function(i) {
    ms_bruce_process_condition_value(labels[[i]], i)
  }, character(1))
  idx <- which(conditions == condition)
  if (!length(idx)) idx <- 1L
  ms_bruce_df_numeric(df, idx[[1L]], "Effect")
}

ms_bruce_process_condition_value <- function(label, index) {
  text <- tolower(as.character(label %||% ""))
  if (grepl("-\\s*sd|low", text)) return("low")
  if (grepl("\\+\\s*sd|high", text)) return("high")
  if (grepl("mean", text)) return("mean")
  paste0("condition_", index)
}

ms_bruce_process_condition_label <- function(label, index) {
  value <- ms_bruce_process_condition_value(label, index)
  switch(value,
    low = "low",
    mean = "mean",
    high = "high",
    trimws(as.character(label %||% value))
  )
}

ms_bruce_process_condition_numeric <- function(label) {
  text <- trimws(as.character(label %||% ""))
  hit <- regmatches(text, regexpr("[-+]?[0-9]*\\.?[0-9]+", text))
  if (!length(hit) || !nzchar(hit)) return(NA_real_)
  ms_safe_numeric(hit[[1L]])
}

ms_bruce_interaction_parts <- function(term) {
  parts <- strsplit(as.character(term %||% ""), ":", fixed = TRUE)[[1L]]
  parts[nzchar(parts)]
}

ms_bruce_find_interaction_term <- function(terms, parts) {
  terms <- as.character(terms %||% character(0))
  target <- sort(as.character(parts %||% character(0)))
  if (length(target) < 2L) return(NA_character_)
  for (term in terms) {
    term_parts <- sort(ms_bruce_interaction_parts(term))
    if (length(term_parts) == length(target) &&
        identical(term_parts, target)) {
      return(term)
    }
  }
  NA_character_
}

ms_bruce_pretty_interaction <- function(term) {
  paste(ms_bruce_interaction_parts(term), collapse = " \u00d7 ")
}

ms_bruce_pretty_ordered_interaction <- function(term, first, second) {
  parts <- ms_bruce_interaction_parts(term)
  first <- as.character(first %||% "")
  second <- as.character(second %||% "")
  if (first %in% parts && second %in% parts) {
    ordered <- c(first, second, setdiff(parts, c(first, second)))
    return(paste(ordered[nzchar(ordered)], collapse = " \u00d7 "))
  }
  ms_bruce_pretty_interaction(term)
}

ms_bruce_clean_process_name <- function(value) {
  text <- trimws(as.character(value %||% ""))
  text <- gsub("^['\"]+|['\"]+$", "", text)
  trimws(text)
}

ms_bruce_df_numeric <- function(df, i, keys) {
  for (key in keys) {
    if (key %in% names(df)) return(ms_safe_numeric(df[[key]][[i]]))
  }
  NA_real_
}

ms_bruce_lm_response <- function(model) {
  vars <- all.vars(stats::formula(model))
  if (length(vars)) vars[[1L]] else NA_character_
}

ms_bruce_lm_terms <- function(model) {
  as.character(attr(stats::terms(model), "term.labels") %||% character(0))
}

ms_bruce_lm_coef <- function(model, term) {
  coefs <- stats::coef(model)
  if (is.null(coefs) || !(term %in% names(coefs))) return(NA_real_)
  ms_safe_numeric(unname(coefs[[term]]))
}

ms_bruce_process_call <- function(.call, fallback) {
  if (!is.null(.call) && length(.call)) {
    call <- trimws(gsub("\\s+", " ", paste(.call, collapse = " ")))
    if (!is.na(call) && nzchar(call)) return(call)
  }
  fallback
}

ms_bruce_first <- function(x, default = NULL) {
  if (is.null(x) || !length(x)) return(default)
  x[[1L]]
}
