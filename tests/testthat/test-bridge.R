# Tests for the R bridge: mellio_payload(), mellio_open(), and helpers.
# Schema: docs/STATS-R-BRIDGE-SCHEMA.md

# ── Internal helpers ─────────────────────────────────────────────────

quiet_mellio_open <- function(x, ..., .call = NULL) {
  withr::local_options(list(lifecycle_verbosity = "quiet"))
  if (is.null(.call)) {
    .call <- paste(deparse(substitute(x), width.cutoff = 500L), collapse = " ")
  }
  mellio_open(x, ..., .call = .call)
}

prepare_lavaan_fit_session <- function() {
  cores <- tryCatch(
    suppressWarnings(parallel::detectCores()),
    error = function(e) NA_integer_
  )
  if (!is.na(cores)) return(invisible(TRUE))

  cache_env <- tryCatch(
    getFromNamespace("lavaan_cache_env", "lavaan"),
    error = function(e) NULL
  )
  options_default <- tryCatch(
    getFromNamespace("lav_options_default", "lavaan"),
    error = function(e) NULL
  )
  if (!is.environment(cache_env) || !is.function(options_default)) {
    return(invisible(FALSE))
  }

  suppressWarnings(options_default())
  if (!exists("opt.default", envir = cache_env, inherits = FALSE) ||
      !exists("opt.check", envir = cache_env, inherits = FALSE)) {
    return(invisible(FALSE))
  }

  opt_default <- get("opt.default", envir = cache_env)
  opt_check <- get("opt.check", envir = cache_env)
  opt_default$ncpus <- 1L
  if (is.null(opt_check$ncpus)) {
    opt_check$ncpus <- list(oklen = c(1L, 1L))
  }
  opt_check$ncpus$nm <- list(
    bounds = c(1, 1),
    first.in = TRUE,
    last.in = TRUE
  )
  assign("opt.default", opt_default, envir = cache_env)
  assign("opt.check", opt_check, envir = cache_env)
  invisible(TRUE)
}

skip_if_lavaan_unusable <- local({
  checked <- FALSE
  ok <- FALSE
  reason <- NULL

  function() {
    skip_if_not_installed("lavaan")
    prepare_lavaan_fit_session()
    if (!checked) {
      checked <<- TRUE
      ok <<- isTRUE(tryCatch({
        env <- new.env(parent = emptyenv())
        suppressWarnings(utils::data("HolzingerSwineford1939",
                                     package = "lavaan",
                                     envir = env))
        model <- "visual =~ x1 + x2 + x3"
        suppressWarnings(suppressMessages(
          lavaan::cfa(model, data = env$HolzingerSwineford1939)
        ))
        TRUE
      }, error = function(e) {
        reason <<- conditionMessage(e)
        FALSE
      }))
    }
    if (!ok) {
      skip(paste("lavaan cannot fit models in this R session:", reason))
    }
  }
})

sem_fixture_param <- function(lhs, op, rhs, estimate = NA_real_,
                              std_error = NA_real_, statistic = NA_real_,
                              p_value = NA_real_, ci_lower = NA_real_,
                              ci_upper = NA_real_, std_estimate = NA_real_,
                              label = NULL, group = NULL,
                              diagram_hidden = FALSE,
                              diagram_hidden_reason = NULL) {
  out <- list(
    lhs = lhs,
    op = op,
    rhs = rhs,
    estimate = estimate,
    std_error = std_error,
    statistic = statistic,
    p_value = p_value,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    std_estimate = std_estimate
  )
  if (!is.null(label)) out$label <- label
  if (!is.null(group)) out$group <- group
  if (isTRUE(diagram_hidden)) out$diagram_hidden <- TRUE
  if (!is.null(diagram_hidden_reason)) {
    out$diagram_hidden_reason <- diagram_hidden_reason
  }
  out
}

sem_fixture_nodes_by_id <- function(diagram) {
  stats::setNames(diagram$nodes, vapply(diagram$nodes, function(node) {
    node$id
  }, character(1)))
}

sem_fixture_edge <- function(diagram, source, target, type) {
  hits <- Filter(function(edge) {
    identical(edge$source, source) &&
      identical(edge$target, target) &&
      identical(edge$type, type)
  }, diagram$edges)
  expect_length(hits, 1L)
  hits[[1]]
}

test_that("lavaan type helpers do not treat every defined parameter as mediation", {
  cfa <- data.frame(
    lhs = c("visual", "visual", "visual"),
    op = c("=~", "=~", "=~"),
    rhs = c("x1", "x2", "x3"),
    stringsAsFactors = FALSE
  )
  expect_equal(mellio:::ms_lavaan_payload_type(cfa), "lavaan_cfa")
  expect_equal(mellio:::ms_lavaan_model_class(cfa), "cfa")

  observed_sem <- data.frame(
    lhs = c("m", "y"),
    op = c("~", "~"),
    rhs = c("x", "m"),
    stringsAsFactors = FALSE
  )
  expect_equal(mellio:::ms_lavaan_payload_type(observed_sem), "lavaan_sem")
  expect_equal(mellio:::ms_lavaan_model_class(observed_sem), "observed_sem")

  latent_sem <- data.frame(
    lhs = c("f1", "f1", "f2", "f2"),
    op = c("=~", "=~", "=~", "~"),
    rhs = c("x1", "x2", "y1", "f1"),
    stringsAsFactors = FALSE
  )
  expect_equal(mellio:::ms_lavaan_payload_type(latent_sem), "lavaan_sem")
  expect_equal(mellio:::ms_lavaan_model_class(latent_sem), "latent_sem")

  contrast_sem <- data.frame(
    lhs = c("m", "y", "contrast_b_paths"),
    op = c("~", "~", ":="),
    rhs = c("x", "m", "b1 - b2"),
    stringsAsFactors = FALSE
  )
  expect_equal(mellio:::ms_lavaan_payload_type(contrast_sem), "lavaan_sem")
  expect_equal(mellio:::ms_lavaan_model_class(contrast_sem), "observed_sem")

  mediation <- data.frame(
    lhs = c("m", "y", "indirect"),
    op = c("~", "~", ":="),
    rhs = c("x", "m", "a*b"),
    stringsAsFactors = FALSE
  )
  expect_equal(mellio:::ms_lavaan_payload_type(mediation), "lavaan_mediation")
  expect_equal(mellio:::ms_lavaan_model_class(mediation), "mediation")

  disconnected <- data.frame(
    lhs = c("m", "y", "indirect"),
    op = c("~", "~", ":="),
    rhs = c("x", "z", "a*b"),
    stringsAsFactors = FALSE
  )
  expect_equal(mellio:::ms_lavaan_payload_type(disconnected), "lavaan_sem")
  expect_equal(mellio:::ms_lavaan_model_class(disconnected), "observed_sem")
})

test_that("lavaan structural diagram helper builds SEM graph data without fitting lavaan", {
  params <- list(
    list(lhs = "realistic_threat", op = "~", rhs = "stance",
         estimate = 0.48, std_error = 0.05, p_value = 0.001,
         ci_lower = 0.38, ci_upper = 0.58, std_estimate = 0.44),
    list(lhs = "symbolic_threat", op = "~", rhs = "stance",
         estimate = 0.17, std_error = 0.07, p_value = 0.02,
         ci_lower = 0.03, ci_upper = 0.31, std_estimate = 0.16),
    list(lhs = "moral_condemnation", op = "~", rhs = "stance",
         estimate = 0.04, std_error = 0.06, p_value = 0.51,
         ci_lower = -0.08, ci_upper = 0.16, std_estimate = 0.04),
    list(lhs = "moral_condemnation", op = "~", rhs = "realistic_threat",
         estimate = 0.47, std_error = 0.06, p_value = 0.001,
         ci_lower = 0.35, ci_upper = 0.59, std_estimate = 0.45),
    list(lhs = "dehumanization", op = "~", rhs = "moral_condemnation",
         estimate = 0.29, std_error = 0.07, p_value = 0.001,
         ci_lower = 0.15, ci_upper = 0.43, std_estimate = 0.28)
  )

  diagram <- mellio:::ms_lavaan_structural_diagram(
    params,
    type = "lavaan_sem",
    type_label = "Structural equation model",
    n = 301
  )
  expect_equal(diagram$type, "structural_path_diagram")
  expect_equal(diagram$model_class, "observed_sem")
  expect_equal(diagram$n, 301)
  expect_length(diagram$nodes, 5L)
  expect_length(diagram$edges, 5L)
  expect_true(all(vapply(diagram$edges, function(edge) {
    edge$type == "structural" && edge$op == "~" &&
      !is.null(edge$estimate) && !is.null(edge$p_value)
  }, logical(1))))
})

test_that("lavaan structural diagram fixture preserves CFA measurement and covariance roles", {
  params <- list(
    sem_fixture_param("visual", "=~", "x1", estimate = 1.00, label = "l1"),
    sem_fixture_param("visual", "=~", "x2", estimate = 0.73, p_value = 0.001,
                      ci_lower = 0.62, ci_upper = 0.84, std_estimate = 0.68,
                      label = "l2"),
    sem_fixture_param("textual", "=~", "x4", estimate = 1.00, label = "l4"),
    sem_fixture_param("textual", "=~", "x5", estimate = 0.91, p_value = 0.001,
                      ci_lower = 0.80, ci_upper = 1.02, std_estimate = 0.77,
                      label = "l5"),
    sem_fixture_param("visual", "~~", "textual", estimate = 0.52,
                      p_value = 0.001, ci_lower = 0.31, ci_upper = 0.73,
                      std_estimate = 0.42)
  )

  diagram <- mellio:::ms_lavaan_structural_diagram(
    params,
    type = "lavaan_cfa",
    type_label = "Confirmatory factor analysis",
    n = 301
  )

  expect_equal(diagram$model_type, "lavaan_cfa")
  expect_equal(diagram$model_class, "cfa")
  expect_equal(diagram$title, "Confirmatory factor analysis")
  expect_equal(diagram$n, 301)
  expect_length(diagram$nodes, 6L)
  expect_length(diagram$edges, 5L)

  nodes <- sem_fixture_nodes_by_id(diagram)
  expect_equal(nodes$visual$role, "latent")
  expect_false(nodes$visual$observed)
  expect_equal(nodes$textual$role, "latent")
  expect_false(nodes$textual$observed)
  expect_equal(nodes$x1$role, "indicator")
  expect_true(nodes$x1$observed)

  loading <- sem_fixture_edge(diagram, "visual", "x2", "measurement")
  expect_equal(loading$op, "=~")
  expect_equal(loading$label, "l2")
  expect_equal(loading$estimate, 0.73)
  expect_equal(loading$std_estimate, 0.68)

  covariance <- sem_fixture_edge(diagram, "visual", "textual", "covariance")
  expect_equal(covariance$op, "~~")
  expect_equal(covariance$estimate, 0.52)
  expect_equal(covariance$std_estimate, 0.42)
})

test_that("lavaan structural diagram fixture preserves latent SEM roles and edge directions", {
  params <- list(
    sem_fixture_param("engagement", "=~", "eng1", estimate = 1.00),
    sem_fixture_param("engagement", "=~", "eng2", estimate = 0.82,
                      p_value = 0.001, std_estimate = 0.74),
    sem_fixture_param("burnout", "=~", "burn1", estimate = 1.00),
    sem_fixture_param("burnout", "=~", "burn2", estimate = 0.88,
                      p_value = 0.001, std_estimate = 0.79),
    sem_fixture_param("burnout", "~", "engagement", estimate = -0.63,
                      std_error = 0.08, statistic = -7.88, p_value = 0.001,
                      ci_lower = -0.79, ci_upper = -0.47, std_estimate = -0.58),
    sem_fixture_param("burnout", "~", "age", estimate = 0.12,
                      std_error = 0.05, statistic = 2.40, p_value = 0.016,
                      ci_lower = 0.02, ci_upper = 0.22, std_estimate = 0.11)
  )

  diagram <- mellio:::ms_lavaan_structural_diagram(
    params,
    type = "lavaan_sem",
    type_label = "Structural equation model",
    n = 214
  )

  expect_equal(diagram$model_class, "latent_sem")
  expect_equal(diagram$n, 214)
  expect_length(diagram$nodes, 7L)
  expect_length(diagram$edges, 6L)

  nodes <- sem_fixture_nodes_by_id(diagram)
  expect_equal(nodes$engagement$role, "latent")
  expect_equal(nodes$burnout$role, "latent")
  expect_equal(nodes$age$role, "exogenous")
  expect_true(nodes$age$observed)

  latent_path <- sem_fixture_edge(diagram, "engagement", "burnout", "structural")
  expect_equal(latent_path$op, "~")
  expect_equal(latent_path$source_label, "engagement")
  expect_equal(latent_path$target_label, "burnout")
  expect_equal(latent_path$estimate, -0.63)
  expect_equal(latent_path$std_estimate, -0.58)

  age_path <- sem_fixture_edge(diagram, "age", "burnout", "structural")
  expect_equal(age_path$estimate, 0.12)
})

test_that("lavaan structural diagram fixture preserves hidden covariate metadata", {
  params <- list(
    sem_fixture_param("moral_condemnation", "~", "realistic_threat",
                      estimate = 0.47, p_value = 0.001,
                      std_estimate = 0.45),
    sem_fixture_param("dehumanization", "~", "moral_condemnation",
                      estimate = 0.29, p_value = 0.001,
                      std_estimate = 0.28),
    sem_fixture_param("moral_condemnation", "~", "ideology",
                      estimate = 0.12, p_value = 0.04,
                      std_estimate = 0.10,
                      diagram_hidden = TRUE,
                      diagram_hidden_reason = "hidden_covariate"),
    sem_fixture_param("dehumanization", "~", "ideology",
                      estimate = 0.08, p_value = 0.09,
                      std_estimate = 0.07,
                      diagram_hidden = TRUE,
                      diagram_hidden_reason = "hidden_covariate")
  )

  diagram <- mellio:::ms_lavaan_structural_diagram(
    params,
    type = "lavaan_sem",
    type_label = "Structural equation model",
    n = 301
  )

  expect_equal(diagram$model_class, "observed_sem")
  expect_length(diagram$nodes, 3L)
  expect_false("ideology" %in% names(sem_fixture_nodes_by_id(diagram)))
  expect_length(diagram$edges, 2L)
  expect_length(diagram$omitted_paths, 2L)
  expect_equal(diagram$omitted_path_count, 2L)
  expect_true(all(vapply(diagram$omitted_paths, function(edge) {
    identical(edge$reason, "hidden_covariate") &&
      identical(edge$source, "ideology") &&
      isTRUE(edge$hidden)
  }, logical(1))))
  expect_length(diagram$hidden_covariates, 1L)
  expect_equal(diagram$hidden_covariates[[1]]$label, "ideology")
  expect_equal(diagram$hidden_covariates[[1]]$path_count, 2L)
  expect_equal(
    sort(unlist(diagram$hidden_covariates[[1]]$targets, use.names = FALSE)),
    c("dehumanization", "moral_condemnation")
  )
})

test_that("lavaan diagram omit helper marks path keys and variable names", {
  params <- list(
    sem_fixture_param("m", "~", "x", estimate = 0.40),
    sem_fixture_param("y", "~", "m", estimate = 0.30),
    sem_fixture_param("m", "~", "ideology", estimate = 0.10),
    sem_fixture_param("y", "~", "ideology", estimate = 0.08)
  )

  marked <- mellio:::ms_lavaan_mark_diagram_omissions(
    params,
    diagram_omit = c("y ~ m", "ideology")
  )

  hidden <- vapply(marked, function(p) isTRUE(p$diagram_hidden), logical(1))
  expect_equal(hidden, c(FALSE, TRUE, TRUE, TRUE))
  expect_equal(marked[[2]]$diagram_hidden_reason, "omitted_path")
  expect_equal(marked[[3]]$diagram_hidden_reason, "hidden_covariate")
  expect_equal(marked[[4]]$diagram_hidden_reason, "hidden_covariate")
})

test_that("lavaan structural diagram fixture declines multi-group diagrams", {
  params <- list(
    sem_fixture_param("visual", "=~", "x1", group = 1),
    sem_fixture_param("visual", "=~", "x2", group = 1),
    sem_fixture_param("visual", "=~", "x1", group = 2),
    sem_fixture_param("visual", "=~", "x2", group = 2)
  )

  expect_null(mellio:::ms_lavaan_structural_diagram(
    params,
    type = "lavaan_cfa",
    type_label = "Confirmatory factor analysis",
    n = 301
  ))
})

test_that("ms_result_id has the rs_ prefix and 8 base36 chars", {
  id <- mellio:::ms_result_id()
  expect_match(id, "^rs_[a-z0-9]{8}$")
  # Two consecutive calls differ
  expect_false(mellio:::ms_result_id() == mellio:::ms_result_id())
})

test_that("ms_now_iso8601 returns a valid ISO 8601 UTC timestamp", {
  t <- mellio:::ms_now_iso8601()
  expect_match(t, "^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}Z$")
})

test_that("ms_safe_numeric coerces non-finite to NA", {
  expect_equal(mellio:::ms_safe_numeric(3.14), 3.14)
  expect_true(is.na(mellio:::ms_safe_numeric(Inf)))
  expect_true(is.na(mellio:::ms_safe_numeric(-Inf)))
  expect_true(is.na(mellio:::ms_safe_numeric(NaN)))
  expect_true(is.na(mellio:::ms_safe_numeric(NULL)))
})

test_that("ms_base64url_encode round-trips", {
  for (s in c("hello", "Mellio", "x", "{\"a\":1,\"b\":2.5}", "")) {
    enc <- mellio:::ms_base64url_encode(s)
    # Convert to standard form for jsonlite to decode
    std <- chartr("-_", "+/", enc)
    pad <- (4 - nchar(std) %% 4) %% 4
    if (pad > 0) std <- paste0(std, strrep("=", pad))
    if (nchar(std) > 0) {
      back <- rawToChar(jsonlite::base64_dec(std))
      expect_equal(back, s)
    }
  }
})

test_that("ms_base64url_encode strips padding and uses URL-safe chars", {
  # 'hello' -> 'aGVsbG8' (no padding)
  expect_equal(mellio:::ms_base64url_encode("hello"), "aGVsbG8")
  # No + or / in any output
  for (s in c("???>>>", "<<???>>>")) {
    enc <- mellio:::ms_base64url_encode(s)
    expect_false(grepl("[+/=]", enc))
  }
})

test_that("editor URL base64 params are compact", {
  raw <- charToRaw(paste(rep("abcdef", 200), collapse = ""))
  enc <- mellio:::mellio_base64_compact(raw)
  expect_false(grepl("[\r\n]", enc))

  tbl <- melliotab(
    data.frame(label = paste0("row ", seq_len(80)), value = seq_len(80)),
    title = "Long URL Table"
  )
  tbl$provenance <- list(
    script = list(file = "/tmp/source-check.R", line = 42L),
    git = list(commit = paste(rep("a", 40), collapse = ""), dirty = TRUE)
  )
  url <- mellio:::build_edit_url(tbl, mode = "table")
  expect_false(grepl("[\r\n]", url))
})

# ── mellio_payload.default ───────────────────────────────────────────────

test_that("mellio_payload keeps guided errors for non-result inputs", {
  expect_error(
    mellio_payload("a string"),
    "not recognized"
  )
  # function / NULL still need guided errors rather than preserved print output.
  expect_error(mellio_payload(NULL), "Nothing to send")
  expect_error(mellio_payload(NULL), "printed output")
  err <- tryCatch(mellio_payload(NULL), error = conditionMessage)
  expect_false(grepl("mediation", err, fixed = TRUE))
})

test_that("mellio_payload.NULL guides printed correlation loops without mediation language", {
  call_txt <- "for (m in motivations) { test <- cor.test(data_clean$SES, data_clean[[m]]); cat(m) }"
  expect_error(mellio_payload(NULL, .call = call_txt), "data.frame")
  expect_error(mellio_payload(NULL, .call = call_txt), "several correlations")
  err <- tryCatch(mellio_payload(NULL, .call = call_txt), error = conditionMessage)
  expect_false(grepl("mediation", err, fixed = TRUE))
})

test_that("mellio_payload.NULL guides printed ANOVA loops toward objects", {
  call_txt <- "for (m in motivations) { fit <- aov(data_clean[[m]] ~ data_clean$SES_band); print(summary(fit)) }"
  expect_error(mellio_payload(NULL, .call = call_txt), "For ANOVA")
  expect_error(mellio_payload(NULL, .call = call_txt), "mellio_open\\(anova\\(fit\\)\\)")
  err <- tryCatch(mellio_payload(NULL, .call = call_txt), error = conditionMessage)
  expect_false(grepl("r = .12", err, fixed = TRUE))
  expect_false(grepl("mediation", err, fixed = TRUE))
})

# ── mellio_payload.htest ─────────────────────────────────────────────────

test_that("Welch t-test extracts to a v0.1 inline payload", {
  p <- mellio_payload(t.test(extra ~ group, data = sleep))

  expect_equal(p$schema_version, "0.1")
  expect_equal(p$card_kind, "inline")
  expect_equal(p$type, "welch_t_test")
  expect_match(p$type_label, "Welch")
  expect_match(p$call, "t\\.test\\(extra ~ group")

  # Statistic block
  expect_equal(p$fields$statistic$name, "t")
  expect_true(is.numeric(p$fields$statistic$value))
  expect_true(is.numeric(p$fields$statistic$df))

  # p-value
  expect_true(p$fields$p_value >= 0 && p$fields$p_value <= 1)

  # Estimate: mean diff with CI
  expect_equal(p$fields$estimate$name, "mean_diff")
  expect_length(p$fields$estimate$ci, 2)

  expect_equal(p$fields$alternative, "two.sided")
  expect_equal(p$fields$conf_level, 0.95)
  expect_equal(p$fields$outcome, "extra")
  expect_equal(p$fields$predictor, "group")
  expect_equal(p$fields$group, "group")
  expect_equal(p$fields$comparison, "1 vs. 2")

  expect_match(p$raw_output, "Welch Two Sample t-test")
})

test_that("Student t-test (var.equal) maps to students_t_test", {
  p <- mellio_payload(t.test(extra ~ group, data = sleep, var.equal = TRUE))
  expect_equal(p$type, "students_t_test")
})

test_that("Paired t-test maps correctly", {
  pre  <- sleep$extra[sleep$group == 1]
  post <- sleep$extra[sleep$group == 2]
  p <- mellio_payload(t.test(pre, post, paired = TRUE))
  expect_equal(p$type, "paired_t_test")
  expect_true(!is.null(p$fields$estimate))
  expect_equal(p$fields$comparison, "pre vs. post")
})

test_that("One-sample t-test maps correctly", {
  p <- mellio_payload(t.test(sleep$extra, mu = 0))
  expect_equal(p$type, "one_sample_t_test")
  expect_equal(p$fields$outcome, "extra")
})

test_that("Pearson correlation extracts r and CI", {
  p <- mellio_payload(cor.test(mtcars$mpg, mtcars$wt))
  expect_equal(p$type, "pearson_correlation")
  expect_equal(p$fields$statistic$name, "t")
  expect_equal(p$fields$estimate$name, "cor")
  expect_length(p$fields$estimate$ci, 2)
  expect_equal(p$fields$comparison, "mpg and wt")
})

test_that("Spearman correlation maps correctly", {
  p <- suppressWarnings(mellio_payload(cor.test(mtcars$mpg, mtcars$wt, method = "spearman")))
  expect_equal(p$type, "spearman_correlation")
  expect_equal(p$fields$p_value_method, "spearman_rank_test")
  expect_equal(p$fields$sample_size, 32)
})

test_that("Kendall correlation maps correctly", {
  p <- suppressWarnings(mellio_payload(cor.test(mtcars$mpg, mtcars$wt, method = "kendall")))
  expect_equal(p$type, "kendall_correlation")
  expect_equal(p$fields$p_value_method, "kendall_rank_test")
  expect_equal(p$fields$sample_size, 32)
})

test_that("Wilcoxon rank-sum maps correctly", {
  p <- suppressWarnings(mellio_payload(wilcox.test(extra ~ group, data = sleep)))
  expect_equal(p$type, "wilcoxon_rank_sum")
  expect_equal(p$fields$statistic$name, "W")
})

test_that("Wilcoxon signed-rank maps correctly", {
  pre  <- c(1, 2, 3, 4, 5)
  post <- c(2, 3, 5, 5, 7)
  p <- suppressWarnings(mellio_payload(wilcox.test(pre, post, paired = TRUE)))
  expect_equal(p$type, "wilcoxon_signed_rank")
  expect_equal(p$fields$statistic$name, "V")
})

test_that("Chi-squared test extracts statistic and df", {
  M <- matrix(c(10, 20, 30, 40), nrow = 2)
  p <- mellio_payload(chisq.test(M))
  expect_equal(p$type, "chi_squared_test")
  expect_true(is.numeric(p$fields$statistic$value))
  expect_true(is.numeric(p$fields$statistic$df))
  expect_equal(p$fields$test_context, "association")
  expect_equal(p$fields$table_type, "chi_square_contingency")
  expect_length(p$fields$rows, 2L)
  expect_equal(p$fields$rows[[1]]$col_1, 10)
})

test_that("Chi-squared goodness-of-fit is tagged separately from association", {
  p <- mellio_payload(chisq.test(c(12, 18, 20)))
  expect_equal(p$type, "chi_squared_test")
  expect_equal(p$fields$test_context, "goodness_of_fit")
  expect_null(p$fields$effect)
  expect_equal(p$fields$table_type, "chi_square_goodness_of_fit")
  expect_length(p$fields$rows, 3L)
  expect_equal(p$fields$rows[[1]]$observed, 12)
  expect_true(is.numeric(p$fields$rows[[1]]$expected))
})

test_that("Fisher exact test extracts odds ratio", {
  M <- matrix(c(10, 20, 30, 40), nrow = 2)
  p <- mellio_payload(fisher.test(M))
  expect_equal(p$type, "fisher_exact_test")
  expect_true(!is.null(p$fields$estimate))
  expect_equal(p$fields$estimate$name, "odds_ratio")
  expect_equal(p$fields$test_context, "association")
})

test_that("Fisher exact omits fields$statistic and surfaces odds_ratio + CI", {
  # Validator on the JS side rejects inline cards whose statistic.value
  # is NA — Fisher's test has no statistic, so we must omit the field
  # entirely. The OR + CI live on fields$estimate.
  p <- mellio_payload(fisher.test(matrix(c(8, 2, 1, 5), nrow = 2)))
  expect_null(p$fields[["statistic"]])
  expect_true(is.finite(p$fields$estimate$value))
  expect_length(p$fields$estimate$ci, 2L)
  expect_true(is.finite(p$fields$p_value))
})

test_that("Vector-input t-test omits fields$groups (no 'mean of x' leak)", {
  # t.test(x, y) names $estimate as c("mean of x", "mean of y") — these
  # are not real group labels and would render as outcome garbage.
  # Formula-call t-tests still emit labels (covered by Welch test above).
  p <- mellio_payload(t.test(c(1, 2, 3, 4, 5), c(1.1, 2.1, 3.1, 4.1, 5.1)))
  expect_equal(p$type, "welch_t_test")
  expect_null(p$fields$groups)
})

test_that("One-sample t-test surfaces null_value (mu) for the paragraph builder", {
  p0  <- mellio_payload(t.test(sleep$extra, mu = 0))
  expect_equal(p0$fields$null_value, 0)
  p25 <- mellio_payload(t.test(sleep$extra, mu = 2.5))
  expect_equal(p25$fields$null_value, 2.5)
})

# ── P4: .data threading (per-group n/M/SD + Cohen's d) ────────────────

test_that("Welch t-test with .data fills per-group n/M/SD aligned to test direction", {
  fit <- t.test(extra ~ group, data = sleep)
  p <- mellio_payload(fit, .data = sleep)

  expect_length(p$fields$groups, 2L)
  labels <- vapply(p$fields$groups, function(g) g$label, character(1))
  expect_setequal(labels, c("1", "2"))

  # Group ordering follows R's $estimate names (groups 1 then 2), so the
  # descriptives align with the mean-difference direction the test
  # reported. Cross-check against aggregate().
  agg <- stats::aggregate(extra ~ group, data = sleep,
                          FUN = function(z) c(n = length(z),
                                              m = mean(z),
                                              s = sd(z)))
  for (i in seq_along(p$fields$groups)) {
    g <- p$fields$groups[[i]]
    j <- which(as.character(agg$group) == g$label)
    expect_equal(g$n,    unname(agg$extra[j, "n"]))
    expect_equal(g$mean, unname(agg$extra[j, "m"]), tolerance = 1e-9)
    expect_equal(g$sd,   unname(agg$extra[j, "s"]), tolerance = 1e-9)
  }
})

test_that("Welch t-test with .data emits Cohen's d (average_group_sd method tag)", {
  fit <- t.test(extra ~ group, data = sleep)
  p <- mellio_payload(fit, .data = sleep)

  es <- p$fields$effect_size
  expect_equal(es$name, "cohens_d")
  expect_equal(es$method, "average_group_sd")
  expect_true(is.numeric(es$value) && is.finite(es$value))

  # Welch d uses the root mean square of the two group SDs.
  g1 <- sleep$extra[sleep$group == 1]
  g2 <- sleep$extra[sleep$group == 2]
  average_group_sd <- sqrt((var(g1) + var(g2)) / 2)
  expect_equal(es$value, (mean(g1) - mean(g2)) / average_group_sd,
               tolerance = 1e-9)
})

test_that("Welch t-test with .data offers observed group means figure", {
  fit <- t.test(extra ~ group, data = sleep)
  p <- mellio_payload(fit, .data = sleep)

  fig <- p$figure_data$adjusted_means
  expect_equal(fig$source, "htest_two_sample")
  expect_equal(fig$mean_kind, "observed")
  expect_equal(fig$outcome, "extra")
  expect_equal(fig$factor$variable, "group")
  expect_equal(fig$ci_method, "group_standard_error")
  expect_length(fig$groups, 2L)
  expect_true(all(vapply(fig$groups, function(g) is.numeric(g$se) && g$se > 0, logical(1))))
  expect_true(any(vapply(p$metadata$available_figures, function(f) {
    identical(f$type, "adjusted_means") && isTRUE(f$default)
  }, logical(1))))
})

test_that("Paired t-test with .data emits 'Differences' descriptives + dz", {
  df <- data.frame(pre = c(1.2, 2.4, 1.8, 3.1, 2.6),
                   post = c(1.8, 3.0, 2.1, 4.0, 3.4))
  fit <- t.test(df$pre, df$post, paired = TRUE)
  p <- mellio_payload(fit, .data = df)

  expect_length(p$fields$groups, 1L)
  expect_equal(p$fields$groups[[1]]$label, "Differences")
  expect_equal(p$fields$groups[[1]]$n, 5L)
  expect_length(p$fields$paired_measures, 2L)
  expect_equal(
    vapply(p$fields$paired_measures, function(row) row$label, character(1)),
    c("pre", "post")
  )
  expect_equal(
    vapply(p$fields$paired_measures, function(row) row$n, numeric(1)),
    c(5, 5)
  )
  diffs <- df$pre - df$post
  expect_equal(p$fields$groups[[1]]$mean, mean(diffs), tolerance = 1e-9)
  expect_equal(p$fields$groups[[1]]$sd,   sd(diffs),   tolerance = 1e-9)

  es <- p$fields$effect_size
  expect_equal(es$name, "cohens_dz")
  expect_equal(es$method, "paired_dz")
  expect_equal(es$value, mean(diffs) / sd(diffs), tolerance = 1e-9)

  fig <- p$figure_data$paired_difference_plot
  expect_equal(fig$source, "paired_t_test")
  expect_equal(fig$first_measure$label, "pre")
  expect_equal(fig$second_measure$label, "post")
  expect_equal(fig$n, 5L)
  expect_length(fig$pairs, 5L)
  expect_equal(fig$pairs[[1]]$first, df$pre[[1]], tolerance = 1e-9)
  expect_equal(fig$pairs[[1]]$second, df$post[[1]], tolerance = 1e-9)
  expect_equal(fig$difference$mean, mean(diffs), tolerance = 1e-9)
  expect_equal(fig$difference$ci_lower, unname(fit$conf.int)[[1]], tolerance = 1e-9)
  expect_true(any(vapply(p$metadata$available_figures, function(f) {
    identical(f$type, "paired_difference_plot") && isTRUE(f$default)
  }, logical(1))))
})

test_that("One-sample t-test with .data emits sample descriptives + d against mu", {
  df <- data.frame(extra = sleep$extra[sleep$group == 1])
  fit <- t.test(df$extra, mu = 0)
  p <- mellio_payload(fit, .data = df)

  expect_length(p$fields$groups, 1L)
  expect_equal(p$fields$groups[[1]]$label, "extra")
  expect_equal(p$fields$groups[[1]]$n, length(df$extra))
  expect_equal(p$fields$groups[[1]]$mean, mean(df$extra), tolerance = 1e-9)
  expect_equal(p$fields$groups[[1]]$sd,   sd(df$extra),   tolerance = 1e-9)

  es <- p$fields$effect_size
  expect_equal(es$name, "cohens_d")
  expect_equal(es$method, "one_sample")
  expect_equal(es$value, mean(df$extra) / sd(df$extra), tolerance = 1e-9)
})

test_that("Paired t-test with .data counts complete pairs", {
  df <- data.frame(
    after = c(10, 13, NA, 15, 18, 21),
    before = c(9, NA, 11, 13, 15, 18)
  )
  fit <- t.test(df$after, df$before, paired = TRUE)
  p <- mellio_payload(fit, .data = df)

  diffs <- df$after - df$before
  diffs <- diffs[is.finite(diffs)]

  expect_equal(p$type, "paired_t_test")
  expect_equal(p$fields$sample_size, length(diffs))
  expect_equal(p$fields$statistic$df, length(diffs) - 1)
  expect_equal(p$fields$groups[[1]]$n, length(diffs))
  expect_length(p$fields$paired_measures, 2L)
  expect_equal(
    vapply(p$fields$paired_measures, function(row) row$n, numeric(1)),
    rep(length(diffs), 2)
  )
  expect_equal(p$fields$paired_measures[[1]]$mean, mean(df$after[is.finite(df$after - df$before)]), tolerance = 1e-9)
  expect_equal(p$fields$paired_measures[[2]]$mean, mean(df$before[is.finite(df$after - df$before)]), tolerance = 1e-9)
  expect_equal(p$fields$groups[[1]]$mean, mean(diffs), tolerance = 1e-9)
  expect_equal(p$fields$groups[[1]]$sd, sd(diffs), tolerance = 1e-9)
  expect_equal(p$fields$effect_size$name, "cohens_dz")
  expect_equal(p$fields$effect_size$value, mean(diffs) / sd(diffs), tolerance = 1e-9)
  expect_length(p$figure_data$paired_difference_plot$pairs, length(diffs))
})

test_that("One-sample t-test with .data keeps null mean and emits mean figure", {
  df <- data.frame(score = c(49, 51, NA, 60, 55, 47, 53))
  fit <- t.test(df$score, mu = 50)
  p <- mellio_payload(fit, .data = df)

  y <- df$score[is.finite(df$score)]

  expect_equal(p$type, "one_sample_t_test")
  expect_equal(p$fields$null_value, 50)
  expect_equal(p$fields$sample_size, length(y))
  expect_equal(p$fields$statistic$df, length(y) - 1)
  expect_equal(p$fields$groups[[1]]$label, "score")
  expect_equal(p$fields$groups[[1]]$n, length(y))
  expect_equal(p$fields$groups[[1]]$mean, mean(y), tolerance = 1e-9)
  expect_equal(p$fields$groups[[1]]$sd, sd(y), tolerance = 1e-9)
  expect_equal(p$fields$effect_size$name, "cohens_d")
  expect_equal(p$fields$effect_size$method, "one_sample")
  expect_equal(p$fields$effect_size$value, (mean(y) - 50) / sd(y), tolerance = 1e-9)
  fig <- p$figure_data$one_sample_mean_plot
  expect_equal(fig$source, "one_sample_t_test")
  expect_equal(fig$sample$label, "score")
  expect_equal(fig$n, length(y))
  expect_equal(fig$mean, mean(y), tolerance = 1e-9)
  expect_equal(fig$sd, sd(y), tolerance = 1e-9)
  expect_equal(fig$null_value, 50)
  expect_equal(fig$ci_lower, unname(fit$conf.int)[[1]], tolerance = 1e-9)
  expect_equal(fig$ci_upper, unname(fit$conf.int)[[2]], tolerance = 1e-9)
  expect_true(any(vapply(p$metadata$available_figures, function(f) {
    identical(f$type, "one_sample_mean_plot") && isTRUE(f$default)
  }, logical(1))))
})

test_that("One-sample t-test auto-enriches from simple caller references", {
  d <- data.frame(score = c(49, 51, NA, 60, 55, 47, 53))
  result <- t.test(d$score, mu = 50)
  p <- mellio_payload(result)

  y <- d$score[is.finite(d$score)]

  expect_equal(p$type, "one_sample_t_test")
  expect_equal(p$fields$groups[[1]]$label, "score")
  expect_equal(p$fields$groups[[1]]$n, length(y))
  expect_equal(p$fields$groups[[1]]$mean, mean(y), tolerance = 1e-9)
  expect_equal(p$fields$groups[[1]]$sd, sd(y), tolerance = 1e-9)
  expect_equal(p$fields$effect_size$name, "cohens_d")
  expect_equal(p$fields$effect_size$value, (mean(y) - 50) / sd(y), tolerance = 1e-9)
  expect_equal(p$figure_data$one_sample_mean_plot$sample$label, "score")
  expect_equal(p$figure_data$one_sample_mean_plot$null_value, 50)
})

test_that("Paired t-test auto-enriches from simple caller references", {
  d <- data.frame(
    after = c(10, 13, NA, 15, 18, 21),
    before = c(9, NA, 11, 13, 15, 18)
  )
  result <- t.test(d$after, d$before, paired = TRUE)
  p <- mellio_payload(result)

  diffs <- d$after - d$before
  diffs <- diffs[is.finite(diffs)]

  expect_equal(p$type, "paired_t_test")
  expect_equal(p$fields$groups[[1]]$label, "Differences")
  expect_equal(p$fields$groups[[1]]$n, length(diffs))
  expect_length(p$fields$paired_measures, 2L)
  expect_equal(
    vapply(p$fields$paired_measures, function(row) row$label, character(1)),
    c("after", "before")
  )
  expect_equal(p$fields$groups[[1]]$mean, mean(diffs), tolerance = 1e-9)
  expect_equal(p$fields$groups[[1]]$sd, sd(diffs), tolerance = 1e-9)
  expect_equal(p$fields$effect_size$name, "cohens_dz")
  expect_equal(p$fields$effect_size$value, mean(diffs) / sd(diffs), tolerance = 1e-9)
})

test_that("Welch formula t-test auto-enriches from one matching data frame", {
  sleep_df <- datasets::sleep
  result <- t.test(extra ~ group, data = sleep_df)
  p <- mellio_payload(result)

  expect_equal(p$type, "welch_t_test")
  expect_length(p$fields$groups, 2L)
  expect_equal(
    vapply(p$fields$groups, function(row) row$label, character(1)),
    c("1", "2")
  )
  expect_equal(
    vapply(p$fields$groups, function(row) row$n, numeric(1)),
    c(10, 10)
  )
  expect_equal(p$fields$groups[[1]]$mean, mean(sleep_df$extra[sleep_df$group == 1]), tolerance = 1e-9)
  expect_equal(p$fields$groups[[2]]$mean, mean(sleep_df$extra[sleep_df$group == 2]), tolerance = 1e-9)
  expect_equal(p$fields$effect_size$name, "cohens_d")
})

test_that("Welch formula t-test chooses the matching data frame among candidates", {
  sleep_a <- datasets::sleep
  sleep_b <- transform(datasets::sleep, extra = extra + 100)
  result <- t.test(extra ~ group, data = sleep_a)
  p <- mellio_payload(result)

  expect_equal(
    vapply(p$fields$groups, function(row) row$n, numeric(1)),
    c(10, 10)
  )
  expect_equal(p$fields$groups[[1]]$mean, mean(sleep_a$extra[sleep_a$group == 1]), tolerance = 1e-9)
  expect_equal(p$fields$groups[[2]]$mean, mean(sleep_a$extra[sleep_a$group == 2]), tolerance = 1e-9)
  expect_equal(p$fields$groups[[1]]$sd, sd(sleep_a$extra[sleep_a$group == 1]), tolerance = 1e-9)
  expect_equal(p$fields$groups[[2]]$sd, sd(sleep_a$extra[sleep_a$group == 2]), tolerance = 1e-9)
  expect_equal(p$fields$effect_size$name, "cohens_d")
})

test_that("mellio_open htest auto-enriches from the user's environment", {
  d <- data.frame(score = c(49, 51, NA, 60, 55, 47, 53))
  result <- t.test(d$score, mu = 50)
  url <- mellio_open(result, browse = FALSE)

  b64 <- sub(".*payload=", "", url)
  b64 <- sub("&.*$", "", b64)
  std <- chartr("-_", "+/", b64)
  pad <- (4 - nchar(std) %% 4) %% 4
  if (pad > 0) std <- paste0(std, strrep("=", pad))
  parsed <- jsonlite::fromJSON(rawToChar(jsonlite::base64_dec(std)),
                               simplifyVector = FALSE)

  y <- d$score[is.finite(d$score)]

  expect_equal(parsed$type, "one_sample_t_test")
  expect_equal(parsed$fields$groups[[1]]$label, "score")
  expect_equal(parsed$fields$groups[[1]]$n, length(y))
  expect_equal(parsed$fields$effect_size$name, "cohens_d")
})

test_that("htest without recoverable data retains group means only", {
  result <- t.test(extra ~ group, data = datasets::sleep)
  p <- mellio_payload(result, .env = new.env(parent = emptyenv()))

  expect_equal(
    vapply(p$fields$groups, function(row) row$mean, numeric(1)),
    unname(result$estimate),
    tolerance = 1e-9
  )
  for (g in p$fields$groups) {
    expect_null(g$n)
    expect_null(g$sd)
  }
  expect_null(p$fields$effect_size)
  expect_null(p$figure_data$adjusted_means)
})

# ── mellio_payload.lm ────────────────────────────────────────────────────

test_that("lm extracts model summary: F, df1/df2, p, R²", {
  m <- lm(mpg ~ wt + cyl, data = mtcars)
  p <- mellio_payload(m)

  expect_equal(p$schema_version, "0.1")
  expect_equal(p$card_kind, "inline")
  expect_equal(p$type, "lm_model_summary")
  expect_equal(p$type_label, "Linear Regression")

  expect_equal(p$fields$statistic$name, "F")
  expect_true(is.numeric(p$fields$statistic$value))
  expect_length(p$fields$statistic$df, 2)
  expect_equal(p$fields$statistic$df[1], 2)
  expect_equal(p$fields$statistic$df[2], 29)

  expect_true(p$fields$p_value < 1e-10)

  expect_equal(p$fields$estimate$name, "R\u00b2")
  expect_true(p$fields$estimate$value > 0 && p$fields$estimate$value < 1)
  expect_true(is.numeric(p$fields$r_squared))
  expect_true(is.numeric(p$fields$adj_r_squared))
  expect_true(is.numeric(p$fields$aic))
  expect_true(is.numeric(p$fields$bic))
  expect_true(is.numeric(p$fields$logLik))
  expect_equal(p$fields$conf_level, 0.95)
  expect_equal(p$fields$coefficient_ci_method, "wald_t")

  expect_equal(p$fields$n, 32)
  expect_match(p$raw_output, "Coefficients")
})

test_that("lm emits a per-coefficient list with B/SE/t/p/CI", {
  m <- lm(mpg ~ wt + cyl, data = mtcars)
  p <- mellio_payload(m)

  coefs <- p$fields$coefficients
  expect_true(is.list(coefs))
  expect_length(coefs, 3L)   # (Intercept) + wt + cyl

  terms <- vapply(coefs, function(c) c$term, character(1))
  expect_setequal(terms, c("(Intercept)", "wt", "cyl"))

  # Every row has the full set of numeric fields
  for (c in coefs) {
    expect_equal(c$estimate_name, "B")
    expect_true(is.numeric(c$estimate))
    expect_true(is.numeric(c$std_error))
    expect_true(is.numeric(c$statistic))
    expect_true(is.numeric(c$p_value))
    expect_true(is.numeric(c$ci_lower))
    expect_true(is.numeric(c$ci_upper))
    expect_equal(c$ci_method, "wald_t")
    expect_true(c$ci_lower < c$ci_upper)
  }

  # Cross-check one row against summary(lm) directly
  s <- summary(m)
  wt_row <- coefs[[which(terms == "wt")]]
  expect_equal(wt_row$estimate,  unname(s$coefficients["wt", "Estimate"]),   tolerance = 1e-6)
  expect_equal(wt_row$std_error, unname(s$coefficients["wt", "Std. Error"]), tolerance = 1e-6)
})

test_that("lm advertises coefficient plot figure data", {
  m <- lm(mpg ~ wt + cyl, data = mtcars)
  p <- mellio_payload(m)

  expect_equal(p$metadata$available_figures[[1]]$type, "coefficient_plot")
  expect_equal(p$metadata$available_figures[[1]]$label, "Coefficient plot")
  expect_true(isTRUE(p$metadata$available_figures[[1]]$default))

  fig <- p$figure_data$coefficient_plot
  expect_equal(fig$estimate_label, "B")
  expect_equal(fig$coefficient_scale, "raw")
  expect_false(fig$default_show_intercept)
  expect_length(fig$coefficients, 3L)
  expect_setequal(
    vapply(fig$coefficients, function(row) row$term, character(1)),
    c("(Intercept)", "wt", "cyl")
  )
  wt <- fig$coefficients[[which(vapply(fig$coefficients, function(row) {
    identical(row$term, "wt")
  }, logical(1)))[1]]]
  expect_true(is.numeric(wt$std_estimate))
  expect_true(is.numeric(wt$std_std_error))
  expect_true(is.numeric(wt$std_ci_lower))
  expect_true(is.numeric(wt$std_ci_upper))
})

test_that("lm advertises cautious two-way interaction plot data", {
  m <- lm(mpg ~ wt * am + hp, data = mtcars)
  p <- mellio_payload(m)

  types <- vapply(p$metadata$available_figures, function(row) row$type, character(1))
  expect_equal(types[[1]], "interaction_plot")
  expect_true(isTRUE(p$metadata$available_figures[[1]]$default))
  expect_true("coefficient_plot" %in% types)

  fig <- p$figure_data$interaction_plot
  expect_equal(fig$interaction_term, "wt:am")
  expect_equal(fig$x$variable, "wt")
  expect_equal(fig$moderator$variable, "am")
  expect_equal(fig$outcome, "mpg")
  expect_equal(fig$scale, "response")
  expect_false(fig$bounded_response)
  expect_length(fig$moderator$levels, 2L)
  expect_length(fig$grid, 160L)
  expect_true(all(vapply(fig$grid[1:4], function(row) {
    is.numeric(row$estimate) && is.numeric(row$ci_lower) && is.numeric(row$ci_upper)
  }, logical(1))))
  expect_length(fig$held_constant, 1L)
  expect_equal(fig$held_constant[[1]]$variable, "hp")

  test <- p$fields$interaction_tests[[1]]
  drop <- drop1(m, test = "F")
  expect_equal(test$term, "wt:am")
  expect_equal(test$statistic$name, "F")
  expect_equal(test$statistic$df, c(1, stats::df.residual(m)))
  expect_equal(test$statistic$value, unname(drop["wt:am", "F value"]), tolerance = 1e-6)
  expect_equal(test$p_value, unname(drop["wt:am", "Pr(>F)"]), tolerance = 1e-6)
  expect_equal(test$effect$name, "eta_sq_partial")
})

test_that("lm advertises categorical-by-categorical interaction plot data", {
  d <- ToothGrowth
  d$dose_f <- factor(d$dose)
  m <- lm(len ~ supp * dose_f, data = d)
  p <- mellio_payload(m)

  types <- vapply(p$metadata$available_figures, function(row) row$type, character(1))
  expect_equal(types[[1]], "interaction_plot")
  expect_true(isTRUE(p$metadata$available_figures[[1]]$default))
  expect_true("coefficient_plot" %in% types)

  fig <- p$figure_data$interaction_plot
  expect_equal(fig$interaction_term, "supp:dose_f")
  expect_equal(fig$interaction_kind, "categorical_by_categorical")
  expect_equal(fig$x$variable, "supp")
  expect_equal(fig$moderator$variable, "dose_f")
  expect_equal(fig$x$type, "categorical")
  expect_equal(fig$moderator$type, "categorical")
  expect_length(fig$x$levels, 2L)
  expect_length(fig$moderator$levels, 3L)
  expect_length(fig$grid, 6L)
  expect_true(all(vapply(fig$grid, function(row) {
    !is.null(row$x_value) && !is.null(row$x_label) &&
      is.numeric(row$estimate) && is.numeric(row$ci_lower) &&
      is.numeric(row$ci_upper)
  }, logical(1))))

  test <- p$fields$interaction_tests[[1]]
  drop <- drop1(m, test = "F")
  expect_equal(test$term, "supp:dose_f")
  expect_equal(test$statistic$name, "F")
  expect_equal(test$statistic$df, c(2, stats::df.residual(m)))
  expect_equal(test$statistic$value, unname(drop["supp:dose_f", "F value"]), tolerance = 1e-6)
  expect_equal(test$p_value, unname(drop["supp:dose_f", "Pr(>F)"]), tolerance = 1e-6)
  expect_equal(test$effect$name, "eta_sq_partial")
})

test_that("lm exposes standardized omnibus model-term tests", {
  d <- ToothGrowth
  d$dose_f <- factor(d$dose)
  m <- lm(len ~ supp * dose_f, data = d)
  p <- mellio_payload(m)

  tests <- p$fields$model_term_tests
  expect_true(is.list(tests))
  term_names <- vapply(tests, function(row) row$term, character(1))
  expect_equal(term_names, c("supp", "dose_f", "supp:dose_f"))

  by_term <- setNames(tests, term_names)
  seq_tbl <- anova(m)
  drop_tbl <- drop1(m, test = "F")

  expect_equal(by_term$supp$method, "anova_f_sequential")
  expect_equal(by_term$supp$statistic$name, "F")
  expect_equal(by_term$supp$statistic$value, unname(seq_tbl["supp", "F value"]), tolerance = 1e-6)
  expect_equal(by_term$supp$ss_type, "type_i_sequential")

  expect_equal(by_term$dose_f$method, "anova_f_sequential")
  expect_equal(by_term$dose_f$statistic$value, unname(seq_tbl["dose_f", "F value"]), tolerance = 1e-6)

  expect_equal(by_term[["supp:dose_f"]]$method, "drop1_f")
  expect_equal(by_term[["supp:dose_f"]]$term_type, "interaction")
  expect_equal(by_term[["supp:dose_f"]]$label, "supp \u00d7 dose f")
  expect_equal(by_term[["supp:dose_f"]]$statistic$df, c(2, stats::df.residual(m)))
  expect_equal(by_term[["supp:dose_f"]]$statistic$value,
               unname(drop_tbl["supp:dose_f", "F value"]), tolerance = 1e-6)
  expect_equal(by_term[["supp:dose_f"]]$test_type, "omnibus")
})

test_that("glm advertises response-scale interaction plot data", {
  set.seed(42)
  n <- 80L
  d <- data.frame(
    y = rbinom(n, 1, 0.5),
    x = rep(seq(-1.5, 1.5, length.out = n / 2L), 2L),
    g = factor(rep(c("control", "treatment"), each = n / 2L))
  )
  m <- glm(y ~ x * g, data = d, family = binomial)
  p <- mellio_payload(m)

  fig <- p$figure_data$interaction_plot
  expect_equal(fig$x$variable, "x")
  expect_equal(fig$moderator$variable, "g")
  expect_true(fig$bounded_response)
  expect_equal(fig$y_label, "Predicted probability")
  expect_equal(fig$model_family, "binomial")
  expect_equal(fig$model_link, "logit")
  estimates <- vapply(fig$grid, function(row) row$estimate, numeric(1))
  expect_true(all(estimates >= 0 & estimates <= 1))
})

test_that("glm exposes standardized chi-square model-term tests", {
  set.seed(2026)
  n <- 160L
  d <- data.frame(
    y = rbinom(n, 1, 0.5),
    g = factor(rep(c("control", "treatment"), each = n / 2L)),
    h = factor(rep(c("easy", "hard"), times = n / 2L)),
    x = rnorm(n)
  )
  m <- glm(y ~ g * h + x, data = d, family = binomial())
  p <- mellio_payload(m)

  tests <- p$fields$model_term_tests
  expect_true(is.list(tests))
  term_names <- vapply(tests, function(row) row$term, character(1))
  expect_equal(term_names, c("g", "h", "x", "g:h"))

  by_term <- setNames(tests, term_names)
  seq_tbl <- anova(m, test = "Chisq")
  drop_tbl <- drop1(m, test = "Chisq")

  expect_equal(by_term$g$method, "anova_chisq_sequential")
  expect_equal(by_term$g$statistic$name, "chi2")
  expect_equal(by_term$g$statistic$df, 1)
  expect_equal(by_term$g$statistic$value, unname(seq_tbl["g", "Deviance"]), tolerance = 1e-6)
  expect_equal(by_term$g$ss_type, "type_i_sequential")

  expect_equal(by_term$x$method, "drop1_chisq")
  expect_equal(by_term$x$statistic$value, unname(drop_tbl["x", "LRT"]), tolerance = 1e-6)

  expect_equal(by_term[["g:h"]]$method, "drop1_chisq")
  expect_equal(by_term[["g:h"]]$term_type, "interaction")
  expect_equal(by_term[["g:h"]]$label, "g \u00d7 h")
  expect_equal(by_term[["g:h"]]$statistic$value, unname(drop_tbl["g:h", "LRT"]), tolerance = 1e-6)
  expect_equal(by_term[["g:h"]]$test_type, "omnibus")
})

test_that("lm advertises continuous-by-continuous simple-slope interaction data", {
  m <- lm(mpg ~ wt * hp + am, data = mtcars)
  p <- mellio_payload(m)

  types <- vapply(p$metadata$available_figures, function(row) row$type, character(1))
  expect_true("interaction_plot" %in% types)

  fig <- p$figure_data$interaction_plot
  expect_equal(fig$interaction_term, "wt:hp")
  expect_equal(fig$interaction_kind, "continuous_by_continuous")
  expect_equal(fig$x$variable, "wt")
  expect_equal(fig$moderator$variable, "hp")
  expect_equal(fig$moderator$type, "numeric")
  expect_length(fig$moderator$levels, 3L)
  expect_equal(
    vapply(fig$moderator$levels, function(level) level$slice, character(1)),
    c("low", "mean", "high")
  )
  expect_true(all(grepl("^(Low|Mean|High) \\(", vapply(
    fig$moderator$levels,
    function(level) level$label,
    character(1)
  ))))
  expect_length(fig$grid, 240L)
  expect_length(unique(vapply(fig$grid, function(row) row$moderator_value, character(1))), 3L)
  expect_equal(fig$moderator_default_preset, "sd")
  expect_length(fig$moderator_value_presets, 3L)
  expect_equal(
    vapply(fig$moderator_value_presets, function(preset) preset$id, character(1)),
    c("sd", "quartiles", "percentiles")
  )
  expect_true(all(vapply(fig$moderator_value_presets, function(preset) {
    length(preset$levels) == 3L && length(preset$grid) == 240L
  }, logical(1))))
  hp <- stats::model.frame(m)$hp
  sd_values <- as.numeric(vapply(fig$moderator_value_presets[[1]]$levels, function(level) {
    level$value
  }, character(1)))
  expect_equal(sd_values, c(mean(hp) - stats::sd(hp), mean(hp), mean(hp) + stats::sd(hp)),
               tolerance = 1e-8)
  expect_equal(
    vapply(fig$moderator_value_presets[[2]]$levels, function(level) level$slice, character(1)),
    c("q1", "median", "q3")
  )
  expect_equal(
    vapply(fig$moderator_value_presets[[3]]$levels, function(level) level$slice, character(1)),
    c("p10", "p50", "p90")
  )
  expect_length(fig$held_constant, 1L)
  expect_equal(fig$held_constant[[1]]$variable, "am")
})

test_that("lm advertises continuous-by-categorical-by-categorical interaction data", {
  m <- lm(mpg ~ wt * factor(am) * factor(vs) + hp, data = mtcars)
  p <- mellio_payload(m)

  types <- vapply(p$metadata$available_figures, function(row) row$type, character(1))
  expect_true("interaction_plot" %in% types)

  fig <- p$figure_data$interaction_plot
  expect_equal(fig$interaction_term, "wt:factor(am):factor(vs)")
  expect_equal(fig$interaction_kind, "continuous_by_categorical_by_categorical")
  expect_equal(fig$x$variable, "wt")
  expect_equal(fig$moderator$variable, "am")
  expect_equal(fig$facet$variable, "vs")
  expect_equal(fig$moderator$type, "categorical")
  expect_equal(fig$facet$type, "categorical")
  expect_length(fig$moderator$levels, 2L)
  expect_length(fig$facet$levels, 2L)
  expect_length(fig$grid, 320L)
  expect_true(all(vapply(fig$grid, function(row) {
    !is.null(row$facet_value) && !is.null(row$facet_label)
  }, logical(1))))
  expect_length(unique(vapply(fig$grid, function(row) row$moderator_value, character(1))), 2L)
  expect_length(unique(vapply(fig$grid, function(row) row$facet_value, character(1))), 2L)
  expect_null(fig$moderator_value_presets)
  expect_length(fig$held_constant, 1L)
  expect_equal(fig$held_constant[[1]]$variable, "hp")
})

test_that("lm advertises continuous-by-continuous-by-categorical interaction data", {
  m <- lm(mpg ~ wt * hp * factor(am) + vs, data = mtcars)
  p <- mellio_payload(m)

  fig <- p$figure_data$interaction_plot
  expect_equal(fig$interaction_term, "wt:hp:factor(am)")
  expect_equal(fig$interaction_kind, "continuous_by_continuous_by_categorical")
  expect_equal(fig$x$variable, "wt")
  expect_equal(fig$moderator$variable, "hp")
  expect_equal(fig$facet$variable, "am")
  expect_equal(fig$moderator$type, "numeric")
  expect_equal(fig$facet$type, "categorical")
  expect_equal(fig$moderator_default_preset, "sd")
  expect_length(fig$moderator_value_presets, 3L)
  expect_true(all(vapply(fig$moderator_value_presets, function(preset) {
    length(preset$levels) == 3L && length(preset$grid) == 480L &&
      all(vapply(preset$grid, function(row) !is.null(row$facet_value), logical(1)))
  }, logical(1))))
  hp <- stats::model.frame(m)$hp
  sd_values <- as.numeric(vapply(fig$moderator_value_presets[[1]]$levels, function(level) {
    level$value
  }, character(1)))
  expect_equal(sd_values, c(mean(hp) - stats::sd(hp), mean(hp), mean(hp) + stats::sd(hp)),
               tolerance = 1e-8)
  expect_length(fig$held_constant, 1L)
  expect_equal(fig$held_constant[[1]]$variable, "vs")
})

test_that("lm advertises continuous-by-continuous-by-continuous interaction data", {
  m <- lm(mpg ~ wt * hp * disp + am, data = mtcars)
  p <- mellio_payload(m)

  fig <- p$figure_data$interaction_plot
  expect_equal(fig$interaction_term, "wt:hp:disp")
  expect_equal(fig$interaction_kind, "continuous_by_continuous_by_continuous")
  expect_equal(fig$x$variable, "wt")
  expect_equal(fig$moderator$variable, "hp")
  expect_equal(fig$facet$variable, "disp")
  expect_equal(fig$moderator$type, "numeric")
  expect_equal(fig$facet$type, "numeric")
  expect_equal(fig$moderator_default_preset, "sd")
  expect_length(fig$moderator$levels, 3L)
  expect_length(fig$facet$levels, 3L)
  expect_length(fig$grid, 720L)
  expect_length(fig$moderator_value_presets, 3L)
  expect_true(all(vapply(fig$moderator_value_presets, function(preset) {
    length(preset$levels) == 3L && length(preset$facet_levels) == 3L &&
      length(preset$grid) == 720L &&
      all(vapply(preset$grid, function(row) !is.null(row$facet_value), logical(1)))
  }, logical(1))))
  hp <- stats::model.frame(m)$hp
  disp <- stats::model.frame(m)$disp
  sd_values <- as.numeric(vapply(fig$moderator_value_presets[[1]]$levels, function(level) {
    level$value
  }, character(1)))
  facet_sd_values <- as.numeric(vapply(fig$moderator_value_presets[[1]]$facet_levels, function(level) {
    level$value
  }, character(1)))
  expect_equal(sd_values, c(mean(hp) - stats::sd(hp), mean(hp), mean(hp) + stats::sd(hp)),
               tolerance = 1e-8)
  expect_equal(facet_sd_values,
               c(mean(disp) - stats::sd(disp), mean(disp), mean(disp) + stats::sd(disp)),
               tolerance = 1e-8)
  expect_length(fig$held_constant, 1L)
  expect_equal(fig$held_constant[[1]]$variable, "am")
})

test_that("lm advertises categorical-by-categorical-by-categorical interaction data", {
  set.seed(42)
  df <- expand.grid(
    a = factor(c("A1", "A2", "A3")),
    b = factor(c("B1", "B2")),
    c = factor(c("C1", "C2")),
    rep = seq_len(3)
  )
  df$z <- seq_len(nrow(df)) / 10
  df$y <- 10 +
    as.numeric(df$a) * 1.2 +
    as.numeric(df$b) * 0.9 +
    as.numeric(df$c) * -0.8 +
    as.numeric(df$a) * as.numeric(df$b) * as.numeric(df$c) * 0.35 +
    df$z * 0.4 +
    stats::rnorm(nrow(df), sd = 0.2)
  m <- lm(y ~ a * b * c + z, data = df)
  p <- mellio_payload(m)

  fig <- p$figure_data$interaction_plot
  expect_equal(fig$interaction_term, "a:b:c")
  expect_equal(fig$interaction_kind, "categorical_by_categorical_by_categorical")
  expect_equal(fig$x$variable, "a")
  expect_equal(fig$moderator$variable, "b")
  expect_equal(fig$facet$variable, "c")
  expect_equal(fig$x$type, "categorical")
  expect_equal(fig$moderator$type, "categorical")
  expect_equal(fig$facet$type, "categorical")
  expect_length(fig$x$levels, 3L)
  expect_length(fig$moderator$levels, 2L)
  expect_length(fig$facet$levels, 2L)
  expect_length(fig$grid, 12L)
  expect_true(all(vapply(fig$grid, function(row) {
    !is.null(row$x_value) && !is.null(row$x_label) &&
      !is.null(row$facet_value) && !is.na(row$estimate)
  }, logical(1))))
  expect_length(fig$held_constant, 1L)
  expect_equal(fig$held_constant[[1]]$variable, "z")
})

test_that("interaction plots stay absent for unsupported interaction shapes", {
  four_way <- mellio_payload(lm(mpg ~ wt * hp * disp * qsec, data = mtcars))
  expect_null(four_way$figure_data$interaction_plot)
})

test_that("unsupported figure payload sections stay absent", {
  p <- mellio_payload(t.test(extra ~ group, data = sleep))
  expect_null(p$figure_data)
  expect_null(p$metadata$available_figures)
})

test_that("lm with a single predictor still emits coefficients", {
  m <- lm(mpg ~ wt, data = mtcars)
  p <- mellio_payload(m)
  expect_length(p$fields$coefficients, 2L)
})

test_that("rank-deficient lm degrades coefficient extraction gracefully", {
  # Collinear predictor — lm drops the redundant column; the remaining
  # coefficients should still extract cleanly.
  m <- suppressWarnings(lm(mpg ~ I(wt + 0) + wt, data = mtcars))
  p <- suppressWarnings(mellio_payload(m))
  expect_true(length(p$fields$coefficients) >= 2L)
})

test_that("lm keeps factor-plus-covariate models as linear regression", {
  df <- data.frame(
    y = c(3.1, 3.4, 3.8, 4.2, 4.0, 4.6, 4.9, 5.1),
    condition = factor(rep(c("control", "treatment"), each = 4)),
    age = c(21, 23, 24, 26, 22, 25, 27, 28)
  )
  p <- mellio_payload(lm(y ~ condition + age, data = df))

  expect_null(p$fields$model_kind)
  expect_equal(p$type_label, "Linear Regression")
  expect_equal(p$fields$outcome, "y")
  expect_equal(p$fields$focal_terms, c("condition", "age"))
  expect_null(p$fields$control_terms)

  roles <- setNames(
    vapply(p$fields$terms, function(t) t$role, character(1)),
    vapply(p$fields$terms, function(t) t$name, character(1))
  )
  expect_equal(roles[["condition"]], "focal")
  expect_equal(roles[["age"]], "focal")

  coefs <- p$fields$coefficients
  term_sources <- vapply(coefs, function(c) {
    if (is.null(c$term_source)) "" else c$term_source
  }, character(1))
  treatment <- coefs[[which(term_sources == "condition")]]
  expect_equal(treatment$term_source, "condition")
})

test_that("lm supports explicit focal and control terms", {
  p <- mellio_payload(lm(mpg ~ wt + cyl + hp, data = mtcars),
                  focal = "wt", controls = c("cyl", "hp"))

  expect_equal(p$fields$model_kind, "controlled_regression")
  expect_equal(p$fields$focal_terms, "wt")
  expect_setequal(p$fields$control_terms, c("cyl", "hp"))
})

test_that("lm(y ~ factor(numeric) + numeric) remains linear regression", {
  # mtcars$cyl is numeric; factor(cyl) appears in model.frame under the
  # column name "factor(cyl)", not "cyl" — a previous bug in
  # ms_model_term_type missed this and labelled the term "other".
  p <- mellio_payload(lm(mpg ~ factor(cyl) + wt, data = mtcars))
  expect_null(p$fields$model_kind)
  expect_equal(p$type_label, "Linear Regression")
})

test_that("lm payload titles regression dispatches as linear regression", {
  p_factor <- mellio_payload(lm(mpg ~ factor(cyl) + wt, data = mtcars))
  expect_equal(p_factor$type_label, "Linear Regression")

  p_reg <- mellio_payload(lm(mpg ~ wt + hp, data = mtcars))
  expect_equal(p_reg$type_label, "Linear Regression")
})

test_that("lm(y ~ numeric + numeric) stays a plain regression", {
  p <- mellio_payload(lm(mpg ~ wt + hp, data = mtcars))
  expect_null(p$fields$model_kind)
  expect_equal(p$fields$focal_terms, c("wt", "hp"))
  expect_null(p$fields$control_terms)
  expect_equal(
    vapply(p$fields$terms, function(term) term$role, character(1)),
    c("focal", "focal")
  )
})

test_that("lm factor casts do not promote the dispatch to ANCOVA", {
  expect_null(
    mellio_payload(lm(mpg ~ as.factor(cyl) + wt, data = mtcars))$fields$model_kind
  )
  expect_null(
    mellio_payload(lm(mpg ~ factor(cyl) + wt, data = mtcars))$fields$model_kind
  )
  expect_null(
    mellio_payload(lm(mpg ~ cut(wt, 2) + hp, data = mtcars))$fields$model_kind
  )
  expect_null(
    mellio_payload(lm(mpg ~ I(cyl > 4) + wt, data = mtcars))$fields$model_kind
  )
  # scale() produces a numeric column → numeric term, no ANCOVA.
  expect_null(
    mellio_payload(lm(mpg ~ scale(wt) + hp, data = mtcars))$fields$model_kind
  )
})

test_that("lavaan CFA produces a structural card", {
  skip_if_lavaan_unusable()
  suppressWarnings(suppressMessages(library(lavaan)))
  HS.model <- "visual =~ x1 + x2 + x3
                textual =~ x4 + x5 + x6
                speed =~ x7 + x8 + x9"
  fit <- cfa(HS.model, data = HolzingerSwineford1939)
  p <- mellio_payload(fit)

  expect_equal(p$card_kind, "structural")
  expect_equal(p$type, "lavaan_cfa")
  expect_match(p$type_label, "Confirmatory factor analysis")
  expect_equal(p$fields$n, 301)

  # Report zone: fit indices
  fi <- p$fields$report_zone$fit_indices
  expect_true(length(fi) >= 6)
  fi_names <- vapply(fi, function(x) x$name, character(1))
  expect_true(all(c("chi2", "CFI", "TLI", "RMSEA", "SRMR", "AIC", "BIC") %in% fi_names))

  # chi2 carries df + p
  chi2 <- fi[[which(fi_names == "chi2")]]
  expect_true(is.numeric(chi2$df))
  expect_true(chi2$p_value >= 0 && chi2$p_value <= 1)

  # RMSEA carries a CI
  rmsea <- fi[[which(fi_names == "RMSEA")]]
  expect_length(rmsea$ci, 2)

  # Inspection zone: parameters
  params <- p$fields$inspection_zone$parameters
  expect_true(length(params) > 0)
  expect_true(all(c("lhs", "op", "rhs", "estimate") %in% names(params[[1]])))

  rel <- p$fields$inspection_zone$reliability
  expect_equal(length(rel), 3L)
  expect_true(all(c("factor", "omega", "ave", "n_indicators") %in% names(rel[[1]])))
  expect_true(all(vapply(rel, function(x) is.numeric(x$omega) && x$omega >= 0 && x$omega <= 1, logical(1))))

  mi <- p$fields$inspection_zone$modification_indices
  expect_true(length(mi) > 0)
  expect_true(length(mi) <= 10L)
  expect_true(all(c("lhs", "op", "rhs", "mi") %in% names(mi[[1]])))

  # P3b: standardized loading range summary in the report zone
  lr <- p$fields$report_zone$loadings_range
  expect_true(is.list(lr))
  expect_true(is.numeric(lr$min) && is.finite(lr$min))
  expect_true(is.numeric(lr$max) && is.finite(lr$max))
  expect_true(lr$min <= lr$max)
  expect_equal(lr$n, 9L)   # 3 factors × 3 indicators

  expect_equal(p$figure_data$structural_path_diagram$model_type, "lavaan_cfa")
  expect_equal(p$metadata$available_figures[[1]]$type, "structural_path_diagram")
  expect_equal(p$metadata$available_figures[[1]]$label, "CFA path diagram")

  diagram <- p$figure_data$structural_path_diagram
  expect_equal(diagram$model_class, "cfa")
  expect_length(diagram$nodes, 12L)
  expect_length(diagram$edges, 12L)
  edge_types <- table(vapply(diagram$edges, function(edge) edge$type, character(1)))
  expect_equal(unname(edge_types[["measurement"]]), 9L)
  expect_equal(unname(edge_types[["covariance"]]), 3L)
  nodes <- sem_fixture_nodes_by_id(diagram)
  expect_equal(nodes$visual$role, "latent")
  expect_false(nodes$visual$observed)
  expect_equal(nodes$x1$role, "indicator")
  expect_true(nodes$x1$observed)
  sem_fixture_edge(diagram, "visual", "textual", "covariance")
})

test_that("lavaan CFA with default ML estimator does not tag scaled fit indices", {
  skip_if_lavaan_unusable()
  suppressWarnings(suppressMessages(library(lavaan)))
  HS.model <- "visual =~ x1 + x2 + x3
                textual =~ x4 + x5 + x6
                speed =~ x7 + x8 + x9"
  fit <- cfa(HS.model, data = HolzingerSwineford1939)
  p <- mellio_payload(fit)

  expect_equal(p$fields$report_zone$estimator, "ML")
  fi <- p$fields$report_zone$fit_indices
  # No scaled flag on any index — ML doesn't produce *.scaled keys.
  for (entry in fi) expect_null(entry$scaled)
})

test_that("lavaan CFA with MLR prefers scaled fit indices and tags estimator", {
  skip_if_lavaan_unusable()
  suppressWarnings(suppressMessages(library(lavaan)))
  HS.model <- "visual =~ x1 + x2 + x3
                textual =~ x4 + x5 + x6
                speed =~ x7 + x8 + x9"
  fit <- cfa(HS.model, data = HolzingerSwineford1939, estimator = "MLR")
  p <- mellio_payload(fit)

  expect_equal(p$fields$report_zone$estimator, "MLR")
  fi <- p$fields$report_zone$fit_indices
  named <- setNames(fi, vapply(fi, function(x) x$name, character(1)))

  # chi2, CFI, TLI, RMSEA should carry scaled=TRUE under MLR.
  for (k in c("chi2", "CFI", "TLI", "RMSEA")) {
    expect_true(isTRUE(named[[k]]$scaled),
                info = paste(k, "should be scaled under MLR"))
  }

  # Scaled chi-square value should equal fitMeasures()$chisq.scaled,
  # which differs from the unscaled chisq under MLR (Yuan-Bentler).
  fm <- lavaan::fitMeasures(fit)
  expect_equal(named$chi2$value, unname(fm[["chisq.scaled"]]), tolerance = 1e-6)
  expect_true(abs(named$chi2$value - unname(fm[["chisq"]])) > 1e-6)

  # RMSEA CI should be the scaled bounds.
  expect_equal(named$RMSEA$ci[1], unname(fm[["rmsea.ci.lower.scaled"]]),
               tolerance = 1e-6)
  expect_equal(named$RMSEA$ci[2], unname(fm[["rmsea.ci.upper.scaled"]]),
               tolerance = 1e-6)
})

test_that("lavaan SEM detected as structural lavaan_sem", {
  skip_if_lavaan_unusable()
  suppressWarnings(suppressMessages(library(lavaan)))
  m <- "ind60 =~ x1 + x2 + x3
        dem60 =~ y1 + y2 + y3 + y4
        dem65 =~ y5 + y6 + y7 + y8
        dem60 ~ ind60
        dem65 ~ ind60 + dem60"
  fit <- sem(m, data = PoliticalDemocracy)
  p <- mellio_payload(fit)
  expect_equal(p$type, "lavaan_sem")
  diagram <- p$figure_data$structural_path_diagram
  expect_equal(diagram$model_type, "lavaan_sem")
  expect_equal(diagram$model_class, "latent_sem")
  expect_equal(p$metadata$available_figures[[1]]$type, "structural_path_diagram")
  nodes <- sem_fixture_nodes_by_id(diagram)
  expect_equal(nodes$ind60$role, "latent")
  expect_false(nodes$ind60$observed)
  expect_equal(nodes$x1$role, "indicator")
  expect_true(nodes$x1$observed)
  sem_fixture_edge(diagram, "ind60", "dem60", "structural")
  sem_fixture_edge(diagram, "dem60", "dem65", "structural")
  expect_true(any(vapply(diagram$edges, function(edge) {
    identical(edge$type, "structural") && identical(edge$op, "~")
  }, logical(1))))
})

test_that("lavaan fitted observed SEM preserves omitted covariates and residual covariance", {
  skip_if_lavaan_unusable()
  suppressWarnings(suppressMessages(library(lavaan)))
  set.seed(7301)
  n <- 180L
  dat <- data.frame(
    x = rnorm(n),
    c1 = rnorm(n)
  )
  dat$m1 <- 0.55 * dat$x + 0.30 * dat$c1 + rnorm(n)
  dat$m2 <- -0.20 * dat$x + 0.35 * dat$c1 +
    0.45 * scale(dat$m1)[, 1] + rnorm(n)
  dat$y <- 0.45 * dat$m1 + 0.25 * dat$m2 + 0.20 * dat$x +
    0.30 * dat$c1 + rnorm(n)

  model <- "m1 ~ x + c1
            m2 ~ x + c1
            y ~ m1 + m2 + x + c1
            m1 ~~ m2"
  fit <- sem(model, data = dat)
  p <- mellio_payload(fit, diagram_omit = c("c1", "y ~ x"))
  diagram <- p$figure_data$structural_path_diagram

  expect_equal(p$type, "lavaan_sem")
  expect_equal(diagram$model_class, "observed_sem")
  expect_length(diagram$nodes, 4L)
  expect_false("c1" %in% names(sem_fixture_nodes_by_id(diagram)))
  expect_length(diagram$edges, 5L)
  edge_types <- table(vapply(diagram$edges, function(edge) edge$type, character(1)))
  expect_equal(unname(edge_types[["structural"]]), 4L)
  expect_equal(unname(edge_types[["covariance"]]), 1L)
  sem_fixture_edge(diagram, "m1", "m2", "covariance")
  expect_false(any(vapply(diagram$edges, function(edge) {
    identical(edge$source, "x") &&
      identical(edge$target, "y") &&
      identical(edge$type, "structural")
  }, logical(1))))

  expect_equal(diagram$omitted_path_count, 5L)
  expect_length(diagram$omitted_paths, 5L)
  expect_length(diagram$hidden_covariates, 1L)
  expect_equal(diagram$hidden_covariates[[1]]$id, "c1")
  expect_equal(diagram$hidden_covariates[[1]]$path_count, 3L)
  expect_equal(
    sort(unlist(diagram$hidden_covariates[[1]]$targets, use.names = FALSE)),
    c("m1", "m2", "y")
  )

  omitted <- diagram$omitted_paths
  expect_true(any(vapply(omitted, function(edge) {
    identical(edge$source, "x") &&
      identical(edge$target, "y") &&
      identical(edge$reason, "omitted_path")
  }, logical(1))))
  expect_true(any(vapply(omitted, function(edge) {
    identical(edge$source, "c1") &&
      identical(edge$target, "m1") &&
      identical(edge$reason, "hidden_covariate")
  }, logical(1))))
  expect_true(any(vapply(omitted, function(edge) {
    identical(edge$source, "x") &&
      identical(edge$target, "c1") &&
      identical(edge$type, "covariance") &&
      identical(edge$reason, "omitted_path")
  }, logical(1))))
})

test_that("lavaan multi-group CFA payload does not auto-create a diagram", {
  skip_if_lavaan_unusable()
  suppressWarnings(suppressMessages(library(lavaan)))
  HS.model <- "visual =~ x1 + x2 + x3
                textual =~ x4 + x5 + x6
                speed =~ x7 + x8 + x9"
  fit <- cfa(HS.model, data = HolzingerSwineford1939, group = "school")
  p <- mellio_payload(fit)

  expect_equal(p$type, "lavaan_cfa")
  expect_match(p$type_label, "2 groups")
  expect_null(p$figure_data$structural_path_diagram)
  expect_null(p$metadata$available_figures)
})

test_that("mellio_payload.numeric produces a descriptive_summary card", {
  x <- c(1, 2, 3, 4, 5, 3, 2, 4, 5, 1, NA, 6)
  p <- mellio_payload(x)
  expect_equal(p$type, "descriptive_summary")
  expect_equal(p$fields$statistic$name, "M")
  expect_equal(p$fields$estimate$name, "SD")
  expect_equal(p$fields$n, 11L)
  expect_equal(p$fields$n_missing, 1L)
  expect_length(p$fields$range, 2)
  expect_true(is.na(p$fields$p_value))
})

test_that("mellio_payload.numeric errors on all-NA input", {
  expect_error(mellio_payload(c(NA_real_, NA_real_)), "No non-missing")
})

test_that("mellio_payload.summaryDefault preserves what summary() gives", {
  s <- summary(c(1, 2, 3, 4, 5, NA))
  p <- mellio_payload(s)
  expect_equal(p$type, "descriptive_summary")
  expect_match(p$type_label, "summary")
  expect_equal(p$fields$statistic$name, "M")
  expect_true(is.numeric(p$fields$statistic$value))
  # No SD field from summary()
  expect_null(p$fields$estimate)
  # Quartiles are preserved
  expect_true(!is.null(p$fields$quartiles))
  expect_equal(p$fields$n_missing, 1)
})

test_that("mellio_payload.character parses captured correlation output", {
  lines <- c(
    "motive_enjoyment       r = +0.019, p = 0.5573",
    "motive_ltb             r = +0.076, p = 0.0180",
    "motive_appearance      r = +0.078, p = 0.0149"
  )

  p <- mellio_payload(lines, title = "SES correlations")
  expect_equal(p$type, "custom_correlation_table")
  expect_equal(p$type_label, "SES correlations")
  expect_equal(p$card_kind, "table")
  expect_equal(p$fields$source, "captured_console")
  expect_equal(length(p$fields$rows), 3L)
  expect_equal(p$fields$rows[[2]]$x, "motive_ltb")
  expect_equal(p$fields$rows[[2]]$r, 0.076)
  expect_match(p$raw_output, "motive_enjoyment")
})

test_that("mellio_open accepts captured correlation output", {
  lines <- c(
    "motive_enjoyment       r = +0.019, p = 0.5573",
    "motive_ltb             r = +0.076, p = 0.0180"
  )

  url <- quiet_mellio_open(lines, title = "SES correlations", browse = FALSE)
  b64 <- sub(".*payload=", "", url)
  b64 <- sub("&.*$", "", b64)
  std <- chartr("-_", "+/", b64)
  pad <- (4 - nchar(std) %% 4) %% 4
  if (pad > 0) std <- paste0(std, strrep("=", pad))
  parsed <- jsonlite::fromJSON(rawToChar(jsonlite::base64_dec(std)),
                               simplifyVector = FALSE)
  expect_equal(parsed$type_label, "SES correlations")
  expect_equal(parsed$fields$source, "captured_console")
})

test_that("mellio_payload.character rejects prose that is not a Stats result", {
  expect_error(mellio_payload("this is just a sentence"), "not recognized")
})

test_that("mellio_payload.character gives ANOVA guidance for captured ANOVA text", {
  lines <- c(
    "ANOVA F-test:",
    "            Df Sum Sq Mean Sq F value Pr(>F)",
    "SES_band     4  2.31    0.58    3.40 0.009",
    "Residuals  966 164.00   0.17"
  )
  expect_error(mellio_payload(lines), "For ANOVA")
  expect_error(mellio_payload(lines), "data.frame")
})

test_that("mellio_payload.table produces a frequency table card", {
  p <- mellio_payload(table(c("a", "b", "a")))

  expect_equal(p$card_kind, "table")
  expect_equal(p$type, "frequency_table")
  expect_equal(p$fields$table_type, "frequency_table")
  expect_equal(length(p$fields$rows), 2L)
  expect_equal(length(p$fields$columns), 2L)
  expect_equal(p$fields$total, 3L)
})

test_that("psych::alpha payload produces a reliability Result Card", {
  skip_if_not_installed("psych")
  set.seed(42)
  latent <- rnorm(80)
  d <- data.frame(
    item1 = latent + rnorm(80, sd = 0.4),
    item2 = latent + rnorm(80, sd = 0.4),
    item3 = latent + rnorm(80, sd = 0.4),
    item4 = latent + rnorm(80, sd = 0.4)
  )
  a <- suppressMessages(suppressWarnings(psych::alpha(d, warnings = FALSE)))

  p <- mellio_payload(a)
  expect_equal(p$type, "cronbach_alpha")
  expect_equal(p$card_kind, "inline")
  expect_equal(p$fields$statistic$name, "\u03b1")
  expect_true(is.numeric(p$fields$statistic$value))
  expect_equal(p$fields$item_count, 4L)
  expect_equal(p$fields$n, 80L)
  pkg_names <- vapply(p$packages, `[[`, character(1), "name")
  expect_true("psych" %in% pkg_names)

  # P3c: 95% CI on \u03b1 via Feldt normal approximation. Bounds bracket
  # the point estimate and conf_level is set.
  expect_length(p$fields$statistic$ci, 2L)
  expect_true(is.numeric(p$fields$statistic$ci[1]) &&
              is.numeric(p$fields$statistic$ci[2]))
  expect_lt(p$fields$statistic$ci[1], p$fields$statistic$value)
  expect_gt(p$fields$statistic$ci[2], p$fields$statistic$value)
  expect_equal(p$fields$conf_level, 0.95)
})

test_that("psych::corr.test payload produces a table Result Card", {
  skip_if_not_installed("psych")
  ct <- suppressWarnings(psych::corr.test(
    x = mtcars[, c("mpg", "wt")],
    y = mtcars[, "hp", drop = FALSE],
    adjust = "none"
  ))

  p <- mellio_payload(ct)
  expect_equal(p$type, "psych_corr_test")
  expect_equal(p$card_kind, "table")
  expect_equal(p$fields$table_type, "correlations")
  expect_equal(length(p$fields$rows), 2L)
  expect_equal(p$fields$rows[[1]]$x, "mpg")
  expect_equal(p$fields$rows[[1]]$y, "hp")
  expect_true(is.numeric(p$fields$rows[[1]]$r))
  expect_true(is.numeric(p$fields$rows[[1]]$p_value))
  expect_equal(p$fields$n, 32L)
  expect_null(p$figure_data)
  expect_null(p$metadata$available_figures)
})

test_that("psych::corr.test symmetric payload includes heatmap figure data", {
  skip_if_not_installed("psych")
  ct <- suppressWarnings(psych::corr.test(
    mtcars[, c("mpg", "wt", "hp")]
  ))

  p <- mellio_payload(ct)
  expect_equal(p$type, "psych_corr_test")
  expect_equal(p$fields$adjust, "holm")
  expect_equal(p$metadata$available_figures[[1]]$type, "correlation_heatmap")
  expect_true(isTRUE(p$metadata$available_figures[[1]]$default))
  expect_equal(p$metadata$available_figures[[2]]$type, "correlation_forest")
  expect_false(isTRUE(p$metadata$available_figures[[2]]$default))

  fig <- p$figure_data$correlation_heatmap
  expect_equal(fig$variables, c("mpg", "wt", "hp"))
  expect_equal(fig$method, "Pearson")
  expect_equal(fig$adjust, "holm")
  expect_equal(fig$n, 32L)
  expect_equal(length(fig$r), 3L)
  expect_equal(length(fig$p), 3L)
  expect_equal(fig$r[[1]][[2]], unname(ct$r[1, 2]), tolerance = 1e-8)
  expect_equal(fig$p[[2]][[1]], unname(ct$p[1, 2]), tolerance = 1e-8)

  forest <- p$figure_data$correlation_forest
  expect_equal(forest$variables, c("mpg", "wt", "hp"))
  expect_equal(forest$method, "Pearson")
  expect_equal(forest$adjust, "holm")
  expect_equal(forest$n, 32L)
  expect_equal(length(forest$pairs), 3L)
  expect_equal(forest$pairs[[1]]$r, unname(ct$r[1, 2]), tolerance = 1e-8)
  expect_true(is.numeric(forest$pairs[[1]]$ci_lower))
  expect_true(is.numeric(forest$pairs[[1]]$ci_upper))
})

test_that("mellio_open routes bare correlation matrices to Stats", {
  mat <- cor(mtcars[, c("mpg", "wt", "hp")])
  url <- quiet_mellio_open(mat, browse = FALSE)
  expect_true(grepl("#stats/payload=", url, fixed = TRUE))

  b64 <- sub(".*payload=", "", url)
  b64 <- sub("&.*$", "", b64)
  std <- chartr("-_", "+/", b64)
  pad <- (4 - nchar(std) %% 4) %% 4
  if (pad > 0) std <- paste0(std, strrep("=", pad))
  parsed <- jsonlite::fromJSON(rawToChar(jsonlite::base64_dec(std)),
                               simplifyVector = FALSE)

  expect_equal(parsed$type, "correlation_matrix")
  expect_equal(parsed$card_kind, "table")
  expect_equal(parsed$metadata$available_figures[[1]]$type, "correlation_heatmap")
  expect_equal(parsed$metadata$available_figures[[2]]$type, "correlation_forest")
})

test_that("psych::corr.test figure data preserves literal correlation methods", {
  skip_if_not_installed("psych")
  spearman <- suppressWarnings(psych::corr.test(
    mtcars[, c("mpg", "wt", "hp")],
    method = "spearman",
    adjust = "none"
  ))
  kendall <- suppressWarnings(psych::corr.test(
    mtcars[, c("mpg", "wt", "hp")],
    method = "kendall",
    adjust = "none"
  ))

  p_spearman <- mellio_payload(spearman)
  p_kendall <- mellio_payload(kendall)

  expect_equal(p_spearman$figure_data$correlation_heatmap$method, "Spearman")
  expect_equal(p_spearman$figure_data$correlation_forest$method, "Spearman")
  expect_null(p_spearman$figure_data$correlation_heatmap$adjust)
  expect_equal(p_kendall$figure_data$correlation_heatmap$method, "Kendall")
  expect_equal(p_kendall$figure_data$correlation_forest$method, "Kendall")
})

test_that("mellio_open preserves small p-values in table payload URLs", {
  skip_if_not_installed("psych")
  ct <- suppressWarnings(psych::corr.test(
    x = mtcars[, c("mpg", "wt")],
    y = mtcars[, "hp", drop = FALSE],
    adjust = "none"
  ))
  url <- quiet_mellio_open(ct, browse = FALSE)
  b64 <- sub(".*payload=", "", url)
  b64 <- sub("&.*$", "", b64)
  std <- chartr("-_", "+/", b64)
  pad <- (4 - nchar(std) %% 4) %% 4
  if (pad > 0) std <- paste0(std, strrep("=", pad))
  parsed <- jsonlite::fromJSON(rawToChar(jsonlite::base64_dec(std)),
                               simplifyVector = FALSE)
  expect_gt(parsed$fields$rows[[1]]$p_value, 0)
  expect_lt(parsed$fields$rows[[1]]$p_value, .001)
})

test_that("correlation-shaped data frames become table Result Cards", {
  df <- data.frame(
    motivation = c("motive_enjoyment", "motive_ltb"),
    r = c(.31, -.22),
    p = c(.014, .091)
  )

  p <- mellio_payload(df)
  expect_equal(p$type, "custom_correlation_table")
  expect_equal(p$card_kind, "table")
  expect_equal(length(p$fields$rows), 2L)
  expect_equal(p$fields$rows[[1]]$x, "motive_enjoyment")
  expect_equal(p$fields$rows[[1]]$r, .31)
})

test_that("correlation-package-shaped data frames become table Result Cards", {
  df <- data.frame(
    Parameter1 = c("motive_enjoyment", "motive_ltb"),
    Parameter2 = c("emotion_resp", "emotion_resp"),
    r = c(.31, -.22),
    CI_low = c(.08, -.46),
    CI_high = c(.51, .03),
    p = c(.014, .091)
  )

  p <- mellio_payload(df)
  expect_equal(p$type, "custom_correlation_table")
  expect_equal(p$fields$rows[[1]]$x, "motive_enjoyment")
  expect_equal(p$fields$rows[[1]]$y, "emotion_resp")
  expect_equal(p$fields$rows[[1]]$ci_lower, .08)
})

test_that("ANOVA-shaped data frames become table Result Cards", {
  df <- data.frame(
    outcome = c("motive_enjoyment", "motive_ltb"),
    term = c("SES_band", "SES_band"),
    df1 = c(4L, 4L),
    df2 = c(966L, 966L),
    F = c(2.31, 3.40),
    p = c(.056, .009)
  )

  p <- mellio_payload(df, title = "ANOVA by SES band")
  expect_equal(p$type, "custom_anova_table")
  expect_equal(p$type_label, "ANOVA by SES band")
  expect_equal(p$card_kind, "table")
  expect_equal(p$fields$table_type, "anova")
  expect_equal(length(p$fields$rows), 2L)
  expect_equal(p$fields$rows[[1]]$outcome, "motive_enjoyment")
  expect_equal(p$fields$rows[[1]]$term, "SES_band")
  expect_equal(p$fields$rows[[1]]$df1, 4L)
  expect_equal(p$fields$rows[[1]]$df2, 966L)
  expect_equal(p$fields$rows[[1]]$f, 2.31)
  expect_equal(p$fields$rows[[1]]$p_value, .056)
})

test_that("summary(aov())-shaped data frames infer residual df", {
  fit <- aov(mpg ~ factor(cyl), data = mtcars)
  df <- as.data.frame(summary(fit)[[1]])

  p <- mellio_payload(df)
  expect_equal(p$type, "custom_anova_table")
  expect_equal(p$fields$table_type, "anova")
  expect_equal(length(p$fields$rows), 1L)
  expect_equal(p$fields$rows[[1]]$term, "factor(cyl)")
  expect_equal(p$fields$rows[[1]]$df1, 2L)
  expect_equal(p$fields$rows[[1]]$df2, 29L)
  expect_true(is.numeric(p$fields$rows[[1]]$f))
  expect_true(is.numeric(p$fields$rows[[1]]$p_value))
})

test_that("summary.aov objects reuse the ANOVA extractor", {
  fit <- aov(mpg ~ factor(cyl), data = mtcars)
  p <- mellio_payload(summary(fit), .call = "summary(fit)", .env = environment())

  expect_equal(p$type, "anova_single_model")
  expect_equal(p$call, "summary(fit)")
  expect_equal(p$fields$term, "factor(cyl)")
  expect_equal(p$fields$outcome, "mpg")
  expect_equal(p$fields$predictor, "cyl")
  expect_equal(p$fields$ss_type_label, "Type I (sequential)")
  expect_equal(p$fields$df_effect, 2L)
  expect_equal(p$fields$df_error, 29L)
  expect_true(is.numeric(p$fields$statistic$value))
  expect_true(is.numeric(p$fields$p_value))
  expect_match(p$raw_output, "factor\\(cyl\\)")
})

test_that("balanced factorial aov payloads disclose sequential SS context", {
  tg <- transform(ToothGrowth, dose = factor(dose))
  fit <- aov(len ~ dose * supp, data = tg)
  p <- mellio_payload(summary(fit), .call = "summary(fit)", .env = environment())

  expect_equal(p$fields$outcome, "len")
  expect_equal(p$fields$ss_type_label, "Type I (sequential)")
  expect_match(p$fields$design_balance_note, "Type I, II, and III")
})

test_that("generic statistic data frames are not mistaken for ANOVA", {
  df <- data.frame(
    term = "x",
    df = 30,
    statistic = 2.1,
    p.value = .044,
    check.names = FALSE
  )

  expect_error(mellio_payload(df), "could not recognize")
})

test_that("coefficient-shaped data frames become table Result Cards", {
  df <- data.frame(
    term = "x",
    estimate = .2,
    std.error = .05,
    statistic = 4,
    conf.low = .1,
    conf.high = .3,
    p.value = .01,
    check.names = FALSE
  )

  p <- mellio_payload(df)
  expect_equal(p$type, "custom_regression_table")
  expect_equal(p$card_kind, "table")
  expect_equal(p$fields$table_type, "coefficients")
  expect_equal(p$fields$rows[[1]]$term, "x")
  expect_equal(p$fields$rows[[1]]$estimate, .2)
  expect_equal(p$fields$rows[[1]]$ci_upper, .3)
})

test_that("descriptive-shaped data frames become table Result Cards", {
  df <- data.frame(
    Variable = c("emotion_resp", "motiv_persist"),
    M = c(4.1, 5.2),
    SD = c(.8, .7),
    n = c(120L, 118L)
  )

  p <- mellio_payload(df)
  expect_equal(p$type, "custom_descriptive_table")
  expect_equal(p$card_kind, "table")
  expect_equal(p$fields$table_type, "descriptives")
  expect_equal(p$fields$rows[[1]]$variable, "emotion_resp")
  expect_equal(p$fields$rows[[1]]$mean, 4.1)
})

test_that("unsupported data frames give short table-shape guidance", {
  expect_error(mellio_payload(data.frame(a = 1:3, b = 4:6)), "could not recognize")
  expect_error(mellio_payload(data.frame(a = 1:3, b = 4:6)), "mellio_open\\(melliotab")
})

test_that("function input gives a guided warning without executing", {
  fn <- function(outcome_var) {
    cor.test(mtcars$mpg, mtcars[[outcome_var]])
  }
  expect_error(mellio_payload(fn), "function, not a result")
  expect_error(mellio_payload(fn), "outcome_var")
  expect_error(mellio_payload(fn), "psych::corr.test")
})

test_that("lavaan focal-path marking surfaces rows in report_zone", {
  skip_if_lavaan_unusable()
  suppressWarnings(suppressMessages(library(lavaan)))
  m <- "ind60 =~ x1 + x2 + x3
        dem60 =~ y1 + y2 + y3 + y4
        dem60 ~ ind60"
  fit <- sem(m, data = PoliticalDemocracy)

  # No focal → no focal_paths field
  p_nofocal <- mellio_payload(fit)
  expect_null(p_nofocal$fields$report_zone$focal_paths)

  # With focal → matching rows surface
  p <- mellio_payload(fit, focal = c("dem60 ~ ind60"))
  fp <- p$fields$report_zone$focal_paths
  expect_equal(length(fp), 1L)
  expect_equal(fp[[1]]$lhs, "dem60")
  expect_equal(fp[[1]]$op,  "~")
  expect_equal(fp[[1]]$rhs, "ind60")

  # Same rows in inspection_zone carry is_focal flag
  insp <- p$fields$inspection_zone$parameters
  matched <- Filter(function(r) isTRUE(r$is_focal), insp)
  expect_equal(length(matched), 1L)
})

test_that("focal matching is whitespace-insensitive", {
  skip_if_lavaan_unusable()
  suppressWarnings(suppressMessages(library(lavaan)))
  m <- "ind60 =~ x1 + x2 + x3
        dem60 =~ y1 + y2
        dem60 ~ ind60"
  fit <- sem(m, data = PoliticalDemocracy)
  p1 <- mellio_payload(fit, focal = "dem60 ~ ind60")
  p2 <- mellio_payload(fit, focal = "dem60~ind60")
  p3 <- mellio_payload(fit, focal = " dem60  ~  ind60 ")
  expect_equal(length(p1$fields$report_zone$focal_paths), 1L)
  expect_equal(length(p2$fields$report_zone$focal_paths), 1L)
  expect_equal(length(p3$fields$report_zone$focal_paths), 1L)
})

test_that("lavaan mediation (:= defined params) is detected", {
  skip_if_lavaan_unusable()
  suppressWarnings(suppressMessages(library(lavaan)))
  set.seed(1)
  df <- data.frame(X = rnorm(50), M = rnorm(50), Y = rnorm(50))
  m <- "M ~ a*X
        Y ~ b*M + cc*X
        indirect := a*b"
  fit <- sem(m, data = df)
  p <- mellio_payload(fit)
  expect_equal(p$type, "lavaan_mediation")

  # The defined parameter should appear in the inspection zone
  ops <- vapply(p$fields$inspection_zone$parameters,
                function(x) x$op, character(1))
  expect_true(":=" %in% ops)
})

test_that("lavaan := rows auto-promote to focal_paths when no explicit focal", {
  # The whole point of declaring `indirect := a*b` is to surface the
  # derived effect; the extractor should auto-promote := rows to
  # focal_paths so the indirect effect gets prominent display in the
  # report zone — without requiring the user to also set
  # attr(fit, "mellio_focal").
  skip_if_lavaan_unusable()
  suppressWarnings(suppressMessages(library(lavaan)))
  set.seed(42)
  n <- 100
  df <- data.frame(X = rnorm(n))
  df$M <- 0.5 * df$X + rnorm(n)
  df$Y <- 0.3 * df$X + 0.4 * df$M + rnorm(n)
  m <- "M ~ a*X
        Y ~ b*M + c*X
        indirect := a*b
        total    := c + (a*b)"
  fit <- sem(m, data = df)
  p <- mellio_payload(fit)

  fp <- p$fields$report_zone$focal_paths
  expect_true(!is.null(fp), info = "focal_paths should be auto-populated")
  expect_equal(length(fp), 2L)

  # Both auto-promoted rows should carry the auto_focal marker so the JS
  # side can distinguish auto vs explicit focal (might want to render
  # them differently later).
  for (row in fp) {
    expect_identical(row$op, ":=")
    expect_true(isTRUE(row$auto_focal),
                info = paste("row", row$lhs, "missing auto_focal marker"))
  }
  fp_names <- vapply(fp, `[[`, character(1), "lhs")
  expect_setequal(fp_names, c("indirect", "total"))
})

test_that("explicit focal= still wins over := auto-promotion", {
  # If the user has marked specific paths via focal=, respect that —
  # don't override with the := auto-rule.
  skip_if_lavaan_unusable()
  suppressWarnings(suppressMessages(library(lavaan)))
  set.seed(7)
  n <- 80
  df <- data.frame(X = rnorm(n))
  df$M <- 0.5 * df$X + rnorm(n)
  df$Y <- 0.3 * df$X + 0.4 * df$M + rnorm(n)
  m <- "M ~ a*X
        Y ~ b*M + c*X
        indirect := a*b"
  fit <- sem(m, data = df)
  # Mark only the Y ~ M path as focal (a regression, not the defined one)
  p <- mellio_payload(fit, focal = "Y ~ M")

  fp <- p$fields$report_zone$focal_paths
  expect_equal(length(fp), 1L)
  expect_identical(fp[[1]]$op, "~")
  expect_identical(fp[[1]]$lhs, "Y")
  expect_identical(fp[[1]]$rhs, "M")
  # Should NOT be marked auto_focal — user picked it explicitly
  expect_null(fp[[1]]$auto_focal)
})

test_that("lavaan parallel mediation payload carries multi-mediator structure", {
  # Parallel-mediator support relies on the lavaan bridge surfacing every
  # regression (M1 ~ X, M2 ~ X, Y ~ M1, Y ~ M2, Y ~ X) and each defined
  # effect (ind_M1, ind_M2, total_indirect, optional contrast) in the
  # inspection-zone parameter list. Topology detection happens JS-side,
  # so this test pins down the shape the JS bridge expects to consume.
  skip_if_lavaan_unusable()
  suppressWarnings(suppressMessages(library(lavaan)))
  set.seed(11)
  n <- 250
  df <- data.frame(X = rnorm(n))
  df$M1 <- 0.5 * df$X + rnorm(n)
  df$M2 <- 0.3 * df$X + rnorm(n)
  df$Y  <- 0.3 * df$M1 + 0.4 * df$M2 + 0.1 * df$X + rnorm(n)

  model <- "
    M1 ~ a1*X
    M2 ~ a2*X
    Y  ~ b1*M1 + b2*M2 + cp*X
    ind_M1         := a1*b1
    ind_M2         := a2*b2
    total_indirect := ind_M1 + ind_M2
    contrast_M1_M2 := ind_M1 - ind_M2
  "
  fit <- sem(model, data = df)
  p <- mellio_payload(fit)

  expect_equal(p$type, "lavaan_mediation")
  params <- p$fields$inspection_zone$parameters

  # Regressions: M1~X, M2~X, Y~M1, Y~M2, Y~X
  regressions <- Filter(function(r) identical(r$op, "~"), params)
  reg_keys <- vapply(regressions, function(r) paste(r$lhs, r$rhs, sep = "~"),
                     character(1))
  expect_true(all(c("M1~X", "M2~X", "Y~M1", "Y~M2", "Y~X") %in% reg_keys))

  # Defined effects: ind_M1, ind_M2, total_indirect, contrast_M1_M2
  defined <- Filter(function(r) identical(r$op, ":="), params)
  def_keys <- vapply(defined, function(r) r$lhs, character(1))
  expect_true(all(c("ind_M1", "ind_M2", "total_indirect",
                    "contrast_M1_M2") %in% def_keys))

  # Each defined effect carries an estimate + CI bounds that the JS-side
  # companion table needs.
  ind_m1 <- Find(function(r) identical(r$lhs, "ind_M1"), defined)
  expect_true(is.numeric(ind_m1$estimate))
  expect_true(is.numeric(ind_m1$ci_lower) && is.numeric(ind_m1$ci_upper))

  # Sample size for the diagram's N anchor.
  expect_equal(p$fields$n, 250)

  # Auto-promoted focal_paths should carry all four defined effects.
  fp <- p$fields$report_zone$focal_paths
  expect_equal(length(fp), 4L)

})

test_that("lavaan serial mediation payload carries chain structure", {
  # Serial-mediator support relies on the lavaan bridge surfacing the
  # mediator-to-mediator path (M2 ~ M1) alongside optional direct paths
  # from X to M2 and M1 to Y. The JS bridge then infers the Model 6
  # topology and lays out X → M1 → M2 → Y.
  skip_if_lavaan_unusable()
  suppressWarnings(suppressMessages(library(lavaan)))
  set.seed(13)
  n <- 250
  df <- data.frame(X = rnorm(n))
  df$M1 <- 0.5 * df$X + rnorm(n)
  df$M2 <- 0.45 * df$M1 + 0.25 * df$X + rnorm(n)
  df$Y  <- 0.4 * df$M2 + 0.18 * df$M1 + 0.1 * df$X + rnorm(n)

  model <- "
    M1 ~ a1*X
    M2 ~ d21*M1 + a2*X
    Y  ~ b2*M2 + b1*M1 + cp*X
    ind_M1          := a1*b1
    ind_M2          := a2*b2
    serial_indirect := a1*d21*b2
    total_indirect  := ind_M1 + ind_M2 + serial_indirect
  "
  fit <- sem(model, data = df)
  p <- mellio_payload(fit)

  expect_equal(p$type, "lavaan_mediation")
  params <- p$fields$inspection_zone$parameters

  regressions <- Filter(function(r) identical(r$op, "~"), params)
  reg_keys <- vapply(regressions, function(r) paste(r$lhs, r$rhs, sep = "~"),
                     character(1))
  expect_true(all(c("M1~X", "M2~M1", "M2~X",
                    "Y~M2", "Y~M1", "Y~X") %in% reg_keys))

  defined <- Filter(function(r) identical(r$op, ":="), params)
  def_keys <- vapply(defined, function(r) r$lhs, character(1))
  expect_true(all(c("ind_M1", "ind_M2", "serial_indirect",
                    "total_indirect") %in% def_keys))

  serial_ind <- Find(function(r) identical(r$lhs, "serial_indirect"), defined)
  expect_true(is.numeric(serial_ind$estimate))
  expect_true(is.numeric(serial_ind$ci_lower) &&
              is.numeric(serial_ind$ci_upper))

  expect_equal(p$fields$n, 250)
  fp <- p$fields$report_zone$focal_paths
  expect_equal(length(fp), 4L)
})

test_that("lavaan moderated mediation payload carries Model 14 scaffold inputs", {
  # Hayes Model 14 support is inferred JS-side from a single mediator plus an
  # outcome-side mediator-by-moderator product term. lavaan users commonly
  # precompute the observed product, so this pins the bridge shape rather than
  # requiring lavaan parser support for `M:W` syntax.
  skip_if_lavaan_unusable()
  suppressWarnings(suppressMessages(library(lavaan)))
  set.seed(17)
  n <- 250
  df <- data.frame(X = rnorm(n), W = rnorm(n))
  df$M <- 0.5 * df$X + rnorm(n)
  df$MW <- df$M * df$W
  df$Y <- 0.35 * df$M + 0.25 * df$W + 0.30 * df$MW +
    0.12 * df$X + rnorm(n)

  model <- "
    M ~ a*X
    Y ~ b*M + w*W + int*MW + cp*X
    ind_low      := a*(b + int*(-1))
    ind_mean     := a*b
    ind_high     := a*(b + int*(1))
    index_modmed := a*int
  "
  fit <- sem(model, data = df)
  p <- mellio_payload(fit)

  expect_equal(p$type, "lavaan_mediation")
  params <- p$fields$inspection_zone$parameters

  regressions <- Filter(function(r) identical(r$op, "~"), params)
  reg_keys <- vapply(regressions, function(r) paste(r$lhs, r$rhs, sep = "~"),
                     character(1))
  expect_true(all(c("M~X", "Y~M", "Y~W", "Y~MW", "Y~X") %in% reg_keys))

  defined <- Filter(function(r) identical(r$op, ":="), params)
  def_keys <- vapply(defined, function(r) r$lhs, character(1))
  expect_true(all(c("ind_low", "ind_mean", "ind_high",
                    "index_modmed") %in% def_keys))

  fp <- p$fields$report_zone$focal_paths
  expect_equal(length(fp), 4L)

  rsq <- p$fields$structural_r_squared
  rsq_vars <- vapply(rsq, function(r) r$variable, character(1))
  expect_true(all(c("M", "Y") %in% rsq_vars))

  obs <- p$fields$observed_variables
  obs_vars <- vapply(obs, function(r) r$variable, character(1))
  expect_true(all(c("M", "W", "MW") %in% obs_vars))
  w_obs <- obs[[match("W", obs_vars)]]
  expect_equal(w_obs$mean, mean(df$W), tolerance = 1e-9)
})

test_that("lavaan moderated mediation payload carries Model 7 scaffold inputs", {
  # Hayes Model 7 support is inferred JS-side from a single mediator plus a
  # mediator-side predictor-by-moderator product term. The index of moderated
  # mediation is the user's defined a3*b parameter, not the Model 14 a*int form.
  skip_if_lavaan_unusable()
  suppressWarnings(suppressMessages(library(lavaan)))
  set.seed(19)
  n <- 250
  df <- data.frame(X = rnorm(n), W = rnorm(n))
  df$XW <- df$X * df$W
  df$M <- 0.50 * df$X + 0.25 * df$W + 0.30 * df$XW + rnorm(n)
  df$Y <- 0.40 * df$M + 0.12 * df$X + rnorm(n)

  model <- "
    M ~ a1*X + a2*W + a3*XW
    Y ~ b*M + cp*X
    ind_low      := (a1 + a3*(-1))*b
    ind_mean     := a1*b
    ind_high     := (a1 + a3*(1))*b
    index_modmed := a3*b
  "
  fit <- sem(model, data = df)
  p <- mellio_payload(fit)

  expect_equal(p$type, "lavaan_mediation")
  params <- p$fields$inspection_zone$parameters

  regressions <- Filter(function(r) identical(r$op, "~"), params)
  reg_keys <- vapply(regressions, function(r) paste(r$lhs, r$rhs, sep = "~"),
                     character(1))
  expect_true(all(c("M~X", "M~W", "M~XW", "Y~M", "Y~X") %in% reg_keys))

  defined <- Filter(function(r) identical(r$op, ":="), params)
  def_keys <- vapply(defined, function(r) r$lhs, character(1))
  expect_true(all(c("ind_low", "ind_mean", "ind_high",
                    "index_modmed") %in% def_keys))

  a3 <- Find(function(r) identical(r$lhs, "M") &&
               identical(r$rhs, "XW"), regressions)
  b <- Find(function(r) identical(r$lhs, "Y") &&
              identical(r$rhs, "M"), regressions)
  idx <- Find(function(r) identical(r$lhs, "index_modmed"), defined)
  expect_equal(idx$estimate, a3$estimate * b$estimate, tolerance = 1e-8)

  fp <- p$fields$report_zone$focal_paths
  expect_equal(length(fp), 4L)

  obs <- p$fields$observed_variables
  obs_vars <- vapply(obs, function(r) r$variable, character(1))
  expect_true(all(c("X", "W", "XW") %in% obs_vars))
})

test_that("lavaan Model 7 payload preserves awkward centered product names", {
  # Regression coverage for the JS-side alias path: the observed product name
  # intentionally drops the centered suffixes used by its components.
  skip_if_lavaan_unusable()
  suppressWarnings(suppressMessages(library(lavaan)))
  set.seed(23)
  n <- 250
  stress <- rnorm(n)
  support <- rnorm(n)
  df <- data.frame(
    Stress_c = stress - mean(stress),
    Support_c = support - mean(support)
  )
  df$Stress_x_Support <- df$Stress_c * df$Support_c
  df$M <- 0.45 * df$Stress_c + 0.20 * df$Support_c +
    0.28 * df$Stress_x_Support + rnorm(n)
  df$Y <- 0.42 * df$M + 0.10 * df$Stress_c + rnorm(n)

  model <- "
    M ~ a1*Stress_c + a2*Support_c + a3*Stress_x_Support
    Y ~ b*M + cp*Stress_c
    ind_low      := (a1 + a3*(-1))*b
    ind_mean     := a1*b
    ind_high     := (a1 + a3*(1))*b
    index_modmed := a3*b
  "
  fit <- sem(model, data = df)
  p <- mellio_payload(fit)

  expect_equal(p$type, "lavaan_mediation")
  params <- p$fields$inspection_zone$parameters

  regressions <- Filter(function(r) identical(r$op, "~"), params)
  reg_keys <- vapply(regressions, function(r) paste(r$lhs, r$rhs, sep = "~"),
                     character(1))
  expect_true(all(c("M~Stress_c", "M~Support_c",
                    "M~Stress_x_Support", "Y~M",
                    "Y~Stress_c") %in% reg_keys))

  defined <- Filter(function(r) identical(r$op, ":="), params)
  idx <- Find(function(r) identical(r$lhs, "index_modmed"), defined)
  a3 <- Find(function(r) identical(r$lhs, "M") &&
               identical(r$rhs, "Stress_x_Support"), regressions)
  b <- Find(function(r) identical(r$lhs, "Y") &&
              identical(r$rhs, "M"), regressions)
  expect_equal(idx$estimate, a3$estimate * b$estimate, tolerance = 1e-8)

  obs <- p$fields$observed_variables
  obs_vars <- vapply(obs, function(r) r$variable, character(1))
  expect_true(all(c("Stress_c", "Support_c",
                    "Stress_x_Support") %in% obs_vars))
})

test_that("lavaan moderated mediation payload carries Model 8 scaffold inputs", {
  # Hayes Model 8 extends Model 7 by using the same X-by-W product on the
  # direct X -> Y path. The JS bridge renders this as Model 7 plus direct-path
  # moderation, not as a generic multi-moderation model.
  skip_if_lavaan_unusable()
  suppressWarnings(suppressMessages(library(lavaan)))
  set.seed(29)
  n <- 250
  df <- data.frame(X = rnorm(n), W = rnorm(n))
  df$XW <- df$X * df$W
  df$M <- 0.50 * df$X + 0.24 * df$W + 0.28 * df$XW + rnorm(n)
  df$Y <- 0.40 * df$M + 0.16 * df$X + 0.18 * df$W +
    0.22 * df$XW + rnorm(n)

  model <- "
    M ~ a1*X + a2*W + a3*XW
    Y ~ b*M + c1*X + c2*W + c3*XW
    ind_low       := (a1 + a3*(-1))*b
    ind_mean      := a1*b
    ind_high      := (a1 + a3*(1))*b
    index_modmed  := a3*b
    direct_low    := c1 + c3*(-1)
    direct_mean   := c1
    direct_high   := c1 + c3*(1)
  "
  fit <- sem(model, data = df)
  p <- mellio_payload(fit)

  expect_equal(p$type, "lavaan_mediation")
  params <- p$fields$inspection_zone$parameters

  regressions <- Filter(function(r) identical(r$op, "~"), params)
  reg_keys <- vapply(regressions, function(r) paste(r$lhs, r$rhs, sep = "~"),
                     character(1))
  expect_true(all(c("M~X", "M~W", "M~XW",
                    "Y~M", "Y~X", "Y~W", "Y~XW") %in% reg_keys))

  defined <- Filter(function(r) identical(r$op, ":="), params)
  def_keys <- vapply(defined, function(r) r$lhs, character(1))
  expect_true(all(c("ind_low", "ind_mean", "ind_high",
                    "index_modmed",
                    "direct_low", "direct_mean", "direct_high") %in% def_keys))

  a3 <- Find(function(r) identical(r$lhs, "M") &&
               identical(r$rhs, "XW"), regressions)
  b <- Find(function(r) identical(r$lhs, "Y") &&
              identical(r$rhs, "M"), regressions)
  c1 <- Find(function(r) identical(r$lhs, "Y") &&
               identical(r$rhs, "X"), regressions)
  c3 <- Find(function(r) identical(r$lhs, "Y") &&
               identical(r$rhs, "XW"), regressions)
  idx <- Find(function(r) identical(r$lhs, "index_modmed"), defined)
  direct_high <- Find(function(r) identical(r$lhs, "direct_high"), defined)
  expect_equal(idx$estimate, a3$estimate * b$estimate, tolerance = 1e-8)
  expect_equal(direct_high$estimate, c1$estimate + c3$estimate,
               tolerance = 1e-8)
})

test_that("lavaan moderated mediation payload carries Model 15 scaffold inputs", {
  # Hayes Model 15 extends Model 14 by also moderating the direct X -> Y path
  # with the same W. The outcome equation therefore has two product terms:
  # M-by-W for the b path and X-by-W for the direct path.
  skip_if_lavaan_unusable()
  suppressWarnings(suppressMessages(library(lavaan)))
  set.seed(31)
  n <- 250
  df <- data.frame(X = rnorm(n), W = rnorm(n))
  df$M <- 0.50 * df$X + rnorm(n)
  df$MW <- df$M * df$W
  df$XW <- df$X * df$W
  df$Y <- 0.40 * df$M + 0.18 * df$W + 0.24 * df$MW +
    0.15 * df$X + 0.20 * df$XW + rnorm(n)

  model <- "
    M ~ a*X
    Y ~ b1*M + b2*W + b3*MW + c1*X + c3*XW
    ind_low       := a*(b1 + b3*(-1))
    ind_mean      := a*b1
    ind_high      := a*(b1 + b3*(1))
    index_modmed  := a*b3
    direct_low    := c1 + c3*(-1)
    direct_mean   := c1
    direct_high   := c1 + c3*(1)
  "
  fit <- sem(model, data = df)
  p <- mellio_payload(fit)

  expect_equal(p$type, "lavaan_mediation")
  params <- p$fields$inspection_zone$parameters

  regressions <- Filter(function(r) identical(r$op, "~"), params)
  reg_keys <- vapply(regressions, function(r) paste(r$lhs, r$rhs, sep = "~"),
                     character(1))
  expect_true(all(c("M~X", "Y~M", "Y~W", "Y~MW",
                    "Y~X", "Y~XW") %in% reg_keys))

  defined <- Filter(function(r) identical(r$op, ":="), params)
  def_keys <- vapply(defined, function(r) r$lhs, character(1))
  expect_true(all(c("ind_low", "ind_mean", "ind_high",
                    "index_modmed",
                    "direct_low", "direct_mean", "direct_high") %in% def_keys))

  a <- Find(function(r) identical(r$lhs, "M") &&
              identical(r$rhs, "X"), regressions)
  b3 <- Find(function(r) identical(r$lhs, "Y") &&
               identical(r$rhs, "MW"), regressions)
  c1 <- Find(function(r) identical(r$lhs, "Y") &&
               identical(r$rhs, "X"), regressions)
  c3 <- Find(function(r) identical(r$lhs, "Y") &&
               identical(r$rhs, "XW"), regressions)
  idx <- Find(function(r) identical(r$lhs, "index_modmed"), defined)
  direct_high <- Find(function(r) identical(r$lhs, "direct_high"), defined)
  expect_equal(idx$estimate, a$estimate * b3$estimate, tolerance = 1e-8)
  expect_equal(direct_high$estimate, c1$estimate + c3$estimate,
               tolerance = 1e-8)
})

test_that("lavaan Model 15 payload preserves awkward centered product names", {
  # Regression coverage for Model 15 aliasing: both product names drop the
  # centered suffixes used by their component variables.
  skip_if_lavaan_unusable()
  suppressWarnings(suppressMessages(library(lavaan)))
  set.seed(37)
  n <- 250
  stress <- rnorm(n)
  support <- rnorm(n)
  df <- data.frame(
    Stress_c = stress - mean(stress),
    Support_c = support - mean(support)
  )
  df$Burnout_c <- 0.50 * df$Stress_c + rnorm(n)
  df$Burnout_x_Support <- df$Burnout_c * df$Support_c
  df$Stress_x_Support <- df$Stress_c * df$Support_c
  df$Turnover <- 0.40 * df$Burnout_c + 0.18 * df$Support_c +
    0.24 * df$Burnout_x_Support + 0.15 * df$Stress_c +
    0.20 * df$Stress_x_Support + rnorm(n)

  model <- "
    Burnout_c ~ a*Stress_c
    Turnover ~ b1*Burnout_c + b2*Support_c + b3*Burnout_x_Support +
      c1*Stress_c + c3*Stress_x_Support
    ind_low       := a*(b1 + b3*(-1))
    ind_mean      := a*b1
    ind_high      := a*(b1 + b3*(1))
    index_modmed  := a*b3
    direct_low    := c1 + c3*(-1)
    direct_mean   := c1
    direct_high   := c1 + c3*(1)
  "
  fit <- sem(model, data = df)
  p <- mellio_payload(fit)

  expect_equal(p$type, "lavaan_mediation")
  params <- p$fields$inspection_zone$parameters

  regressions <- Filter(function(r) identical(r$op, "~"), params)
  reg_keys <- vapply(regressions, function(r) paste(r$lhs, r$rhs, sep = "~"),
                     character(1))
  expect_true(all(c("Burnout_c~Stress_c",
                    "Turnover~Burnout_c",
                    "Turnover~Support_c",
                    "Turnover~Burnout_x_Support",
                    "Turnover~Stress_c",
                    "Turnover~Stress_x_Support") %in% reg_keys))

  defined <- Filter(function(r) identical(r$op, ":="), params)
  idx <- Find(function(r) identical(r$lhs, "index_modmed"), defined)
  a <- Find(function(r) identical(r$lhs, "Burnout_c") &&
              identical(r$rhs, "Stress_c"), regressions)
  b3 <- Find(function(r) identical(r$lhs, "Turnover") &&
               identical(r$rhs, "Burnout_x_Support"), regressions)
  expect_equal(idx$estimate, a$estimate * b3$estimate, tolerance = 1e-8)

  obs <- p$fields$observed_variables
  obs_vars <- vapply(obs, function(r) r$variable, character(1))
  expect_true(all(c("Stress_c", "Support_c", "Burnout_c",
                    "Burnout_x_Support",
                    "Stress_x_Support") %in% obs_vars))
})

test_that("mediation::mediate produces a table Result Card", {
  skip_if_not_installed("mediation")
  set.seed(1)
  d <- data.frame(
    x = rnorm(80),
    m = rnorm(80),
    y = rnorm(80)
  )
  d$m <- 0.4 * d$x + d$m
  d$y <- 0.2 * d$x + 0.5 * d$m + d$y

  res <- mediation::mediate(
    lm(m ~ x, data = d),
    lm(y ~ x + m, data = d),
    treat = "x",
    mediator = "m",
    boot = FALSE,
    sims = 20
  )

  p <- mellio_payload(res)
  expect_equal(p$type, "mediation_mediate")
  expect_equal(p$card_kind, "table")
  expect_equal(p$fields$table_type, "mediation")
  expect_equal(length(p$fields$rows), 4L)
  expect_equal(p$fields$rows[[1]]$effect, "ACME")
  expect_true(is.numeric(p$fields$rows[[1]]$estimate))
  expect_true(is.numeric(p$fields$rows[[1]]$p_value))
  expect_equal(p$fields$treatment, "x")
  expect_equal(p$fields$mediator, "m")
  expect_equal(p$fields$outcome, "y")
  expect_match(p$fields$note, "proportion mediated")
  pkg_names <- vapply(p$packages, `[[`, character(1), "name")
  expect_true("mediation" %in% pkg_names)
})

test_that("mediation labels can preserve original variable names", {
  skip_if_not_installed("mediation")
  set.seed(2)
  d <- data.frame(
    treat_var = rnorm(60),
    mediator_var = rnorm(60),
    outcome_var = rnorm(60)
  )
  d$mediator_var <- 0.4 * d$treat_var + d$mediator_var
  d$outcome_var <- 0.2 * d$treat_var + 0.5 * d$mediator_var + d$outcome_var

  res <- mediation::mediate(
    lm(mediator_var ~ treat_var, data = d),
    lm(outcome_var ~ treat_var + mediator_var, data = d),
    treat = "treat_var",
    mediator = "mediator_var",
    boot = FALSE,
    sims = 20
  )
  attr(res, "mellio_labels") <- list(
    treatment = "ex.days.current",
    mediator = "with_others",
    outcome = "emotion_resp"
  )

  p <- mellio_payload(res)
  expect_equal(p$fields$treatment, "ex.days.current")
  expect_equal(p$fields$mediator, "with_others")
  expect_equal(p$fields$outcome, "emotion_resp")
})

test_that("mediation::mediate extracts a, b, c' path coefficients", {
  # Mediation diagrams need the underlying regression coefficients (a, b,
  # c') in addition to the high-level ACME/ADE/Total decomposition. Verify
  # the extractor pulls them from model.m / model.y and that the values
  # match what coef() returns on those fits.
  skip_if_not_installed("mediation")
  set.seed(123)
  n <- 200
  d <- data.frame(x = rnorm(n))
  d$m <- 0.5 * d$x + rnorm(n)
  d$y <- 0.3 * d$x + 0.4 * d$m + rnorm(n)

  model_m <- lm(m ~ x, data = d)
  model_y <- lm(y ~ x + m, data = d)
  res <- mediation::mediate(
    model_m, model_y, treat = "x", mediator = "m",
    boot = FALSE, sims = 20
  )
  p <- mellio_payload(res)

  expect_true(!is.null(p$fields$paths), info = "fields$paths should be populated")
  expect_equal(length(p$fields$paths), 3L)

  paths_by_id <- stats::setNames(p$fields$paths,
                                 vapply(p$fields$paths, `[[`, character(1), "path"))
  expect_named(paths_by_id, c("a", "b", "c_prime"), ignore.order = TRUE)

  # Path estimates should match coef() on the underlying models exactly.
  expect_equal(paths_by_id$a$estimate,       unname(coef(model_m)["x"]), tolerance = 1e-10)
  expect_equal(paths_by_id$b$estimate,       unname(coef(model_y)["m"]), tolerance = 1e-10)
  expect_equal(paths_by_id$c_prime$estimate, unname(coef(model_y)["x"]), tolerance = 1e-10)

  # SE / p / CI should all be present and numeric for an lm-based fit.
  for (p_id in c("a", "b", "c_prime")) {
    row <- paths_by_id[[p_id]]
    expect_true(is.numeric(row$std_error) && !is.na(row$std_error),
                info = paste(p_id, "std_error missing"))
    expect_true(is.numeric(row$p_value) && !is.na(row$p_value),
                info = paste(p_id, "p_value missing"))
    expect_true(is.numeric(row$ci_lower) && !is.na(row$ci_lower),
                info = paste(p_id, "ci_lower missing"))
    expect_true(is.numeric(row$ci_upper) && !is.na(row$ci_upper),
                info = paste(p_id, "ci_upper missing"))
  }

  # Labels should describe the path direction using the variable names.
  # Skip the arrow-character check (Unicode round-tripping varies between
  # R versions); just verify both variable names appear in the label.
  expect_true(grepl("x", paths_by_id$a$label) &&
              grepl("m", paths_by_id$a$label))
  expect_true(grepl("m", paths_by_id$b$label) &&
              grepl("y", paths_by_id$b$label))
  expect_true(grepl("x", paths_by_id$c_prime$label) &&
              grepl("y", paths_by_id$c_prime$label))
  expect_match(paths_by_id$c_prime$label, "direct")
})

test_that("processR medSummary renders as a Stats table payload", {
  x <- data.frame(
    Effect = c("indirect", "direct", "total", "prop.mediated"),
    equation = c("a*b", "c", "c + a*b", "indirect/total"),
    est = c(0.18, 0.22, 0.40, 0.45),
    ci.lower = c(0.08, 0.05, 0.21, 0.21),
    ci.upper = c(0.30, 0.39, 0.58, 1.12),
    pvalue = c(0.004, 0.020, 0.001, 0.004),
    stringsAsFactors = FALSE
  )
  attr(x, "boot.ci.type") <- "perc"
  attr(x, "se") <- "bootstrap"
  class(x) <- c("medSummary", "data.frame")

  p <- mellio_payload(x, .call = "processR::medSummary(fit)")
  expect_equal(p$type, "processr_mediation_summary")
  expect_equal(p$card_kind, "table")
  expect_equal(p$fields$table_type, "processr_mediation_summary")
  expect_equal(length(p$fields$rows), 4L)
  expect_equal(p$fields$rows[[1]]$effect, "indirect")
  expect_equal(p$fields$rows[[1]]$ci_lower, 0.08)
  expect_equal(p$fields$columns[[3]]$label, "Estimate")
  expect_match(p$fields$note, "table only")
  expect_match(p$fields$note, "proportion mediated is a ratio")
  pkg_names <- vapply(p$packages, `[[`, character(1), "name")
  expect_true("processR" %in% pkg_names)
  expect_true(!is.null(getS3method("mellio_open_dispatch", "medSummary", optional = TRUE)))
})

test_that("processR modmedSummary becomes conditional effect rows", {
  x <- data.frame(
    values = c(-1, 0, 1),
    indirect = c(0.28, 0.18, 0.08),
    lower = c(0.15, 0.09, 0.01),
    upper = c(0.42, 0.29, 0.17),
    indirectp = c(0.001, 0.004, 0.040),
    direct = c(0.52, 0.26, 0.01),
    lowerd = c(0.31, 0.08, -0.18),
    upperd = c(0.75, 0.46, 0.20),
    directp = c(0.001, 0.010, 0.920),
    label = c("M1", "M1", "M1"),
    stringsAsFactors = FALSE
  )
  attr(x, "mod") <- "Support_c"
  attr(x, "indirect") <- "(a1+a3*W)*b"
  attr(x, "direct") <- "c1+c3*W"
  attr(x, "boot.ci.type") <- "perc"
  class(x) <- c("modmedSummary", "data.frame")

  p <- mellio_payload(x, .call = "processR::modmedSummary(fit)")
  expect_equal(p$type, "processr_moderated_mediation_summary")
  expect_equal(p$fields$columns[[1]]$label, "Support (W)")
  expect_false(any(vapply(p$fields$columns, function(col) {
    identical(col$key, "label")
  }, logical(1))))
  expect_equal(length(p$fields$rows), 6L)
  expect_equal(p$fields$rows[[1]]$effect, "Conditional indirect effect")
  expect_equal(p$fields$rows[[2]]$effect, "Conditional direct effect")
  expect_equal(p$fields$rows[[6]]$p_value, 0.920)
  expect_match(p$fields$note, "processR moderated-mediation summary")
  expect_true(!is.null(getS3method("mellio_open_dispatch", "modmedSummary", optional = TRUE)))
})

test_that("bruceR PROCESS Model 4 list renders as a simple mediation payload", {
  set.seed(404)
  n <- 90
  dat <- data.frame(
    Stress = stats::rnorm(n),
    Tenure = stats::rnorm(n)
  )
  dat$Burnout <- 0.65 * dat$Stress - 0.10 * dat$Tenure + stats::rnorm(n, sd = 0.50)
  dat$TurnoverIntent <- -0.05 * dat$Stress + 0.50 * dat$Burnout +
    0.02 * dat$Tenure + stats::rnorm(n, sd = 0.60)

  model_m <- stats::lm(Burnout ~ Tenure + Stress, data = dat)
  model_y <- stats::lm(TurnoverIntent ~ Tenure + Stress + Burnout, data = dat)
  a <- unname(stats::coef(model_m)[["Stress"]])
  b <- unname(stats::coef(model_y)[["Burnout"]])
  c_prime <- unname(stats::coef(model_y)[["Stress"]])
  indirect <- a * b
  total <- indirect + c_prime

  med <- data.frame(
    Effect = c(indirect, c_prime, total),
    "S.E." = c(0.043, 0.049, 0.043),
    LLCI = c(indirect - 0.08, c_prime - 0.10, total - 0.08),
    ULCI = c(indirect + 0.08, c_prime + 0.10, total + 0.08),
    p = c(0.000, 0.336, 0.000),
    z = c(7.55, -0.94, 6.43),
    pval = c(4.4e-14, 0.346, 1.3e-10),
    check.names = FALSE
  )
  rownames(med) <- c("Indirect (ab)", "Direct (c')", "Total (c)")
  x <- list(
    process.id = 4,
    process.type = "Simple Mediation",
    model.m = list(model.m.1 = model_m),
    model.y = model_y,
    results = list(list(mediation = med))
  )

  p <- mellio_payload(x, .call = "bruceR::PROCESS(data, y, x, meds, model = 4)")
  expect_equal(p$type, "mediation_mediate")
  expect_equal(p$card_kind, "table")
  expect_equal(p$fields$table_type, "mediation")
  expect_equal(p$fields$source, "bruceR::PROCESS")
  expect_equal(p$fields$process_id, 4L)
  expect_equal(p$fields$treatment, "Stress")
  expect_equal(p$fields$mediator, "Burnout")
  expect_equal(p$fields$outcome, "TurnoverIntent")
  expect_equal(as.character(p$fields$controls), "Tenure")
  expect_equal(p$fields$n, n)
  expect_true(isTRUE(p$fields$boot))
  expect_equal(p$fields$ci_type, "perc")
  expect_equal(p$fields$conf_level, 0.95)
  expect_equal(length(p$fields$rows), 3L)
  expect_equal(p$fields$rows[[1]]$effect, "Indirect effect")
  expect_equal(p$fields$rows[[1]]$raw_effect, "Indirect (ab)")
  expect_true(is.na(p$fields$rows[[1]]$p_value))
  expect_true(is.na(p$fields$rows[[1]]$statistic))
  expect_equal(p$fields$rows[[1]]$statistic_type, "bootstrap_ci")
  expect_equal(p$fields$columns[[6]]$label, "Bootstrap 95% CI")
  expect_match(p$fields$note, "path coefficients")
  expect_match(p$fields$note, "z and p are omitted")
  expect_match(p$fields$note, "N = 90", fixed = TRUE)
  expect_match(p$fields$effect_note, "mediation-effect intervals")
  expect_match(p$fields$effect_note, "z and p are omitted")
  expect_match(p$fields$effect_note, "N = 90", fixed = TRUE)
  expect_match(p$fields$path_note, "lm models")
  expect_match(p$fields$path_note, "N = 90", fixed = TRUE)
  expect_match(p$fields$path_note, "R\u00b2 by regression", fixed = TRUE)
  expect_false(grepl("R^2", p$fields$path_note, fixed = TRUE))
  expect_equal(length(p$fields$path_rows), 5L)
  expect_equal(vapply(p$fields$path_rows[1:3], `[[`, character(1), "path"),
               c("a", "b", "c_prime"))
  expect_true(any(vapply(p$fields$path_rows, function(row) {
    identical(row$parameter, "Tenure \u2192 Burnout")
  }, logical(1))))
  expect_true(any(vapply(p$fields$path_rows, function(row) {
    identical(row$parameter, "Tenure \u2192 TurnoverIntent")
  }, logical(1))))
  expect_equal(p$fields$path_columns[[4]]$label, "t")
  paths <- stats::setNames(p$fields$paths, vapply(p$fields$paths, `[[`, character(1), "path"))
  expect_true(all(c("a", "b", "c_prime") %in% names(paths)))
  expect_equal(paths$a$estimate, a)
  expect_equal(paths$b$estimate, b)
  expect_equal(paths$c_prime$estimate, c_prime)
  pkg_names <- vapply(p$packages, `[[`, character(1), "name")
  expect_true("bruceR" %in% pkg_names)
})

test_that("bruceR PROCESS Model 4 list renders parallel mediation", {
  set.seed(414)
  n <- 92
  dat <- data.frame(
    Stress = stats::rnorm(n),
    Tenure = stats::rnorm(n)
  )
  dat$Burnout <- 0.58 * dat$Stress - 0.08 * dat$Tenure +
    stats::rnorm(n, sd = 0.60)
  dat$Exhaustion <- 0.32 * dat$Stress - 0.06 * dat$Tenure +
    stats::rnorm(n, sd = 0.70)
  dat$TurnoverIntent <- 0.06 * dat$Stress + 0.42 * dat$Burnout +
    0.35 * dat$Exhaustion - 0.04 * dat$Tenure +
    stats::rnorm(n, sd = 0.75)

  model_m1 <- stats::lm(Burnout ~ Tenure + Stress, data = dat)
  model_m2 <- stats::lm(Exhaustion ~ Tenure + Stress, data = dat)
  model_y <- stats::lm(TurnoverIntent ~ Tenure + Stress + Burnout + Exhaustion,
                       data = dat)
  a1 <- unname(stats::coef(model_m1)[["Stress"]])
  a2 <- unname(stats::coef(model_m2)[["Stress"]])
  b1 <- unname(stats::coef(model_y)[["Burnout"]])
  b2 <- unname(stats::coef(model_y)[["Exhaustion"]])
  c_prime <- unname(stats::coef(model_y)[["Stress"]])
  ind_1 <- a1 * b1
  ind_2 <- a2 * b2
  total_indirect <- ind_1 + ind_2

  med1 <- data.frame(
    Effect = c(ind_1, c_prime),
    "S.E." = c(0.030, 0.050),
    LLCI = c(ind_1 - 0.05, c_prime - 0.10),
    ULCI = c(ind_1 + 0.05, c_prime + 0.10),
    p = c(0, 0.44),
    z = c(2.8, 0.7),
    pval = c(0.005, 0.480),
    check.names = FALSE
  )
  rownames(med1) <- c("Indirect (ab)", "Direct (c')")
  med2 <- data.frame(
    Effect = c(ind_2, c_prime),
    "S.E." = c(0.020, 0.050),
    LLCI = c(0.002, c_prime - 0.10),
    ULCI = c(ind_2 + 0.04, c_prime + 0.10),
    p = c(0, 0.44),
    z = c(1.6, 0.7),
    pval = c(0.103, 0.480),
    check.names = FALSE
  )
  rownames(med2) <- c("Indirect (ab)", "Direct (c')")

  x <- list(
    process.id = 4,
    process.type = "Parallel Multiple Mediation (2 meds)",
    model.m = list(model.m.1 = model_m1, model.m.2 = model_m2),
    model.y = model_y,
    results = list(list(mediation = med1), list(mediation = med2)),
    nsim = 1000
  )

  p <- mellio_payload(
    x,
    .call = "bruceR::PROCESS(data, y, x, meds, med.type = 'parallel')"
  )
  expect_equal(p$type, "mediation_mediate")
  expect_equal(p$type_label, "Parallel mediation analysis (bruceR PROCESS Model 4)")
  expect_equal(p$card_kind, "table")
  expect_equal(p$fields$source, "bruceR::PROCESS")
  expect_equal(p$fields$process_id, 4L)
  expect_equal(p$fields$topology, "parallel")
  expect_equal(p$fields$hayes_model, "4")
  expect_equal(p$fields$treatment, "Stress")
  expect_equal(p$fields$mediator, "Burnout")
  expect_equal(as.character(p$fields$mediators), c("Burnout", "Exhaustion"))
  expect_equal(p$fields$outcome, "TurnoverIntent")
  expect_equal(as.character(p$fields$controls), "Tenure")
  expect_equal(p$fields$n, n)
  expect_true(isTRUE(p$fields$boot))
  expect_equal(p$fields$sims, 1000L)

  expect_equal(length(p$fields$rows), 5L)
  expect_equal(vapply(p$fields$rows, `[[`, character(1), "effect"),
               c("Indirect effect via Burnout",
                 "Indirect effect via Exhaustion",
                 "Total indirect effect",
                 "Direct effect",
                 "Total effect"))
  expect_equal(p$fields$rows[[1]]$estimate, ind_1, tolerance = 1e-8)
  expect_equal(p$fields$rows[[2]]$estimate, ind_2, tolerance = 1e-8)
  expect_true(is.na(p$fields$rows[[2]]$p_value))
  expect_true(is.na(p$fields$rows[[2]]$statistic))
  expect_equal(p$fields$rows[[2]]$statistic_type, "bootstrap_ci")
  expect_equal(p$fields$rows[[2]]$ci_lower, 0.002)
  expect_equal(p$fields$rows[[3]]$estimate, total_indirect, tolerance = 1e-8)
  expect_true(isTRUE(p$fields$rows[[3]]$synthesised))
  expect_equal(p$fields$rows[[3]]$ci_note, "not returned")
  expect_equal(p$fields$rows[[4]]$estimate, c_prime, tolerance = 1e-8)
  expect_equal(p$fields$rows[[5]]$estimate, c_prime + total_indirect,
               tolerance = 1e-8)
  expect_match(p$fields$note, "specific indirect-effect intervals")
  expect_match(p$fields$note, "z and p are omitted")
  expect_match(p$fields$note, "total indirect and total effects are computed")
  expect_match(p$fields$effect_note, "Parallel mediation effects")
  expect_match(p$fields$effect_note, "z and p are omitted")
  expect_match(p$fields$path_note, "R\u00b2 by regression", fixed = TRUE)

  paths <- stats::setNames(p$fields$path_rows,
                           vapply(p$fields$path_rows, `[[`, character(1), "path"))
  expect_true(all(c("a1", "a2", "b1", "b2", "c_prime",
                    "control_m1_Tenure", "control_m2_Tenure",
                    "control_y_Tenure") %in% names(paths)))
  expect_equal(paths$a1$estimate, a1)
  expect_equal(paths$a2$estimate, a2)
  expect_equal(paths$b1$estimate, b1)
  expect_equal(paths$b2$estimate, b2)
  expect_equal(paths$c_prime$estimate, c_prime)
  expect_equal(paths$a2$parameter,
               "a2: Stress \u2192 Exhaustion")
  main_paths <- stats::setNames(p$fields$paths,
                                vapply(p$fields$paths, `[[`, character(1), "path"))
  expect_true(all(c("a1", "a2", "b1", "b2", "c_prime") %in% names(main_paths)))
  pkg_names <- vapply(p$packages, `[[`, character(1), "name")
  expect_true(all(c("bruceR", "mediation") %in% pkg_names))
})

test_that("bruceR PROCESS Model 6 list renders serial mediation", {
  set.seed(606)
  n <- 95
  dat <- data.frame(
    Stress = stats::rnorm(n),
    Tenure = stats::rnorm(n)
  )
  dat$Burnout <- 0.55 * dat$Stress - 0.08 * dat$Tenure +
    stats::rnorm(n, sd = 0.60)
  dat$Exhaustion <- 0.20 * dat$Stress + 0.45 * dat$Burnout -
    0.07 * dat$Tenure + stats::rnorm(n, sd = 0.65)
  dat$TurnoverIntent <- 0.08 * dat$Stress + 0.25 * dat$Burnout +
    0.50 * dat$Exhaustion - 0.05 * dat$Tenure +
    stats::rnorm(n, sd = 0.70)

  model_m1 <- stats::lm(Burnout ~ Tenure + Stress, data = dat)
  model_m2 <- stats::lm(Exhaustion ~ Tenure + Stress + Burnout, data = dat)
  model_y <- stats::lm(TurnoverIntent ~ Tenure + Stress + Burnout + Exhaustion,
                       data = dat)
  a1 <- unname(stats::coef(model_m1)[["Stress"]])
  a2 <- unname(stats::coef(model_m2)[["Stress"]])
  d21 <- unname(stats::coef(model_m2)[["Burnout"]])
  b1 <- unname(stats::coef(model_y)[["Burnout"]])
  b2 <- unname(stats::coef(model_y)[["Exhaustion"]])
  c_prime <- unname(stats::coef(model_y)[["Stress"]])
  ind_m1 <- a1 * b1
  ind_m2 <- a2 * b2
  serial <- a1 * d21 * b2
  total_indirect <- ind_m1 + ind_m2 + serial
  total <- c_prime + total_indirect

  lav <- data.frame(
    Estimate = c(total_indirect, ind_m1, ind_m2, serial, c_prime, total),
    "S.E." = c(0.07, 0.04, 0.04, 0.05, 0.06, 0.08),
    z = c(5.1, 3.0, 2.4, 3.3, 1.2, 4.8),
    pval = c(0.000001, 0.003, 0.016, 0.001, 0.230, 0.000002),
    BootLLCI = c(total_indirect - 0.10, ind_m1 - 0.07, ind_m2 - 0.06,
                 serial - 0.07, c_prime - 0.09, total - 0.11),
    BootULCI = c(total_indirect + 0.10, ind_m1 + 0.07, ind_m2 + 0.06,
                 serial + 0.07, c_prime + 0.09, total + 0.11),
    Beta = NA_real_,
    check.names = FALSE
  )
  rownames(lav) <- c(" Indirect_All", " Ind_X_M1_Y", " Ind_X_M2_Y",
                     " Ind_X_M1_M2_Y", " Direct", " Total")

  x <- list(
    process.id = 6,
    process.type = "Serial Multiple Mediation (2 meds)",
    model.m = list(model.m.1 = model_m1, model.m.2 = model_m2),
    model.y = model_y,
    results = list(
      list(
        lavaan.syntax = paste(
          "Burnout ~ Tenure + a1*Stress",
          "Exhaustion ~ Tenure + a2*Stress + d12*Burnout",
          "TurnoverIntent ~ Tenure + c.*Stress + b1*Burnout + b2*Exhaustion",
          "Indirect_All := a1*b1 + a2*b2 + a1*d12*b2",
          sep = "\n"
        ),
        lavaan.mediation = lav
      )
    ),
    nsim = 1000
  )

  p <- mellio_payload(
    x,
    .call = "bruceR::PROCESS(data, y, x, meds, med.type = 'serial')"
  )
  expect_equal(p$type, "mediation_mediate")
  expect_equal(p$type_label, "Serial mediation analysis (bruceR PROCESS Model 6)")
  expect_equal(p$card_kind, "table")
  expect_equal(p$fields$source, "bruceR::PROCESS")
  expect_equal(p$fields$process_id, 6L)
  expect_equal(p$fields$topology, "serial")
  expect_equal(p$fields$hayes_model, "6")
  expect_equal(p$fields$treatment, "Stress")
  expect_equal(p$fields$mediator, "Burnout")
  expect_equal(as.character(p$fields$mediators), c("Burnout", "Exhaustion"))
  expect_equal(p$fields$outcome, "TurnoverIntent")
  expect_equal(as.character(p$fields$controls), "Tenure")
  expect_equal(p$fields$n, n)
  expect_true(isTRUE(p$fields$boot))
  expect_equal(p$fields$sims, 1000L)

  expect_equal(length(p$fields$rows), 6L)
  expect_equal(vapply(p$fields$rows, `[[`, character(1), "effect"),
               c("Indirect effect via Burnout",
                 "Indirect effect via Exhaustion",
                 "Serial indirect effect",
                 "Total indirect effect",
                 "Direct effect",
                 "Total effect"))
  expect_equal(p$fields$rows[[1]]$raw_effect, "Ind_X_M1_Y")
  expect_equal(p$fields$rows[[3]]$estimate, serial, tolerance = 1e-8)
  expect_equal(p$fields$rows[[4]]$estimate, total_indirect, tolerance = 1e-8)
  expect_equal(p$fields$rows[[5]]$estimate, c_prime, tolerance = 1e-8)
  expect_equal(p$fields$rows[[1]]$ci_lower, ind_m1 - 0.07)
  expect_match(p$fields$note, "serial mediation-effect intervals")
  expect_match(p$fields$effect_note, "Serial mediation effects")
  expect_match(p$fields$path_note, "R\u00b2 by regression", fixed = TRUE)

  expect_equal(length(p$fields$serial_edges), 1L)
  expect_equal(p$fields$serial_edges[[1]]$from, "Burnout")
  expect_equal(p$fields$serial_edges[[1]]$to, "Exhaustion")

  paths <- stats::setNames(p$fields$path_rows,
                           vapply(p$fields$path_rows, `[[`, character(1), "path"))
  expect_true(all(c("a1", "a2", "d21", "b1", "b2", "c_prime",
                    "control_m1_Tenure", "control_m2_Tenure",
                    "control_y_Tenure") %in% names(paths)))
  expect_equal(paths$a1$estimate, a1)
  expect_equal(paths$a2$estimate, a2)
  expect_equal(paths$d21$estimate, d21)
  expect_equal(paths$b1$estimate, b1)
  expect_equal(paths$b2$estimate, b2)
  expect_equal(paths$c_prime$estimate, c_prime)
  expect_equal(paths$d21$parameter,
               "d21: Burnout \u2192 Exhaustion")

  main_paths <- stats::setNames(p$fields$paths,
                                vapply(p$fields$paths, `[[`, character(1), "path"))
  expect_true(all(c("a1", "a2", "d21", "b1", "b2", "c_prime") %in% names(main_paths)))
  pkg_names <- vapply(p$packages, `[[`, character(1), "name")
  expect_true(all(c("bruceR", "lavaan") %in% pkg_names))
})

test_that("bruceR PROCESS Model 7 list renders moderated mediation", {
  set.seed(707)
  n <- 100
  dat <- data.frame(
    Stress = stats::rnorm(n),
    Support = stats::rnorm(n),
    Tenure = stats::rnorm(n)
  )
  dat$Burnout <- 0.55 * dat$Stress - 0.20 * dat$Support -
    0.25 * dat$Stress * dat$Support - 0.08 * dat$Tenure +
    stats::rnorm(n, sd = 0.70)
  dat$TurnoverIntent <- 0.45 * dat$Burnout + 0.12 * dat$Stress +
    0.03 * dat$Support - 0.04 * dat$Tenure +
    stats::rnorm(n, sd = 0.80)

  model_m <- stats::lm(Burnout ~ Tenure + Stress * Support, data = dat)
  model_y <- stats::lm(TurnoverIntent ~ Tenure + Stress + Support + Burnout,
                       data = dat)
  a1 <- unname(stats::coef(model_m)[["Stress"]])
  a3 <- unname(stats::coef(model_m)[["Stress:Support"]])
  b <- unname(stats::coef(model_y)[["Burnout"]])
  low <- (a1 + a3 * -1) * b
  mean_eff <- a1 * b
  high <- (a1 + a3 * 1) * b

  med <- data.frame(
    Effect = c(low, mean_eff, high),
    "S.E." = c(0.09, 0.06, 0.05),
    LLCI = c(low - 0.10, mean_eff - 0.08, high - 0.07),
    ULCI = c(low + 0.10, mean_eff + 0.08, high + 0.07),
    p = c(0, 0, 0.03),
    z = c(4.8, 4.1, 2.2),
    pval = c(0.000002, 0.00004, 0.028),
    check.names = FALSE
  )
  med <- cbind(
    setNames(data.frame(c("-1.000 (- SD)", "0.000 (Mean)", "1.000 (+ SD) "),
                        check.names = FALSE),
             "\"Support\"    "),
    med
  )
  x <- list(
    process.id = 7,
    process.type = "Moderated Mediation",
    model.m = list(model.m.1 = model_m),
    model.y = model_y,
    results = list(
      list(conditional = data.frame(Support = c(-1, 0, 1))),
      list(mediation = med)
    )
  )

  p <- mellio_payload(
    x,
    .call = "bruceR::PROCESS(data, y, x, meds, mods, mod.path = 'x-m')"
  )
  expect_equal(p$type, "mediation_mediate")
  expect_equal(p$type_label, "Moderated mediation analysis (bruceR PROCESS Model 7)")
  expect_equal(p$fields$process_id, 7L)
  expect_equal(p$fields$topology, "moderated_mediation")
  expect_equal(p$fields$hayes_model, "7")
  expect_equal(p$fields$treatment, "Stress")
  expect_equal(p$fields$moderator, "Support")
  expect_equal(p$fields$moderation$moderated_path, "a")
  expect_equal(p$fields$moderation$interaction_term, "Stress:Support")
  expect_equal(p$fields$mediator, "Burnout")
  expect_equal(p$fields$outcome, "TurnoverIntent")
  expect_equal(as.character(p$fields$controls), "Tenure")
  expect_equal(p$fields$n, n)

  expect_equal(length(p$fields$rows), 4L)
  expect_equal(vapply(p$fields$rows[1:3], `[[`, character(1), "condition"),
               c("low", "mean", "high"))
  expect_equal(p$fields$rows[[1]]$effect,
               "Conditional indirect effect (low Support)")
  expect_equal(p$fields$rows[[4]]$effect,
               "Index of moderated mediation")
  expect_equal(p$fields$rows[[4]]$estimate, a3 * b, tolerance = 1e-8)
  expect_true(is.na(p$fields$rows[[4]]$ci_lower))
  expect_match(p$fields$effect_note, "point estimate")

  paths <- stats::setNames(p$fields$path_rows,
                           vapply(p$fields$path_rows, `[[`, character(1), "path"))
  expect_true(all(c("a1", "a2", "a3", "b", "c_prime", "c2",
                    "control_m_Tenure", "control_y_Tenure") %in% names(paths)))
  expect_equal(paths$a1$estimate, a1)
  expect_equal(paths$a3$estimate, a3)
  expect_equal(paths$b$estimate, b)
  expect_equal(paths$a3$parameter,
               "a3: Stress \u00d7 Support \u2192 Burnout")
  expect_true(is.null(paths$a3[["std_estimate"]]))
  expect_match(p$fields$path_note, "\u03b2 is not reported for interaction paths",
               fixed = TRUE)
  expect_match(p$fields$path_note, "R\u00b2 by regression", fixed = TRUE)
})

test_that("bruceR PROCESS Model 8 list renders conditional indirect and direct effects", {
  set.seed(808)
  n <- 100
  dat <- data.frame(
    Stress = stats::rnorm(n),
    Support = stats::rnorm(n),
    Tenure = stats::rnorm(n)
  )
  dat$Burnout <- 0.55 * dat$Stress - 0.20 * dat$Support -
    0.25 * dat$Stress * dat$Support - 0.08 * dat$Tenure +
    stats::rnorm(n, sd = 0.70)
  dat$TurnoverIntent <- 0.45 * dat$Burnout + 0.12 * dat$Stress -
    0.18 * dat$Stress * dat$Support + 0.03 * dat$Support -
    0.04 * dat$Tenure + stats::rnorm(n, sd = 0.80)

  model_m <- stats::lm(Burnout ~ Tenure + Stress * Support, data = dat)
  model_y <- stats::lm(TurnoverIntent ~ Tenure + Stress * Support + Burnout,
                       data = dat)
  a1 <- unname(stats::coef(model_m)[["Stress"]])
  a3 <- unname(stats::coef(model_m)[["Stress:Support"]])
  b <- unname(stats::coef(model_y)[["Burnout"]])
  c1 <- unname(stats::coef(model_y)[["Stress"]])
  c3 <- unname(stats::coef(model_y)[["Stress:Support"]])
  low_ind <- (a1 + a3 * -1) * b
  mean_ind <- a1 * b
  high_ind <- (a1 + a3 * 1) * b
  low_direct <- c1 + c3 * -1
  mean_direct <- c1
  high_direct <- c1 + c3 * 1
  condition_values <- c("-1.000 (- SD)", "0.000 (Mean)", "1.000 (+ SD) ")

  direct_slopes <- data.frame(
    Effect = c(low_direct, mean_direct, high_direct),
    "S.E." = c(0.08, 0.07, 0.09),
    LLCI = c(low_direct - 0.10, mean_direct - 0.08, high_direct - 0.11),
    ULCI = c(low_direct + 0.10, mean_direct + 0.08, high_direct + 0.11),
    t = c(3.1, 1.8, -0.4),
    pval = c(0.002, 0.074, 0.690),
    check.names = FALSE
  )
  direct_slopes <- cbind(
    setNames(data.frame(condition_values, check.names = FALSE),
             "\"Support\"    "),
    direct_slopes
  )
  a_slopes <- data.frame(
    Effect = c(a1 + a3 * -1, a1, a1 + a3 * 1),
    "S.E." = c(0.07, 0.05, 0.08),
    LLCI = c(a1 + a3 * -1 - 0.09, a1 - 0.07, a1 + a3 * 1 - 0.10),
    ULCI = c(a1 + a3 * -1 + 0.09, a1 + 0.07, a1 + a3 * 1 + 0.10),
    t = c(8.0, 9.2, 5.1),
    pval = c(0, 0, 0.00001),
    check.names = FALSE
  )
  a_slopes <- cbind(
    setNames(data.frame(condition_values, check.names = FALSE),
             "\"Support\"    "),
    a_slopes
  )
  med <- data.frame(
    Effect = c(low_ind, mean_ind, high_ind),
    "S.E." = c(0.09, 0.06, 0.05),
    LLCI = c(low_ind - 0.10, mean_ind - 0.08, high_ind - 0.07),
    ULCI = c(low_ind + 0.10, mean_ind + 0.08, high_ind + 0.07),
    p = c(0, 0, 0.03),
    z = c(4.8, 4.1, 2.2),
    pval = c(0.000002, 0.00004, 0.028),
    check.names = FALSE
  )
  med <- cbind(
    setNames(data.frame(condition_values, check.names = FALSE),
             "\"Support\"    "),
    med
  )
  x <- list(
    process.id = 8,
    process.type = "Moderated Mediation",
    model.m = list(model.m.1 = model_m),
    model.y = model_y,
    results = list(
      list(
        conditional = data.frame(Support = c(-1, 0, 1)),
        simple.slopes = direct_slopes
      ),
      list(
        conditional = data.frame(Support = c(-1, 0, 1)),
        simple.slopes = a_slopes
      ),
      list(mediation = med)
    )
  )

  p <- mellio_payload(
    x,
    .call = "bruceR::PROCESS(data, y, x, meds, mods, mod.path = c('x-m', 'x-y'))"
  )
  expect_equal(p$type, "mediation_mediate")
  expect_equal(p$type_label, "Moderated mediation analysis (bruceR PROCESS Model 8)")
  expect_equal(p$fields$process_id, 8L)
  expect_equal(p$fields$topology, "moderated_mediation")
  expect_equal(p$fields$hayes_model, "8")
  expect_equal(p$fields$moderation$moderated_path, "a")
  expect_equal(p$fields$moderation$direct_moderated_path, "X_to_Y")
  expect_equal(p$fields$moderation$interaction_term, "Stress:Support")
  expect_equal(p$fields$moderation$direct_interaction_term, "Stress:Support")
  expect_equal(as.character(p$fields$controls), "Tenure")
  expect_equal(p$fields$n, n)

  expect_equal(length(p$fields$rows), 7L)
  expect_equal(vapply(p$fields$rows[1:3], `[[`, character(1), "condition"),
               c("low", "mean", "high"))
  expect_equal(p$fields$rows[[4]]$effect,
               "Index of moderated mediation")
  expect_equal(p$fields$rows[[4]]$estimate, a3 * b, tolerance = 1e-8)
  expect_equal(vapply(p$fields$rows[5:7], `[[`, character(1), "condition"),
               c("low", "mean", "high"))
  expect_equal(p$fields$rows[[5]]$effect,
               "Conditional direct effect (low Support)")
  expect_equal(p$fields$rows[[5]]$estimate, low_direct, tolerance = 1e-8)
  expect_equal(p$fields$rows[[5]]$statistic_type, "t")
  expect_equal(p$fields$columns[[4]]$label, "test")
  expect_match(p$fields$effect_note, "conditional direct intervals")
  expect_match(p$fields$effect_note, "point estimate")

  paths <- stats::setNames(p$fields$path_rows,
                           vapply(p$fields$path_rows, `[[`, character(1), "path"))
  expect_true(all(c("a1", "a2", "a3", "b", "c1", "c2", "c3",
                    "control_m_Tenure", "control_y_Tenure") %in% names(paths)))
  expect_equal(paths$a1$estimate, a1)
  expect_equal(paths$a3$estimate, a3)
  expect_equal(paths$b$estimate, b)
  expect_equal(paths$c1$estimate, c1)
  expect_equal(paths$c3$estimate, c3)
  expect_equal(paths$c3$parameter,
               "c3: Stress \u00d7 Support \u2192 TurnoverIntent")
  expect_true(is.null(paths$a3[["std_estimate"]]))
  expect_true(is.null(paths$c3[["std_estimate"]]))
  expect_match(p$fields$path_note, "\u03b2 is not reported for interaction paths",
               fixed = TRUE)
})

test_that("bruceR PROCESS Model 14 list renders b-path moderated mediation", {
  set.seed(1414)
  n <- 100
  dat <- data.frame(
    Stress = stats::rnorm(n),
    Support = stats::rnorm(n),
    Tenure = stats::rnorm(n)
  )
  dat$Burnout <- 0.55 * dat$Stress + 0.08 * dat$Support -
    0.08 * dat$Tenure + stats::rnorm(n, sd = 0.70)
  dat$TurnoverIntent <- 0.20 * dat$Stress + 0.50 * dat$Burnout -
    0.25 * dat$Burnout * dat$Support + 0.05 * dat$Support -
    0.04 * dat$Tenure + stats::rnorm(n, sd = 0.80)

  model_m <- stats::lm(Burnout ~ Tenure + Stress + Support, data = dat)
  model_y <- stats::lm(TurnoverIntent ~ Tenure + Stress +
                         Burnout * Support, data = dat)
  a <- unname(stats::coef(model_m)[["Stress"]])
  b1 <- unname(stats::coef(model_y)[["Burnout"]])
  b3 <- unname(stats::coef(model_y)[["Burnout:Support"]])
  low <- a * (b1 + b3 * -1)
  mean_eff <- a * b1
  high <- a * (b1 + b3 * 1)

  simple_slopes <- data.frame(
    Effect = c(b1 + b3 * -1, b1, b1 + b3 * 1),
    "S.E." = c(0.09, 0.07, 0.08),
    LLCI = c(b1 + b3 * -1 - 0.10, b1 - 0.08, b1 + b3 * 1 - 0.09),
    ULCI = c(b1 + b3 * -1 + 0.10, b1 + 0.08, b1 + b3 * 1 + 0.09),
    t = c(5.2, 4.8, 2.9),
    pval = c(0.000002, 0.00001, 0.004),
    check.names = FALSE
  )
  condition_values <- c("-1.000 (- SD)", "0.000 (Mean)", "1.000 (+ SD) ")
  simple_slopes <- cbind(
    setNames(data.frame(condition_values, check.names = FALSE),
             "\"Support\"    "),
    simple_slopes
  )

  med <- data.frame(
    Effect = c(low, mean_eff, high),
    "S.E." = c(0.08, 0.06, 0.05),
    LLCI = c(low - 0.10, mean_eff - 0.08, high - 0.07),
    ULCI = c(low + 0.10, mean_eff + 0.08, high + 0.07),
    p = c(0, 0, 0.02),
    z = c(5.1, 4.7, 2.4),
    pval = c(0.000001, 0.000003, 0.016),
    check.names = FALSE
  )
  med <- cbind(
    setNames(data.frame(condition_values, check.names = FALSE),
             "\"Support\"    "),
    med
  )
  x <- list(
    process.id = 14,
    process.type = "Moderated Mediation",
    model.m = list(model.m.1 = model_m),
    model.y = model_y,
    results = list(
      list(
        conditional = data.frame(Support = c(-1, 0, 1)),
        simple.slopes = simple_slopes
      ),
      list(mediation = med)
    )
  )

  p <- mellio_payload(
    x,
    .call = "bruceR::PROCESS(data, y, x, meds, mods, mod.path = 'm-y')"
  )
  expect_equal(p$type, "mediation_mediate")
  expect_equal(p$type_label, "Moderated mediation analysis (bruceR PROCESS Model 14)")
  expect_equal(p$fields$process_id, 14L)
  expect_equal(p$fields$topology, "moderated_mediation")
  expect_equal(p$fields$hayes_model, "14")
  expect_equal(p$fields$treatment, "Stress")
  expect_equal(p$fields$moderator, "Support")
  expect_equal(p$fields$moderation$moderated_path, "b")
  expect_equal(p$fields$moderation$path, "M_to_Y")
  expect_equal(p$fields$moderation$interaction_term, "Burnout:Support")
  expect_equal(p$fields$mediator, "Burnout")
  expect_equal(p$fields$outcome, "TurnoverIntent")
  expect_equal(as.character(p$fields$controls), "Tenure")
  expect_equal(p$fields$n, n)

  expect_equal(length(p$fields$rows), 4L)
  expect_equal(vapply(p$fields$rows[1:3], `[[`, character(1), "condition"),
               c("low", "mean", "high"))
  expect_equal(p$fields$rows[[1]]$effect,
               "Conditional indirect effect (low Support)")
  expect_equal(p$fields$rows[[4]]$effect,
               "Index of moderated mediation")
  expect_equal(p$fields$rows[[4]]$estimate, a * b3, tolerance = 1e-8)
  expect_true(is.na(p$fields$rows[[4]]$ci_lower))
  expect_match(p$fields$effect_note, "a \u00d7 b3", fixed = TRUE)

  paths <- stats::setNames(p$fields$path_rows,
                           vapply(p$fields$path_rows, `[[`, character(1), "path"))
  expect_true(all(c("a", "b1", "b2", "b3", "c_prime",
                    "control_m_Support", "control_m_Tenure",
                    "control_y_Tenure") %in% names(paths)))
  expect_equal(paths$a$estimate, a)
  expect_equal(paths$b1$estimate, b1)
  expect_equal(paths$b3$estimate, b3)
  expect_equal(paths$b3$parameter,
               "b3: Burnout \u00d7 Support \u2192 TurnoverIntent")
  expect_true(is.null(paths$b3[["std_estimate"]]))
  expect_match(p$fields$path_note, "\u03b2 is not reported for interaction paths",
               fixed = TRUE)
})

test_that("bruceR PROCESS Model 15 list renders b-path and direct moderated mediation", {
  set.seed(1515)
  n <- 100
  dat <- data.frame(
    Stress = stats::rnorm(n),
    Support = stats::rnorm(n),
    Tenure = stats::rnorm(n)
  )
  dat$Burnout <- 0.55 * dat$Stress + 0.08 * dat$Support -
    0.08 * dat$Tenure + stats::rnorm(n, sd = 0.70)
  dat$TurnoverIntent <- 0.20 * dat$Stress + 0.50 * dat$Burnout -
    0.25 * dat$Burnout * dat$Support -
    0.18 * dat$Stress * dat$Support + 0.05 * dat$Support -
    0.04 * dat$Tenure + stats::rnorm(n, sd = 0.80)

  model_m <- stats::lm(Burnout ~ Tenure + Stress + Support, data = dat)
  model_y <- stats::lm(TurnoverIntent ~ Tenure + Stress * Support +
                         Burnout * Support, data = dat)
  coefs <- stats::coef(model_y)
  a <- unname(stats::coef(model_m)[["Stress"]])
  b1 <- unname(coefs[["Burnout"]])
  b_interaction <- grep("^(Burnout:Support|Support:Burnout)$",
                        names(coefs), value = TRUE)[[1L]]
  direct_interaction <- grep("^(Stress:Support|Support:Stress)$",
                             names(coefs), value = TRUE)[[1L]]
  b3 <- unname(coefs[[b_interaction]])
  c1 <- unname(coefs[["Stress"]])
  c3 <- unname(coefs[[direct_interaction]])
  low_ind <- a * (b1 + b3 * -1)
  mean_ind <- a * b1
  high_ind <- a * (b1 + b3 * 1)
  low_direct <- c1 + c3 * -1
  mean_direct <- c1
  high_direct <- c1 + c3 * 1
  condition_values <- c("-1.000 (- SD)", "0.000 (Mean)", "1.000 (+ SD) ")

  direct_slopes <- data.frame(
    Effect = c(low_direct, mean_direct, high_direct),
    "S.E." = c(0.08, 0.07, 0.09),
    LLCI = c(low_direct - 0.10, mean_direct - 0.08, high_direct - 0.11),
    ULCI = c(low_direct + 0.10, mean_direct + 0.08, high_direct + 0.11),
    t = c(3.1, 1.8, -0.4),
    pval = c(0.002, 0.074, 0.690),
    check.names = FALSE
  )
  direct_slopes <- cbind(
    setNames(data.frame(condition_values, check.names = FALSE),
             "\"Support\"    "),
    direct_slopes
  )
  b_slopes <- data.frame(
    Effect = c(b1 + b3 * -1, b1, b1 + b3 * 1),
    "S.E." = c(0.09, 0.07, 0.08),
    LLCI = c(b1 + b3 * -1 - 0.10, b1 - 0.08, b1 + b3 * 1 - 0.09),
    ULCI = c(b1 + b3 * -1 + 0.10, b1 + 0.08, b1 + b3 * 1 + 0.09),
    t = c(5.2, 4.8, 2.9),
    pval = c(0.000002, 0.00001, 0.004),
    check.names = FALSE
  )
  b_slopes <- cbind(
    setNames(data.frame(condition_values, check.names = FALSE),
             "\"Support\"    "),
    b_slopes
  )
  med <- data.frame(
    Effect = c(low_ind, mean_ind, high_ind),
    "S.E." = c(0.08, 0.06, 0.05),
    LLCI = c(low_ind - 0.10, mean_ind - 0.08, high_ind - 0.07),
    ULCI = c(low_ind + 0.10, mean_ind + 0.08, high_ind + 0.07),
    p = c(0, 0, 0.02),
    z = c(5.1, 4.7, 2.4),
    pval = c(0.000001, 0.000003, 0.016),
    check.names = FALSE
  )
  med <- cbind(
    setNames(data.frame(condition_values, check.names = FALSE),
             "\"Support\"    "),
    med
  )
  x <- list(
    process.id = 15,
    process.type = "Moderated Mediation",
    model.m = list(model.m.1 = model_m),
    model.y = model_y,
    results = list(
      list(
        conditional = data.frame(Support = c(-1, 0, 1)),
        simple.slopes = direct_slopes
      ),
      list(
        conditional = data.frame(Support = c(-1, 0, 1)),
        simple.slopes = b_slopes
      ),
      list(mediation = med)
    )
  )

  p <- mellio_payload(
    x,
    .call = "bruceR::PROCESS(data, y, x, meds, mods, mod.path = c('m-y', 'x-y'))"
  )
  expect_equal(p$type, "mediation_mediate")
  expect_equal(p$type_label, "Moderated mediation analysis (bruceR PROCESS Model 15)")
  expect_equal(p$fields$process_id, 15L)
  expect_equal(p$fields$topology, "moderated_mediation")
  expect_equal(p$fields$hayes_model, "15")
  expect_equal(p$fields$treatment, "Stress")
  expect_equal(p$fields$moderator, "Support")
  expect_equal(p$fields$moderation$moderated_path, "b")
  expect_equal(p$fields$moderation$path, "M_to_Y")
  expect_equal(p$fields$moderation$interaction_term, b_interaction)
  expect_equal(p$fields$moderation$direct_moderated_path, "X_to_Y")
  expect_equal(p$fields$moderation$direct_interaction_term,
               direct_interaction)
  expect_equal(as.character(p$fields$controls), "Tenure")
  expect_equal(p$fields$n, n)

  expect_equal(length(p$fields$rows), 7L)
  expect_equal(vapply(p$fields$rows[1:3], `[[`, character(1), "condition"),
               c("low", "mean", "high"))
  expect_equal(p$fields$rows[[4]]$effect,
               "Index of moderated mediation")
  expect_equal(p$fields$rows[[4]]$estimate, a * b3, tolerance = 1e-8)
  expect_equal(vapply(p$fields$rows[5:7], `[[`, character(1), "condition"),
               c("low", "mean", "high"))
  expect_equal(p$fields$rows[[5]]$effect,
               "Conditional direct effect (low Support)")
  expect_equal(p$fields$rows[[5]]$estimate, low_direct, tolerance = 1e-8)
  expect_equal(p$fields$rows[[5]]$statistic_type, "t")
  expect_equal(p$fields$columns[[4]]$label, "test")
  expect_match(p$fields$effect_note, "a \u00d7 b3", fixed = TRUE)
  expect_match(p$fields$effect_note, "conditional direct intervals")

  paths <- stats::setNames(p$fields$path_rows,
                           vapply(p$fields$path_rows, `[[`, character(1), "path"))
  expect_true(all(c("a", "b1", "b2", "b3", "c1", "c3",
                    "control_m_Support", "control_m_Tenure",
                    "control_y_Tenure") %in% names(paths)))
  expect_equal(paths$a$estimate, a)
  expect_equal(paths$b1$estimate, b1)
  expect_equal(paths$b3$estimate, b3)
  expect_equal(paths$c1$estimate, c1)
  expect_equal(paths$c3$estimate, c3)
  expect_equal(paths$b3$parameter,
               "b3: Burnout \u00d7 Support \u2192 TurnoverIntent")
  expect_equal(paths$c3$parameter,
               "c3: Stress \u00d7 Support \u2192 TurnoverIntent")
  expect_true(is.null(paths$b3[["std_estimate"]]))
  expect_true(is.null(paths$c3[["std_estimate"]]))
  expect_match(p$fields$path_note, "\u03b2 is not reported for interaction paths",
               fixed = TRUE)
})

test_that("unsupported bruceR PROCESS lists fail visibly", {
  x <- list(
    process.id = 59,
    process.type = "Moderated Mediation",
    results = list(raw = "not yet supported")
  )

  p <- mellio_payload(x, .call = "bruceR::PROCESS(..., model = 59)")
  expect_equal(p$type, "unsupported")
  expect_equal(p$card_kind, "unsupported")
  expect_match(p$fields$message, "bruceR::PROCESS")
  expect_match(p$fields$message, "Models 4, 6, 7, 8, 14, and 15")
  expect_equal(p$fields$process_id, 59)
})

test_that("glmer payload — model-summary with family in type_label", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("broom.mixed")
  suppressMessages(library(lme4))
  m <- glmer(cbind(incidence, size - incidence) ~ period + (1 | herd),
             data = cbpp, family = binomial)
  p <- mellio_payload(m)
  expect_equal(p$type, "glmer_model_summary")
  expect_match(p$type_label, "Generalized mixed model \\(binomial; logit\\)")
  expect_null(p$fields[["statistic"]])
  expect_true(is.na(p$fields$p_value))
  expect_equal(p$fields$n, 56L)
  expect_equal(p$fields$conf_level, 0.95)
  expect_equal(p$fields$coefficient_ci_method, "wald")
  expect_equal(p$fields$coefficient_ci_scale, "link")
  expect_equal(p$fields$model_family, "binomial")
  expect_equal(p$fields$model_link, "logit")
  expect_equal(p$fields$outcome, "cbind(incidence, size - incidence)")
  expect_equal(p$fields$focal_terms, "period")
  coefs <- p$fields$coefficients
  expect_true(any(vapply(coefs, function(row) identical(row$term, "period2"), logical(1))))
  period2 <- coefs[[which(vapply(coefs, function(row) identical(row$term, "period2"), logical(1)))[1]]]
  expect_equal(period2$statistic_label, "z")
  expect_equal(period2$ci_method, "wald")
  expect_true(is.numeric(period2$p_value))
  expect_true(period2$p_value < 0.01)
  expect_equal(p$fields$groups[[1]]$label, "herd")
  expect_equal(p$fields$groups[[1]]$n, 15L)
  expect_true(length(p$fields$random_effects) > 0L)

  tests <- p$fields$model_term_tests
  expect_true(is.list(tests))
  expect_length(tests, 1L)
  drop_tbl <- drop1(m, test = "Chisq")
  expect_equal(tests[[1]]$term, "period")
  expect_equal(tests[[1]]$label, "period")
  expect_equal(tests[[1]]$method, "drop1_chisq")
  expect_equal(tests[[1]]$test_type, "omnibus")
  expect_equal(tests[[1]]$term_type, "main")
  expect_equal(tests[[1]]$statistic$name, "chi2")
  expect_equal(tests[[1]]$statistic$df, 3)
  expect_equal(tests[[1]]$statistic$value, unname(drop_tbl["period", "LRT"]), tolerance = 1e-6)
  expect_equal(tests[[1]]$p_value, unname(drop_tbl["period", "Pr(Chi)"]), tolerance = 1e-6)

  fig <- p$figure_data$coefficient_plot
  expect_equal(fig$coefficient_scale, "odds_ratio")
  expect_equal(fig$estimate_label, "OR")
  expect_equal(fig$model_family, "binomial")
  expect_equal(fig$model_link, "logit")
  expect_equal(fig$groups[[1]]$label, "herd")
  period2_fig <- fig$coefficients[[which(vapply(fig$coefficients, function(row) {
    identical(row$term, "period2")
  }, logical(1)))[1]]]
  expect_equal(period2_fig$estimate_name, "OR")
  expect_equal(period2_fig$log_estimate, period2$estimate, tolerance = 1e-8)
  expect_equal(period2_fig$estimate, exp(period2$estimate), tolerance = 1e-8)
})

test_that("glmer payload exposes continuous main-effect predicted probability figures", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("broom.mixed")
  skip_if_not_installed("emmeans")
  suppressMessages(library(lme4))

  set.seed(2032)
  n_id <- 80
  obs_per_id <- 4
  d <- data.frame(
    id = factor(rep(seq_len(n_id), each = obs_per_id)),
    group = factor(rep(rep(c("Control", "Treatment"), each = n_id / 2), each = obs_per_id))
  )
  d$stress <- rnorm(nrow(d), 0, 1)
  subject_effect <- rnorm(n_id, 0, 0.8)
  eta <- -0.5 +
    rep(subject_effect, each = obs_per_id) +
    1.1 * d$stress +
    ifelse(d$group == "Treatment", 0.8, 0)
  d$success <- stats::rbinom(nrow(d), size = 1, prob = stats::plogis(eta))

  m <- glmer(success ~ stress + group + (1 | id), data = d, family = binomial)
  p <- mellio_payload(m)

  expect_equal(p$type, "glmer_model_summary")
  expect_equal(p$fields$model_kind, "controlled_glmm")
  types <- vapply(p$metadata$available_figures, function(row) row$type, character(1))
  expect_equal(types[[1]], "adjusted_means")
  expect_equal(p$metadata$available_figures[[1]]$label, "Predicted probabilities")
  expect_true(isTRUE(p$metadata$available_figures[[1]]$default))
  expect_true("interaction_plot" %in% types)
  expect_true("coefficient_plot" %in% types)

  means <- p$figure_data$adjusted_means
  expect_equal(means$source, "glmer_emmeans")
  expect_equal(means$mean_kind, "predicted_probability")
  expect_equal(means$factor$variable, "group")
  expect_equal(means$outcome, "success")
  expect_equal(means$subject$variable, "id")
  expect_equal(means$subject$n, n_id)
  expect_equal(means$y_label, "Predicted probability")
  expect_equal(means$model_family, "binomial")
  expect_equal(means$model_link, "logit")
  expect_true(isTRUE(means$bounded_response))
  expect_length(means$groups, 2L)
  expect_true(any(vapply(means$covariates, function(row) {
    identical(row$variable, "stress")
  }, logical(1))))
  expect_true(all(vapply(means$groups, function(row) {
    is.numeric(row$mean) &&
      row$mean >= 0 &&
      row$mean <= 1 &&
      is.numeric(row$ci_lower) &&
      is.numeric(row$ci_upper) &&
      row$ci_lower < row$ci_upper
  }, logical(1))))

  fig <- p$figure_data$interaction_plot
  expect_equal(fig$source, "glmer_emmeans")
  expect_equal(fig$mean_kind, "predicted_probability")
  expect_equal(fig$interaction_kind, "continuous_main_effect")
  expect_equal(fig$interaction_term, "stress")
  expect_equal(fig$x$variable, "stress")
  expect_equal(fig$x$type, "numeric")
  expect_equal(fig$moderator$variable, "estimate")
  expect_equal(fig$outcome, "success")
  expect_equal(fig$subject$variable, "id")
  expect_equal(fig$subject$n, n_id)
  expect_equal(fig$y_label, "Predicted probability")
  expect_equal(fig$model_family, "binomial")
  expect_equal(fig$model_link, "logit")
  expect_true(isTRUE(fig$bounded_response))
  expect_length(fig$grid, 80L)
  expect_true(any(vapply(fig$marginalized_terms, function(row) {
    identical(row$variable, "group")
  }, logical(1))))
  expect_true(all(vapply(fig$grid, function(row) {
    is.numeric(row$x) &&
      is.numeric(row$estimate) &&
      row$estimate >= 0 &&
      row$estimate <= 1 &&
      is.numeric(row$se) &&
      is.numeric(row$ci_lower) &&
      is.numeric(row$ci_upper) &&
      row$ci_lower < row$ci_upper
  }, logical(1))))

  pkg_names <- vapply(p$packages, `[[`, character(1), "name")
  expect_true("lme4" %in% pkg_names)
  expect_true("emmeans" %in% pkg_names)
})

test_that("glmer payload exposes continuous interaction predicted probability figures", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("broom.mixed")
  skip_if_not_installed("emmeans")
  suppressMessages(library(lme4))

  set.seed(2033)
  n_id <- 80
  obs_per_id <- 5
  d <- data.frame(
    id = factor(rep(seq_len(n_id), each = obs_per_id)),
    group = factor(rep(rep(c("Control", "Treatment"), each = n_id / 2), each = obs_per_id))
  )
  d$stress <- rnorm(nrow(d), 0, 1)
  subject_effect <- rnorm(n_id, 0, 0.7)
  eta <- -0.5 +
    rep(subject_effect, each = obs_per_id) +
    0.9 * d$stress +
    ifelse(d$group == "Treatment", 0.6, 0) +
    ifelse(d$group == "Treatment", 0.8 * d$stress, 0)
  d$success <- stats::rbinom(nrow(d), size = 1, prob = stats::plogis(eta))

  m <- glmer(
    success ~ stress * group + (1 | id),
    data = d,
    family = binomial,
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  )
  p <- mellio_payload(m)

  expect_equal(p$type, "glmer_model_summary")
  expect_null(p$figure_data$adjusted_means)
  types <- vapply(p$metadata$available_figures, function(row) row$type, character(1))
  expect_equal(types[[1]], "interaction_plot")
  expect_equal(p$metadata$available_figures[[1]]$label, "Interaction plot")
  expect_true(isTRUE(p$metadata$available_figures[[1]]$default))
  expect_true("coefficient_plot" %in% types)

  fig <- p$figure_data$interaction_plot
  expect_equal(fig$source, "glmer_emmeans")
  expect_equal(fig$mean_kind, "predicted_probability")
  expect_equal(fig$interaction_kind, "continuous_by_categorical")
  expect_equal(fig$interaction_term, "stress:group")
  expect_equal(fig$variables, c("stress", "group"))
  expect_equal(fig$x$variable, "stress")
  expect_equal(fig$x$type, "numeric")
  expect_equal(fig$moderator$variable, "group")
  expect_equal(fig$moderator$type, "categorical")
  expect_equal(fig$outcome, "success")
  expect_equal(fig$subject$variable, "id")
  expect_equal(fig$subject$n, n_id)
  expect_equal(fig$y_label, "Predicted probability")
  expect_equal(fig$model_family, "binomial")
  expect_equal(fig$model_link, "logit")
  expect_true(isTRUE(fig$bounded_response))
  expect_length(fig$grid, 160L)
  expect_true(all(vapply(fig$grid, function(row) {
    is.numeric(row$x) &&
      is.numeric(row$estimate) &&
      row$estimate >= 0 &&
      row$estimate <= 1 &&
      is.numeric(row$se) &&
      is.numeric(row$ci_lower) &&
      is.numeric(row$ci_upper) &&
      row$ci_lower < row$ci_upper
  }, logical(1))))

  pkg_names <- vapply(p$packages, `[[`, character(1), "name")
  expect_true("emmeans" %in% pkg_names)
})

test_that("glmer payload exposes categorical interaction predicted probability figures", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("broom.mixed")
  skip_if_not_installed("emmeans")
  suppressMessages(library(lme4))

  set.seed(2034)
  n_id <- 80
  obs_per_id <- 4
  d <- data.frame(
    id = factor(rep(seq_len(n_id), each = obs_per_id)),
    group = factor(rep(rep(c("Control", "Treatment"), each = n_id / 2), each = obs_per_id)),
    condition = factor(rep(rep(c("Neutral", "Cue"), each = obs_per_id / 2), times = n_id))
  )
  d$stress <- rnorm(nrow(d), 0, 1)
  subject_effect <- rnorm(n_id, 0, 0.7)
  eta <- -0.7 +
    rep(subject_effect, each = obs_per_id) +
    0.5 * d$stress +
    ifelse(d$group == "Treatment", 0.5, 0) +
    ifelse(d$condition == "Cue", 0.6, 0) +
    ifelse(d$group == "Treatment" & d$condition == "Cue", 0.8, 0)
  d$success <- stats::rbinom(nrow(d), size = 1, prob = stats::plogis(eta))

  m <- glmer(
    success ~ group * condition + stress + (1 | id),
    data = d,
    family = binomial,
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  )
  p <- mellio_payload(m)

  expect_equal(p$type, "glmer_model_summary")
  expect_null(p$figure_data$adjusted_means)
  types <- vapply(p$metadata$available_figures, function(row) row$type, character(1))
  expect_equal(types[[1]], "interaction_plot")
  expect_equal(p$metadata$available_figures[[1]]$label, "Interaction plot")
  expect_true(isTRUE(p$metadata$available_figures[[1]]$default))
  expect_true("coefficient_plot" %in% types)

  fig <- p$figure_data$interaction_plot
  expect_equal(fig$source, "glmer_emmeans")
  expect_equal(fig$mean_kind, "predicted_probability")
  expect_equal(fig$interaction_kind, "categorical_by_categorical")
  expect_equal(fig$interaction_term, "group:condition")
  expect_equal(fig$variables, c("condition", "group"))
  expect_equal(fig$x$variable, "condition")
  expect_equal(fig$x$type, "categorical")
  expect_equal(fig$moderator$variable, "group")
  expect_equal(fig$moderator$type, "categorical")
  expect_equal(fig$outcome, "success")
  expect_equal(fig$subject$variable, "id")
  expect_equal(fig$subject$n, n_id)
  expect_equal(fig$y_label, "Predicted probability")
  expect_equal(fig$model_family, "binomial")
  expect_equal(fig$model_link, "logit")
  expect_true(isTRUE(fig$bounded_response))
  expect_length(fig$grid, 4L)
  expect_equal(sort(unique(vapply(fig$grid, `[[`, character(1), "x_value"))),
               c("Cue", "Neutral"))
  expect_equal(sort(unique(vapply(fig$grid, `[[`, character(1), "moderator_value"))),
               c("Control", "Treatment"))
  expect_true(any(vapply(fig$held_constant, function(row) {
    identical(row$variable, "stress")
  }, logical(1))))
  expect_true(all(vapply(fig$grid, function(row) {
    is.numeric(row$x) &&
      is.character(row$x_value) &&
      is.numeric(row$estimate) &&
      row$estimate >= 0 &&
      row$estimate <= 1 &&
      is.numeric(row$se) &&
      is.numeric(row$ci_lower) &&
      is.numeric(row$ci_upper) &&
      row$ci_lower < row$ci_upper
  }, logical(1))))

  pkg_names <- vapply(p$packages, `[[`, character(1), "name")
  expect_true("emmeans" %in% pkg_names)
})

test_that("glm payload routes through broom generic instead of lm", {
  m <- glm(am ~ wt + hp, data = mtcars, family = binomial)
  p <- mellio_payload(m)
  expect_equal(p$type, "logistic_regression")
  expect_match(p$type_label, "Logistic regression")
  expect_equal(p$fields$model_family, "binomial")
  expect_equal(p$fields$model_link, "logit")
  expect_equal(p$fields$outcome, "am")
  expect_equal(p$fields$focal_terms, c("wt", "hp"))
  expect_true(is.numeric(p$fields$aic))
  expect_true(is.numeric(p$fields$bic))
  expect_true(is.numeric(p$fields$logLik))
  expect_equal(p$fields$conf_level, 0.95)
  expect_equal(p$fields$coefficient_ci_method, "profile_likelihood")
  expect_equal(p$fields$coefficient_ci_scale, "link")
  expect_true(is.numeric(p$fields$null_deviance))
  expect_true(is.numeric(p$fields$residual_deviance))
  coefs <- p$fields$coefficients
  expect_true(any(vapply(coefs, function(row) identical(row$term, "wt"), logical(1))))
  wt <- coefs[[which(vapply(coefs, function(row) identical(row$term, "wt"), logical(1)))[1]]]
  expect_equal(wt$statistic_label, "z")
  expect_true(is.numeric(wt$p_value))
  expect_equal(wt$ci_method, "profile_likelihood")
  expect_null(p$fields[["coefficient_scale"]])

  fig <- p$figure_data$coefficient_plot
  expect_equal(fig$coefficient_scale, "odds_ratio")
  expect_equal(fig$estimate_label, "OR")
  expect_equal(fig$model_family, "binomial")
  expect_equal(fig$model_link, "logit")
  wt_fig <- fig$coefficients[[which(vapply(fig$coefficients, function(row) {
    identical(row$term, "wt")
  }, logical(1)))[1]]]
  expect_equal(wt_fig$estimate_name, "OR")
  expect_equal(wt_fig$log_estimate, wt$estimate, tolerance = 1e-8)
  expect_equal(wt_fig$estimate, exp(wt$estimate), tolerance = 1e-8)
  if (!is.null(wt$ci_lower)) {
    expect_equal(wt_fig$ci_lower, exp(wt$ci_lower), tolerance = 1e-8)
    expect_equal(wt_fig$ci_upper, exp(wt$ci_upper), tolerance = 1e-8)
  }
})

test_that("glm payload can exponentiate logistic coefficients", {
  m <- glm(am ~ wt + hp, data = mtcars, family = binomial)
  p <- mellio_payload(m, exponentiate = TRUE, conf.int = FALSE)
  wt <- p$fields$coefficients[[which(vapply(p$fields$coefficients, function(row) {
    identical(row$term, "wt")
  }, logical(1)))[1]]]
  expect_equal(wt$estimate_name, "OR")
  expect_equal(p$fields$coefficient_scale, "odds_ratio")
  expect_true(wt$estimate > 0)

  wt_fig <- p$figure_data$coefficient_plot$coefficients[[which(vapply(p$figure_data$coefficient_plot$coefficients, function(row) {
    identical(row$term, "wt")
  }, logical(1)))[1]]]
  expect_equal(p$figure_data$coefficient_plot$coefficient_scale, "odds_ratio")
  expect_equal(wt_fig$estimate, wt$estimate, tolerance = 1e-8)
  expect_null(wt_fig$log_estimate)
})

test_that("glm payload surfaces boundary and separation diagnostics", {
  d <- data.frame(
    y = c(rep(0, 10), rep(1, 10)),
    x = c(rep(0, 10), rep(1, 10))
  )
  m <- suppressWarnings(glm(y ~ x, data = d, family = binomial()))
  p <- mellio_payload(m)
  warnings <- p$fields$model_warnings
  expect_true(is.list(warnings))
  types <- vapply(warnings, function(row) if (is.null(row$type)) "" else row$type, character(1))
  expect_true(any(types %in% c("boundary_fit", "separation_or_boundary")))
  messages <- vapply(warnings, function(row) if (is.null(row$message)) "" else row$message, character(1))
  expect_true(any(grepl("boundary|probabilities", messages, ignore.case = TRUE)))
  expect_null(p$figure_data$coefficient_plot)
  figure_types <- vapply(p$metadata$available_figures, function(row) row$type, character(1))
  expect_false("coefficient_plot" %in% figure_types)
})

test_that("gam payload uses a generalized-additive model type", {
  skip_if_not_installed("mgcv")

  fit <- mgcv::gam(mpg ~ s(wt) + hp, data = mtcars)
  p <- mellio_payload(fit)

  expect_equal(p$type, "generalized_additive_model")
  expect_match(p$type_label, "Generalized additive model")
  expect_equal(p$fields$outcome, "mpg")
  expect_true(any(vapply(p$fields$coefficients, function(row) {
    identical(row$term, "hp")
  }, logical(1))))
})

test_that("coxph payload includes hazard-ratio table-note metadata", {
  skip_if_not_installed("survival")

  fit <- survival::coxph(survival::Surv(time, status) ~ age + sex,
                         data = survival::lung)
  p <- mellio_payload(fit)

  expect_equal(p$type, "cox_proportional_hazards")
  expect_equal(p$fields$coefficient_scale, "hazard_ratio")
  expect_equal(p$fields$coefficient_ci_method, "wald")
  expect_equal(p$fields$coefficient_ci_scale, "hazard_ratio")
  expect_equal(p$fields$coefficient_p_value_method, "wald_z")
  expect_equal(p$fields$conf_level, 0.95)
  expect_true(is.numeric(p$fields$aic))
  expect_true(is.numeric(p$fields$bic))
  expect_true(is.numeric(p$fields$logLik))
  expect_true(p$fields$events > 0)

  sex <- p$fields$coefficients[[which(vapply(p$fields$coefficients, function(row) {
    identical(row$term, "sex")
  }, logical(1)))[1]]]
  expect_equal(sex$estimate_name, "HR")
  expect_equal(sex$ci_method, "wald")
  expect_true(is.numeric(sex$ci_lower))
  expect_true(is.numeric(sex$ci_upper))
})

test_that("geeglm payload includes cluster and robust-Wald table-note metadata", {
  skip_if_not_installed("geepack")

  data(dietox, package = "geepack")
  fit <- geepack::geeglm(Weight ~ Time, id = Pig, data = dietox,
                         corstr = "exchangeable")
  p <- mellio_payload(fit)

  expect_equal(p$type, "gee_model_summary")
  expect_equal(p$fields$model_family, "gaussian")
  expect_equal(p$fields$model_link, "identity")
  expect_equal(p$fields$id_variable, "Pig")
  expect_equal(p$fields$n_clusters, length(unique(dietox$Pig)))
  expect_equal(p$fields$correlation_structure, "exchangeable")
  expect_equal(p$fields$coefficient_ci_method, "wald_robust")
  expect_equal(p$fields$coefficient_ci_scale, "response")
  expect_equal(p$fields$coefficient_p_value_method, "robust_wald")
  expect_equal(p$fields$std_error_method, "robust_sandwich")
  expect_equal(p$fields$statistic_label, "Wald \u03c7\u00b2")
  expect_equal(p$fields$working_correlation_parameter_name, "alpha")
  expect_true(is.numeric(p$fields$working_correlation_parameter))
  expect_equal(p$fields$conf_level, 0.95)

  time <- p$fields$coefficients[[which(vapply(p$fields$coefficients, function(row) {
    identical(row$term, "Time")
  }, logical(1)))[1]]]
  expect_equal(time$estimate_name, "B")
  expect_equal(time$ci_method, "wald_robust")
  expect_equal(time$statistic_label, "Wald \u03c7\u00b2")
  expect_equal(time$df, 1)
  expect_true(is.numeric(time$ci_lower))
  expect_true(is.numeric(time$ci_upper))
  expect_true(time$ci_lower < time$ci_upper)
})

test_that("polr payload includes ordinal table-note metadata", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("broom")

  set.seed(20260525)
  n <- 120
  x <- stats::rnorm(n)
  latent <- 0.8 * x + stats::rnorm(n)
  y <- ordered(
    cut(latent,
        breaks = stats::quantile(latent, probs = c(0, 0.33, 0.66, 1)),
        include.lowest = TRUE,
        labels = c("low", "medium", "high")),
    levels = c("low", "medium", "high")
  )
  fit <- MASS::polr(y ~ x, data = data.frame(y = y, x = x), Hess = TRUE)
  p <- suppressWarnings(suppressMessages(mellio_payload(fit)))

  expect_equal(p$type, "ordinal_regression")
  expect_equal(p$fields$model_family, "ordinal")
  expect_equal(p$fields$model_link, "logit")
  expect_equal(p$fields$coefficient_scale, "proportional_odds")
  expect_equal(p$fields$coefficient_ci_scale, "proportional_odds")
  expect_equal(p$fields$coefficient_p_value_method, "wald_z")
  expect_equal(p$fields$conf_level, 0.95)
  expect_equal(p$fields$coefficient_ci_method, "profile_likelihood")
  expect_equal(unlist(p$fields$outcome_levels, use.names = FALSE),
               c("low", "medium", "high"))
  expect_equal(p$fields$threshold_count, 2L)

  x_row <- p$fields$coefficients[[which(vapply(p$fields$coefficients, function(row) {
    identical(row$term, "x")
  }, logical(1)))[1]]]
  expect_equal(x_row$estimate_name, "OR")
  expect_equal(x_row$ci_method, "profile_likelihood")
  expect_true(is.numeric(x_row$ci_lower))
  expect_true(is.numeric(x_row$ci_upper))
})

test_that("clm payload includes ordinal profile-likelihood CI metadata", {
  skip_if_not_installed("ordinal")

  set.seed(20260525)
  n <- 100
  x <- stats::rnorm(n)
  latent <- 0.7 * x + stats::rnorm(n)
  y <- ordered(
    cut(latent,
        breaks = stats::quantile(latent, probs = c(0, 0.33, 0.66, 1)),
        include.lowest = TRUE,
        labels = c("low", "medium", "high")),
    levels = c("low", "medium", "high")
  )
  fit <- ordinal::clm(y ~ x, data = data.frame(y = y, x = x))
  p <- suppressWarnings(suppressMessages(mellio_payload(fit)))

  expect_equal(p$type, "ordinal_regression")
  expect_equal(p$fields$coefficient_scale, "proportional_odds")
  expect_equal(p$fields$coefficient_ci_method, "profile_likelihood")
  expect_equal(p$fields$coefficient_ci_scale, "proportional_odds")
  expect_equal(p$fields$coefficient_p_value_method, "wald_z")
  expect_equal(p$fields$conf_level, 0.95)
  expect_true(is.numeric(p$fields$bic))
  expect_equal(p$fields$threshold_count, 2L)

  x_row <- p$fields$coefficients[[which(vapply(p$fields$coefficients, function(row) {
    identical(row$term, "x")
  }, logical(1)))[1]]]
  expect_equal(x_row$ci_method, "profile_likelihood")
  expect_true(is.numeric(x_row$ci_lower))
  expect_true(is.numeric(x_row$ci_upper))
})

test_that("multinom payload includes table-note metadata and Wald CIs", {
  skip_if_not_installed("nnet")

  fit <- nnet::multinom(Species ~ Sepal.Length + Sepal.Width,
                        data = iris,
                        trace = FALSE)
  p <- mellio_payload(fit)

  expect_equal(p$type, "multinomial_logistic_regression")
  expect_equal(p$card_kind, "table")
  expect_equal(p$fields$table_type, "multinomial_coefficients")
  expect_equal(p$fields$model_family, "multinomial logistic")
  expect_equal(p$fields$model_link, "logit")
  expect_equal(p$fields$coefficient_scale, "multinomial_logit")
  expect_equal(p$fields$coefficient_ci_method, "wald")
  expect_equal(p$fields$coefficient_ci_scale, "link")
  expect_equal(p$fields$coefficient_p_value_method, "wald_z")
  expect_equal(p$fields$conf_level, 0.95)
  expect_equal(unlist(p$fields$outcome_levels, use.names = FALSE),
               levels(iris$Species))
  expect_equal(p$fields$reference_level, "setosa")
  expect_true(is.numeric(p$fields$bic))
  expect_true(is.numeric(p$fields$logLik))

  expect_true(any(vapply(p$fields$columns, function(col) {
    identical(col$key, "ci")
  }, logical(1))))
  row <- p$fields$rows[[which(vapply(p$fields$rows, function(row) {
    identical(row$comparison, "versicolor") && identical(row$term, "Sepal.Length")
  }, logical(1)))[1]]]
  expect_equal(row$estimate_name, "B")
  expect_equal(row$ci_method, "wald")
  expect_true(is.numeric(row$ci_lower))
  expect_true(is.numeric(row$ci_upper))
})

test_that("lme payload uses the generic mixed-model type", {
  skip_if_not_installed("nlme")
  skip_if_not_installed("broom.mixed")

  fit <- nlme::lme(mpg ~ wt, random = ~1 | cyl, data = mtcars)
  p <- mellio_payload(fit)

  expect_equal(p$type, "mixed_model_summary")
  expect_match(p$type_label, "Linear mixed model")
  expect_equal(p$fields$outcome, "mpg")
  expect_equal(p$fields$focal_terms, "wt")
  expect_true(any(vapply(p$fields$coefficients, function(row) {
    identical(row$term, "wt")
  }, logical(1))))
})

test_that("lmer payload — model summary includes fixed and random effects", {
  skip_if_not_installed("lme4")
  suppressMessages(library(lme4))
  m <- lmer(Reaction ~ Days + (Days | Subject), data = sleepstudy)

  p <- mellio_payload(m)
  expect_equal(p$type, "lmer_model_summary")
  expect_match(p$type_label, "Linear mixed model \\((REML|ML)\\)")
  expect_equal(p$fields$statistic$name, "AIC")
  expect_true(is.numeric(p$fields$statistic$value))
  expect_true(is.na(p$fields$p_value))
  expect_equal(p$fields$estimate$name, "logLik")
  expect_true(is.numeric(p$fields$estimate$value))
  expect_equal(p$fields$bic, ms_safe_numeric(BIC(m)), tolerance = 1e-6)
  expect_equal(p$fields$logLik, ms_safe_numeric(as.numeric(logLik(m))), tolerance = 1e-6)
  expect_equal(p$fields$n, 180L)
  expect_equal(p$fields$model_fit, "REML")
  expect_equal(p$fields$conf_level, 0.95)
  has_lmer_test <- requireNamespace("lmerTest", quietly = TRUE)
  if (has_lmer_test) {
    expect_equal(p$fields$coefficient_ci_method, "satterthwaite_t")
    expect_equal(p$fields$coefficient_p_value_method, "satterthwaite_t")
    expect_equal(p$fields$fixed_effect_df_method, "satterthwaite")
    expect_equal(p$fields$model_term_test_method, "satterthwaite_f")
    expect_equal(p$fields$model_term_ss_type, "type_iii")
  } else {
    expect_equal(p$fields$coefficient_ci_method, "wald")
    expect_null(p$fields$coefficient_p_value_method)
    expect_null(p$fields$fixed_effect_df_method)
  }
  expect_equal(p$fields$outcome, "Reaction")
  expect_equal(p$fields$focal_terms, "Days")
  if (requireNamespace("performance", quietly = TRUE)) {
    r2_fun <- getExportedValue("performance", "r2")
    r2_df <- as.data.frame(r2_fun(m))
    expect_equal(p$fields$r2_marginal, r2_df$R2_marginal[[1]], tolerance = 1e-6)
    expect_equal(p$fields$r2_conditional, r2_df$R2_conditional[[1]], tolerance = 1e-6)
    expect_equal(p$fields$r2_method, "nakagawa")
  } else {
    expect_null(p$fields$r2_marginal)
    expect_null(p$fields$r2_conditional)
  }

  tests <- p$fields$model_term_tests
  expect_true(is.list(tests))
  expect_length(tests, 1L)
  expect_equal(tests[[1]]$term, "Days")
  expect_equal(tests[[1]]$label, "Days")
  expect_equal(tests[[1]]$test_type, "omnibus")
  expect_equal(tests[[1]]$term_type, "main")
  expect_equal(tests[[1]]$test_scope, "fixed_effect")
  if (has_lmer_test) {
    lmer_test_fit <- lmerTest::as_lmerModLmerTest(m)
    satt_tbl <- as.data.frame(stats::anova(lmer_test_fit, type = 3, ddf = "Satterthwaite"))
    expect_equal(tests[[1]]$method, "anova_f_satterthwaite")
    expect_equal(tests[[1]]$model_fit, "REML")
    expect_equal(tests[[1]]$ddf_method, "satterthwaite")
    expect_equal(tests[[1]]$ss_type, "type_iii")
    expect_equal(tests[[1]]$statistic$name, "F")
    expect_equal(unclass(tests[[1]]$statistic$df),
                 c(unname(satt_tbl["Days", "NumDF"]), unname(satt_tbl["Days", "DenDF"])),
                 tolerance = 1e-6)
    expect_equal(tests[[1]]$statistic$value, unname(satt_tbl["Days", "F value"]), tolerance = 1e-6)
    expect_equal(tests[[1]]$p_value, unname(satt_tbl["Days", "Pr(>F)"]), tolerance = 1e-6)
  } else {
    drop_tbl <- drop1(refitML(m), test = "Chisq")
    expect_equal(tests[[1]]$method, "drop1_chisq_ml_refit")
    expect_equal(tests[[1]]$model_fit, "ML")
    expect_equal(tests[[1]]$refit, "ML")
    expect_equal(tests[[1]]$refit_from, "REML")
    expect_equal(tests[[1]]$statistic$name, "chi2")
    expect_equal(tests[[1]]$statistic$df, 1)
    expect_equal(tests[[1]]$statistic$value, unname(drop_tbl["Days", "LRT"]), tolerance = 1e-6)
    expect_equal(tests[[1]]$p_value, unname(drop_tbl["Days", "Pr(Chi)"]), tolerance = 1e-6)
  }

  coefs <- p$fields$coefficients
  expect_length(coefs, 2L)
  days <- coefs[[which(vapply(coefs, `[[`, character(1), "term") == "Days")]]
  expect_equal(days$estimate, 10.467, tolerance = 0.001)
  expect_equal(days$std_error, 1.546, tolerance = 0.001)
  expect_equal(days$statistic, 6.771, tolerance = 0.001)
  expect_equal(days$statistic_label, "t")
  expect_equal(days$estimate_name, "B")
  if (has_lmer_test) {
    expect_equal(days$df, unname(satt_tbl["Days", "DenDF"]), tolerance = 1e-6)
    expect_equal(days$p_value, unname(satt_tbl["Days", "Pr(>F)"]), tolerance = 1e-6)
    expect_equal(days$p_value_method, "satterthwaite_t")
    expect_equal(days$ci_method, "satterthwaite_t")
    expect_equal(days$ci_lower, days$estimate - qt(0.975, days$df) * days$std_error, tolerance = 1e-6)
    expect_equal(days$ci_upper, days$estimate + qt(0.975, days$df) * days$std_error, tolerance = 1e-6)
  } else {
    expect_equal(days$ci_method, "wald")
    expect_equal(days$ci_lower, days$estimate - qnorm(0.975) * days$std_error, tolerance = 1e-6)
    expect_equal(days$ci_upper, days$estimate + qnorm(0.975) * days$std_error, tolerance = 1e-6)
  }

  fig <- p$figure_data$coefficient_plot
  expect_equal(fig$coefficient_scale, "raw")
  expect_equal(fig$estimate_label, "B")
  expect_equal(fig$model_fit, "REML")
  expect_equal(fig$groups[[1]]$label, "Subject")
  expect_true(any(vapply(fig$coefficients, function(row) {
    identical(row$term, "Days") && identical(row$estimate, days$estimate)
  }, logical(1))))

  expect_equal(p$fields$groups[[1]]$label, "Subject")
  expect_equal(p$fields$groups[[1]]$n, 18L)
  expect_equal(p$fields$random_terms, "(Days | Subject)")

  random <- p$fields$random_effects
  expect_true(any(vapply(random, function(row) {
    identical(row$group, "Subject") && identical(row$name, "Days")
  }, logical(1))))
  days_random <- random[[which(vapply(random, function(row) {
    identical(row$group, "Subject") && identical(row$name, "Days")
  }, logical(1)))[1]]]
  expect_equal(days_random$variance, 35.07, tolerance = 0.01)
  expect_equal(days_random$std_dev, 5.922, tolerance = 0.001)
  expect_equal(days_random$corr, 0.066, tolerance = 0.01)

  # lme4 should appear in packages
  pkg_names <- vapply(p$packages, `[[`, character(1), "name")
  expect_true("lme4" %in% pkg_names)
  if (has_lmer_test) {
    expect_true("lmerTest" %in% pkg_names)
  }
  if (requireNamespace("performance", quietly = TRUE)) {
    expect_true("performance" %in% pkg_names)
  }

  # Data hash should be present (lmer carries @frame)
  skip_if_not_installed("digest")
  expect_match(p$provenance$data$hash, "^[0-9a-f]{40}$")
  expect_equal(p$provenance$data$n, 180L)
})

test_that("lmer payload surfaces singular-fit diagnostics", {
  skip_if_not_installed("lme4")
  suppressMessages(library(lme4))
  m <- suppressMessages(suppressWarnings(
    lmer(mpg ~ wt + (wt | cyl), data = mtcars)
  ))
  p <- mellio_payload(m)
  warnings <- p$fields$model_warnings
  expect_true(is.list(warnings))
  types <- vapply(warnings, function(row) if (is.null(row$type)) "" else row$type, character(1))
  expect_true("singular_fit" %in% types)
  singular <- warnings[[which(types == "singular_fit")[[1]]]]
  expect_equal(singular$severity, "warning")
  expect_match(singular$message, "singular-fit")
})

test_that("lmer payload exposes emmeans interaction profile figures", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("emmeans")
  suppressMessages(library(lme4))

  set.seed(2028)
  n_id <- 40
  d <- data.frame(
    id = factor(rep(seq_len(n_id), each = 3)),
    time = factor(
      rep(c("Pre", "Post", "Follow-up"), times = n_id),
      levels = c("Pre", "Post", "Follow-up")
    )
  )
  d$group <- factor(rep(rep(c("Control", "Treatment"), each = n_id / 2), each = 3))
  subject_effect <- rnorm(n_id, 0, 5)
  d$score <- 50 + rep(subject_effect, each = 3) +
    ifelse(d$group == "Treatment", 3, 0) +
    ifelse(d$time == "Post", 4, ifelse(d$time == "Follow-up", 6, 0)) +
    ifelse(d$group == "Treatment" & d$time == "Follow-up", 4, 0) +
    rnorm(nrow(d), 0, 4)

  m <- lmer(score ~ time * group + (1 | id), data = d)
  p <- mellio_payload(m)

  expect_equal(p$type, "lmer_model_summary")
  types <- vapply(p$metadata$available_figures, function(row) row$type, character(1))
  expect_equal(types[[1]], "interaction_plot")
  expect_true(isTRUE(p$metadata$available_figures[[1]]$default))
  expect_true("coefficient_plot" %in% types)

  fig <- p$figure_data$interaction_plot
  expect_equal(fig$source, "lmer_emmeans")
  expect_equal(fig$mean_kind, "estimated_marginal")
  expect_equal(fig$interaction_kind, "categorical_by_categorical")
  expect_equal(fig$interaction_term, "time:group")
  expect_equal(fig$x$variable, "time")
  expect_equal(fig$moderator$variable, "group")
  expect_equal(fig$outcome, "score")
  expect_equal(fig$subject$variable, "id")
  expect_equal(fig$subject$n, n_id)
  expect_equal(fig$ci_method, "emmeans")
  expect_length(fig$grid, 6L)
  expect_true(all(vapply(fig$grid, function(row) {
    is.numeric(row$estimate) &&
      is.numeric(row$se) &&
      is.numeric(row$ci_lower) &&
      is.numeric(row$ci_upper) &&
      row$ci_lower < row$ci_upper
  }, logical(1))))

  pkg_names <- vapply(p$packages, `[[`, character(1), "name")
  expect_true("lme4" %in% pkg_names)
  expect_true("emmeans" %in% pkg_names)
})

test_that("lmer payload exposes main-effect estimated marginal means figures", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("emmeans")
  suppressMessages(library(lme4))

  set.seed(2029)
  n_id <- 40
  d <- data.frame(
    id = factor(rep(seq_len(n_id), each = 3)),
    time = factor(
      rep(c("Pre", "Post", "Follow-up"), times = n_id),
      levels = c("Pre", "Post", "Follow-up")
    )
  )
  d$group <- factor(rep(rep(c("Control", "Treatment"), each = n_id / 2), each = 3))
  subject_effect <- rnorm(n_id, 0, 5)
  d$score <- 50 + rep(subject_effect, each = 3) +
    ifelse(d$group == "Treatment", 3, 0) +
    ifelse(d$time == "Post", 4, ifelse(d$time == "Follow-up", 6, 0)) +
    rnorm(nrow(d), 0, 4)

  m <- lmer(score ~ time + group + (1 | id), data = d)
  p <- mellio_payload(m)

  expect_equal(p$type, "lmer_model_summary")
  types <- vapply(p$metadata$available_figures, function(row) row$type, character(1))
  expect_equal(types[[1]], "adjusted_means")
  expect_true(isTRUE(p$metadata$available_figures[[1]]$default))
  expect_true("coefficient_plot" %in% types)

  fig <- p$figure_data$adjusted_means
  expect_equal(fig$source, "lmer_emmeans")
  expect_equal(fig$mean_kind, "estimated_marginal")
  expect_equal(fig$factor$variable, "time")
  expect_true(isTRUE(fig$connect_levels))
  expect_equal(fig$outcome, "score")
  expect_equal(fig$subject$variable, "id")
  expect_equal(fig$subject$n, n_id)
  expect_equal(fig$ci_method, "emmeans")
  expect_length(fig$groups, 3L)
  expect_true(any(vapply(fig$marginalized_terms, function(row) {
    identical(row$variable, "group")
  }, logical(1))))
  expect_true(all(vapply(fig$groups, function(row) {
    is.numeric(row$mean) &&
      is.numeric(row$se) &&
      is.numeric(row$ci_lower) &&
      is.numeric(row$ci_upper) &&
      row$ci_lower < row$ci_upper
  }, logical(1))))

  pkg_names <- vapply(p$packages, `[[`, character(1), "name")
  expect_true("emmeans" %in% pkg_names)
})

test_that("lmer payload exposes continuous interaction estimated effect figures", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("emmeans")
  suppressMessages(library(lme4))

  set.seed(2030)
  n_id <- 60
  obs_per_id <- 4
  d <- data.frame(
    id = factor(rep(seq_len(n_id), each = obs_per_id)),
    group = factor(rep(rep(c("Control", "Treatment"), each = n_id / 2), each = obs_per_id))
  )
  d$stress <- rnorm(nrow(d), 0, 1)
  subject_effect <- rnorm(n_id, 0, 4)
  d$score <- 50 + rep(subject_effect, each = obs_per_id) +
    3 * d$stress +
    ifelse(d$group == "Treatment", 4, 0) +
    ifelse(d$group == "Treatment", 2 * d$stress, 0) +
    rnorm(nrow(d), 0, 4)

  m <- lmer(score ~ stress * group + (1 | id), data = d)
  p <- mellio_payload(m)

  expect_equal(p$type, "lmer_model_summary")
  types <- vapply(p$metadata$available_figures, function(row) row$type, character(1))
  expect_equal(types[[1]], "interaction_plot")
  expect_true(isTRUE(p$metadata$available_figures[[1]]$default))
  expect_true("coefficient_plot" %in% types)

  fig <- p$figure_data$interaction_plot
  expect_equal(fig$source, "lmer_emmeans")
  expect_equal(fig$mean_kind, "estimated_marginal")
  expect_equal(fig$interaction_kind, "continuous_by_categorical")
  expect_equal(fig$interaction_term, "stress:group")
  expect_equal(fig$x$variable, "stress")
  expect_equal(fig$x$type, "numeric")
  expect_equal(fig$moderator$variable, "group")
  expect_equal(fig$moderator$type, "categorical")
  expect_equal(fig$outcome, "score")
  expect_equal(fig$subject$variable, "id")
  expect_equal(fig$subject$n, n_id)
  expect_length(fig$grid, 160L)
  expect_true(all(vapply(fig$grid, function(row) {
    is.numeric(row$x) &&
      is.numeric(row$estimate) &&
      is.numeric(row$se) &&
      is.numeric(row$ci_lower) &&
      is.numeric(row$ci_upper) &&
      row$ci_lower < row$ci_upper
  }, logical(1))))

  pkg_names <- vapply(p$packages, `[[`, character(1), "name")
  expect_true("emmeans" %in% pkg_names)
})

test_that("lmer payload exposes continuous main-effect estimated effect figures", {
  skip_if_not_installed("lme4")
  skip_if_not_installed("emmeans")
  suppressMessages(library(lme4))

  set.seed(2031)
  n_id <- 60
  obs_per_id <- 4
  d <- data.frame(
    id = factor(rep(seq_len(n_id), each = obs_per_id)),
    group = factor(rep(rep(c("Control", "Treatment"), each = n_id / 2), each = obs_per_id))
  )
  d$stress <- rnorm(nrow(d), 0, 1)
  subject_effect <- rnorm(n_id, 0, 4)
  d$score <- 50 + rep(subject_effect, each = obs_per_id) +
    3 * d$stress +
    ifelse(d$group == "Treatment", 4, 0) +
    rnorm(nrow(d), 0, 4)

  m <- lmer(score ~ stress + group + (1 | id), data = d)
  p <- mellio_payload(m)

  expect_equal(p$type, "lmer_model_summary")
  types <- vapply(p$metadata$available_figures, function(row) row$type, character(1))
  expect_equal(types[[1]], "adjusted_means")
  expect_true("interaction_plot" %in% types)
  expect_true(isTRUE(p$metadata$available_figures[[1]]$default))

  means <- p$figure_data$adjusted_means
  expect_equal(means$source, "lmer_emmeans")
  expect_equal(means$mean_kind, "estimated_marginal")
  expect_equal(means$factor$variable, "group")
  expect_equal(means$outcome, "score")
  expect_true(any(vapply(means$covariates, function(row) {
    identical(row$variable, "stress")
  }, logical(1))))
  expect_length(means$groups, 2L)

  fig <- p$figure_data$interaction_plot
  expect_equal(fig$source, "lmer_emmeans")
  expect_equal(fig$mean_kind, "estimated_marginal")
  expect_equal(fig$interaction_kind, "continuous_main_effect")
  expect_equal(fig$interaction_term, "stress")
  expect_equal(fig$x$variable, "stress")
  expect_equal(fig$x$type, "numeric")
  expect_equal(fig$moderator$variable, "estimate")
  expect_equal(fig$outcome, "score")
  expect_equal(fig$subject$variable, "id")
  expect_equal(fig$subject$n, n_id)
  expect_length(fig$grid, 80L)
  expect_true(any(vapply(fig$marginalized_terms, function(row) {
    identical(row$variable, "group")
  }, logical(1))))
  expect_true(all(vapply(fig$grid, function(row) {
    is.numeric(row$x) &&
      is.numeric(row$estimate) &&
      is.numeric(row$se) &&
      is.numeric(row$ci_lower) &&
      is.numeric(row$ci_upper) &&
      row$ci_lower < row$ci_upper
  }, logical(1))))

  pkg_names <- vapply(p$packages, `[[`, character(1), "name")
  expect_true("emmeans" %in% pkg_names)
})

test_that("Single-predictor lm uses linear regression type_label", {
  p <- mellio_payload(lm(mpg ~ wt, data = mtcars))
  expect_equal(p$type_label, "Linear Regression")
})

test_that("lm preserves its own call when no .call is passed", {
  m <- lm(mpg ~ wt, data = mtcars)
  p <- mellio_payload(m)
  expect_match(p$call, "lm\\(formula = mpg ~ wt")
})

test_that("payload includes minimal provenance and packages by default", {
  p <- mellio_payload(t.test(extra ~ group, data = sleep))

  expect_true(!is.null(p$provenance))
  expect_equal(p$provenance$r_version, R.version.string)
  expect_true(nzchar(p$provenance$platform))
  expect_true(nzchar(p$provenance$mellio_version))
  expect_null(p$provenance$working_dir)
  expect_null(p$provenance$sender)
  expect_null(p$provenance$git)
  expect_null(p$provenance$script)

  expect_true(is.list(p$packages))
  pkg_names <- vapply(p$packages, `[[`, character(1), "name")
  expect_true(all(c("R", "mellio", "stats") %in% pkg_names))
})

test_that("ms_script_source filters out bridge-* internals", {
  src <- mellio:::ms_script_source()
  expect_true(is.null(src) ||
              !grepl("/bridge-", src$file %||% ""))
})

test_that("ms_rstudio_selection_line tolerates RStudio cursor shapes", {
  expect_equal(mellio:::ms_rstudio_selection_line(list(row = 12L, column = 3L)), 12L)
  expect_equal(mellio:::ms_rstudio_selection_line(c(row = 13L, column = 4L)), 13L)
  expect_equal(mellio:::ms_rstudio_selection_line(c(14L, 5L)), 14L)
  expect_true(is.na(mellio:::ms_rstudio_selection_line(NULL)))
  expect_true(is.na(mellio:::ms_rstudio_selection_line(character())))
})

test_that("ms_script_provenance returns NULL or a list with file", {
  prov <- mellio:::ms_script_provenance()
  expect_true(is.null(prov) || is.list(prov))
  if (!is.null(prov)) {
    expect_true(!is.null(prov$file) && nzchar(prov$file))
    if (!is.null(prov$hash)) {
      expect_match(prov$hash, "^[0-9a-f]{40}$")
    }
  }
})

test_that("lm payload includes data hash + n", {
  skip_if_not_installed("digest")
  p <- mellio_payload(lm(mpg ~ wt + cyl, data = mtcars))
  expect_true(!is.null(p$provenance$data))
  expect_match(p$provenance$data$hash, "^[0-9a-f]{40}$")
  expect_equal(p$provenance$data$n, 32)
})

test_that("data hash is reproducible across calls", {
  skip_if_not_installed("digest")
  h1 <- mellio_payload(lm(mpg ~ wt, data = mtcars))$provenance$data$hash
  h2 <- mellio_payload(lm(mpg ~ wt, data = mtcars))$provenance$data$hash
  expect_equal(h1, h2)
})

test_that("data hash changes when data changes", {
  skip_if_not_installed("digest")
  h1 <- mellio_payload(lm(mpg ~ wt, data = mtcars))$provenance$data$hash
  h2 <- mellio_payload(lm(mpg ~ wt, data = mtcars[1:20, ]))$provenance$data$hash
  expect_false(h1 == h2)
})

test_that("htest payload has no data hash (no model.frame available)", {
  p <- mellio_payload(t.test(extra ~ group, data = sleep))
  expect_null(p$provenance$data)
})

test_that("full provenance captures sender identity", {
  p <- withr::with_options(
    list(mellio.provenance = "full"),
    mellio_payload(t.test(extra ~ group, data = sleep))
  )
  expect_true(!is.null(p$provenance$sender))
  expect_true(nzchar(p$provenance$sender$user) || nzchar(p$provenance$sender$host))
})

test_that("full provenance captures git state when in a git repo", {
  skip_if_not(mellio:::ms_git_available(), "git not on PATH")
  skip_if_not(mellio:::ms_in_git_repo(),   "not inside a git working tree")

  p <- withr::with_options(
    list(mellio.provenance = "full"),
    mellio_payload(t.test(extra ~ group, data = sleep))
  )
  expect_true(!is.null(p$provenance$git))
  expect_match(p$provenance$git$commit, "^[0-9a-f]{40}$")
  expect_true(nzchar(p$provenance$git$branch))
  expect_true(is.logical(p$provenance$git$dirty))
})

test_that("ms_git_state returns NULL outside a git repo", {
  withr::with_dir(tempdir(), {
    # tempdir() is virtually never inside a git work tree
    res <- mellio:::ms_git_state()
    # Either NULL (not in repo) or has fields (if user has set up tempdir as git)
    expect_true(is.null(res) || is.list(res))
  })
})

test_that("anova on a single lm picks last non-Residuals term", {
  m <- lm(mpg ~ wt, data = mtcars)
  p <- mellio_payload(anova(m))
  expect_equal(p$type, "anova_single_model")
  expect_match(p$type_label, "focal term: wt")
  expect_equal(p$fields$statistic$name, "F")
  expect_length(p$fields$statistic$df, 2)
  expect_equal(p$fields$statistic$df[1], 1)   # term df
  expect_equal(p$fields$statistic$df[2], 30)  # residual df
  expect_equal(p$fields$term, "wt")
  expect_equal(p$fields$df_effect, 1)
  expect_equal(p$fields$df_error, 30)
  expect_equal(p$fields$ss_type, "type_i_sequential")
  expect_equal(p$fields$ss_type_label, "Type I (sequential)")
  expect_match(p$type_label, "Type I")
  expect_true(is.numeric(p$fields$sum_sq))
  expect_true(is.numeric(p$fields$mean_sq))
  expect_true(is.numeric(p$fields$residual_sum_sq))
  expect_true(is.numeric(p$fields$residual_mean_sq))
  expect_true(p$fields$p_value < 1e-9)
  expect_null(p$metadata$available_figures)
  expect_null(p$figure_data)
})

test_that("car::Anova tables preserve requested sums-of-squares type", {
  skip_if_not_installed("car")

  m <- lm(mpg ~ factor(cyl) + wt + hp, data = mtcars)
  p1 <- mellio_payload(anova(m), .call = "anova(m)")
  expect_equal(p1$fields$ss_type, "type_i_sequential")
  expect_equal(p1$fields$term, "hp")
  expect_null(p1$fields$model_kind)
  expect_equal(
    vapply(p1$fields$all_terms, function(row) row$is_focal, logical(1)),
    c(FALSE, FALSE, TRUE)
  )

  p2 <- mellio_payload(car::Anova(m, type = 2))
  expect_equal(p2$fields$ss_type, "type_ii")
  expect_equal(p2$fields$ss_type_label, "Type II")
  expect_match(p2$type_label, "Type II")
  expect_match(p2$type_label, "ANCOVA")
  expect_equal(p2$fields$model_kind, "ancova")
  expect_equal(p2$fields$term, "factor(cyl)")
  expect_equal(p2$fields$focal_terms, "factor(cyl)")
  expect_equal(p2$fields$control_terms, c("wt", "hp"))
  cyl <- p2$fields$all_terms[[1]]
  expect_equal(cyl$term, "factor(cyl)")
  expect_true(cyl$is_focal)
  expect_equal(cyl$f, 2.8776, tolerance = 1e-4)
  hp <- p2$fields$all_terms[[3]]
  expect_equal(hp$term, "hp")
  expect_false(hp$is_focal)

  old <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(old), add = TRUE)
  m3 <- lm(mpg ~ factor(cyl) + wt + hp, data = mtcars)
  p3 <- mellio_payload(car::Anova(m3, type = 3))
  expect_equal(p3$fields$ss_type, "type_iii")
  expect_equal(p3$fields$ss_type_label, "Type III")
  expect_match(p3$fields$ss_type_note, "contrast coding")
  expect_match(p3$type_label, "Type III")
  expect_match(p3$type_label, "ANCOVA")
  expect_equal(p3$fields$term, "factor(cyl)")
  expect_equal(p3$fields$focal_terms, "factor(cyl)")
  expect_equal(p3$fields$control_terms, c("wt", "hp"))
  expect_false("(Intercept)" %in% vapply(p3$fields$all_terms, function(row) row$term, character(1)))
})

test_that("car::Anova factorial tables preserve Type II and Type III labels", {
  skip_if_not_installed("car")

  set.seed(20260526)
  df <- expand.grid(
    a = factor(c("Control", "Treatment")),
    b = factor(c("Low", "Medium", "High")),
    rep = seq_len(8)
  )
  df$y <- rnorm(nrow(df)) +
    0.5 * (df$a == "Treatment") +
    0.3 * (df$b == "High") +
    0.7 * (df$a == "Treatment" & df$b == "High")
  m <- lm(y ~ a * b, data = df)

  p2 <- mellio_payload(car::Anova(m, type = 2))
  expect_equal(p2$fields$ss_type, "type_ii")
  expect_equal(p2$fields$ss_type_label, "Type II")
  expect_match(p2$type_label, "Type II")
  expect_null(p2$fields$model_kind)
  expect_equal(
    vapply(p2$fields$all_terms, function(row) row$term, character(1)),
    c("a", "b", "a:b")
  )
  expect_equal(p2$fields$all_terms[[3]]$term_type, "interaction")

  old <- options(contrasts = c("contr.sum", "contr.poly"))
  on.exit(options(old), add = TRUE)
  m3 <- lm(y ~ a * b, data = df)
  p3 <- mellio_payload(car::Anova(m3, type = 3))
  expect_equal(p3$fields$ss_type, "type_iii")
  expect_equal(p3$fields$ss_type_label, "Type III")
  expect_match(p3$fields$ss_type_note, "contrast coding")
  expect_match(p3$type_label, "Type III")
  expect_null(p3$fields$model_kind)
  expect_equal(
    vapply(p3$fields$all_terms, function(row) row$term, character(1)),
    c("a", "b", "a:b")
  )
})

test_that("direct aov objects route to ANOVA payloads", {
  m <- aov(mpg ~ factor(cyl), data = mtcars)
  p <- mellio_payload(m)
  expect_equal(p$type, "anova_single_model")
  expect_match(p$type_label, "focal term: factor\\(cyl\\)")
  expect_equal(p$fields$statistic$name, "F")
  expect_equal(p$fields$outcome, "mpg")
  expect_equal(p$fields$predictor, "cyl")
  expect_true(is.numeric(p$fields$p_value))
  expect_equal(p$fields$all_terms[[1]]$effect$name, "eta_sq_partial")
  expect_true(is.numeric(p$fields$all_terms[[1]]$eta_sq_partial))

  expect_equal(p$metadata$available_figures[[1]]$type, "adjusted_means")
  expect_equal(p$metadata$available_figures[[1]]$label, "Means plot")
  expect_true(isTRUE(p$metadata$available_figures[[1]]$default))

  fig <- p$figure_data$adjusted_means
  expect_equal(fig$source, "aov_one_way")
  expect_equal(fig$mean_kind, "estimated")
  expect_equal(fig$factor$variable, "cyl")
  expect_equal(fig$factor$term, "factor(cyl)")
  expect_equal(fig$outcome, "mpg")
  expect_equal(fig$ci_method, "residual_mean_square")
  expect_equal(fig$residual_df, stats::df.residual(m))
  expect_length(fig$groups, 3L)
  expect_equal(
    vapply(fig$groups, function(row) row$level, character(1)),
    c("4", "6", "8")
  )
  expected_means <- tapply(stats::model.frame(m)$mpg,
                           stats::model.frame(m)[["factor(cyl)"]],
                           mean)
  expect_equal(
    vapply(fig$groups, function(row) row$mean, numeric(1)),
    as.numeric(expected_means[c("4", "6", "8")]),
    tolerance = 1e-8
  )
  expect_true(all(vapply(fig$groups, function(row) {
    is.numeric(row$se) &&
      is.numeric(row$ci_lower) &&
      is.numeric(row$ci_upper) &&
      row$ci_lower < row$ci_upper
  }, logical(1))))
})

test_that("direct aov preserves labelled factor levels in means figures", {
  d <- mtcars
  d$transmission <- factor(
    d$am,
    levels = c(0, 1),
    labels = c("Automatic", "Manual")
  )
  m <- aov(mpg ~ transmission, data = d)
  p <- mellio_payload(m)

  expect_equal(p$type, "anova_single_model")
  expect_equal(p$fields$predictor, "transmission")
  expect_equal(p$figure_data$adjusted_means$factor$variable, "transmission")
  expect_equal(
    vapply(p$figure_data$adjusted_means$groups, function(row) row$level, character(1)),
    c("Automatic", "Manual")
  )
})

test_that("repeated-measures aovlist exposes within-subject means figure", {
  set.seed(2026)
  n <- 36
  subject_effect <- rnorm(n, 0, 6)
  d <- data.frame(
    id = factor(rep(seq_len(n), each = 3)),
    time = factor(
      rep(c("Pre", "Post", "Follow-up"), times = n),
      levels = c("Pre", "Post", "Follow-up")
    )
  )
  d$score <- 50 + rep(subject_effect, each = 3) +
    ifelse(d$time == "Post", 5, ifelse(d$time == "Follow-up", 8, 0)) +
    rnorm(nrow(d), 0, 4)

  fit <- aov(score ~ time + Error(id / time), data = d)
  p <- mellio_payload(fit)

  expect_equal(p$card_kind, "table")
  expect_equal(p$type_label, "Repeated-measures ANOVA table")
  expect_equal(p$fields$model_kind, "repeated_measures_anova")
  expect_equal(p$fields$outcome, "score")
  expect_equal(p$fields$subject, "id")
  expect_equal(p$fields$subject_variable, "id")
  expect_equal(p$fields$n, n)
  expect_equal(p$fields$within_terms, "time")
  expect_equal(p$fields$factors[[1]]$variable, "time")
  expect_equal(p$fields$factors[[1]]$role, "within")
  expect_equal(
    vapply(p$fields$factors[[1]]$levels, function(row) row$value, character(1)),
    c("Pre", "Post", "Follow-up")
  )
  expect_match(p$fields$table_note, "All listed factors are within-subjects")
  expect_match(p$fields$table_note, "Partial η²")
  expect_match(p$fields$table_note, "Degrees of freedom")
  expect_false(grepl("df_effect|df_error|SS_effect|SS_error", p$fields$table_note))
  expect_match(p$fields$sphericity_note, "Sphericity tests")
  expect_true("eta_sq_partial" %in% vapply(p$fields$columns, function(row) row$key, character(1)))
  eta_col <- p$fields$columns[[which(vapply(
    p$fields$columns, function(row) row$key, character(1)
  ) == "eta_sq_partial")]]
  expect_equal(eta_col$format, "bounded")
  expect_equal(eta_col$label, "partial η²")
  expect_true("meansq" %in% vapply(p$fields$columns, function(row) row$key, character(1)))
  types <- vapply(p$metadata$available_figures, function(row) row$type, character(1))
  expect_true("adjusted_means" %in% types)

  rows <- p$fields$rows
  idx <- which(vapply(rows, function(row) identical(row$term, "time"), logical(1)))
  expect_length(idx, 1L)
  time_row <- rows[[idx]]
  expect_equal(time_row$df1, 2)
  expect_equal(time_row$df2, 70)
  expect_equal(time_row$f, time_row$statistic)
  expect_equal(time_row$stratum, "time × subjects")
  expect_equal(time_row$stratum_raw, "id:time")
  expect_equal(time_row$effect$name, "eta_sq_partial")
  expect_true(is.numeric(time_row$eta_sq_partial))

  fig <- p$figure_data$adjusted_means
  expect_equal(fig$source, "repeated_measures_aov")
  expect_equal(fig$mean_kind, "within_subject")
  expect_equal(fig$factor$variable, "time")
  expect_equal(fig$subject$variable, "id")
  expect_equal(fig$subject$n, n)
  expect_equal(fig$ci_method, "within_subject_morey")
  expect_equal(
    vapply(fig$groups, function(row) row$level, character(1)),
    c("Pre", "Post", "Follow-up")
  )
  expect_true(all(vapply(fig$groups, function(row) {
    is.numeric(row$se) &&
      is.numeric(row$ci_lower) &&
      is.numeric(row$ci_upper) &&
      row$ci_lower < row$ci_upper
  }, logical(1))))
  expect_equal(fig$focal_test$statistic$name, "F")
  expect_equal(as.numeric(fig$focal_test$statistic$df), c(2, 70))
  expect_equal(fig$focal_test$effect$name, "eta_sq_partial")
})

test_that("multi-factor repeated-measures aovlist exposes within-subject interaction figure", {
  set.seed(2027)
  n <- 24
  d <- expand.grid(
    id = factor(seq_len(n)),
    time = factor(c("Pre", "Post"), levels = c("Pre", "Post")),
    dose = factor(c("Low", "High"), levels = c("Low", "High"))
  )
  subject_effect <- rnorm(n, 0, 5)
  d$score <- 50 + subject_effect[as.integer(d$id)] +
    ifelse(d$time == "Post", 4, 0) +
    ifelse(d$dose == "High", 2, 0) +
    ifelse(d$time == "Post" & d$dose == "High", 3, 0) +
    rnorm(nrow(d), 0, 3)

  fit <- aov(score ~ time * dose + Error(id / (time * dose)), data = d)
  p <- mellio_payload(fit)

  expect_equal(p$card_kind, "table")
  expect_equal(p$type_label, "Repeated-measures ANOVA table")
  expect_equal(p$fields$model_kind, "repeated_measures_anova")
  expect_equal(p$fields$outcome, "score")
  expect_equal(p$fields$subject, "id")
  expect_equal(p$fields$n, n)
  expect_equal(p$fields$within_terms, c("time", "dose", "time:dose"))
  expect_match(p$fields$sphericity_note, "not applicable")
  expect_match(p$fields$table_note, "all within-subject factors have two levels")
  expect_equal(p$fields$columns[[1]]$label, "Error stratum")
  expect_setequal(
    vapply(p$fields$factors, function(row) row$variable, character(1)),
    c("time", "dose")
  )

  terms <- vapply(p$fields$rows, function(row) row$term, character(1))
  expect_true(all(c("time", "dose", "time:dose", "Residuals") %in% terms))
  strata <- vapply(p$fields$rows, function(row) row$stratum, character(1))
  expect_false("Between subjects" %in% strata)
  expect_true(all(c("time × subjects", "dose × subjects", "time × dose × subjects") %in% strata))
  effect_rows <- Filter(function(row) !is.null(row$f) && !is.na(row$f), p$fields$rows)
  expect_true(all(vapply(effect_rows, function(row) {
    is.numeric(row$eta_sq_partial) &&
      identical(row$effect$name, "eta_sq_partial")
  }, logical(1))))
  types <- vapply(p$metadata$available_figures, function(row) row$type, character(1))
  expect_false("adjusted_means" %in% types)
  expect_true("interaction_plot" %in% types)

  fig <- p$figure_data$interaction_plot
  expect_equal(fig$source, "repeated_measures_aov")
  expect_equal(fig$mean_kind, "observed_cell_means")
  expect_equal(fig$interaction_kind, "categorical_by_categorical")
  expect_equal(fig$interaction_term, "time:dose")
  expect_equal(fig$x$variable, "time")
  expect_equal(fig$moderator$variable, "dose")
  expect_equal(fig$subject$variable, "id")
  expect_equal(fig$subject$n, n)
  expect_equal(fig$ci_method, "within_subject_morey")
  expect_length(fig$grid, 4L)
  expect_true(all(vapply(fig$grid, function(row) {
    is.numeric(row$estimate) &&
      is.numeric(row$ci_lower) &&
      is.numeric(row$ci_upper) &&
      row$ci_lower < row$ci_upper
  }, logical(1))))
  expect_equal(fig$focal_test$effect$name, "eta_sq_partial")
})

test_that("repeated-measures aovlist uses natural level order for report metadata", {
  set.seed(2031)
  n <- 12
  d <- expand.grid(
    id = factor(seq_len(n)),
    time = factor(c("Pre", "Post")),
    dose = factor(c("Low", "High"))
  )
  expect_equal(levels(d$time), c("Post", "Pre"))
  expect_equal(levels(d$dose), c("High", "Low"))

  subject_effect <- rnorm(n, 0, 4)
  d$score <- 40 + subject_effect[as.integer(d$id)] +
    ifelse(d$time == "Post", 3, 0) +
    ifelse(d$dose == "High", 2, 0) +
    rnorm(nrow(d), 0, 2)

  fit <- aov(score ~ time * dose + Error(id / (time * dose)), data = d)
  p <- mellio_payload(fit)
  factors <- setNames(
    p$fields$factors,
    vapply(p$fields$factors, function(row) row$variable, character(1))
  )

  expect_equal(
    vapply(factors$time$levels, function(row) row$value, character(1)),
    c("Pre", "Post")
  )
  expect_equal(
    vapply(factors$dose$levels, function(row) row$value, character(1)),
    c("Low", "High")
  )
})

test_that("afex_aov exposes repeated-measures ANOVA, sphericity, and means figure", {
  skip_if_not_installed("afex")

  set.seed(2026)
  n <- 36
  subject_effect <- rnorm(n, 0, 6)
  d <- data.frame(
    id = factor(rep(seq_len(n), each = 3)),
    time = factor(
      rep(c("Pre", "Post", "Follow-up"), times = n),
      levels = c("Pre", "Post", "Follow-up")
    )
  )
  d$score <- 50 + rep(subject_effect, each = 3) +
    ifelse(d$time == "Post", 5, ifelse(d$time == "Follow-up", 8, 0)) +
    rnorm(nrow(d), 0, 4)

  fit <- afex::aov_ez(id = "id", dv = "score", within = "time", data = d)
  p <- mellio_payload(fit)

  expect_equal(p$type, "custom_anova_table")
  expect_equal(p$card_kind, "table")
  expect_match(p$type_label, "Repeated-measures ANOVA table")
  expect_equal(p$fields$source, "R afex_aov")
  expect_equal(p$fields$model_kind, "repeated_measures_anova")
  expect_equal(p$fields$outcome, "score")
  expect_equal(p$fields$subject, "id")
  expect_equal(p$fields$n, n)
  expect_equal(p$fields$within_terms, "time")
  expect_equal(p$fields$ss_type_label, "Type 3")
  expect_equal(p$fields$effect_size, "ges")
  expect_equal(
    p$fields$columns[[which(vapply(p$fields$columns, function(row) row$key, character(1)) == "ges")]]$label,
    "generalized η²"
  )

  rows <- p$fields$rows
  idx <- which(vapply(rows, function(row) identical(row$term, "time"), logical(1)))
  expect_length(idx, 1L)
  time_row <- rows[[idx]]
  expect_true(is.numeric(time_row$df1))
  expect_true(is.numeric(time_row$df2))
  expect_true(is.numeric(time_row$f))
  expect_true(is.numeric(time_row$p_value))
  expect_true(is.numeric(time_row$ges))
  expect_equal(time_row$effect$name, "eta_sq_generalized")

  expect_length(p$fields$sphericity_tests, 1L)
  expect_equal(p$fields$sphericity_tests[[1]]$term, "time")
  expect_true(is.numeric(p$fields$sphericity_tests[[1]]$w))
  expect_true(is.numeric(p$fields$sphericity_tests[[1]]$p_value))
  expect_length(p$fields$sphericity_corrections, 1L)
  expect_true(is.numeric(p$fields$sphericity_corrections[[1]]$gg_epsilon))

  types <- vapply(p$metadata$available_figures, function(row) row$type, character(1))
  expect_true("adjusted_means" %in% types)
  fig <- p$figure_data$adjusted_means
  expect_equal(fig$source, "afex_aov")
  expect_equal(fig$mean_kind, "within_subject")
  expect_equal(fig$factor$variable, "time")
  expect_equal(fig$subject$n, n)
  expect_equal(fig$ci_method, "within_subject_morey")
  expect_equal(fig$focal_test$effect$name, "eta_sq_generalized")
})

test_that("afex_aov exposes mixed ANOVA interaction profile figures", {
  skip_if_not_installed("afex")

  set.seed(2027)
  n_per_group <- 18
  ids <- factor(seq_len(n_per_group * 2))
  group_by_id <- factor(rep(c("Control", "Treatment"), each = n_per_group))
  d <- data.frame(
    id = factor(rep(ids, each = 3)),
    time = factor(
      rep(c("Pre", "Post", "Follow-up"), times = length(ids)),
      levels = c("Pre", "Post", "Follow-up")
    )
  )
  d$group <- factor(rep(group_by_id, each = 3))
  subject_effect <- rnorm(n_per_group * 2, 0, 5)
  d$score <- 50 + rep(subject_effect, each = 3) +
    ifelse(d$group == "Treatment", 3, 0) +
    ifelse(d$time == "Post", 4, ifelse(d$time == "Follow-up", 6, 0)) +
    ifelse(d$group == "Treatment" & d$time == "Follow-up", 4, 0) +
    rnorm(nrow(d), 0, 4)

  fit <- afex::aov_ez(id = "id", dv = "score", within = "time", between = "group", data = d)
  p <- mellio_payload(fit)

  expect_equal(p$type, "custom_anova_table")
  expect_match(p$type_label, "Mixed ANOVA table")
  expect_equal(p$fields$source, "R afex_aov")
  expect_equal(p$fields$model_kind, "mixed_anova")
  expect_equal(p$fields$outcome, "score")
  expect_equal(p$fields$subject, "id")
  expect_equal(p$fields$n, n_per_group * 2)
  expect_equal(p$fields$within_terms, "time")
  expect_equal(p$fields$between_terms, "group")
  expect_equal(
    vapply(p$fields$factors, function(row) paste(row$variable, row$role, sep = ":"), character(1)),
    c("time:within", "group:between")
  )

  terms <- vapply(p$fields$rows, function(row) row$term, character(1))
  expect_true(all(c("group", "time", "group:time") %in% terms))
  expect_true("time" %in% vapply(p$fields$sphericity_tests, function(row) row$term, character(1)))
  expect_true("group:time" %in% vapply(p$fields$sphericity_tests, function(row) row$term, character(1)))

  types <- vapply(p$metadata$available_figures, function(row) row$type, character(1))
  expect_true("interaction_plot" %in% types)
  fig <- p$figure_data$interaction_plot
  expect_equal(fig$source, "afex_aov")
  expect_equal(fig$mean_kind, "observed_cell_means")
  expect_equal(fig$interaction_kind, "categorical_by_categorical")
  expect_equal(fig$interaction_term, "group:time")
  expect_equal(fig$x$variable, "time")
  expect_equal(fig$moderator$variable, "group")
  expect_equal(fig$subject$n, n_per_group * 2)
  expect_length(fig$grid, 6L)
  expect_true(all(vapply(fig$grid, function(row) {
    is.numeric(row$estimate) &&
      is.numeric(row$ci_lower) &&
      is.numeric(row$ci_upper) &&
      row$ci_lower < row$ci_upper
  }, logical(1))))
  expect_equal(fig$focal_test$effect$name, "eta_sq_generalized")
})

test_that("lm factor-plus-covariate objects expose adjusted means from emmeans", {
  skip_if_not_installed("emmeans")

  m <- lm(mpg ~ factor(cyl) + wt, data = mtcars)
  p <- mellio_payload(m)

  expect_equal(p$type, "lm_model_summary")
  expect_equal(p$type_label, "Linear Regression")
  expect_null(p$fields$model_kind)
  types <- vapply(p$metadata$available_figures, function(row) row$type, character(1))
  expect_equal(types[[1]], "adjusted_means")
  expect_true("coefficient_plot" %in% types)

  fig <- p$figure_data$adjusted_means
  expect_equal(fig$source, "ancova_emmeans")
  expect_equal(fig$mean_kind, "estimated_marginal")
  expect_equal(fig$factor$term, "factor(cyl)")
  expect_equal(fig$factor$variable, "cyl")
  expect_equal(fig$outcome, "mpg")
  expect_equal(fig$ci_method, "emmeans")
  expect_length(fig$covariates, 1L)
  expect_equal(fig$covariates[[1]]$variable, "wt")
  expect_equal(fig$covariates[[1]]$rule, "sample_mean")
  expect_equal(fig$covariates[[1]]$value, mean(mtcars$wt), tolerance = 1e-8)

  emm <- as.data.frame(summary(
    emmeans::emmeans(m, ~ factor(cyl), at = list(wt = mean(mtcars$wt))),
    infer = c(TRUE, FALSE)
  ))
  expect_equal(
    vapply(fig$groups, function(row) row$level, character(1)),
    as.character(emm$cyl)
  )
  expect_equal(
    vapply(fig$groups, function(row) row$mean, numeric(1)),
    emm$emmean,
    tolerance = 1e-8
  )
  expect_equal(
    vapply(fig$groups, function(row) row$se, numeric(1)),
    emm$SE,
    tolerance = 1e-8
  )
  expect_equal(fig$focal_test$statistic$name, "F")
  expect_equal(as.numeric(fig$focal_test$statistic$df), c(2, 28))
  expect_true(is.numeric(fig$focal_test$p_value))
  expect_equal(fig$focal_test$effect$name, "eta_sq_partial")
})

test_that("ANCOVA supports multiple numeric covariates", {
  skip_if_not_installed("emmeans")

  m <- lm(mpg ~ factor(cyl) + wt + hp, data = mtcars)
  p <- mellio_payload(m)
  fig <- p$figure_data$adjusted_means

  expect_equal(fig$source, "ancova_emmeans")
  expect_equal(
    vapply(fig$covariates, function(row) row$variable, character(1)),
    c("wt", "hp")
  )
  expect_equal(
    vapply(fig$covariates, function(row) row$value, numeric(1)),
    c(mean(mtcars$wt), mean(mtcars$hp)),
    tolerance = 1e-8
  )
})

test_that("ANCOVA aov objects use the categorical factor as focal term", {
  skip_if_not_installed("emmeans")

  m <- aov(mpg ~ factor(cyl) + wt, data = mtcars)
  p <- mellio_payload(m)

  expect_equal(p$type, "anova_single_model")
  expect_match(p$type_label, "ANCOVA")
  expect_equal(p$fields$model_kind, "ancova")
  expect_equal(p$fields$term, "factor(cyl)")
  expect_equal(p$fields$focal_terms, "factor(cyl)")
  expect_equal(p$fields$control_terms, "wt")
  expect_equal(p$metadata$available_figures[[1]]$type, "adjusted_means")
  expect_equal(p$figure_data$adjusted_means$source, "ancova_emmeans")
})

test_that("group by covariate interaction does not emit adjusted means", {
  m <- lm(mpg ~ factor(cyl) * wt, data = mtcars)
  p <- mellio_payload(m)

  types <- vapply(p$metadata$available_figures, function(row) row$type, character(1))
  expect_true("interaction_plot" %in% types)
  expect_null(p$figure_data$adjusted_means)
})

test_that("two-way ANOVA populates fields$all_terms additively", {
  set.seed(42)
  df <- data.frame(
    y = rnorm(120),
    a = factor(sample(c("x", "y", "z"), 120, replace = TRUE)),
    b = factor(sample(c("p", "q"), 120, replace = TRUE))
  )
  df$y <- df$y + 0.6 * (df$a == "y") + 1.2 * (df$b == "q")
  p <- mellio_payload(aov(y ~ a * b, data = df))

  # Existing inline-shape contract is unchanged (card stays inline,
  # focal-term scalars present, partial η² subscript preserved).
  expect_equal(p$card_kind, "inline")
  expect_equal(p$type, "anova_single_model")
  expect_match(p$type_label, "terms: a, b, a:b", fixed = TRUE)
  expect_equal(p$fields$statistic$name, "F")
  expect_equal(p$fields$term, "a:b")        # focal = last non-Residuals row
  expect_equal(p$fields$effect$name, "eta_sq_partial")

  # Additive: all_terms carries every non-Residuals term, tags type +
  # focal flag. Residuals row is NOT in this array — it lives in the
  # field-level residual_* scalars to avoid F(—) artifacts in renderers
  # that don't filter by row.f.
  expect_true(is.list(p$fields$all_terms))
  expect_length(p$fields$all_terms, 3L)
  term_names <- vapply(p$fields$all_terms, function(r) r$term, character(1))
  expect_setequal(term_names, c("a", "b", "a:b"))

  types <- setNames(
    vapply(p$fields$all_terms, function(r) r$term_type, character(1)),
    term_names
  )
  expect_equal(types[["a"]],   "main")
  expect_equal(types[["b"]],   "main")
  expect_equal(types[["a:b"]], "interaction")

  focal_flags <- setNames(
    vapply(p$fields$all_terms, function(r) r$is_focal, logical(1)),
    term_names
  )
  expect_false(focal_flags[["a"]])
  expect_false(focal_flags[["b"]])
  expect_true(focal_flags[["a:b"]])
  expect_equal(p$metadata$available_figures[[1]]$type, "interaction_plot")
  expect_true(isTRUE(p$metadata$available_figures[[1]]$default))
  expect_equal(p$fields$ss_type_label, "Type I (sequential)")
  expect_equal(p$figure_data$interaction_plot$interaction_kind, "categorical_by_categorical")
  expect_equal(p$figure_data$interaction_plot$source, "aov_interaction")
  expect_equal(p$figure_data$interaction_plot$mean_kind, "estimated_marginal")
  expect_equal(p$figure_data$interaction_plot$y_label, "Estimated marginal mean y")
  expect_equal(p$figure_data$interaction_plot$x$variable, "a")
  expect_equal(p$figure_data$interaction_plot$moderator$variable, "b")
  expect_length(p$figure_data$interaction_plot$grid, 6L)
})

test_that("additive multi-term ANOVA labels all terms and keeps residual scalars", {
  set.seed(20260526)
  df <- data.frame(
    y = rnorm(90),
    a = factor(rep(c("Control", "Treatment", "Placebo"), each = 30)),
    b = factor(rep(c("Low", "High"), length.out = 90))
  )
  p <- mellio_payload(aov(y ~ a + b, data = df))

  expect_equal(p$type, "anova_single_model")
  expect_match(p$type_label, "terms: a, b", fixed = TRUE)
  expect_equal(
    vapply(p$fields$all_terms, function(row) row$term, character(1)),
    c("a", "b")
  )
  expect_equal(
    vapply(p$fields$all_terms, function(row) row$term_type, character(1)),
    c("main", "main")
  )
  expect_true(is.numeric(p$fields$df_error))
  expect_true(is.numeric(p$fields$residual_sum_sq))
  expect_true(is.numeric(p$fields$residual_mean_sq))
  expect_equal(p$metadata$available_figures[[1]]$type, "adjusted_means")
  expect_equal(p$fields$term, "b")
  expect_equal(p$figure_data$adjusted_means$factor$term, "b")
  expect_equal(
    vapply(p$figure_data$adjusted_means$marginalized_terms, function(row) row$term, character(1)),
    "a"
  )

  p_focal <- mellio_payload(aov(y ~ a + b, data = df), focal = "a")
  expect_equal(p_focal$fields$term, "a")
  expect_equal(p_focal$figure_data$adjusted_means$factor$term, "a")
})

test_that("three-term additive ANOVA keeps each main-effect row", {
  set.seed(20260526)
  df <- expand.grid(
    a = factor(c("A1", "A2")),
    b = factor(c("B1", "B2", "B3")),
    c = factor(c("C1", "C2")),
    rep = seq_len(5)
  )
  df$y <- rnorm(nrow(df)) +
    0.4 * (df$a == "A2") +
    0.2 * (df$b == "B3") -
    0.3 * (df$c == "C2")
  p <- mellio_payload(aov(y ~ a + b + c, data = df))

  expect_equal(p$type, "anova_single_model")
  expect_match(p$type_label, "terms: a, b, c", fixed = TRUE)
  expect_equal(
    vapply(p$fields$all_terms, function(row) row$term, character(1)),
    c("a", "b", "c")
  )
  expect_equal(
    vapply(p$fields$all_terms, function(row) row$term_type, character(1)),
    c("main", "main", "main")
  )
  expect_equal(p$metadata$available_figures[[1]]$type, "adjusted_means")
  expect_equal(p$fields$term, "c")
  expect_equal(p$figure_data$adjusted_means$factor$term, "c")
  expect_equal(
    vapply(p$figure_data$adjusted_means$marginalized_terms, function(row) row$term, character(1)),
    c("a", "b")
  )
})

test_that("additive factorial lm objects expose main-effect estimated marginal means", {
  skip_if_not_installed("emmeans")

  d <- ToothGrowth
  d$dose_f <- factor(d$dose)
  m <- lm(len ~ supp + dose_f, data = d)
  p <- mellio_payload(m)

  types <- vapply(p$metadata$available_figures, function(row) row$type, character(1))
  expect_equal(types[[1]], "adjusted_means")
  expect_true("coefficient_plot" %in% types)

  fig <- p$figure_data$adjusted_means
  expect_equal(fig$source, "factorial_emmeans")
  expect_equal(fig$mean_kind, "estimated_marginal")
  expect_equal(fig$factor$variable, "supp")
  expect_equal(fig$factor$term, "supp")
  expect_equal(fig$outcome, "len")
  expect_equal(fig$ci_method, "emmeans")
  expect_length(fig$marginalized_terms, 1L)
  expect_equal(fig$marginalized_terms[[1]]$variable, "dose_f")
  expect_equal(fig$marginalized_terms[[1]]$rule, "estimated_marginal")

  emm <- as.data.frame(summary(
    emmeans::emmeans(m, ~ supp),
    infer = c(TRUE, FALSE)
  ))
  expect_equal(
    vapply(fig$groups, function(row) row$level, character(1)),
    as.character(emm$supp)
  )
  expect_equal(
    vapply(fig$groups, function(row) row$mean, numeric(1)),
    emm$emmean,
    tolerance = 1e-8
  )
})

test_that("anova supports explicit ANCOVA focal and control roles", {
  df <- data.frame(
    y = c(3.1, 3.4, 3.8, 4.2, 4.0, 4.6, 4.9, 5.1),
    condition = factor(rep(c("control", "treatment"), each = 4)),
    age = c(21, 23, 24, 26, 22, 25, 27, 28)
  )
  m <- lm(y ~ condition + age, data = df)
  p <- mellio_payload(anova(m), focal = "condition", controls = "age")

  expect_equal(p$type, "anova_single_model")
  expect_equal(p$fields$model_kind, "ancova")
  expect_equal(p$fields$term, "condition")
  expect_equal(p$fields$focal_terms, "condition")
  expect_equal(p$fields$control_terms, "age")
})

test_that("anova model comparison emits the last-row F", {
  m1 <- lm(mpg ~ wt, data = mtcars)
  m2 <- lm(mpg ~ wt + cyl, data = mtcars)
  p <- mellio_payload(anova(m1, m2))
  expect_equal(p$type, "anova_model_comparison")
  expect_equal(p$type_label, "Model comparison (2 models)")
  expect_equal(p$fields$statistic$name, "F")
  expect_length(p$fields$statistic$df, 2)
  expect_equal(p$fields$model_count, 2)
  expect_equal(p$fields$comparison, "model 1 vs. model 2")
  expect_equal(p$fields$df_change, 1)
  expect_equal(p$fields$residual_df, 29)
  expect_equal(as.numeric(p$fields$statistic$df), c(1, 29))
  expect_true(is.numeric(p$fields$sum_of_sq))
  expect_true(is.numeric(p$fields$residual_sum_sq))
  expect_true(is.numeric(p$fields$statistic$value))
  expect_true(p$fields$p_value < 0.01)
})

test_that("mellio_compare emits hierarchical regression comparison cards", {
  m1 <- lm(mpg ~ wt + hp, data = mtcars)
  m2 <- lm(mpg ~ wt + hp + cyl, data = mtcars)
  p <- mellio_compare(m1, m2)

  expect_equal(p$type, "hierarchical_regression_comparison")
  expect_equal(p$fields$model_kind, "hierarchical_regression")
  expect_equal(p$fields$model_count, 2)
  expect_equal(p$fields$outcome, "mpg")
  expect_length(p$fields$models, 2)
  expect_length(p$fields$comparisons, 1)
  expect_equal(as.character(p$fields$comparisons[[1]]$added_terms), "cyl")
  expect_equal(p$fields$df_change, 1)
  expect_equal(p$fields$residual_df, 28)
  expect_equal(p$fields$comparisons[[1]]$df_change, 1)
  expect_equal(p$fields$comparisons[[1]]$residual_df, 28)
  expect_equal(p$fields$statistic$name, "F")
  expect_length(p$fields$statistic$df, 2)
  expect_equal(as.numeric(p$fields$statistic$df), c(1, 28))
  expect_true(is.numeric(p$fields$r_squared))
  expect_true(is.numeric(p$fields$adj_r_squared))
  expect_true(is.numeric(p$fields$r_squared_change))
  expect_equal(
    p$fields$r_squared_change,
    summary(m2)$r.squared - summary(m1)$r.squared,
    tolerance = 1e-10
  )
  expect_true(is.numeric(p$fields$p_change))
})

test_that("mellio_open accepts prebuilt comparison payloads", {
  m1 <- lm(mpg ~ wt + hp, data = mtcars)
  m2 <- lm(mpg ~ wt + hp + cyl, data = mtcars)
  url <- quiet_mellio_open(mellio_compare(m1, m2), browse = FALSE)
  expect_match(url, "#stats/payload=")
})

test_that("mellio_open on lm round-trips through the URL", {
  url <- quiet_mellio_open(lm(mpg ~ wt + cyl, data = mtcars), browse = FALSE)
  b64 <- sub(".*payload=", "", url)
  b64 <- sub("&.*$", "", b64)
  std <- chartr("-_", "+/", b64)
  pad <- (4 - nchar(std) %% 4) %% 4
  if (pad > 0) std <- paste0(std, strrep("=", pad))
  parsed <- jsonlite::fromJSON(rawToChar(jsonlite::base64_dec(std)),
                               simplifyVector = FALSE)
  expect_equal(parsed$type, "lm_model_summary")
  expect_match(parsed$call, "lm\\(mpg ~ wt \\+ cyl, data = mtcars\\)")
})

# ── JSON serialisation ───────────────────────────────────────────────

test_that("mellio_to_json produces valid JSON matching the schema", {
  p <- mellio_payload(t.test(extra ~ group, data = sleep))
  json <- mellio_to_json(p)
  parsed <- jsonlite::fromJSON(json, simplifyVector = FALSE)

  expect_equal(parsed$schema_version, "0.1")
  expect_equal(parsed$card_kind, "inline")
  expect_equal(parsed$type, "welch_t_test")
  expect_true("statistic" %in% names(parsed$fields))
  expect_true("p_value"   %in% names(parsed$fields))
})

test_that("CI is serialised as a JSON array, not a scalar", {
  p <- mellio_payload(t.test(extra ~ group, data = sleep))
  json <- mellio_to_json(p)
  # Look for the CI as an array literal
  expect_match(json, "\"ci\":\\[-?[0-9.]+,-?[0-9.]+\\]")
})

# ── mellio_open ──────────────────────────────────────────────────────────

test_that("mellio_open returns a URL with payload in the hash", {
  withr::with_options(
    list(mellio.editor_url = "https://example.com"),
    {
      url <- quiet_mellio_open(t.test(extra ~ group, data = sleep), browse = FALSE)
      expect_match(url, "^https://example.com/#stats/payload=")
    }
  )
})

test_that("options('mellio.editor_url') overrides the base URL", {
  withr::with_options(
    list(mellio.editor_url = "https://example.com"),
    {
      url <- quiet_mellio_open(t.test(extra ~ group, data = sleep), browse = FALSE)
      expect_match(url, "^https://example.com/#stats/payload=")
    }
  )
})

test_that("MELLIO_URL env var does not override the public default", {
  withr::with_options(list(mellio.editor_url = NULL), {
    withr::with_envvar(c(MELLIO_URL = "https://env-example.com"), {
      url <- quiet_mellio_open(t.test(extra ~ group, data = sleep), browse = FALSE)
      expect_no_match(url, "^https://env-example.com/#stats/payload=")
      expect_match(url, "^https://www.mellioapp.com/#stats/payload=")
    })
  })
})

test_that("default base URL is the production Mellio app", {
  withr::with_options(list(mellio.editor_url = NULL), {
    withr::with_envvar(c(MELLIO_URL = ""), {
      url <- quiet_mellio_open(t.test(extra ~ group, data = sleep), browse = FALSE)
      expect_match(url, "^https://")
      expect_false(grepl("127.0.0.1|localhost", url))
    })
  })
})

test_that("mellio_open captures the user's original call", {
  url <- quiet_mellio_open(cor.test(mtcars$mpg, mtcars$wt), browse = FALSE)
  b64 <- sub(".*payload=", "", url)
  b64 <- sub("&.*$", "", b64)
  std <- chartr("-_", "+/", b64)
  pad <- (4 - nchar(std) %% 4) %% 4
  if (pad > 0) std <- paste0(std, strrep("=", pad))
  parsed <- jsonlite::fromJSON(rawToChar(jsonlite::base64_dec(std)),
                               simplifyVector = FALSE)
  expect_equal(parsed$call, "cor.test(mtcars$mpg, mtcars$wt)")
})

test_that("two-sample t-test emits group level labels in fields$groups", {
  # sleep$group is a factor with levels "1" / "2"
  p <- mellio_payload(t.test(extra ~ group, data = sleep))
  expect_true(is.list(p$fields$groups))
  expect_length(p$fields$groups, 2)
  labels <- vapply(p$fields$groups, function(g) g$label, character(1))
  expect_setequal(labels, c("1", "2"))
})

test_that("two-sample t-test preserves descriptive group level labels", {
  df <- data.frame(
    motive_enjoyment = c(3.2, 3.8, 4.1, 2.9, 2.7, 3.0),
    current_member_clean = factor(
      c("current members", "current members", "current members",
        "non-members", "non-members", "non-members")
    )
  )
  p <- mellio_payload(t.test(motive_enjoyment ~ current_member_clean, data = df))
  labels <- vapply(p$fields$groups, function(g) g$label, character(1))
  expect_setequal(labels, c("current members", "non-members"))
})

test_that("Student t-test (var.equal) also emits group labels", {
  p <- mellio_payload(t.test(extra ~ group, data = sleep, var.equal = TRUE))
  expect_length(p$fields$groups, 2)
})

test_that("paired t-test with non-reference inputs does not auto-enrich", {
  p <- mellio_payload(t.test(c(1, 2, 4, 7), c(2, 3, 5, 9), paired = TRUE))
  expect_null(p$fields$groups)
  expect_null(p$figure_data)
  expect_null(p$metadata$available_figures)
})

test_that("one-sample t-test with non-reference input does not auto-enrich", {
  p <- mellio_payload(t.test(c(1, 2, 3, 4), mu = 0))
  expect_null(p$fields$groups)
  expect_equal(p$figure_data$one_sample_mean_plot$source, "one_sample_t_test")
  expect_equal(p$figure_data$one_sample_mean_plot$null_value, 0)
})

test_that("chi-square test emits sample_size and Cramer's V", {
  M <- matrix(c(10, 20, 30, 40), nrow = 2)
  p <- mellio_payload(chisq.test(M))
  expect_equal(p$fields$sample_size, 100)
  expect_equal(p$fields$effect$name, "cramers_v")
  expect_true(is.numeric(p$fields$effect$value))
  # V is bounded in [0, 1] for 2x2; check sanity
  expect_true(p$fields$effect$value >= 0 && p$fields$effect$value <= 1)
})

test_that("chi-square goodness-of-fit (1D) omits Cramer's V", {
  # Single-vector chisq has no contingency table dim → V undefined
  p <- mellio_payload(chisq.test(c(5, 8, 9, 10, 12)))
  expect_null(p$fields$effect)
})

test_that("anova on a single lm emits partial eta-squared", {
  m <- lm(mpg ~ wt, data = mtcars)
  p <- mellio_payload(anova(m))
  expect_equal(p$fields$effect$name, "eta_sq_partial")
  expect_equal(p$fields$outcome, "mpg")
  expect_true(is.numeric(p$fields$effect$value))
  # Computed as sum_sq / (sum_sq + residual_sum_sq) — bounded [0, 1]
  expect_true(p$fields$effect$value > 0 && p$fields$effect$value < 1)
  expect_equal(
    p$fields$effect$value,
    p$fields$sum_sq / (p$fields$sum_sq + p$fields$residual_sum_sq)
  )
})

test_that("anova nested lm calls surface formula variables when available", {
  m <- lm(mpg ~ wt, data = mtcars)
  p <- mellio_payload(anova(m), .call = "anova(lm(mpg ~ wt, data = mtcars))")
  expect_equal(p$fields$outcome, "mpg")
  expect_equal(p$fields$predictor, "wt")
})
