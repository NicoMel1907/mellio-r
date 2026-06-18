# Column type auto-detection
# Ported from apa-generator.js lines 81-251

# ── Header pattern sets ────────────────────────────────

.p_value_headers <- tolower(c(

  "p", "Pr(>|t|)", "Pr(>|z|)", "Pr(>F)", "Pr(>Chi)",
  "p-value", "p value", "p.value", "pvalue",
  "p (adjusted)", "p (raw)", "p adj", "p.adj", "p_adj",
  "padj", "p.adjust", "p.adjusted", "adjusted p", "adj p",
  "Sig.", "Sig", "P", "P Value", "P.Value",
  "Pr(>|F|)", "Pr..t..", "Pr..z..", "Pr..F..",
  "Significance", "Sig. (2-tailed)", "Sig. (1-tailed)",
  "Asymp. Sig.", "Asymp. Sig. (2-tailed)", "Asymp. Sig. (2-sided)",
  "Asymptotic Significance (2-sided)",
  "Exact Sig.", "Exact Sig. (2-tailed)", "Exact Sig. (1-tailed)",
  "Exact Sig. (2-sided)", "Exact Sig. (1-sided)",
  "Sig. F Change", "Point Probability", "Approx. Sig.",
  "Levene's Test for Equality of Variances Sig.",
  "t-test for Equality of Means Sig. (2-tailed)",
  "Variables in the Equation Sig.",
  "Pr(>Chisq)"
))

.estimate_headers <- tolower(c(
  "Estimate", "Coefficient", "B", "b", "Beta", "\u03B2",
  "Std. Error", "Std.Error", "SE", "SEM", "std.error", "Std.Err",
  "M", "SD", "Mean", "Std. Dev",
  "Lower CI", "Upper CI", "CI Lower", "CI Upper",
  "lower.CL", "upper.CL", "Lower.CL", "Upper.CL",
  "lower.CI", "upper.CI", "Lower.CI", "Upper.CI",
  "asymp.LCL", "asymp.UCL", "lwr", "upr",
  "conf.low", "conf.high", "2.5 %", "97.5 %", "2.5%", "97.5%",
  "Odds Ratio", "Hazard Ratio", "OR", "HR",
  "r", "R", "R2", "R\u00B2",
  "Adjusted R2", "Adjusted R\u00B2", "Adj. R2", "Adj. R\u00B2",
  "Delta R2", "Delta R\u00B2", "\u0394R2", "\u0394R\u00B2",
  "Cohen's d", "d", "\u03B7\u00B2", "\u03B7p\u00B2",
  "emmean", "response", "estimate", "prob", "rate",
  "cn_mean", "marginal.mean",
  "Sum of Squares", "Mean Square", "Mean Squares",
  "Type III Sum of Squares", "SS", "MS",
  "Unstandardized Coefficients B", "Unstandardized Coefficients Std. Error",
  "Standardized Coefficients Beta",
  "Std. Deviation", "Std. Error Mean", "Std. Error of the Estimate",
  "Mean Difference", "Std. Error Difference",
  "R Square", "Adjusted R Square", "R Square Change",
  "Pearson Correlation", "Correlation Coefficient",
  "Cronbach's Alpha", "Cronbach's Alpha if Item Deleted",
  "Corrected Item-Total Correlation", "Squared Multiple Correlation",
  "Partial Eta Squared", "Eta Squared", "Omega Squared",
  "Pillai's Trace", "Wilks' Lambda", "Hotelling's Trace",
  "Roy's Largest Root",
  "Variables in the Equation B", "Variables in the Equation Exp(B)",
  "95.0% Confidence Interval for B Lower Bound",
  "95.0% Confidence Interval for B Upper Bound",
  "95% Confidence Interval Lower Bound", "95% Confidence Interval Upper Bound",
  "95% Confidence Interval of the Difference Lower",
  "95% Confidence Interval of the Difference Upper",
  "95% C.I.for EXP(B) Lower", "95% C.I.for EXP(B) Upper",
  "t-test for Equality of Means Mean Difference",
  "t-test for Equality of Means Std. Error Difference",
  "Initial Eigenvalues Total",
  "Extraction Sums of Squared Loadings Total",
  "Rotation Sums of Squared Loadings Total",
  "est", "std.all", "std.lv", "est.std"
))

.t_stat_headers <- tolower(c(
  "t", "t value", "z", "z value", "F", "F value",
  "t.ratio", "z.ratio", "F.ratio",
  "Wald", "Chi.Sq", "Chi Sq", "LR Chisq", "LR.Chisq", "Chisq", "Chisq diff",
  "W", "U", "H", "statistic",
  "Levene's Test for Equality of Variances F",
  "t-test for Equality of Means t",
  "Variables in the Equation Wald",
  "F Change", "\u0394F",
  "Approx. Chi-Square",
  "Mann-Whitney U", "Wilcoxon W", "Kruskal-Wallis H",
  "Mauchly's W"
))

.keep_integer_headers <- tolower(c(
  "n", "N", "df", "df1", "df2", "k", "Count", "Freq", "Frequency",
  "Df", "DF", "nobs", "n.obs", "Sample", "Sample Size",
  "N of Items",
  "Missing", "Missing n", "N missing", "NA's",
  "Variables in the Equation df", "Residual df",
  "t-test for Equality of Means df",
  "Df diff"
))

.text_label_headers <- tolower(c(
  ".y.", "y", "outcome", "variable", "term", "predictor",
  "group", "group1", "group2", "condition", "level",
  "contrast", "comparison", "method", "p.adj.signif",
  "significance"
))

.leading_zero_headers <- tolower(c(
  "p", "p.value", "pvalue", "p-value", "p value", "P", "P Value", "P.Value",
  "p (adjusted)", "p (raw)", "p adj", "p.adj", "p_adj",
  "padj", "p.adjust", "p.adjusted",
  "Sig.", "Sig",
  "Pr(>|t|)", "Pr(>|z|)", "Pr(>F)", "Pr(>Chi)", "Pr(>|F|)",
  "Pr..t..", "Pr..z..", "Pr..F..",
  "r", "R", "R2", "R\u00B2",
  "Adjusted R2", "Adjusted R\u00B2", "Adj. R2", "Adj. R\u00B2",
  "Delta R2", "Delta R\u00B2", "\u0394R2", "\u0394R\u00B2",
  "\u03B2", "Beta", "beta",
  "\u03B7\u00B2", "\u03B7p\u00B2",
  "Cohen's d", "d",
  "Standardized Coefficients Beta",
  "Pearson Correlation", "Correlation Coefficient",
  "Corrected Item-Total Correlation",
  "Partial Eta Squared", "Eta Squared", "Omega Squared",
  "R Square", "Adjusted R Square", "R Square Change",
  "Cronbach's Alpha", "Cronbach's Alpha if Item Deleted",
  "Squared Multiple Correlation"
))

.keep_leading_zero_headers <- tolower(c(
  "SE", "SEM", "Std. Error", "Std.Error", "std.error", "Std.Err",
  "M", "Mean", "SD", "Std. Dev",
  "t", "t value", "t.ratio",
  "z", "z value", "z.ratio",
  "F", "F value", "F.ratio",
  "B", "b", "Estimate", "Coefficient",
  "df", "df1", "df2", "N", "n",
  "Missing", "Missing n", "N missing", "NA's",
  "Lower CI", "Upper CI", "CI Lower", "CI Upper",
  "lower.CL", "upper.CL", "Lower.CL", "Upper.CL",
  "lower.CI", "upper.CI", "Lower.CI", "Upper.CI",
  "conf.low", "conf.high", "lwr", "upr",
  "asymp.LCL", "asymp.UCL",
  "2.5 %", "97.5 %", "2.5%", "97.5%",
  "W", "U", "H", "Wald", "Chi.Sq",
  "SS", "MS", "OR", "HR",
  "emmean", "response", "estimate", "prob", "rate",
  "Unstandardized Coefficients B", "Unstandardized Coefficients Std. Error",
  "Std. Deviation", "Std. Error Mean", "Std. Error of the Estimate",
  "Mean Difference", "Std. Error Difference",
  "Sum of Squares", "Mean Square", "Type III Sum of Squares",
  "F Change",
  "Levene's Test for Equality of Variances F",
  "t-test for Equality of Means t", "t-test for Equality of Means df",
  "t-test for Equality of Means Mean Difference",
  "t-test for Equality of Means Std. Error Difference",
  "Variables in the Equation B", "Variables in the Equation S.E.",
  "Variables in the Equation Wald", "Variables in the Equation df",
  "Variables in the Equation Exp(B)",
  "Approx. Chi-Square",
  "Mann-Whitney U", "Wilcoxon W", "Kruskal-Wallis H", "Mauchly's W",
  "Hypothesis df", "Error df",
  "Initial Eigenvalues Total",
  "Extraction Sums of Squared Loadings Total",
  "Rotation Sums of Squared Loadings Total",
  "Hotelling's Trace",
  "95.0% Confidence Interval for B Lower Bound",
  "95.0% Confidence Interval for B Upper Bound",
  "95% Confidence Interval Lower Bound", "95% Confidence Interval Upper Bound",
  "95% Confidence Interval of the Difference Lower",
  "95% Confidence Interval of the Difference Upper",
  "95% C.I.for EXP(B) Lower", "95% C.I.for EXP(B) Upper",
  "Scale Mean if Item Deleted", "Scale Variance if Item Deleted",
  "N of Items", "Std. Error of Skewness", "Std. Error of Kurtosis",
  "est", "std.all", "std.lv", "est.std",
  "Chisq", "Chisq diff", "Df diff"
))

# Headers that should be italicized in APA style
.apa_italic_headers <- c(
  "n", "N", "M", "B", "SE", "t", "z", "F", "p", "r", "R",
  "df", "df1", "df2", "k", "d", "g",
  "Df", "W", "U", "H", "V"
)

# Primary estimate headers (for sig star target detection)
.primary_estimate_headers <- tolower(c(
  "Estimate", "Coefficient", "B", "b", "Beta", "\u03B2",
  "estimate",
  "r", "R", "R2", "R\u00B2",
  "Odds Ratio", "OR", "HR", "Hazard Ratio",
  "Cohen's d", "d", "\u03B7\u00B2", "\u03B7p\u00B2",
  "M", "Mean", "emmean", "response", "prob", "rate",
  "Unstandardized Coefficients B",
  "Standardized Coefficients Beta",
  "Variables in the Equation B",
  "Pearson Correlation", "Correlation Coefficient"
))

# ── Detection functions ────────────────────────────────

#' Match a header against a pattern set (case-insensitive)
#' @keywords internal
match_header <- function(header, pattern_set) {
  tolower(trimws(header)) %in% pattern_set
}

#' Normalize a header for compact pattern checks
#' @keywords internal
compact_header <- function(header) {
  h <- tolower(gsub("\\s+", " ", trimws(as.character(header))))
  gsub("[._\\-\\s]", "", h, perl = TRUE)
}

#' Check whether a header denotes a whole-number count column
#' @keywords internal
is_count_header <- function(header) {
  if (match_header(header, .keep_integer_headers)) return(TRUE)
  h <- compact_header(header)
  grepl("^n\\d+$", h, perl = TRUE) ||
    grepl("^npairs?$", h, perl = TRUE) ||
    grepl("^pairs?n$", h, perl = TRUE)
}

#' Check whether a header denotes a label column
#' @keywords internal
is_text_label_header <- function(header) {
  if (match_header(header, .text_label_headers)) return(TRUE)
  h <- compact_header(header)
  grepl("^(group|level|condition)\\d+$", h, perl = TRUE)
}

#' Check whether a header denotes a p-value column
#' @keywords internal
is_p_value_header <- function(header) {
  h <- tolower(gsub("\\s+", " ", trimws(as.character(header))))
  if (!nzchar(h)) return(FALSE)
  if (h %in% .p_value_headers) return(TRUE)
  if (grepl("^p\\s*\\([^)]+\\)$", h, perl = TRUE)) return(TRUE)
  if (grepl("^p[._\\-\\s]*(adj\\.?|adjust(ed)?|raw|unadjusted|uncorrected|corrected|fdr|holm|bonferroni|tukey|sidak|scheffe|dunnett|mvt)$", h, perl = TRUE)) {
    return(TRUE)
  }
  if (grepl("^(adj\\.?|adjusted|raw|unadjusted|uncorrected|corrected|fdr|holm|bonferroni|tukey|sidak|scheffe|dunnett|mvt)[._\\-\\s]*p([._\\-\\s]*(value|val))?$", h, perl = TRUE)) {
    return(TRUE)
  }
  FALSE
}

#' Detect column types from header names
#'
#' @param headers Character vector of column names
#' @return Character vector: "stub", "pvalue", "estimate", "statistic", "integer", or "default"
#' @keywords internal
detect_column_types <- function(headers) {
  vapply(seq_along(headers), function(j) {
    h <- headers[j]
    if (j == 1 && !is_p_value_header(h) &&
        !is_text_label_header(h) &&
        !match_header(h, .estimate_headers) &&
        !match_header(h, .t_stat_headers) &&
        !is_count_header(h)) {
      return("stub")
    }
    if (is_p_value_header(h)) return("pvalue")
    if (is_text_label_header(h)) return("stub")
    if (match_header(h, .estimate_headers)) return("estimate")
    if (match_header(h, .t_stat_headers)) return("statistic")
    if (is_count_header(h)) return("integer")
    "default"
  }, character(1))
}

#' Detect which columns should have leading zeros removed
#'
#' @param headers Character vector of column names
#' @return Logical vector
#' @keywords internal
detect_leading_zero_cols <- function(headers) {
  vapply(headers, function(h) {
    lh <- tolower(trimws(h))
    if (is_count_header(lh)) return(FALSE)
    if (lh %in% .keep_leading_zero_headers) return(FALSE)
    if (lh %in% .leading_zero_headers || is_p_value_header(lh)) return(TRUE)
    FALSE
  }, logical(1), USE.NAMES = FALSE)
}

#' Detect p-value column index
#'
#' @param headers Character vector of column names
#' @return Integer index (1-based) or NULL
#' @keywords internal
detect_p_value_col <- function(headers) {
  for (j in seq_along(headers)) {
    if (is_p_value_header(headers[j])) return(j)
  }
  NULL
}

#' Detect the primary estimate column (target for significance stars)
#'
#' @param headers Character vector of column names
#' @param p_col Integer index of the p-value column
#' @return Integer index (1-based)
#' @keywords internal
detect_estimate_col <- function(headers, p_col) {
  # First pass: primary estimates
  for (j in seq_along(headers)) {
    if (j == p_col) next
    if (match_header(headers[j], .primary_estimate_headers)) return(j)
  }
  # Second pass: test statistics
  test_stat <- tolower(c(
    "F", "t", "z", "W", "U", "H", "\u03C7\u00B2",
    "t.ratio", "z.ratio", "F.ratio", "statistic"
  ))
  for (j in seq_along(headers)) {
    if (j == p_col) next
    if (tolower(trimws(headers[j])) %in% test_stat) return(j)
  }
  # Third pass: secondary estimates (SE, SD, etc.)
  secondary <- tolower(c("SS", "MS", "SE", "SEM", "SD",
                          "Std. Error", "Std.Error", "std.error", "Std.Err"))
  for (j in seq_along(headers)) {
    if (j == p_col) next
    if (tolower(trimws(headers[j])) %in% secondary) return(j)
  }
  # Fallback: column before p
  if (!is.null(p_col) && p_col > 1) return(p_col - 1L)
  1L
}

#' Detect CI columns
#'
#' @param headers Character vector of column names
#' @return Logical vector
#' @keywords internal
detect_ci_cols <- function(headers) {
  ci_patterns <- c("ci", ".cl", "confidence", "ll, ul", "[ll",
                    "95%", "lwr", "upr", "conf.low", "conf.high",
                    "lcl", "ucl")
  vapply(headers, function(h) {
    lh <- tolower(trimws(h))
    any(vapply(ci_patterns, function(p) grepl(p, lh, fixed = TRUE), logical(1)))
  }, logical(1), USE.NAMES = FALSE)
}

#' Check if a header is a statistical symbol that should be italicized
#'
#' @param header Column header string
#' @return Logical
#' @keywords internal
is_stat_symbol <- function(header) {
  trimmed <- trimws(header)
  if (!nzchar(trimmed)) return(FALSE)
  if (trimmed %in% .apa_italic_headers) return(TRUE)
  if (grepl("^n[\\s._\\-]*(\\d+|pairs?)$", trimmed, ignore.case = TRUE, perl = TRUE)) return(TRUE)
  if (is_p_value_header(trimmed) && grepl("^p($|[\\s._\\-(]|adj|adjust|raw|unadjusted|uncorrected|corrected|fdr|holm|bonferroni|tukey|sidak|scheffe|dunnett|mvt|value)", trimmed, ignore.case = TRUE, perl = TRUE)) return(TRUE)
  if (grepl("^\u0394?R[\u00B2\u00B3 2]?$", trimmed)) return(TRUE)
  FALSE
}
