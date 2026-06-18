# SPSS header simplification and APA renaming
# Ported from apa-generator.js lines 254-411 and parser.js lines 765-801

# ── SPSS verbose header -> short APA symbol mapping ────

.spss_header_map <- c(
  # Significance (all variants)
  "sig." = "p", "sig" = "p", "significance" = "p",
  "sig. (2-tailed)" = "p", "sig. (1-tailed)" = "p",
  "asymp. sig." = "p", "asymp. sig. (2-tailed)" = "p",
  "asymp. sig. (2-sided)" = "p",
  "asymptotic significance (2-sided)" = "p",
  "exact sig." = "p",
  "exact sig. (2-tailed)" = "p", "exact sig. (1-tailed)" = "p",
  "exact sig. (2-sided)" = "p", "exact sig. (1-sided)" = "p",
  "sig. f change" = "p", "point probability" = "p",

  # Regression (compound headers)
  "unstandardized coefficients b" = "B",
  "unstandardized coefficients std. error" = "SE",
  "standardized coefficients beta" = "\u03B2",

  # ANOVA
  "type iii sum of squares" = "SS",
  "sum of squares" = "SS", "mean square" = "MS",

  # Model summary
  "r square" = "R\u00B2", "adjusted r square" = "Adj. R\u00B2",
  "r square change" = "\u0394R\u00B2", "f change" = "\u0394F",
  "std. error of the estimate" = "SE",

  # Descriptives
  "std. deviation" = "SD", "std. error" = "SE",
  "std. error mean" = "SEM",
  "mean difference" = "MD", "std. error difference" = "SED",
  "std. error of skewness" = "SE(Skew.)",
  "std. error of kurtosis" = "SE(Kurt.)",

  # T-test (compound headers)
  "levene's test for equality of variances f" = "F",
  "levene's test for equality of variances sig." = "p",
  "t-test for equality of means t" = "t",
  "t-test for equality of means df" = "df",
  "t-test for equality of means sig. (2-tailed)" = "p",
  "t-test for equality of means mean difference" = "MD",
  "t-test for equality of means std. error difference" = "SED",

  # Correlations
  "pearson correlation" = "r",
  "correlation coefficient" = "r",
  "spearman's rho" = "\u03C1",
  "kendall's tau_b" = "\u03C4b",

  # Chi-Square
  "approx. chi-square" = "\u03C7\u00B2",
  "approx. sig." = "p",
  "minimum expected count" = "Exp. Count",

  # Reliability
  "cronbach's alpha" = "\u03B1",
  "cronbach's alpha if item deleted" = "\u03B1 if Deleted",
  "cronbach's alpha based on standardized items" = "\u03B1 (std.)",
  "corrected item-total correlation" = "r(IT)",
  "squared multiple correlation" = "R\u00B2",
  "scale mean if item deleted" = "M if Deleted",
  "scale variance if item deleted" = "Var. if Deleted",
  "n of items" = "k",

  # Factor analysis (compound headers)
  "initial eigenvalues total" = "Eigenvalue",
  "initial eigenvalues % of variance" = "% Var.",
  "initial eigenvalues cumulative %" = "Cum. %",
  "extraction sums of squared loadings total" = "Eigenvalue",
  "extraction sums of squared loadings % of variance" = "% Var.",
  "extraction sums of squared loadings cumulative %" = "Cum. %",
  "rotation sums of squared loadings total" = "Eigenvalue",
  "rotation sums of squared loadings % of variance" = "% Var.",
  "rotation sums of squared loadings cumulative %" = "Cum. %",

  # Logistic regression (compound headers)
  "variables in the equation b" = "B",
  "variables in the equation s.e." = "SE",
  "variables in the equation wald" = "Wald",
  "variables in the equation df" = "df",
  "variables in the equation sig." = "p",
  "variables in the equation exp(b)" = "OR",
  "95% c.i.for exp(b) lower" = "LL",
  "95% c.i.for exp(b) upper" = "UL",

  # Effect sizes
  "partial eta squared" = "\u03B7p\u00B2",
  "eta squared" = "\u03B7\u00B2",
  "omega squared" = "\u03C9\u00B2",
  "cohen's d" = "d",

  # MANOVA / Multivariate
  "pillai's trace" = "V",
  "wilks' lambda" = "\u039B",
  "hotelling's trace" = "T\u00B2",
  "roy's largest root" = "\u03B8",
  "hypothesis df" = "df\u2081", "error df" = "df\u2082",

  # Nonparametric
  "mann-whitney u" = "U", "wilcoxon w" = "W",
  "kruskal-wallis h" = "H",
  "mauchly's w" = "W",

  # CI compound headers
  "95.0% confidence interval for b lower bound" = "LL",
  "95.0% confidence interval for b upper bound" = "UL",
  "95% confidence interval lower bound" = "LL",
  "95% confidence interval upper bound" = "UL",
  "95% confidence interval of the difference lower" = "LL",
  "95% confidence interval of the difference upper" = "UL",

  # Collinearity
  "collinearity statistics tolerance" = "Tolerance",
  "collinearity statistics vif" = "VIF",

  # Suffix-level mappings
  "beta" = "\u03B2", "s.e." = "SE",
  "lower bound" = "LL", "upper bound" = "UL",
  "lower" = "LL", "upper" = "UL",
  "total" = "Eigenvalue",
  "% of variance" = "% Var.", "cumulative %" = "Cum. %",
  "exp(b)" = "OR"
)

# SPSS compound prefixes (for prefix stripping)
.spss_compound_prefixes <- tolower(c(
  "Unstandardized Coefficients", "Standardized Coefficients",
  "Collinearity Statistics", "Bootstrap for",
  "Levene's Test for Equality of Variances",
  "t-test for Equality of Means",
  "Variables in the Equation",
  "95% C.I.for EXP(B)",
  "95.0% Confidence Interval for B", "95% Confidence Interval for B",
  "95.0% Confidence Interval of the Difference",
  "95% Confidence Interval of the Difference",
  "95.0% Confidence Interval", "95% Confidence Interval",
  "Initial Eigenvalues", "Extraction Sums of Squared Loadings",
  "Rotation Sums of Squared Loadings",
  "Correlations"
))

# R-specific column name -> APA abbreviation
.apa_header_map <- c(
  "Pr(>|t|)" = "p", "Pr(>|z|)" = "p", "Pr(>F)" = "p", "Pr(>Chi)" = "p",
  "Pr(>|F|)" = "p", "Pr..t.." = "p", "Pr..z.." = "p", "Pr..F.." = "p",
  "p-value" = "p", "p value" = "p", "P Value" = "p", "P-Value" = "p",
  "p.value" = "p", "pvalue" = "p", "P.Value" = "p",
  "p.adjust" = "p (adjusted)", "p.adj" = "p (adjusted)",
  "p_adj" = "p (adjusted)", "padj" = "p (adjusted)",
  "Sig." = "p", "Sig" = "p",
  "t value" = "t", "z value" = "z", "F value" = "F",
  "t.ratio" = "t", "z.ratio" = "z", "F.ratio" = "F",
  "Std. Error" = "SE", "Std.Error" = "SE", "Std Error" = "SE",
  "std.error" = "SE", "Std.Err" = "SE",
  "Estimate" = "B", "Coefficient" = "B",
  "Sum Sq" = "SS", "Sum of Sq" = "SS", "Mean Sq" = "MS",
  "RSS" = "RSS", "Res.Df" = "Residual df", "Res. Df" = "Residual df",
  "Num DF" = "df1", "Den DF" = "df2",
  "Resid. Dev" = "Residual Deviance", "Resid. Df" = "Residual df",
  "Odds Ratio" = "OR", "Hazard Ratio" = "HR",
  "LR Chisq" = "\u03C7\u00B2", "LR.Chisq" = "\u03C7\u00B2",
  "Chi.Sq" = "\u03C7\u00B2", "Chi Sq" = "\u03C7\u00B2",
  "lower.CL" = "Lower CI", "upper.CL" = "Upper CI",
  "Lower.CL" = "Lower CI", "Upper.CL" = "Upper CI",
  "lower.CI" = "Lower CI", "upper.CI" = "Upper CI",
  "Lower.CI" = "Lower CI", "Upper.CI" = "Upper CI",
  "conf.low" = "Lower CI", "conf.high" = "Upper CI",
  "lwr" = "Lower CI", "upr" = "Upper CI",
  "asymp.LCL" = "Lower CI", "asymp.UCL" = "Upper CI",
  "emmean" = "M", "marginal.mean" = "M",
  "statistic" = "t"
)

#' Simplify an SPSS verbose header to APA abbreviation
#'
#' @param header Character string
#' @return Simplified header string
#' @keywords internal
simplify_spss_header <- function(header) {
  trimmed <- trimws(header)
  if (!nzchar(trimmed)) return(header)

  # Direct lookup
  mapped <- .spss_header_map[tolower(trimmed)]
  if (!is.na(mapped)) return(unname(mapped))

  # Compound prefix stripping
  lower <- tolower(trimmed)
  for (prefix in .spss_compound_prefixes) {
    if (nchar(lower) > nchar(prefix) && startsWith(lower, paste0(prefix, " "))) {
      suffix <- trimws(substring(trimmed, nchar(prefix) + 2))
      suffix_mapped <- .spss_header_map[tolower(suffix)]
      if (!is.na(suffix_mapped)) return(unname(suffix_mapped))
      return(suffix)
    }
  }

  header
}

#' Rename R/broom column names to APA standard abbreviations
#'
#' @param headers Character vector of column names
#' @return Character vector with renamed headers
#' @keywords internal
rename_to_apa <- function(headers) {
  vapply(headers, function(h) {
    trimmed <- trimws(h)
    if (!nzchar(trimmed)) return(h)
    if (trimmed %in% names(.apa_header_map)) {
      return(unname(.apa_header_map[trimmed]))
    }
    h
  }, character(1), USE.NAMES = FALSE)
}
