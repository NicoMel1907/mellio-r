# mellio

R tools for sending statistical results, model comparisons, tables, and figures
to the Mellio web app and creating publication-ready tables in R.

Public beta: Mellio for R supports selected models, tests, model comparisons,
tables, and figures while the workflow is validated with real users. Please
review outputs before publication.

## Installation

```r
remotes::install_github("NicoMel1907/mellio-r", upgrade = "never")
```

## Quick Start

```r
library(mellio)

# Send statistical results to Mellio
mellio_open(t.test(extra ~ group, data = sleep))

model <- lm(mpg ~ wt + hp, data = mtcars)
mellio_open(model)

# Create a publication-ready table in R
tab <- melliotab(model, title = "Predictors of Fuel Efficiency")
tab

# Compare models side by side
model_1 <- lm(mpg ~ wt, data = mtcars)
model_2 <- lm(mpg ~ wt + hp + disp, data = mtcars)
melliotab(model_1, model_2, labels = c("Step 1", "Step 2"))

# Open a table in Mellio
mellio_open(tab)
```

`mellio_open()` also accepts tabular data and supported plot/image inputs:

```r
# Tables workspace
mellio_open(mtcars[1:10, 1:4])

# Figures workspace
library(ggplot2)
p <- ggplot(mtcars, aes(wt, mpg)) + geom_point()
mellio_open(p, title = "Weight vs. Fuel Efficiency")
mellio_open("figure.png", title = "Experimental Setup")
```

## Tables

### What It Does

- Detects column types such as p-values, estimates, test statistics, and integers
- Formats p-values as `< .001` when appropriate
- Removes leading zeros for bounded statistics such as r, p, and beta in APA style
- Rounds numeric columns to your decimal setting
- Italicizes statistical symbols in headers such as B, t, p, F, and df
- Generates table notes with model fit information where available
- Supports APA 7th and IEEE table styling

### Supported Inputs

| Object | Function | What it produces |
|---|---|---|
| `data.frame` | `melliotab(df)` | Formatted table with auto-detected columns |
| `lm` | `melliotab(model)` | Coefficients, confidence intervals, model note |
| `glm` | `melliotab(model)` | Coefficients with optional odds ratios |
| `aov` | `melliotab(model)` | ANOVA table with effect sizes |
| `htest` | `melliotab(t.test(...))` | Test statistic, df, p-value, confidence interval |
| `matrix` | `melliotab(cor_matrix)` | Correlation matrix |
| `lavaan` | `melliotab(fit)` | CFA/SEM loadings, paths, fit indices |
| `lm`, `glm`, ... | `melliotab(m1, m2, ...)` | Side-by-side model comparison |

### Create

```r
melliotab(x, style = "apa7", title = NULL, number = NULL,
          note = NULL, decimals = 2, p_decimals = 3)
```

### Modify

All modifiers return the modified object, so you can chain them with `|>`:

```r
melliotab(model, title = "Regression Results") |>
  mt_sig_stars() |>
  mt_note("Robustness checks are reported in the supplement.")
```

The same modifiers also work with ordinary assignment:

```r
model <- lm(mpg ~ wt + hp + factor(cyl), data = mtcars)

tab <- melliotab(
  model,
  title = "Predictors of Fuel Efficiency",
  number = 1,
  decimals = 2,
  p_decimals = 3
)
tab <- mt_sig_stars(tab)
tab <- mt_note(
  tab,
  "Note. Estimates are unstandardized regression coefficients. Cylinder group is included as a categorical predictor."
)
tab
```

| Function | Description |
|---|---|
| `mt_title(x, title)` | Set table title |
| `mt_number(x, number)` | Set table number |
| `mt_note(x, note)` | Set table note |
| `mt_source(x, source)` | Set source text |
| `mt_set_style(x, style)` | Change citation style |
| `mt_decimals(x, decimals, p_decimals)` | Set decimal places |
| `mt_sig_stars(x)` | Add significance stars |
| `mt_spanner(x, label, columns)` | Group columns under a header |
| `mt_indent(x, rows, level)` | Indent row labels |
| `mt_section_title(x, label, before)` | Insert a section title row |
| `mt_diagonal(x, mode, triangle)` | Control correlation matrix display |
| `mt_format_ci(x)` | Merge CI columns into `[low, high]` |
| `mt_remove_leading_zeros(x)` | Toggle leading zero removal |
| `mt_simplify_headers(x)` | Shorten verbose SPSS-style headers |

## Model Comparison

```r
model_1 <- lm(mpg ~ wt, data = mtcars)
model_2 <- lm(mpg ~ wt + hp + disp, data = mtcars)

melliotab(
  model_1, model_2,
  title = "Hierarchical Regression: Fuel Efficiency",
  labels = c("Step 1", "Step 2"),
  dep.var.labels = "Miles per Gallon"
)
```

`mt_compare(model_1, model_2, ...)` remains available as the explicit
table-comparison helper. For a richer Mellio Stats card with model-level
R-squared, adjusted R-squared, delta R-squared, and F-change, use:

```r
mellio_open(mellio_compare(model_1, model_2, labels = c("Step 1", "Step 2")))
```

## Mellio Web Handoff

Use `mellio_open()` to send supported R objects to the Mellio web app:

```r
mellio_open(t.test(score ~ group, data = my_data))
mellio_open(lm(mpg ~ wt + cyl, data = mtcars))
mellio_open(model_1, model_2, labels = c("Step 1", "Step 2"))
mellio_open(melliotab(model, title = "My Table", number = 1))
mellio_open(my_data)
mellio_open(p, title = "My Plot")
```

By default `mellio_open()` opens `https://www.mellioapp.com`. Advanced
users can override the destination with `options("mellio.editor_url")`, but
should only point it at a trusted Mellio deployment.

```r
options(mellio.editor_url = "https://www.mellioapp.com")
```

The R payload is encoded in the URL fragment. URL fragments are not sent as
HTTP requests to the server, but the full URL can still be visible to the
browser, the opened web app, browser history, extensions, and anyone the URL
is shared with.

Mellio payloads include R/package-version metadata and data fingerprints
where available. Local machine details such as user name, host name,
working directory, git state, and script path are opt-in:

```r
options(mellio.provenance = "full")
```

To omit provenance metadata:

```r
options(mellio.provenance = FALSE)
```

Use `citation("mellio")` for the package citation.
