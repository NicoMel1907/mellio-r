# mellio

Create publication-ready statistical tables in R and send supported results,
tables, and figures to the Mellio web app.

Mellio for R has two main workflows:

- `melliotab()` creates formatted tables locally in R.
- `mellio_open()` opens supported R objects in Mellio for editing, organizing,
  and reporting.

Review statistical output before publication, especially for complex models or
objects created by optional packages.

## Installation

```r
remotes::install_github("NicoMel1907/mellio-r")
```

## Quick Start

```r
library(mellio)

model <- lm(mpg ~ wt + hp, data = mtcars)

# Create a table in R
tab <- melliotab(model, title = "Predictors of Fuel Efficiency")
tab

# Open the model result or finished table in Mellio
mellio_open(model)
mellio_open(tab)

# Compare models side by side
model_1 <- lm(mpg ~ wt, data = mtcars)
model_2 <- lm(mpg ~ wt + hp + disp, data = mtcars)

melliotab(model_1, model_2, labels = c("Step 1", "Step 2"))
mellio_open(model_1, model_2, labels = c("Step 1", "Step 2"))
```

`mellio_open()` also accepts tabular data and supported plot/image inputs:

```r
mellio_open(mtcars[1:10, 1:4])

library(ggplot2)
p <- ggplot(mtcars, aes(wt, mpg)) + geom_point()
mellio_open(p, title = "Weight vs. Fuel Efficiency")
mellio_open("figure.png", title = "Experimental Setup")
```

## Tables In R

`melliotab()` formats common statistical outputs for manuscript-style tables.
It detects p-values, estimates, confidence intervals, test statistics, integer
counts, and bounded statistics such as correlations.

| Object | Example | Table output |
|---|---|---|
| `data.frame` | `melliotab(df)` | Formatted data table |
| `lm` | `melliotab(model)` | Coefficients, confidence intervals, model note |
| `glm` | `melliotab(model)` | Coefficients, with optional odds ratios |
| `aov` | `melliotab(model)` | ANOVA table with effect sizes |
| `htest` | `melliotab(t.test(...))` | Test statistic, df, p-value, confidence interval |
| `matrix` | `melliotab(cor_matrix)` | Correlation matrix |
| `lavaan` | `melliotab(fit, section = "loadings")` | CFA/SEM fit, loadings, paths, effects, reliability |
| `psych::fa` | `melliotab(efa_fit)` | EFA loadings, variance, fit indices |
| multiple models | `melliotab(m1, m2)` | Side-by-side model comparison |

## Choosing Table Sections

Some statistical objects can produce more than one useful table. When there is
no single safe default, Mellio asks you to choose a section.

For SEM/CFA models:

```r
melliotab(sem_fit, section = "loadings")
melliotab(sem_fit, section = "fit")
melliotab(sem_fit, section = "paths")
melliotab(sem_fit, section = "reliability")
```

EFA defaults to factor loadings, with optional alternatives:

```r
melliotab(efa_fit)
melliotab(efa_fit, what = "variance")
melliotab(efa_fit, what = "fit")
```

## Common Table Options

Core table settings are regular `melliotab()` arguments:

```r
model <- lm(mpg ~ wt + hp + factor(cyl), data = mtcars)

melliotab(
  model,
  style = "ieee",
  title = "Predictors of Fuel Efficiency",
  number = 1,
  note = "Estimates are unstandardized regression coefficients.",
  decimals = 2,
  p_decimals = 3
)
```

Table modifiers can be added with the base R pipe:

```r
melliotab(
  model,
  style = "ieee",
  title = "Predictors of Fuel Efficiency",
  number = 1,
  note = "Estimates are unstandardized regression coefficients.",
  decimals = 2,
  p_decimals = 3
) |>
  mt_sig_stars(remove_p = FALSE) |>
  mt_spanner("95% CI", columns = c("Lower CI", "Upper CI"))
```

Common customization tasks:

| Task | Use |
|---|---|
| Change citation style | `style = "apa7"` or `style = "ieee"` |
| Set title, number, or note | `title =`, `number =`, `note =` |
| Control rounding | `decimals = 2`, `p_decimals = 3` |
| Add significance stars | `mt_sig_stars(remove_p = FALSE)` |
| Group related columns | `mt_spanner("95% CI", columns = c("Lower CI", "Upper CI"))` |

Significance stars are never added by default. Use `mt_sig_stars()` only when
that convention is appropriate for your manuscript, course, or journal.

## Correlation Tables

Correlation matrices can be shown as full matrices or as lower/upper triangles:

```r
cor_tab <- melliotab(cor(mtcars[, c("mpg", "wt", "hp")]))

cor_tab |>
  mt_diagonal(mode = "dash", triangle = "lower")
```

`mode` controls the diagonal cells: `"dash"`, `"one"`, or `"blank"`.
`triangle` controls which half of the matrix is shown: `"all"`, `"lower"`, or
`"upper"`.

## Advanced Layout

Use these helpers when you need more control over a table's structure:

```r
melliotab(model, title = "Predictors of Fuel Efficiency") |>
  mt_section_title("Cylinder terms", before = 4) |>
  mt_indent(rows = 4:5, level = 1)
```

`before` inserts the section title before a row number. `after` can be used
instead when it is more natural to place the section title after a row.
`level` controls indentation depth.

## Copy Or Save Tables

The default RStudio workflow is to preview the table and copy it for writing
tools. When you need a file, use the manual output helpers:

```r
mt_copy(tab)
mt_save(tab, "regression-table.html")
mt_save(tab, "regression-table.tex")
mt_save(tab, "regression-table.md")
```

`mt_copy()` uses the system clipboard on supported desktop platforms. On other
systems, use `mt_save()`.

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

For a richer Mellio Stats card with model-level R-squared, adjusted R-squared,
delta R-squared, and F-change:

```r
mellio_open(mellio_compare(model_1, model_2, labels = c("Step 1", "Step 2")))
```

## Mellio Web Handoff

Use `mellio_open()` when you explicitly want to open an object in the Mellio web
app:

```r
mellio_open(t.test(score ~ group, data = my_data))
mellio_open(lm(mpg ~ wt + cyl, data = mtcars))
mellio_open(model_1, model_2, labels = c("Step 1", "Step 2"))
mellio_open(melliotab(model, title = "My Table", number = 1))
mellio_open(my_data)
mellio_open(p, title = "My Plot")
```

By default `mellio_open()` opens `https://www.mellioapp.com`. Advanced users can
override the destination with `options("mellio.editor_url")`, but should only
point it at a trusted Mellio deployment.

```r
options(mellio.editor_url = "https://www.mellioapp.com")
```

The handoff data is encoded in the URL fragment. URL fragments are not sent as
HTTP requests to the server, but the full URL can still be visible to the
browser, the opened web app, browser history, extensions, and anyone the URL is
shared with.

Mellio includes R/package-version metadata and data fingerprints where
available. Local machine details such as user name, host name, working
directory, git state, and script path are opt-in:

```r
options(mellio.provenance = "full")
```

To omit provenance metadata:

```r
options(mellio.provenance = FALSE)
```

Use `citation("mellio")` for the package citation.
