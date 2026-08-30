# SEDesign: experiment design and contrasts object (S7 class)

`SEDesign` enforces the relationship between individual samples, design
groups, and groups involved in statistical contrasts.

## Usage

``` r
SEDesign(
  design = matrix(nrow = 0, ncol = 0),
  contrasts = matrix(nrow = 0, ncol = 0),
  samples = character(0),
  factors = character(0)
)
```

## Arguments

- design:

  `matrix` with rownames as samples, colnames as design groups,
  containing 0/1 sample-to-group association values.

- contrasts:

  `matrix` with rownames matching `colnames(design)`, and colnames as
  contrast names, containing contrast coefficients.

- samples:

  `character` vector of sample identifiers, typically equal to
  `rownames(design)`.

- factors:

  `character` vector of labels for the underlying experimental factors,
  equivalent to `colnames(design_df)`, i.e. the columns produced by
  splitting each design group name (`groups(object)`) on `"_"`. When
  empty, or when its length does not match the number of underlying
  factors, it is auto-populated with `"factor1"`, `"factor2"`, etc.
  Editing
  [`factors()`](https://jmw86069.github.io/jamses/reference/factors.md)
  only renames these labels (used primarily by
  [`print()`](https://rdrr.io/r/base/print.html)); it has no effect on
  `design`, `contrasts`, or
  [`groups()`](https://jmw86069.github.io/jamses/reference/groups.md).

## Details

In addition to the four constructor arguments below, `SEDesign` objects
carry two internal-use properties, accessible via `@`:

- `design_df`: `data.frame` with one row per design group
  (`groups(object)`, i.e. `colnames(design)`), derived by splitting each
  group name on `"_"` into one column per underlying experimental
  factor. Column names default to `"factor1"`, `"factor2"`, etc., and
  can be customized via `factors()<-`.

- `contrasts_df`: `data.frame`, cached result of
  [`contrasts_to_factors()`](https://jmw86069.github.io/jamses/reference/contrasts_to_factors.md),
  refreshed whenever `design` or `contrasts` are updated.

## See also

Other jam experiment design:
[`[,SEDesign-method`](https://jmw86069.github.io/jamses/reference/sub-SEDesign-method.md),
[`check_sedesign()`](https://jmw86069.github.io/jamses/reference/check_sedesign.md),
[`contrast2comp()`](https://jmw86069.github.io/jamses/reference/contrast2comp.md),
[`contrast_colors_by_group()`](https://jmw86069.github.io/jamses/reference/contrast_colors_by_group.md),
[`contrast_names_to_sedesign()`](https://jmw86069.github.io/jamses/reference/contrast_names_to_sedesign.md),
[`contrastnames()`](https://jmw86069.github.io/jamses/reference/contrastnames.md),
[`contrasts()`](https://jmw86069.github.io/jamses/reference/contrasts.md),
[`contrasts<-()`](https://jmw86069.github.io/jamses/reference/contrasts-set.md),
[`contrasts_to_factors()`](https://jmw86069.github.io/jamses/reference/contrasts_to_factors.md),
[`contrasts_to_venn_setlists()`](https://jmw86069.github.io/jamses/reference/contrasts_to_venn_setlists.md),
[`design,SEDesign-method`](https://jmw86069.github.io/jamses/reference/design.md),
[`draw_oneway_contrast()`](https://jmw86069.github.io/jamses/reference/draw_oneway_contrast.md),
[`draw_twoway_contrast()`](https://jmw86069.github.io/jamses/reference/draw_twoway_contrast.md),
[`factors()`](https://jmw86069.github.io/jamses/reference/factors.md),
[`filter_contrast_names()`](https://jmw86069.github.io/jamses/reference/filter_contrast_names.md),
[`groups()`](https://jmw86069.github.io/jamses/reference/groups.md),
[`groups_to_sedesign()`](https://jmw86069.github.io/jamses/reference/groups_to_sedesign.md),
[`plot.SEDesign()`](https://jmw86069.github.io/jamses/reference/plot.SEDesign.md),
[`plot_sedesign()`](https://jmw86069.github.io/jamses/reference/plot_sedesign.md),
[`print,SEDesign-method`](https://jmw86069.github.io/jamses/reference/print-SEDesign-method.md),
[`samples()`](https://jmw86069.github.io/jamses/reference/samples.md),
[`sedesign_to_factors()`](https://jmw86069.github.io/jamses/reference/sedesign_to_factors.md),
[`sort_contrasts()`](https://jmw86069.github.io/jamses/reference/sort_contrasts.md),
[`validate_sedesign()`](https://jmw86069.github.io/jamses/reference/validate_sedesign.md)
