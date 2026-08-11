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
  Editing [`factors()`](factors.md) only renames these labels (used
  primarily by [`print()`](https://rdrr.io/r/base/print.html)); it has
  no effect on `design`, `contrasts`, or [`groups()`](groups.md).

## Details

In addition to the four constructor arguments below, `SEDesign` objects
carry two internal-use properties, accessible via `@`:

- `design_df`: `data.frame` with one row per design group
  (`groups(object)`, i.e. `colnames(design)`), derived by splitting each
  group name on `"_"` into one column per underlying experimental
  factor. Column names default to `"factor1"`, `"factor2"`, etc., and
  can be customized via `factors()<-`.

- `contrasts_df`: `data.frame`, cached result of
  [`contrasts_to_factors()`](contrasts_to_factors.md), refreshed
  whenever `design` or `contrasts` are updated.

## See also

Other jam experiment design:
[`[,SEDesign-method`](sub-SEDesign-method.md),
[`check_sedesign()`](check_sedesign.md),
[`contrast2comp()`](contrast2comp.md),
[`contrast_colors_by_group()`](contrast_colors_by_group.md),
[`contrast_names_to_sedesign()`](contrast_names_to_sedesign.md),
[`contrastnames()`](contrastnames.md), [`contrasts()`](contrasts.md),
[`contrasts<-()`](contrasts-set.md),
[`contrasts_to_factors()`](contrasts_to_factors.md),
[`contrasts_to_venn_setlists()`](contrasts_to_venn_setlists.md),
[`design,SEDesign-method`](design-SEDesign-method.md),
`design<-,SEDesign,matrix-method`,
[`draw_oneway_contrast()`](draw_oneway_contrast.md),
[`draw_twoway_contrast()`](draw_twoway_contrast.md),
[`factors()`](factors.md),
[`filter_contrast_names()`](filter_contrast_names.md),
[`groups()`](groups.md),
[`groups_to_sedesign()`](groups_to_sedesign.md),
[`plot.SEDesign()`](plot.SEDesign.md),
[`plot_sedesign()`](plot_sedesign.md),
[`print.SEDesign()`](print.SEDesign.md), [`samples()`](samples.md),
[`sedesign_to_factors()`](sedesign_to_factors.md),
[`sort_contrasts()`](sort_contrasts.md),
[`validate_sedesign()`](validate_sedesign.md)
