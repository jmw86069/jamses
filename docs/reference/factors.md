# Experimental factor labels for SEDesign objects

`factors()` returns the labels of the underlying experimental factors,
equivalent to `colnames(x@design_df)`, i.e. the columns produced by
splitting each design group name (`groups(object)`) on `"_"`.
`factors()<-` renames these labels; it has no effect on `design`,
`contrasts`, [`groups()`](groups.md), or `contrasts_df` – it is used
only for internal bookkeeping and for the
[`print()`](https://rdrr.io/r/base/print.html) summary.

## Usage

``` r
factors(object, ...)

factors(object, ...) <- value
```

## Arguments

- object:

  `SEDesign` object

- ...:

  additional arguments are ignored

- value:

  `character` vector, length equal to `ncol(x@design_df)` (the number of
  underlying experimental factors)

## See also

Other jam experiment design: [`SEDesign()`](SEDesign.md),
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
[`filter_contrast_names()`](filter_contrast_names.md),
[`groups()`](groups.md),
[`groups_to_sedesign()`](groups_to_sedesign.md),
[`plot.SEDesign()`](plot.SEDesign.md),
[`plot_sedesign()`](plot_sedesign.md),
[`print.SEDesign()`](print.SEDesign.md), [`samples()`](samples.md),
[`sedesign_to_factors()`](sedesign_to_factors.md),
[`sort_contrasts()`](sort_contrasts.md),
[`validate_sedesign()`](validate_sedesign.md)
