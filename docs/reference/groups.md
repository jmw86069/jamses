# Design group names for SEDesign objects

`groups()` returns the design group names, equivalent to
`colnames(design(object))`. `groups()<-` renames design groups, updating
`colnames(design)` and `rownames(contrasts)` together.

## Usage

``` r
groups(object, ...)

groups(object, ...) <- value
```

## Arguments

- object:

  `SEDesign` object

- ...:

  additional arguments are ignored

- value:

  `character` vector, length equal to `ncol(design(object))`

## Details

In principle, design groups are derived identifiers and should not need
to be renamed directly (prefer re-running
[`groups_to_sedesign()`](groups_to_sedesign.md) with updated inputs).
`groups()<-` is provided for convenience, but use it with care:
`design_df` (and therefore [`factors()`](factors.md) labels derived from
group names) is rebuilt from the new group names, which resets any
customized `factors()<-` labels if the number of underlying factors
changes.

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
[`factors()`](factors.md),
[`filter_contrast_names()`](filter_contrast_names.md),
[`groups_to_sedesign()`](groups_to_sedesign.md),
[`plot.SEDesign()`](plot.SEDesign.md),
[`plot_sedesign()`](plot_sedesign.md),
[`print.SEDesign()`](print.SEDesign.md), [`samples()`](samples.md),
[`sedesign_to_factors()`](sedesign_to_factors.md),
[`sort_contrasts()`](sort_contrasts.md),
[`validate_sedesign()`](validate_sedesign.md)
