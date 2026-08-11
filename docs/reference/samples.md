# Sample identifiers for SEDesign objects

`samples()` returns sample identifiers, equivalent to
`rownames(design(object))`. `samples()<-` renames samples, and updates
`rownames(design(object))` to match. It accepts a `character` vector
`value` with length equal to `length(samples(object))`.

## Usage

``` r
samples(object, ...)

samples(object, ...) <- value
```

## Arguments

- object:

  `SEDesign` object

- ...:

  additional arguments are ignored

- value:

  `character` vector, length equal to `length(samples(object))`

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
[`groups()`](groups.md),
[`groups_to_sedesign()`](groups_to_sedesign.md),
[`plot.SEDesign()`](plot.SEDesign.md),
[`plot_sedesign()`](plot_sedesign.md),
[`print.SEDesign()`](print.SEDesign.md),
[`sedesign_to_factors()`](sedesign_to_factors.md),
[`sort_contrasts()`](sort_contrasts.md),
[`validate_sedesign()`](validate_sedesign.md)
