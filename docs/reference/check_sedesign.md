# Check SEDesign object

Check whether a SEDesign object is valid

## Usage

``` r
check_sedesign(object)
```

## Arguments

- object:

  `SEDesign` object

## Details

This function checks whether an `SEDesign` object is valid:

- if `samples` is provided, and if `design` is provided, `samples` must
  match `rownames(design)`.

- if `design` is provided, and if `contrasts` is provided,
  `colnames(design)` must match `rownames(contrasts)`.

- if `contrasts` is provided, `design` must also be provided.

Note that `samples` can be a subset of `rownames(design)`, in which case
the `design` will also be subset accordingly.

Similarly, `colnames(design)` can be a subset of `rownames(contrasts)`,
which would force `contrasts` to be subset accordingly.

Typically the order of `samples` should match the order of
`rownames(design)` but this is not required. Downstream methods should
confirm this order.

Typically the order of `colnames(design)` should match the order of
`rownames(contrast)` but this is not required. Downstream methods should
confirm this order.

## See also

Other jam experiment design: [`SEDesign()`](SEDesign.md),
[`[,SEDesign-method`](sub-SEDesign-method.md),
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
