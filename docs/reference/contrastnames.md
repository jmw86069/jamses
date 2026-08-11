# Contrast names for SEDesign objects

`contrastnames()` and `contrast_names()` are equivalent, returning
`colnames(contrasts(object))`. `contrastnames()<-` renames contrast
columns directly, accepting a `character` vector `value` with length
equal to `ncol(contrasts(object))`. `contrast_names()<-` instead
rebuilds the contrasts matrix using
[`limma::makeContrasts()`](https://rdrr.io/pkg/limma/man/makeContrasts.html),
accepting a `character` vector `value` of contrast names (must not
contain duplicated values).

## Usage

``` r
contrastnames(object, ...)

contrastnames(object, ...) <- value

contrast_names(object, ...)

contrast_names(object, ...) <- value
```

## Arguments

- object:

  `SEDesign` object

- ...:

  additional arguments are ignored

- value:

  `character` vector: contrast column names (for `contrastnames<-`) or
  contrast names to rebuild via
  [`limma::makeContrasts()`](https://rdrr.io/pkg/limma/man/makeContrasts.html)
  (for `contrast_names<-`)

## See also

Other jam experiment design: [`SEDesign()`](SEDesign.md),
[`[,SEDesign-method`](sub-SEDesign-method.md),
[`check_sedesign()`](check_sedesign.md),
[`contrast2comp()`](contrast2comp.md),
[`contrast_colors_by_group()`](contrast_colors_by_group.md),
[`contrast_names_to_sedesign()`](contrast_names_to_sedesign.md),
[`contrasts()`](contrasts.md), [`contrasts<-()`](contrasts-set.md),
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
