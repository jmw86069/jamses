# Subset a SEDesign object by samples and/or design groups

Subset a SEDesign object by samples and/or design groups

## Usage

``` r
## S7 method for class <SEDesign>
x[i = NULL, j = NULL, k = NULL, ...]
```

## Arguments

- x:

  `SEDesign` object

- i, j:

  optional `character` or `integer` vectors used to subset samples (`i`)
  and/or design groups (`j`). See
  [`validate_sedesign()`](https://jmw86069.github.io/jamses/reference/validate_sedesign.md).

- ...:

  additional arguments are ignored

## See also

Other jam experiment design:
[`SEDesign()`](https://jmw86069.github.io/jamses/reference/SEDesign.md),
[`check_sedesign()`](https://jmw86069.github.io/jamses/reference/check_sedesign.md),
[`contrast2comp()`](https://jmw86069.github.io/jamses/reference/contrast2comp.md),
[`contrast_colors_by_group()`](https://jmw86069.github.io/jamses/reference/contrast_colors_by_group.md),
[`contrast_names_to_sedesign()`](https://jmw86069.github.io/jamses/reference/contrast_names_to_sedesign.md),
[`contrastnames()`](https://jmw86069.github.io/jamses/reference/contrastnames.md),
[`contrasts()`](https://jmw86069.github.io/jamses/reference/contrasts.md),
[`contrasts<-()`](https://jmw86069.github.io/jamses/reference/contrasts-set.md),
[`contrasts_to_factors()`](https://jmw86069.github.io/jamses/reference/contrasts_to_factors.md),
[`contrasts_to_venn_setlists()`](https://jmw86069.github.io/jamses/reference/contrasts_to_venn_setlists.md),
[`design,SEDesign-method`](https://jmw86069.github.io/jamses/reference/design-SEDesign-method.md),
`design<-,SEDesign,matrix-method`,
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
