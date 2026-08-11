# Print / show method for SEDesign objects

S7 objects bypass the classic S4 `show` generic for console
auto-printing; a plain S3 `print.SEDesign` method is used instead, since
S7 objects retain "SEDesign" in their S3
[`class()`](https://rdrr.io/r/base/class.html) vector.

## Usage

``` r
# S3 method for class 'SEDesign'
print(x, ...)
```

## Arguments

- x:

  `SEDesign` object

- ...:

  additional arguments are ignored

## Details

The summary includes:

- the number of samples, groups, and contrasts

- Each factor name, with factor levels.

Factor levels are printed in order, so that the first entry is
prioritized as the control in subsequent contrasts.

For example: 'level1', 'level2', 'level3' should always produce
contrasts 'level2-level1', 'level3-level1', 'level3-level2'. It should
not produce a contrast 'level2-level3'.

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
[`plot_sedesign()`](plot_sedesign.md), [`samples()`](samples.md),
[`sedesign_to_factors()`](sedesign_to_factors.md),
[`sort_contrasts()`](sort_contrasts.md),
[`validate_sedesign()`](validate_sedesign.md)
