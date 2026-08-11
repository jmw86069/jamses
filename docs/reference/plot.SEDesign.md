# Plot method for SEDesign objects

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method for
`SEDesign` objects, a thin wrapper around
[`plot_sedesign()`](plot_sedesign.md). S7 objects retain their class in
the S3 [`class()`](https://rdrr.io/r/base/class.html) vector (see
[`print.SEDesign()`](print.SEDesign.md)), so this S3 method is
dispatched by the base
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) generic for
`SEDesign` objects.

## Usage

``` r
# S3 method for class 'SEDesign'
plot(x, ...)
```

## Arguments

- x:

  `SEDesign` object as returned by
  [`groups_to_sedesign()`](groups_to_sedesign.md).

- ...:

  additional arguments passed to [`plot_sedesign()`](plot_sedesign.md).

## Value

see [`plot_sedesign()`](plot_sedesign.md).

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
[`plot_sedesign()`](plot_sedesign.md),
[`print.SEDesign()`](print.SEDesign.md), [`samples()`](samples.md),
[`sedesign_to_factors()`](sedesign_to_factors.md),
[`sort_contrasts()`](sort_contrasts.md),
[`validate_sedesign()`](validate_sedesign.md)

## Examples

``` r
isamples_1 <- paste0(
   rep(c("DMSO", "Etop", "DMSO", "Etop"), each=6),
   "_",
   rep(c("NF", "Flag"), each=12),
   "_",
   rep(c("WT", "KO", "WT", "KO", "WT", "D955N", "WT", "D955N"), each=3),
   "_",
   LETTERS[1:3])
idf <- data.frame(jamba::rbindList(strsplit(isamples_1, "_")))[,1:3]
rownames(idf) <- isamples_1;
sedesign_1 <- groups_to_sedesign(idf)

plot(sedesign_1)

```
