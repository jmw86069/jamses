# Convert contrast names to sedesign object

Convert contrast names to sedesign object

## Usage

``` r
contrast_names_to_sedesign(
  contrast_names,
  factor_sep = "_",
  contrast_sep = "-",
  ...
)
```

## Arguments

- contrast_names:

  `character` vector of contrast names where factors are separated by
  `factor_sep` and contrasts are separated by `contrast_sep`.

- factor_sep, contrast_sep:

  `character` strings used as delimiters between factors and contrasts,
  respectively.

- ...:

  additional arguments are ignored.

## Value

`sedesign` object, or if the input `contrast_names` contain mixed number
of factors per contrast, the output is split into a `list` of `sedesign`
objects based upon the number of factors.

## Details

This function is a convenience function intended only to convert a
vector of contrast names into `sedesign` format for use in other
functions. It assumes only one sample replicate per group for this
purpose.

One utility of this function is to convert two-way contrast names into a
contrast matrix, to test whether the contrast defined is equivalent for
the two names.

## See also

Other jam experiment design: [`SEDesign()`](SEDesign.md),
[`[,SEDesign-method`](sub-SEDesign-method.md),
[`check_sedesign()`](check_sedesign.md),
[`contrast2comp()`](contrast2comp.md),
[`contrast_colors_by_group()`](contrast_colors_by_group.md),
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

## Examples

``` r
contrast_names_3fac <- c(
   "(CellA_Treated-CellA_Control)-(CellB_Treated-CellB_Control)",
   "(CellA_Treated-CellB_Treated)-(CellA_Control-CellB_Control)")
contrast_names_to_sedesign(contrast_names_3fac)
#> <SEDesign>
#>  @ design      : num [1:4, 1:4] 1 0 0 0 0 1 0 0 0 0 ...
#>  .. - attr(*, "dimnames")=List of 2
#>  ..  ..$ : chr [1:4] "CellA_Treated" "CellA_Control" "CellB_Treated" "CellB_Control"
#>  ..  ..$ : chr [1:4] "CellA_Treated" "CellA_Control" "CellB_Treated" "CellB_Control"
#>  @ contrasts   : num [1:4, 1:2] 1 -1 -1 1 1 -1 -1 1
#>  .. - attr(*, "dimnames")=List of 2
#>  ..  ..$ Levels   : chr [1:4] "CellA_Treated" "CellA_Control" "CellB_Treated" "CellB_Control"
#>  ..  ..$ Contrasts: chr [1:2] "(CellA_Treated-CellA_Control)-(CellB_Treated-CellB_Control)" "(CellA_Treated-CellB_Treated)-(CellA_Control-CellB_Control)"
#>  @ samples     : chr [1:4] "CellA_Treated" "CellA_Control" "CellB_Treated" ...
#>  @ factors     : chr [1:2] "factor1" "factor2"
#>  @ design_df   :'data.frame':    4 obs. of  2 variables:
#>  .. $ factor1: chr  "CellA" "CellA" "CellB" "CellB"
#>  .. $ factor2: chr  "Treated" "Control" "Treated" "Control"
#>  @ contrasts_df:'data.frame':    2 obs. of  2 variables:
#>  .. $ factor1: chr  "CellA-CellB" "CellA-CellB"
#>  .. $ factor2: chr  "Treated-Control" "Treated-Control"

contrast_names_3way <- c(
   paste0("((CellA_Treated_Mut-CellB_Treated_Mut)-",
      "(CellA_Control_Mut-CellB_Control_Mut))-",
      "((CellA_Treated_WT-CellB_Treated_WT)-",
      "(CellA_Control_WT-CellB_Control_WT))"),
   paste0("((CellA_Treated_Mut-CellA_Control_Mut)-",
      "(CellB_Treated_Mut-CellB_Control_Mut))-",
      "((CellA_Treated_WT-CellA_Control_WT)-",
      "(CellB_Treated_WT-CellB_Control_WT))"))
contrast_names_to_sedesign(contrast_names_3way)
#> <SEDesign>
#>  @ design      : num [1:8, 1:8] 1 0 0 0 0 0 0 0 0 0 ...
#>  .. - attr(*, "dimnames")=List of 2
#>  ..  ..$ : chr [1:8] "CellA_Treated_Mut" "CellB_Treated_Mut" "CellA_Control_Mut" "CellB_Control_Mut" ...
#>  ..  ..$ : chr [1:8] "CellA_Treated_Mut" "CellA_Treated_WT" "CellA_Control_Mut" "CellA_Control_WT" ...
#>  @ contrasts   : num [1:8, 1:2] 1 -1 -1 1 -1 1 1 -1 1 -1 ...
#>  .. - attr(*, "dimnames")=List of 2
#>  ..  ..$ Levels   : chr [1:8] "CellA_Treated_Mut" "CellA_Treated_WT" "CellA_Control_Mut" "CellA_Control_WT" ...
#>  ..  ..$ Contrasts: chr [1:2] "((CellA_Treated_Mut-CellB_Treated_Mut)-(CellA_Control_Mut-CellB_Control_Mut))-((CellA_Treated_WT-CellB_Treated_"| __truncated__ "((CellA_Treated_Mut-CellA_Control_Mut)-(CellB_Treated_Mut-CellB_Control_Mut))-((CellA_Treated_WT-CellA_Control_"| __truncated__
#>  @ samples     : chr [1:8] "CellA_Treated_Mut" "CellB_Treated_Mut" "CellA_Control_Mut" ...
#>  @ factors     : chr [1:3] "factor1" "factor2" "factor3"
#>  @ design_df   :'data.frame':    8 obs. of  3 variables:
#>  .. $ factor1: chr  "CellA" "CellA" "CellA" "CellA" ...
#>  .. $ factor2: chr  "Treated" "Treated" "Control" "Control" ...
#>  .. $ factor3: chr  "Mut" "WT" "Mut" "WT" ...
#>  @ contrasts_df:'data.frame':    2 obs. of  3 variables:
#>  .. $ factor1: chr  "CellA-CellB" "CellA-CellB"
#>  .. $ factor2: chr  "Treated-Control" "Treated-Control"
#>  .. $ factor3: chr  "Mut-WT" "Mut-WT"
```
