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

Other jam experiment design:
[`SEDesign()`](https://jmw86069.github.io/jamses/reference/SEDesign.md),
[`[,SEDesign-method`](https://jmw86069.github.io/jamses/reference/sub-SEDesign-method.md),
[`check_sedesign()`](https://jmw86069.github.io/jamses/reference/check_sedesign.md),
[`contrast2comp()`](https://jmw86069.github.io/jamses/reference/contrast2comp.md),
[`contrast_colors_by_group()`](https://jmw86069.github.io/jamses/reference/contrast_colors_by_group.md),
[`contrastnames()`](https://jmw86069.github.io/jamses/reference/contrastnames.md),
[`contrasts()`](https://jmw86069.github.io/jamses/reference/contrasts.md),
[`contrasts<-()`](https://jmw86069.github.io/jamses/reference/contrasts-set.md),
[`contrasts_to_factors()`](https://jmw86069.github.io/jamses/reference/contrasts_to_factors.md),
[`contrasts_to_venn_setlists()`](https://jmw86069.github.io/jamses/reference/contrasts_to_venn_setlists.md),
[`design,SEDesign-method`](https://jmw86069.github.io/jamses/reference/design.md),
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

## Examples

``` r
contrast_names_3fac <- c(
   "(CellA_Treated-CellA_Control)-(CellB_Treated-CellB_Control)",
   "(CellA_Treated-CellB_Treated)-(CellA_Control-CellB_Control)")
contrast_names_to_sedesign(contrast_names_3fac)
#> <SEDesign> 4 samples, 4 groups, 2 contrasts
#> factors:
#>   - factor1: CellA, CellB
#>   - factor2: Treated, Control

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
#> <SEDesign> 8 samples, 8 groups, 2 contrasts
#> factors:
#>   - factor1: CellA, CellB
#>   - factor2: Treated, Control
#>   - factor3: Mut, WT
```
