# Convert contrasts to data.frame of design factors

Convert contrasts to data.frame of design factors, each factor in a
column, with factor levels or with factor level contrasts.

## Usage

``` r
contrasts_to_factors(
  contrast_names = NULL,
  sedesign = NULL,
  factor_names = NULL,
  factor_sep = "_",
  rowname = c("contrast", "comp"),
  verbose = FALSE,
  ...
)
```

## Arguments

- contrast_names:

  `character` vector of contrast names, or `SEDesign` object, which is
  equivalent to argument 'sedesign'.

- sedesign:

  `SEDesign` object, used when `contrast_names` is not supplied.

- factor_names:

  `character` vector of colnames to use for the resulting `data.frame`,
  typically the name of each experimental factor.

- factor_sep:

  `character` string delimited between factors used in each group name.
  It is passed to
  [`contrast2comp()`](https://jmw86069.github.io/jamses/reference/contrast2comp.md)
  as `contrast_factor_delim`.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are passed to
  [`contrast2comp()`](https://jmw86069.github.io/jamses/reference/contrast2comp.md)
  as relevant.

## Value

`data.frame` with factors in each column, and values indicating either
individual factor levels, or comparison of two factor levels.

## Details

This function is intended to summarize contrasts by factor columns, with
factor levels or with factor level contrasts, to make it easier to
review which factor or factors are being compared in each contrast.

**Note:** It assumes all contrasts are "proper", which is defined as a
contrast where only one factor is compared at a time for each contrast.

A proper contrast:

- `'A_B-C_B'`

  - This contrast compares 'factor1': `(A-C)`, both with 'factor2': `B`

An improper contrast:

- `'A_B-C_D'`

  - This contrast compares 'factor1': `(A-C)`, while 'factor2' is also
    changing: `(B-D)`.

  - However, it is "improper" because the factors are not independently
    controlled, instead the factors are confounded into one contrast.

  - The "proper" form, for the purpose of the assumptions in this
    function, is defined as a "two-way" style contrast, a "diff of
    diffs":

  `(A_B-C_B)-(A_D-B_D)`

### Todo

- Currently, "improper" contrasts are returned in the first column, with
  no additional information.

- In future, "improper" contrasts should be flagged, or handled in a way
  that preserves the factor levels represented while still representing
  the original contrast.

## See also

Other jam experiment design:
[`SEDesign()`](https://jmw86069.github.io/jamses/reference/SEDesign.md),
[`[,SEDesign-method`](https://jmw86069.github.io/jamses/reference/sub-SEDesign-method.md),
[`check_sedesign()`](https://jmw86069.github.io/jamses/reference/check_sedesign.md),
[`contrast2comp()`](https://jmw86069.github.io/jamses/reference/contrast2comp.md),
[`contrast_colors_by_group()`](https://jmw86069.github.io/jamses/reference/contrast_colors_by_group.md),
[`contrast_names_to_sedesign()`](https://jmw86069.github.io/jamses/reference/contrast_names_to_sedesign.md),
[`contrastnames()`](https://jmw86069.github.io/jamses/reference/contrastnames.md),
[`contrasts()`](https://jmw86069.github.io/jamses/reference/contrasts.md),
[`contrasts<-()`](https://jmw86069.github.io/jamses/reference/contrasts-set.md),
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

## Examples

``` r
group_names <- paste0(
   rep(c("UL3", "dH1A", "dH1B"), each=5), "_",
   c("Veh", "DEX", "PMA", "SF", "Ins"))
sedesign <- groups_to_sedesign(group_names)
contrasts_to_factors(sedesign)
#>                                           factor1 factor2
#> dH1A_Veh-UL3_Veh                         dH1A-UL3     Veh
#> dH1B_Veh-UL3_Veh                         dH1B-UL3     Veh
#> dH1B_Veh-dH1A_Veh                       dH1B-dH1A     Veh
#> dH1A_DEX-UL3_DEX                         dH1A-UL3     DEX
#> dH1B_DEX-UL3_DEX                         dH1B-UL3     DEX
#> dH1B_DEX-dH1A_DEX                       dH1B-dH1A     DEX
#> dH1A_PMA-UL3_PMA                         dH1A-UL3     PMA
#> dH1B_PMA-UL3_PMA                         dH1B-UL3     PMA
#> dH1B_PMA-dH1A_PMA                       dH1B-dH1A     PMA
#> dH1A_SF-UL3_SF                           dH1A-UL3      SF
#> dH1B_SF-UL3_SF                           dH1B-UL3      SF
#> dH1B_SF-dH1A_SF                         dH1B-dH1A      SF
#> dH1A_Ins-UL3_Ins                         dH1A-UL3     Ins
#> dH1B_Ins-UL3_Ins                         dH1B-UL3     Ins
#> dH1B_Ins-dH1A_Ins                       dH1B-dH1A     Ins
#> UL3_DEX-UL3_Veh                               UL3 DEX-Veh
#> UL3_PMA-UL3_Veh                               UL3 PMA-Veh
#> UL3_SF-UL3_Veh                                UL3  SF-Veh
#> UL3_Ins-UL3_Veh                               UL3 Ins-Veh
#> UL3_PMA-UL3_DEX                               UL3 PMA-DEX
#> UL3_SF-UL3_DEX                                UL3  SF-DEX
#> UL3_Ins-UL3_DEX                               UL3 Ins-DEX
#> UL3_SF-UL3_PMA                                UL3  SF-PMA
#> UL3_Ins-UL3_PMA                               UL3 Ins-PMA
#> UL3_Ins-UL3_SF                                UL3  Ins-SF
#> dH1A_DEX-dH1A_Veh                            dH1A DEX-Veh
#> dH1A_PMA-dH1A_Veh                            dH1A PMA-Veh
#> dH1A_SF-dH1A_Veh                             dH1A  SF-Veh
#> dH1A_Ins-dH1A_Veh                            dH1A Ins-Veh
#> dH1A_PMA-dH1A_DEX                            dH1A PMA-DEX
#> dH1A_SF-dH1A_DEX                             dH1A  SF-DEX
#> dH1A_Ins-dH1A_DEX                            dH1A Ins-DEX
#> dH1A_SF-dH1A_PMA                             dH1A  SF-PMA
#> dH1A_Ins-dH1A_PMA                            dH1A Ins-PMA
#> dH1A_Ins-dH1A_SF                             dH1A  Ins-SF
#> dH1B_DEX-dH1B_Veh                            dH1B DEX-Veh
#> dH1B_PMA-dH1B_Veh                            dH1B PMA-Veh
#> dH1B_SF-dH1B_Veh                             dH1B  SF-Veh
#> dH1B_Ins-dH1B_Veh                            dH1B Ins-Veh
#> dH1B_PMA-dH1B_DEX                            dH1B PMA-DEX
#> dH1B_SF-dH1B_DEX                             dH1B  SF-DEX
#> dH1B_Ins-dH1B_DEX                            dH1B Ins-DEX
#> dH1B_SF-dH1B_PMA                             dH1B  SF-PMA
#> dH1B_Ins-dH1B_PMA                            dH1B Ins-PMA
#> dH1B_Ins-dH1B_SF                             dH1B  Ins-SF
#> (dH1A_DEX-UL3_DEX)-(dH1A_Veh-UL3_Veh)    dH1A-UL3 DEX-Veh
#> (dH1A_PMA-UL3_PMA)-(dH1A_Veh-UL3_Veh)    dH1A-UL3 PMA-Veh
#> (dH1A_SF-UL3_SF)-(dH1A_Veh-UL3_Veh)      dH1A-UL3  SF-Veh
#> (dH1A_Ins-UL3_Ins)-(dH1A_Veh-UL3_Veh)    dH1A-UL3 Ins-Veh
#> (dH1A_PMA-UL3_PMA)-(dH1A_DEX-UL3_DEX)    dH1A-UL3 PMA-DEX
#> (dH1A_SF-UL3_SF)-(dH1A_DEX-UL3_DEX)      dH1A-UL3  SF-DEX
#> (dH1A_Ins-UL3_Ins)-(dH1A_DEX-UL3_DEX)    dH1A-UL3 Ins-DEX
#> (dH1A_SF-UL3_SF)-(dH1A_PMA-UL3_PMA)      dH1A-UL3  SF-PMA
#> (dH1A_Ins-UL3_Ins)-(dH1A_PMA-UL3_PMA)    dH1A-UL3 Ins-PMA
#> (dH1A_Ins-UL3_Ins)-(dH1A_SF-UL3_SF)      dH1A-UL3  Ins-SF
#> (dH1B_DEX-UL3_DEX)-(dH1B_Veh-UL3_Veh)    dH1B-UL3 DEX-Veh
#> (dH1B_PMA-UL3_PMA)-(dH1B_Veh-UL3_Veh)    dH1B-UL3 PMA-Veh
#> (dH1B_SF-UL3_SF)-(dH1B_Veh-UL3_Veh)      dH1B-UL3  SF-Veh
#> (dH1B_Ins-UL3_Ins)-(dH1B_Veh-UL3_Veh)    dH1B-UL3 Ins-Veh
#> (dH1B_PMA-UL3_PMA)-(dH1B_DEX-UL3_DEX)    dH1B-UL3 PMA-DEX
#> (dH1B_SF-UL3_SF)-(dH1B_DEX-UL3_DEX)      dH1B-UL3  SF-DEX
#> (dH1B_Ins-UL3_Ins)-(dH1B_DEX-UL3_DEX)    dH1B-UL3 Ins-DEX
#> (dH1B_SF-UL3_SF)-(dH1B_PMA-UL3_PMA)      dH1B-UL3  SF-PMA
#> (dH1B_Ins-UL3_Ins)-(dH1B_PMA-UL3_PMA)    dH1B-UL3 Ins-PMA
#> (dH1B_Ins-UL3_Ins)-(dH1B_SF-UL3_SF)      dH1B-UL3  Ins-SF
#> (dH1B_DEX-dH1A_DEX)-(dH1B_Veh-dH1A_Veh) dH1B-dH1A DEX-Veh
#> (dH1B_PMA-dH1A_PMA)-(dH1B_Veh-dH1A_Veh) dH1B-dH1A PMA-Veh
#> (dH1B_SF-dH1A_SF)-(dH1B_Veh-dH1A_Veh)   dH1B-dH1A  SF-Veh
#> (dH1B_Ins-dH1A_Ins)-(dH1B_Veh-dH1A_Veh) dH1B-dH1A Ins-Veh
#> (dH1B_PMA-dH1A_PMA)-(dH1B_DEX-dH1A_DEX) dH1B-dH1A PMA-DEX
#> (dH1B_SF-dH1A_SF)-(dH1B_DEX-dH1A_DEX)   dH1B-dH1A  SF-DEX
#> (dH1B_Ins-dH1A_Ins)-(dH1B_DEX-dH1A_DEX) dH1B-dH1A Ins-DEX
#> (dH1B_SF-dH1A_SF)-(dH1B_PMA-dH1A_PMA)   dH1B-dH1A  SF-PMA
#> (dH1B_Ins-dH1A_Ins)-(dH1B_PMA-dH1A_PMA) dH1B-dH1A Ins-PMA
#> (dH1B_Ins-dH1A_Ins)-(dH1B_SF-dH1A_SF)   dH1B-dH1A  Ins-SF
```
