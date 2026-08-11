# Sort contrasts by factor and level

Sort contrasts by factor and level, given contrast names, or `SEDesign`,
or `data.frame` of factors.

## Usage

``` r
sort_contrasts(x, factor_order = NULL, ...)
```

## Arguments

- x:

  one of:

  - `character` vector of contrast names

  - `SEDesign`

  - `data.frame` with factors

- factor_order:

  `integer`, default NULL uses all factors, or specify the factor
  ordering to use for sorting.

- ...:

  additional arguments are passed to internal functions. For example
  'factor_names' is passed to
  [`contrasts_to_factors()`](contrasts_to_factors.md).

## Value

`data.frame` with factors in each column, and factor levels or factor
level contrasts as column values. The
[`rownames()`](https://rdrr.io/r/base/colnames.html) contain contrast
names.

## Sort order

- Contrasts are sorted by depth: one-factor comparisons, two-factor
  comparisons, etc.

- Contrasts are sorted for factor level contrasts in each factor in
  `factor_order`. If `factor_order=1` then contrasts appear first in the
  first factor column.

- Contrasts are then sorted by that factor column using the observed
  contrast order.

- Contrasts are then sorted by other factor columns as they appear in
  `factor_order`, using the observed level order in each factor column.

## Sort order, described another way:

- All one-way contrasts will appear at the top.

  - The first one-way contrasts will be comparisons using the first
    value in `factor_order`, sorted by the comparison, then sorted by
    each remaining column in `factor_order`.

  - The next set of one-way contrasts will be the second value in
    `factor_order`, sorted by that comparison, then sorted by each
    remaining column in `factor_order`.

- All two-way contrasts will appear at the end, since each two-way
  contrast involves two factors.

  - They will be sorted to involve contrasts in the same order as as
    `factor_order`.

  - The first two-way contrasts will involve comparisons involving the
    first two values in `factor_order`.

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
[`print.SEDesign()`](print.SEDesign.md), [`samples()`](samples.md),
[`sedesign_to_factors()`](sedesign_to_factors.md),
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
# simple data.frame with group information
idf <- data.frame(jamba::rbindList(strsplit(isamples_1, "_")))[,1:3]
rownames(idf) <- isamples_1;
colnames(idf) <- c("Treatment", "Flag", "Genotype")
# convert to sedesign
sedesign <- groups_to_sedesign(idf)
sort_contrasts(sedesign)
#>                                                               Treatment    Flag
#> Etop_NF_WT-DMSO_NF_WT                                         Etop-DMSO      NF
#> Etop_NF_KO-DMSO_NF_KO                                         Etop-DMSO      NF
#> Etop_Flag_WT-DMSO_Flag_WT                                     Etop-DMSO    Flag
#> Etop_Flag_D955N-DMSO_Flag_D955N                               Etop-DMSO    Flag
#> DMSO_Flag_WT-DMSO_NF_WT                                            DMSO Flag-NF
#> Etop_Flag_WT-Etop_NF_WT                                            Etop Flag-NF
#> DMSO_NF_KO-DMSO_NF_WT                                              DMSO      NF
#> Etop_NF_KO-Etop_NF_WT                                              Etop      NF
#> DMSO_Flag_D955N-DMSO_Flag_WT                                       DMSO    Flag
#> Etop_Flag_D955N-Etop_Flag_WT                                       Etop    Flag
#> (Etop_Flag_WT-DMSO_Flag_WT)-(Etop_NF_WT-DMSO_NF_WT)           Etop-DMSO Flag-NF
#> (Etop_NF_KO-DMSO_NF_KO)-(Etop_NF_WT-DMSO_NF_WT)               Etop-DMSO      NF
#> (Etop_Flag_D955N-DMSO_Flag_D955N)-(Etop_Flag_WT-DMSO_Flag_WT) Etop-DMSO    Flag
#>                                                               Genotype
#> Etop_NF_WT-DMSO_NF_WT                                               WT
#> Etop_NF_KO-DMSO_NF_KO                                               KO
#> Etop_Flag_WT-DMSO_Flag_WT                                           WT
#> Etop_Flag_D955N-DMSO_Flag_D955N                                  D955N
#> DMSO_Flag_WT-DMSO_NF_WT                                             WT
#> Etop_Flag_WT-Etop_NF_WT                                             WT
#> DMSO_NF_KO-DMSO_NF_WT                                            KO-WT
#> Etop_NF_KO-Etop_NF_WT                                            KO-WT
#> DMSO_Flag_D955N-DMSO_Flag_WT                                  D955N-WT
#> Etop_Flag_D955N-Etop_Flag_WT                                  D955N-WT
#> (Etop_Flag_WT-DMSO_Flag_WT)-(Etop_NF_WT-DMSO_NF_WT)                 WT
#> (Etop_NF_KO-DMSO_NF_KO)-(Etop_NF_WT-DMSO_NF_WT)                  KO-WT
#> (Etop_Flag_D955N-DMSO_Flag_D955N)-(Etop_Flag_WT-DMSO_Flag_WT) D955N-WT
```
