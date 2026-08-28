# Convert contrast names to Venn setlists for visual comparison

Convert contrast names to Venn setlists for visual comparison

## Usage

``` r
contrasts_to_venn_setlists(
  contrast_names = NULL,
  sestats = NULL,
  sedesign = NULL,
  include_multifactor = TRUE,
  include_singlefactor = TRUE,
  factor_names = NULL,
  contrast_style = c("contrast", "comp", "factors"),
  max_venn_size = 4,
  verbose = FALSE,
  ...
)
```

## Arguments

- contrast_names:

  `character` vector of contrast names, used as the priority input when
  supplied.

- sestats:

  `SEStats` or `list` output from
  [`se_contrast_stats()`](https://jmw86069.github.io/jamses/reference/se_contrast_stats.md),
  used only when `contrast_names` is not supplied.

- sedesign:

  `SEDesign` object, used only when neither `contrast_names` nor
  `sestats` are supplied.

- include_multifactor:

  `logical` indicating whether to include twoway contrasts in the Venn
  diagram logic. Currently the logic includes comparisons only across
  compatible twoway comparisons, it does not (yet) include Venn diagrams
  that include oneway and corresponding twoway comparisons together.
  Note: One of `include_multifactor` and `include_singlefactor` must be
  true.

- include_singlefactor:

  `logical` indicating whether to include single-factor contrasts in the
  Venn diagram logic. Note: One of `include_multifactor` and
  `include_singlefactor` must be true.

- contrast_style:

  `character` string indicating how to return the resulting Venn set
  lists:

  - `"contrast"` (default) returns a `list` with contrast names

  - `"comp"` returns a `list` with comp names from
    [`contrast2comp()`](https://jmw86069.github.io/jamses/reference/contrast2comp.md)

  - `"factors"` returns a `data.frame` from
    [`contrasts_to_factors()`](https://jmw86069.github.io/jamses/reference/contrasts_to_factors.md)
    with one column per design factor. In this case it may be useful to
    pass `factor_names` in order to assign column names to the factor
    columns.

- max_venn_size:

  `numeric` maximum number of groups to include per Venn diagram. When
  the number of contrasts to be included in a Venn set contains one
  extra contrast, the last two sets will be adjusted to accomodate the
  extra set:

  - when `max_venn_size=2` the last set will contain 3 members;

  - when `max_venn_size=3` (or higher) the last set will contain 2
    members, and the previous set will contain `(max_venn_size - 1)`
    members.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are passed to
  [`contrasts_to_factors()`](https://jmw86069.github.io/jamses/reference/contrasts_to_factors.md).

## Value

`list` with contrast names suitable for use in Venn diagrams.

## Details

This function is still under active development to be improved, feedback
is welcomed.

The motivation is to take a set of contrast names, and return reasonable
subsets of contrasts suitable for visual comparison using Venn diagrams.
Ultimately, the process is analogous to defining contrasts themselves:
keep experimental factors fixed while varying one factor at a time. The
difference is that experimental "factors" may themselves involve a
comparison.

The process is currently being tested for two-factor design scenarios,
and will be extended to handle higher factor designs in future.

### The process

- Contrasts are converted to `data.frame` with
  [`contrasts_to_factors()`](https://jmw86069.github.io/jamses/reference/contrasts_to_factors.md)

- Each factor column is iterated to produce sets of contrasts as
  follows:

  - The factor `data.frame` is subset for rows with a comparison in the
    factor column.

  - When `include_multifactor=FALSE` (not default) then data is filtered
    to remove rows with comparisons in any other factor columns.

  - When `include_singlefactor=FALSE` (not default) then data is
    filtered to remove rows with single values in any other factor
    columns.

  - The remaining rows are iteratively split using values in the other
    factor columns.

  - Remaining rows are also iteratively split by the depth of the
    contrast, oneway comparisons, and twoway comparisons.

  - If any subset contains more than `max_venn_size` rows, it is first
    split by the control factor level in the factor comparison, then it
    is split by the depth of the comparison.

  - In all cases, the resulting sets are split into subsets with size
    `max_venn_size`.

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
[`contrasts_to_factors()`](https://jmw86069.github.io/jamses/reference/contrasts_to_factors.md),
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

# by default it returns contrast names
venn_setlists <- contrasts_to_venn_setlists(sedesign=sedesign,
   include_multifactor=FALSE,
   factor_names=c("Genotype", "Treatment"))
jamba::sdim(venn_setlists)
#>                                       rows     class
#> dH1A-UL3,dH1B-UL3,dH1B-dH1A : Veh        3 character
#> dH1A-UL3,dH1B-UL3,dH1B-dH1A : DEX        3 character
#> dH1A-UL3,dH1B-UL3,dH1B-dH1A : PMA        3 character
#> dH1A-UL3,dH1B-UL3,dH1B-dH1A : SF         3 character
#> dH1A-UL3,dH1B-UL3,dH1B-dH1A : Ins        3 character
#> dH1A-UL3 : Veh,DEX,PMA                   3 character
#> dH1A-UL3 : SF,Ins                        2 character
#> dH1B-UL3 : Veh,DEX,PMA                   3 character
#> dH1B-UL3 : SF,Ins                        2 character
#> dH1B-dH1A : Veh,DEX,PMA                  3 character
#> dH1B-dH1A : SF,Ins                       2 character
#> UL3 : DEX-Veh,PMA-Veh,SF-Veh,Ins-Veh     4 character
#> UL3 : PMA-DEX,SF-DEX,Ins-DEX             3 character
#> UL3 : SF-PMA,Ins-PMA                     2 character
#> dH1A : DEX-Veh,PMA-Veh,SF-Veh,Ins-Veh    4 character
#> dH1A : PMA-DEX,SF-DEX,Ins-DEX            3 character
#> dH1A : SF-PMA,Ins-PMA                    2 character
#> dH1B : DEX-Veh,PMA-Veh,SF-Veh,Ins-Veh    4 character
#> dH1B : PMA-DEX,SF-DEX,Ins-DEX            3 character
#> dH1B : SF-PMA,Ins-PMA                    2 character
#> UL3,dH1A,dH1B : DEX-Veh                  3 character
#> UL3,dH1A,dH1B : Ins-DEX                  3 character
#> UL3,dH1A,dH1B : Ins-PMA                  3 character
#> UL3,dH1A,dH1B : Ins-SF                   3 character
#> UL3,dH1A,dH1B : Ins-Veh                  3 character
#> UL3,dH1A,dH1B : PMA-DEX                  3 character
#> UL3,dH1A,dH1B : PMA-Veh                  3 character
#> UL3,dH1A,dH1B : SF-DEX                   3 character
#> UL3,dH1A,dH1B : SF-PMA                   3 character
#> UL3,dH1A,dH1B : SF-Veh                   3 character

# plot the contrasts included in one particular Venn setlist
withr::with_par(list("mfrow"=c(2, 3)), {
for (n in 1:5) {
setest <- sedesign[, , comp2contrast(venn_setlists[[n]])]
plot_sedesign(setest, contrast_style="none")
}
})


# subset some groups to simplify
x <- jamba::unvigrep("Ins|SF", groups(sedesign));
sedesignsub <- validate_sedesign(sedesign, groups=x)
venn_set_comps <- contrasts_to_venn_setlists(sedesign=sedesignsub,
   contrast_style="comp",
   include_multifactor=FALSE)
withr::with_par(list("mfrow"=c(2, 3)), {
for (n in 1:6 + 6) {
setest <- sedesignsub[, , comp2contrast(venn_set_comps[[n]])]
plot_sedesign(setest, contrast_style="none")
}
})

```
