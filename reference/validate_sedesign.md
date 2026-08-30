# Validate SEDesign object contents

Validate SEDesign object contents

## Usage

``` r
validate_sedesign(
  object,
  min_reps = 1,
  samples = NULL,
  groups = NULL,
  contrasts = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- object:

  `SEDesign` object

- min_reps:

  `integer` indicating the minimum required replicate samples per design
  group to be used during analysis. Any design groups with fewer
  replicates will be removed from the design matrix, and subsequently
  will be removed from the contrasts matrix.

- samples, groups, contrasts:

  `character` vectors used to subset the samples, groups, or contrasts.

- verbose:

  `logical` indicating whether to print verbose output, where
  `verbose=TRUE` will print messages at the end of operations, and
  `verbose=2` will also print messages during operations.

- ...:

  additional arguments are ignored.

## Value

`SEDesign` object after validation updates have been applied.

## Details

This function validates and enforces constraints on `SEDesign` objects:

- `samples` must match `rownames(design)`

- `colnames(design)` must match `rownames(contrasts)`

If `samples` does not exist, and `rownames(design)` does exist, then
`samples` will be defined as `rownames(design)`.

If `design` and `samples` are provided, but `rownames(design)` is empty,
it must be the same length as `samples`. then `rownames(design)` will be
defined as `samples`.

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
[`sort_contrasts()`](https://jmw86069.github.io/jamses/reference/sort_contrasts.md)

## Examples

``` r
factors2 <- rep(c("one", "two", "three", "four"), each=3)
factors2 <- factor(factors2,
   levels=unique(factors2))
names(factors2) <- paste0("sample", seq_along(factors2))
factors2
#>  sample1  sample2  sample3  sample4  sample5  sample6  sample7  sample8 
#>      one      one      one      two      two      two    three    three 
#>  sample9 sample10 sample11 sample12 
#>    three     four     four     four 
#> Levels: one two three four

mm2 <- model.matrix(~0 + factors2)
rownames(mm2) <- names(factors2)
colnames(mm2) <- levels(factors2);
mm2
#>          one two three four
#> sample1    1   0     0    0
#> sample2    1   0     0    0
#> sample3    1   0     0    0
#> sample4    0   1     0    0
#> sample5    0   1     0    0
#> sample6    0   1     0    0
#> sample7    0   0     1    0
#> sample8    0   0     1    0
#> sample9    0   0     1    0
#> sample10   0   0     0    1
#> sample11   0   0     0    1
#> sample12   0   0     0    1
#> attr(,"assign")
#> [1] 1 1 1 1
#> attr(,"contrasts")
#> attr(,"contrasts")$factors2
#> [1] "contr.treatment"
#> 

icontrastnames <- c("two-one",
   "four-three",
   "(four-three)-(two-one)");
icon <- c(-1, 1, 0, 0,
   0, 0, -1, 1,
   1, -1, -1, 1)
icontrasts2 <- matrix(icon,
   ncol=3,
   dimnames=list(levels(factors2),
      icontrastnames))
icontrasts2
#>       two-one four-three (four-three)-(two-one)
#> one        -1          0                      1
#> two         1          0                     -1
#> three       0         -1                     -1
#> four        0          1                      1

condes2 <- SEDesign(
   design=mm2,
   contrasts=icontrasts2)
condes2
#> <SEDesign> 12 samples, 4 groups, 3 contrasts
#> factors:
#>   - factor1: one, two, three, four

# now subset samples
validate_sedesign(condes2,
   samples=paste0("sample", 2:12))
#> <SEDesign> 11 samples, 4 groups, 3 contrasts
#> factors:
#>   - factor1: one, two, three, four

# now subset enough samples to remove one group
validate_sedesign(condes2,
   samples=paste0("sample", 4:12))
#> <SEDesign> 9 samples, 3 groups, 1 contrasts
#> factors:
#>   - factor1: two, three, four

validate_sedesign(condes2, groups=c("one", "two"))
#> <SEDesign> 12 samples, 2 groups, 1 contrasts
#> factors:
#>   - factor1: one, two

condes2[, c("one", "two"), ]
#> <SEDesign> 12 samples, 2 groups, 1 contrasts
#> factors:
#>   - factor1: one, two
```
