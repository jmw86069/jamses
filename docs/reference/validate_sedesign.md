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
[`sort_contrasts()`](sort_contrasts.md)

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
#> <SEDesign>
#>  @ design      : num [1:12, 1:4] 1 1 1 0 0 0 0 0 0 0 ...
#>  .. - attr(*, "dimnames")=List of 2
#>  ..  ..$ : chr [1:12] "sample1" "sample2" "sample3" "sample4" ...
#>  ..  ..$ : chr [1:4] "one" "two" "three" "four"
#>  .. - attr(*, "assign")= int [1:4] 1 1 1 1
#>  .. - attr(*, "contrasts")=List of 1
#>  ..  ..$ factors2: chr "contr.treatment"
#>  @ contrasts   : num [1:4, 1:3] -1 1 0 0 0 0 -1 1 1 -1 ...
#>  .. - attr(*, "dimnames")=List of 2
#>  ..  ..$ : chr [1:4] "one" "two" "three" "four"
#>  ..  ..$ : chr [1:3] "two-one" "four-three" "(four-three)-(two-one)"
#>  @ samples     : chr(0) 
#>  @ factors     : chr "factor1"
#>  @ design_df   :'data.frame':    4 obs. of  1 variable:
#>  .. $ factor1: chr  "one" "two" "three" "four"
#>  @ contrasts_df:'data.frame':    3 obs. of  1 variable:
#>  .. $ factor1: chr  "two-one" "four-three" "(four-three)-(two-one)"

# now subset samples
validate_sedesign(condes2,
   samples=paste0("sample", 2:12))
#> <SEDesign>
#>  @ design      : num [1:11, 1:4] 1 1 0 0 0 0 0 0 0 0 ...
#>  .. - attr(*, "dimnames")=List of 2
#>  ..  ..$ : chr [1:11] "sample2" "sample3" "sample4" "sample5" ...
#>  ..  ..$ : chr [1:4] "one" "two" "three" "four"
#>  @ contrasts   : num [1:4, 1:3] -1 1 0 0 0 0 -1 1 1 -1 ...
#>  .. - attr(*, "dimnames")=List of 2
#>  ..  ..$ : chr [1:4] "one" "two" "three" "four"
#>  ..  ..$ : chr [1:3] "two-one" "four-three" "(four-three)-(two-one)"
#>  @ samples     : chr [1:11] "sample2" "sample3" "sample4" "sample5" "sample6" ...
#>  @ factors     : chr "factor1"
#>  @ design_df   :'data.frame':    4 obs. of  1 variable:
#>  .. $ factor1: chr  "one" "two" "three" "four"
#>  @ contrasts_df:'data.frame':    3 obs. of  1 variable:
#>  .. $ factor1: chr  "two-one" "four-three" "(four-three)-(two-one)"

# now subset enough samples to remove one group
validate_sedesign(condes2,
   samples=paste0("sample", 4:12))
#> <SEDesign>
#>  @ design      : num [1:9, 1:3] 1 1 1 0 0 0 0 0 0 0 ...
#>  .. - attr(*, "dimnames")=List of 2
#>  ..  ..$ : chr [1:9] "sample4" "sample5" "sample6" "sample7" ...
#>  ..  ..$ : chr [1:3] "two" "three" "four"
#>  @ contrasts   : num [1:3, 1] 0 -1 1
#>  .. - attr(*, "dimnames")=List of 2
#>  ..  ..$ : chr [1:3] "two" "three" "four"
#>  ..  ..$ : chr "four-three"
#>  @ samples     : chr [1:9] "sample4" "sample5" "sample6" "sample7" "sample8" "sample9" ...
#>  @ factors     : chr "factor1"
#>  @ design_df   :'data.frame':    3 obs. of  1 variable:
#>  .. $ factor1: chr  "two" "three" "four"
#>  @ contrasts_df:'data.frame':    1 obs. of  1 variable:
#>  .. $ factor1: chr "four-three"

validate_sedesign(condes2, groups=c("one", "two"))
#> <SEDesign>
#>  @ design      : num [1:12, 1:2] 1 1 1 0 0 0 0 0 0 0 ...
#>  .. - attr(*, "dimnames")=List of 2
#>  ..  ..$ : chr [1:12] "sample1" "sample2" "sample3" "sample4" ...
#>  ..  ..$ : chr [1:2] "one" "two"
#>  @ contrasts   : num [1:2, 1] -1 1
#>  .. - attr(*, "dimnames")=List of 2
#>  ..  ..$ : chr [1:2] "one" "two"
#>  ..  ..$ : chr "two-one"
#>  @ samples     : chr [1:12] "sample1" "sample2" "sample3" "sample4" "sample5" ...
#>  @ factors     : chr "factor1"
#>  @ design_df   :'data.frame':    2 obs. of  1 variable:
#>  .. $ factor1: chr  "one" "two"
#>  @ contrasts_df:'data.frame':    1 obs. of  1 variable:
#>  .. $ factor1: chr "two-one"

condes2[, c("one", "two"), ]
#> <SEDesign>
#>  @ design      : num [1:12, 1:2] 1 1 1 0 0 0 0 0 0 0 ...
#>  .. - attr(*, "dimnames")=List of 2
#>  ..  ..$ : chr [1:12] "sample1" "sample2" "sample3" "sample4" ...
#>  ..  ..$ : chr [1:2] "one" "two"
#>  @ contrasts   : num [1:2, 1] -1 1
#>  .. - attr(*, "dimnames")=List of 2
#>  ..  ..$ : chr [1:2] "one" "two"
#>  ..  ..$ : chr "two-one"
#>  @ samples     : chr [1:12] "sample1" "sample2" "sample3" "sample4" "sample5" ...
#>  @ factors     : chr "factor1"
#>  @ design_df   :'data.frame':    2 obs. of  1 variable:
#>  .. $ factor1: chr  "one" "two"
#>  @ contrasts_df:'data.frame':    1 obs. of  1 variable:
#>  .. $ factor1: chr "two-one"
```
