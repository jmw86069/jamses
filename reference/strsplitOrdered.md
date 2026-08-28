# Split the elements of an ordered factor vector

Split the elements of an ordered factor vector

## Usage

``` r
strsplitOrdered(
  x,
  split = "_",
  fixed = FALSE,
  perl = FALSE,
  useBytes = FALSE,
  sortFunc = jamba::mixedSort,
  keepOrder = TRUE,
  ...
)
```

## Arguments

- x:

  character or factor vector.

- split:

  character split value sent to
  [`base::strsplit()`](https://rdrr.io/r/base/strsplit.html).

- fixed, perl, useBytes:

  additional arguments sent to
  [`base::split()`](https://rdrr.io/r/base/split.html).

- sortFunc:

  function used to sort character values when the input `x` is a
  character vector. The default
  [`jamba::mixedSort()`](https://jmw86069.github.io/jamba/reference/mixedSort.html)
  applies alphanumeric sort.

- keepOrder:

  logical indicating whether to keep the order of values in the input
  data, for example with character input the values will be ordered by
  the first appearance of each term.

- ...:

  additional arguments are ignored.

## Value

list of factor vectors, where each factor shares the same global factor
levels based upon the input data.

## Details

This function performs
[`base::strsplit()`](https://rdrr.io/r/base/strsplit.html) while trying
to maintain the order of factor levels in the output, based upon the
order of factor levels in the input data.

## See also

Other jamses utilities:
[`choose_annotation_colnames()`](https://jmw86069.github.io/jamses/reference/choose_annotation_colnames.md),
[`combine_sestats()`](https://jmw86069.github.io/jamses/reference/combine_sestats.md),
[`contrast2comp_dev()`](https://jmw86069.github.io/jamses/reference/contrast2comp_dev.md),
[`fold_to_log2fold()`](https://jmw86069.github.io/jamses/reference/fold_to_log2fold.md),
[`intercalate()`](https://jmw86069.github.io/jamses/reference/intercalate.md),
[`list2im_opt()`](https://jmw86069.github.io/jamses/reference/list2im_opt.md),
[`list_to_sestats()`](https://jmw86069.github.io/jamses/reference/list_to_sestats.md),
[`log2fold_to_fold()`](https://jmw86069.github.io/jamses/reference/log2fold_to_fold.md),
[`make_block_arrow_polygon()`](https://jmw86069.github.io/jamses/reference/make_block_arrow_polygon.md),
[`mark_stat_hits()`](https://jmw86069.github.io/jamses/reference/mark_stat_hits.md),
[`matrix_normalize()`](https://jmw86069.github.io/jamses/reference/matrix_normalize.md),
[`merge_statdf_all_test()`](https://jmw86069.github.io/jamses/reference/merge_statdf_all_test.md),
[`point_handedness()`](https://jmw86069.github.io/jamses/reference/point_handedness.md),
[`point_slope_intercept()`](https://jmw86069.github.io/jamses/reference/point_slope_intercept.md),
[`shortest_unique_abbreviation()`](https://jmw86069.github.io/jamses/reference/shortest_unique_abbreviation.md),
[`shrinkDataFrame()`](https://jmw86069.github.io/jamses/reference/shrinkDataFrame.md),
[`shrink_df()`](https://jmw86069.github.io/jamses/reference/shrink_df.md),
[`shrink_matrix()`](https://jmw86069.github.io/jamses/reference/shrink_matrix.md),
[`sort_samples()`](https://jmw86069.github.io/jamses/reference/sort_samples.md),
[`sub_split_vector()`](https://jmw86069.github.io/jamses/reference/sub_split_vector.md),
[`update_function_params()`](https://jmw86069.github.io/jamses/reference/update_function_params.md),
[`update_list_elements()`](https://jmw86069.github.io/jamses/reference/update_list_elements.md)

## Examples

``` r
# first define a vector of sample groups
iGroups <- jamba::nameVector(paste(rep(c("WT", "KO"), each=6),
   rep(c("Control", "Treated"), each=3),
   sep="_"));
iGroups <- factor(iGroups, levels=unique(iGroups));
iGroups;
#> WT_Control_v1 WT_Control_v2 WT_Control_v3 WT_Treated_v1 WT_Treated_v2 
#>    WT_Control    WT_Control    WT_Control    WT_Treated    WT_Treated 
#> WT_Treated_v3 KO_Control_v1 KO_Control_v2 KO_Control_v3 KO_Treated_v1 
#>    WT_Treated    KO_Control    KO_Control    KO_Control    KO_Treated 
#> KO_Treated_v2 KO_Treated_v3 
#>    KO_Treated    KO_Treated 
#> Levels: WT_Control WT_Treated KO_Control KO_Treated
strsplitOrdered(iGroups, "_");
#> [[1]]
#> [1] WT      Control
#> Levels: WT KO Control Treated
#> 
#> [[2]]
#> [1] WT      Control
#> Levels: WT KO Control Treated
#> 
#> [[3]]
#> [1] WT      Control
#> Levels: WT KO Control Treated
#> 
#> [[4]]
#> [1] WT      Treated
#> Levels: WT KO Control Treated
#> 
#> [[5]]
#> [1] WT      Treated
#> Levels: WT KO Control Treated
#> 
#> [[6]]
#> [1] WT      Treated
#> Levels: WT KO Control Treated
#> 
#> [[7]]
#> [1] KO      Control
#> Levels: WT KO Control Treated
#> 
#> [[8]]
#> [1] KO      Control
#> Levels: WT KO Control Treated
#> 
#> [[9]]
#> [1] KO      Control
#> Levels: WT KO Control Treated
#> 
#> [[10]]
#> [1] KO      Treated
#> Levels: WT KO Control Treated
#> 
#> [[11]]
#> [1] KO      Treated
#> Levels: WT KO Control Treated
#> 
#> [[12]]
#> [1] KO      Treated
#> Levels: WT KO Control Treated
#> 
```
