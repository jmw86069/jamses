# Sub-split a split vector (Internal)

Sub-split a split vector (Internal)

## Usage

``` r
sub_split_vector(icolsplit, max_size = 4, ...)
```

## Arguments

- icolsplit:

  `numeric` vector used to split another vector

- max_size:

  `numeric` max number of identical values permitted in `icolsplit`.

- ...:

  additional arguments are ignored

## Value

`factor` vector suitable to use for splitting a vector into subsets, in
order.

## Details

This function is used to support
[`contrasts_to_venn_setlists()`](https://jmw86069.github.io/jamses/reference/contrasts_to_venn_setlists.md)
in order to limit the number of Venn sets included by subgroup. This
function adds another level of splitting, while keeping the elements in
the original order after splitting.

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
[`strsplitOrdered()`](https://jmw86069.github.io/jamses/reference/strsplitOrdered.md),
[`update_function_params()`](https://jmw86069.github.io/jamses/reference/update_function_params.md),
[`update_list_elements()`](https://jmw86069.github.io/jamses/reference/update_list_elements.md)

## Examples

``` r
icolsplit <- rep(c(1, 2, 3), c(6, 5, 4))
newsplit <- sub_split_vector(icolsplit)
split(icolsplit, newsplit)
#> $`1.1`
#> [1] 1 1 1 1
#> 
#> $`1.2`
#> [1] 1 1
#> 
#> $`2.1`
#> [1] 2 2 2
#> 
#> $`2.2`
#> [1] 2 2
#> 
#> $`3`
#> [1] 3 3 3 3
#> 

newsplit3 <- sub_split_vector(icolsplit, max_size=3)
split(icolsplit, newsplit3)
#> $`1.1`
#> [1] 1 1 1
#> 
#> $`1.2`
#> [1] 1 1 1
#> 
#> $`2.1`
#> [1] 2 2 2
#> 
#> $`2.2`
#> [1] 2 2
#> 
#> $`3.1`
#> [1] 3 3
#> 
#> $`3.2`
#> [1] 3 3
#> 
```
