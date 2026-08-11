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
[`contrasts_to_venn_setlists()`](contrasts_to_venn_setlists.md) in order
to limit the number of Venn sets included by subgroup. This function
adds another level of splitting, while keeping the elements in the
original order after splitting.

## See also

Other jamses utilities:
[`choose_annotation_colnames()`](choose_annotation_colnames.md),
[`combine_sestats()`](combine_sestats.md),
[`contrast2comp_dev()`](contrast2comp_dev.md),
[`fold_to_log2fold()`](fold_to_log2fold.md),
[`intercalate()`](intercalate.md), [`list2im_opt()`](list2im_opt.md),
[`list_to_sestats()`](list_to_sestats.md),
[`log2fold_to_fold()`](log2fold_to_fold.md),
[`make_block_arrow_polygon()`](make_block_arrow_polygon.md),
[`mark_stat_hits()`](mark_stat_hits.md),
[`matrix_normalize()`](matrix_normalize.md),
[`merge_statdf_all_test()`](merge_statdf_all_test.md),
[`point_handedness()`](point_handedness.md),
[`point_slope_intercept()`](point_slope_intercept.md),
[`shortest_unique_abbreviation()`](shortest_unique_abbreviation.md),
[`shrinkDataFrame()`](shrinkDataFrame.md),
[`shrink_df()`](shrink_df.md), [`shrink_matrix()`](shrink_matrix.md),
[`sort_samples()`](sort_samples.md),
[`strsplitOrdered()`](strsplitOrdered.md),
[`update_function_params()`](update_function_params.md),
[`update_list_elements()`](update_list_elements.md)

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
