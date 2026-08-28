# Update a subset of list elements

Update a subset of list elements

## Usage

``` r
update_list_elements(
  source_list,
  update_list,
  list_layer_num = 1,
  verbose = FALSE,
  ...
)
```

## Details

This function is intended to help update a nested `source_list`, a
subset of whose values should be replaced with entries in `update_list`,
leaving any original entries in `source_list` which were not defined in
`update_list`.

This function may be useful when manipulating lattice or ggplot2
graphical parameters, which are often stored in a nested list structure.

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
[`sub_split_vector()`](https://jmw86069.github.io/jamses/reference/sub_split_vector.md),
[`update_function_params()`](https://jmw86069.github.io/jamses/reference/update_function_params.md)
