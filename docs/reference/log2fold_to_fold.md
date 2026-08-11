# Convert log2 fold change to signed fold change

Convert log2 fold change to signed fold change

## Usage

``` r
log2fold_to_fold(x, ...)
```

## Arguments

- x:

  `numeric` vector

## Value

`numeric` vector representing signed fold change values.

## Details

This function takes log2 fold change values as input, and returns normal
space fold change values that retain the positive and negative sign, and
the magnitude.

For example:

- `log2 fold change = 2` becomes `fold change = 4`.

- `log2 fold change = -2` becomes `fold change = -4`.

This function therefore differs from similar functions that convert log2
fold change into a ratio. Instead, `log2fold_to_fold()` specifically
retains the magnitude of negative changes.

## See also

Other jamses utilities:
[`choose_annotation_colnames()`](choose_annotation_colnames.md),
[`combine_sestats()`](combine_sestats.md),
[`contrast2comp_dev()`](contrast2comp_dev.md),
[`fold_to_log2fold()`](fold_to_log2fold.md),
[`intercalate()`](intercalate.md), [`list2im_opt()`](list2im_opt.md),
[`list_to_sestats()`](list_to_sestats.md),
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
[`sub_split_vector()`](sub_split_vector.md),
[`update_function_params()`](update_function_params.md),
[`update_list_elements()`](update_list_elements.md)

## Examples

``` r
x <- c(-3, -2, -1, 0, 1, 2, 3);
fc <- log2fold_to_fold(x);
fc;
#> [1] -8 -4 -2  1  2  4  8

fold_to_log2fold(fc);
#> [1] -3 -2 -1  0  1  2  3
```
