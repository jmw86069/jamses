# Convert normal signed fold change to log2 fold change

Convert normal signed fold change to log2 fold change

## Usage

``` r
fold_to_log2fold(x, ...)
```

## Arguments

- x:

  `numeric` vector

## Value

`numeric` vector representing log2 fold change values.

## Details

This function takes fold change values as input, and returns log2 fold
change values.

This function recognizes two forms of input:

- ratio, which includes values between 0 and 1, but no negative values;

- fold change, as from [`log2fold_to_fold()`](log2fold_to_fold.md) which
  includes no values between 0 and 1, but may include negative values.

For example, for ratio input:

- `ratio = 4` becomes `log2 fold change = 2`.

- `ratio = 0.25` becomes `log2 fold change = -2`.

For example, for fold change input:

- `fold change = 4` becomes `log2 fold change = 2`.

- `fold change = -4` becomes `log2 fold change = -2`.

## See also

Other jamses utilities:
[`choose_annotation_colnames()`](choose_annotation_colnames.md),
[`combine_sestats()`](combine_sestats.md),
[`contrast2comp_dev()`](contrast2comp_dev.md),
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
