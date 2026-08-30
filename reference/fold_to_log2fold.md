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

- fold change, as from
  [`log2fold_to_fold()`](https://jmw86069.github.io/jamses/reference/log2fold_to_fold.md)
  which includes no values between 0 and 1, but may include negative
  values.

For example, for ratio input:

- `ratio = 4` becomes `log2 fold change = 2`.

- `ratio = 0.25` becomes `log2 fold change = -2`.

For example, for fold change input:

- `fold change = 4` becomes `log2 fold change = 2`.

- `fold change = -4` becomes `log2 fold change = -2`.

## See also

Other jamses utilities:
[`choose_annotation_colnames()`](https://jmw86069.github.io/jamses/reference/choose_annotation_colnames.md),
[`combine_sestats()`](https://jmw86069.github.io/jamses/reference/combine_sestats.md),
[`contrast2comp_dev()`](https://jmw86069.github.io/jamses/reference/contrast2comp_dev.md),
[`intercalate()`](https://jmw86069.github.io/jamses/reference/intercalate.md),
[`list2im_opt()`](https://jmw86069.github.io/jamses/reference/list2im_opt.md),
[`list2im_value_internal()`](https://jmw86069.github.io/jamses/reference/list2im_value_internal.md),
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
[`update_function_params()`](https://jmw86069.github.io/jamses/reference/update_function_params.md),
[`update_list_elements()`](https://jmw86069.github.io/jamses/reference/update_list_elements.md)

## Examples

``` r
x <- c(-3, -2, -1, 0, 1, 2, 3);
fc <- log2fold_to_fold(x);
fc;
#> [1] -8 -4 -2  1  2  4  8

fold_to_log2fold(fc);
#> [1] -3 -2 -1  0  1  2  3
```
