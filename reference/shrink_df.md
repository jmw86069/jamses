# Shrink data.frame by row groups

Shrink data.frame by row groups

## Usage

``` r
shrink_df(
  x,
  by,
  string_func = function(x) jamba::cPasteSU(x, na.rm = TRUE),
  num_func = function(x) mean(x, na.rm = TRUE),
  add_string_cols = NULL,
  num_to_string_func = as.character,
  keep_na_groups = TRUE,
  include_num_reps = FALSE,
  extra_funcs = NULL,
  do_test = FALSE,
  use_new_method = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- by:

  `character` vector of one or more `colnames(df)`, used to define the
  row grouping.

- string_func:

  `function` used for `character` and other non-numeric column types.
  For efficiency, `string_func` by default is applied to the entire
  column, with `list` input, expecting `vector` output. It is not
  applied using `data.table`.

- num_func:

  `function` used for `numeric` column types. This function is applied
  using `data.table` and should expect a `vector` input, and provide a
  single atomic value output.

- extra_funcs:

  `list`, default `NULL`, containing `function` objects. The list names
  should match `colnames(x)`, in order to apply a function to a specific
  column in `x`. These functions will therefore override the default
  functions defined by `string_func` and `num_func`. Only one function
  is applied per column.

- do_test:

  `logical`, default FALSE, indicating whether to perform an internal
  test with internally-generated argument values.

- use_new_method:

  `logical` default FALSE, whether to call newer tidy/data.table methods
  (TRUE), or call
  [`shrinkDataFrame()`](https://jmw86069.github.io/jamses/reference/shrinkDataFrame.md)
  (FALSE). Currently
  [`shrinkDataFrame()`](https://jmw86069.github.io/jamses/reference/shrinkDataFrame.md)
  is remarkably faster. More research necessary.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are ignored.

- df:

  `data.frame` or compatible input class.

## Details

This function is currently a wrapper for
[`shrinkDataFrame()`](https://jmw86069.github.io/jamses/reference/shrinkDataFrame.md),
it was formerly a simplified version of
[`shrinkDataFrame()`](https://jmw86069.github.io/jamses/reference/shrinkDataFrame.md)
which is intended to use more modern methods from the R package
`data.table`.

The general idea is to collapse `numeric` columns using `num_func`, and
collapse `character` and all other columns using `string_func`.

Any exceptions, where a different function should be applied, are passed
via argument `extra_funcs` which is a `list` of functions named by
values in `colnames(df)`.

## See also

Other jamses utilities:
[`choose_annotation_colnames()`](https://jmw86069.github.io/jamses/reference/choose_annotation_colnames.md),
[`combine_sestats()`](https://jmw86069.github.io/jamses/reference/combine_sestats.md),
[`contrast2comp_dev()`](https://jmw86069.github.io/jamses/reference/contrast2comp_dev.md),
[`fold_to_log2fold()`](https://jmw86069.github.io/jamses/reference/fold_to_log2fold.md),
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
[`shrink_matrix()`](https://jmw86069.github.io/jamses/reference/shrink_matrix.md),
[`sort_samples()`](https://jmw86069.github.io/jamses/reference/sort_samples.md),
[`strsplitOrdered()`](https://jmw86069.github.io/jamses/reference/strsplitOrdered.md),
[`sub_split_vector()`](https://jmw86069.github.io/jamses/reference/sub_split_vector.md),
[`update_function_params()`](https://jmw86069.github.io/jamses/reference/update_function_params.md),
[`update_list_elements()`](https://jmw86069.github.io/jamses/reference/update_list_elements.md)
