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
  (TRUE), or call [`shrinkDataFrame()`](shrinkDataFrame.md) (FALSE).
  Currently [`shrinkDataFrame()`](shrinkDataFrame.md) is remarkably
  faster. More research necessary.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are ignored.

- df:

  `data.frame` or compatible input class.

## Details

This function is currently a wrapper for
[`shrinkDataFrame()`](shrinkDataFrame.md), it was formerly a simplified
version of [`shrinkDataFrame()`](shrinkDataFrame.md) which is intended
to use more modern methods from the R package `data.table`.

The general idea is to collapse `numeric` columns using `num_func`, and
collapse `character` and all other columns using `string_func`.

Any exceptions, where a different function should be applied, are passed
via argument `extra_funcs` which is a `list` of functions named by
values in `colnames(df)`.

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
[`shrink_matrix()`](shrink_matrix.md),
[`sort_samples()`](sort_samples.md),
[`strsplitOrdered()`](strsplitOrdered.md),
[`sub_split_vector()`](sub_split_vector.md),
[`update_function_params()`](update_function_params.md),
[`update_list_elements()`](update_list_elements.md)
