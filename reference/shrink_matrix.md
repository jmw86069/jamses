# Shrink a numeric matrix by groups of rows

Shrink a numeric matrix across groups of rows, by applying a summary
function.

## Usage

``` r
shrink_matrix(
  x,
  groupBy,
  shrink_func = function(x) {
     mean(x, na.rm = TRUE)
 },
  return_class = c("data.frame", "matrix"),
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  `numeric` matrix

- groupBy:

  `character` or `factor` vector of group labels, whose length equals
  `nrow(x)`. These values will become rownames in the output data.

- shrink_func:

  `function` that takes vector input and returns **single value**
  output. The vector class can be checked, in order to call a function
  on numeric or character data separately, as needed.

- return_class:

  `character` string indicating the return data type.

  - `"data.frame"` returns a `data.frame` whose first column contains
    entries from `groups`.

  - `"matrix"` returns a numeric matrix whose rownames are entries from
    `groups`.

- verbose:

  logical indicating whether to print verbose output.

## Value

`data.frame` or `matrix` based upon argument `return_class`.

## Details

This function is mainly a wrapper to use the amazingly fast `data.table`
package, with the ability to provide a custom function to shrink row
values.

The default function uses `mean(x, na.rm=TRUE)` so that `NA` values are
ignored where possible.

This function applies the same `shrink_func` to all columns, and it
optimal for `numeric` values. For more control over which function to
apply to specific columns, see
[`shrinkDataFrame()`](https://jmw86069.github.io/jamses/reference/shrinkDataFrame.md).

Trivia: This function is identical to `splicejam::shrinkDataFrame()`
except that the default `shrink_func` includes `na.rm=TRUE` and no
longer calls the [`.Internal()`](https://rdrr.io/r/base/Internal.html)
function, since that is not permitted by CRAN package guidelines.

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
[`sort_samples()`](https://jmw86069.github.io/jamses/reference/sort_samples.md),
[`strsplitOrdered()`](https://jmw86069.github.io/jamses/reference/strsplitOrdered.md),
[`sub_split_vector()`](https://jmw86069.github.io/jamses/reference/sub_split_vector.md),
[`update_function_params()`](https://jmw86069.github.io/jamses/reference/update_function_params.md),
[`update_list_elements()`](https://jmw86069.github.io/jamses/reference/update_list_elements.md)
