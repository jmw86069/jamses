# Shrink data.frame by row groups

Shrink data.frame by row groups

## Usage

``` r
shrinkDataFrame(
  x,
  groupBy,
  na.rm = TRUE,
  string_func = function(x) jamba::cPasteSU(x, na.rm = TRUE),
  num_func = function(x) {
     mean(x, na.rm = TRUE)
 },
  add_string_cols = NULL,
  num_to_string_func = as.character,
  keep_na_groups = TRUE,
  include_num_reps = FALSE,
  collapse_method = 2,
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  `data.frame` (or equivalent)

- groupBy:

  `character` vector with one of the following:

  - one or more columns in `colnames(x)`. The values in these columns
    will define the row groups used.

  - `character` or `factor` with length equal to `nrow(x)`. These values
    will define the row groups used.

- string_func:

  `function`, default uses
  [`jamba::cPasteSU()`](https://jmw86069.github.io/jamba/reference/cPaste.html),
  used for `character` or `factor` columns. Note that string columns are
  handled differently than `numeric` columns by applying vectorized
  operations across the complete set of rows in one step, rather than
  calling `data.table` on each subgroup.

- num_func:

  `function`, default `function(x)mean(x, na.rm=TRUE)`, used for
  `numeric` columns. Note that this function is applied to each row
  group by `data.table`, and is typically very efficient for `numeric`
  values.

- add_string_cols:

  `character` with optional `numeric` columns that should be handled as
  if they were `character` columns. Default `NULL`.

- num_to_string_func:

  `function` used for `add_string_cols` when converting `numeric`
  columns to `character`. Default
  [`as.character()`](https://rdrr.io/r/base/character.html) retains the
  full `numeric` value, however it may be useful to use something like
  `function(x)signif(x, digits=3)` to limit the output to only three
  significant digits, or `function(x)format(x, digits=3)`.

- keep_na_groups:

  `logical`, default TRUE, whether to convert `NA` values in row groups
  to `""` so they are retained in the output.

  - You may want to use `keep_na_groups=FALSE` when there are a large
    number of un-annotated rows that should not be aggregated together.
    This situation may occur if converting a probe to a gene symbol,
    where a subset of probes cannot be converted to a gene symbol and
    instead receive `NA`.

- include_num_reps:

  `logical` indicating whether to add a column `"num_reps"` to the
  output, with the `integer` number of rows in each row group.

- collapse_method:

  `integer` default 2, indicating the internal collapse method used.
  Experimental.

  - `1` collapses each `numeric` column independently.

  - `2` collapses each set of `numeric` columns that use the same
    numeric shrink function. When all `numeric` columns use the same
    shrink function, they are all calculated in a single step, which is
    typically much faster.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are ignored.

## Details

Purpose is to shrink a `data.frame` to have one row per row grouping.
The row grouping can use a single column of identifiers, or multiple
columns. The challenge is to apply a relevant function to each column,
expecting there will be columns with `numeric`, `character`, or `factor`
types.

The default behavior:

- `numeric` columns are summarized with `mean(x, na.rm=TRUE)`, so that
  NA values are ignored when there are non-NA values present.

- `character` columns are combined using unique, sorted `character`
  strings.

  - This step uses
    [`jamba::cPasteSU()`](https://jmw86069.github.io/jamba/reference/cPaste.html)
    where the `S` activates sorting using
    [`jamba::mixedSort()`](https://jmw86069.github.io/jamba/reference/mixedSort.html),
    and `U` calls [`unique()`](https://rdrr.io/r/base/unique.html).

  - To retain all values, remove the `U` and call
    [`jamba::cPasteS()`](https://jmw86069.github.io/jamba/reference/cPaste.html)

  - To skip the sort, remove the `S` and call
    [`jamba::cPasteU()`](https://jmw86069.github.io/jamba/reference/cPaste.html)

  - To keep all values, and skip sorting, call
    [`jamba::cPaste()`](https://jmw86069.github.io/jamba/reference/cPaste.html)

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
[`shrink_df()`](https://jmw86069.github.io/jamses/reference/shrink_df.md),
[`shrink_matrix()`](https://jmw86069.github.io/jamses/reference/shrink_matrix.md),
[`sort_samples()`](https://jmw86069.github.io/jamses/reference/sort_samples.md),
[`strsplitOrdered()`](https://jmw86069.github.io/jamses/reference/strsplitOrdered.md),
[`sub_split_vector()`](https://jmw86069.github.io/jamses/reference/sub_split_vector.md),
[`update_function_params()`](https://jmw86069.github.io/jamses/reference/update_function_params.md),
[`update_list_elements()`](https://jmw86069.github.io/jamses/reference/update_list_elements.md)

## Examples

``` r
testdf <- data.frame(check.names=FALSE,
   SYMBOL=rep(c("ACTB", "GAPDH", "PPIA"), c(2, 3, 1)),
   `logFC B-A`=c(1.4, 1.4, 2.3, NA, 2.5, 5.1),
   probe=paste0("probe", 1:6))
shrink_df(testdf, by="SYMBOL")
#>       SYMBOL logFC B-A                probe
#> ACTB    ACTB       1.4        probe1,probe2
#> GAPDH  GAPDH       2.4 probe3,probe4,probe5
#> PPIA    PPIA       5.1               probe6

shrink_df(testdf, by="SYMBOL", num_func=mean)
#>       SYMBOL logFC B-A                probe
#> ACTB    ACTB       1.4        probe1,probe2
#> GAPDH  GAPDH        NA probe3,probe4,probe5
#> PPIA    PPIA       5.1               probe6

shrink_df(testdf, by="SYMBOL", add_string_cols="logFC B-A")
#>       SYMBOL logFC B-A                probe
#> ACTB    ACTB       1.4        probe1,probe2
#> GAPDH  GAPDH   2.3,2.5 probe3,probe4,probe5
#> PPIA    PPIA       5.1               probe6

testdftall <- do.call(rbind, lapply(1:10000, function(i){
   idf <- testdf;
   idf$SYMBOL <- paste0(idf$SYMBOL, "_", i);
   idf;
}))
shrunk_tall <- shrink_df(testdftall,
   by="SYMBOL")
head(shrunk_tall, 6)
#>          SYMBOL logFC B-A                probe
#> ACTB_1   ACTB_1       1.4        probe1,probe2
#> GAPDH_1 GAPDH_1       2.4 probe3,probe4,probe5
#> PPIA_1   PPIA_1       5.1               probe6
#> ACTB_2   ACTB_2       1.4        probe1,probe2
#> GAPDH_2 GAPDH_2       2.4 probe3,probe4,probe5
#> PPIA_2   PPIA_2       5.1               probe6

shrunk_tall2 <- jamses::shrinkDataFrame(testdftall,
   groupBy="SYMBOL")
head(shrunk_tall2, 6)
#>          SYMBOL logFC B-A                probe
#> ACTB_1   ACTB_1       1.4        probe1,probe2
#> GAPDH_1 GAPDH_1       2.4 probe3,probe4,probe5
#> PPIA_1   PPIA_1       5.1               probe6
#> ACTB_2   ACTB_2       1.4        probe1,probe2
#> GAPDH_2 GAPDH_2       2.4 probe3,probe4,probe5
#> PPIA_2   PPIA_2       5.1               probe6
```
