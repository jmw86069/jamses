# Merge stats data.frame from two sestats results

Merge stats `data.frame` rows from two sestats results, specifically
when one object has "all rows" and is intended only to produce the
non-P-value statistics such as group mean, max group mean (mgm), fold,
and logFC.

## Usage

``` r
merge_statdf_all_test(
  sestats_all,
  sestats_test,
  exclude_from_all = "P.val|^hit ",
  ...
)
```

## Arguments

- sestats_all:

  `SEStats` object or equivalent legacy `list` format, as returned by
  [`se_contrast_stats()`](https://jmw86069.github.io/jamses/reference/se_contrast_stats.md)
  when using argument `define_hits=FALSE`. These results should
  therefore only contain the group mean and log fold change values, for
  all rows in the input data, regardless if they were tested in
  `sestts_test`.

- sestats_test:

  `SEStats` object or equivalent legacy `list` format, as returned by
  [`se_contrast_stats()`](https://jmw86069.github.io/jamses/reference/se_contrast_stats.md)
  when using argument `define_hits=TRUE`. These results should therefore
  contain the 'hit' column and with 'hit_array' data populated. It
  should only contain rows that were tested, which is assumed to be a
  subset of rows contained in `sestats_all`.

- exclude_from_all:

  `character` pattern used with
  [`jamba::unvigrep()`](https://jmw86069.github.io/jamba/reference/unvigrep.html)
  to exclude colnames from each `data.frame` in the 'stats_dfs' property
  of `sestats_all`. The default `'P.val|^hit '` will exclude P-value and
  adjusted P-value columns if they exist (they should not, but it can be
  tolerated).

- ...:

  additional arguments are ignored.

## Value

`list` object in the same format as `sestats_test`, with only the
elements 'stats_df' and 'stats_dfs' replaced with the expanded complete
set of data.

## Details

- `sestats_all` should have all rows, not just those tested. It will
  specifically exclude stats columns with `'P.val'` or `'^hit '` in the
  name, using argument `exclude_from_all`.

- `sestats_test` should have only the subset of rows tested.

This function may be used internally within jamses to automate the job
of testing a subset of rows, but reporting summary values for all rows -
including those rows not formally tested.

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
