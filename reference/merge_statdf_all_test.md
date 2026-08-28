# Merge stats data.frame from two sestats results

Merge stats data.frame from two sestats results, specifically when one
object has "all rows" and is intended only to produce the non-P-value
statistics such as group mean, max group mean (mgm), fold, and logFC.

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

  `list` with sestats data from
  [`se_contrast_stats()`](https://jmw86069.github.io/jamses/reference/se_contrast_stats.md)
  which used "all rows" of the source data, along with
  `define_hits=FALSE`.

- sestats_test:

  `list` with sestats data from
  [`se_contrast_stats()`](https://jmw86069.github.io/jamses/reference/se_contrast_stats.md)
  which used a subset of rows, usually "detected rows" for the
  statistical test.

- exclude_from_all:

  `character` pattern used with
  [`jamba::unvigrep()`](https://jmw86069.github.io/jamba/reference/unvigrep.html)
  to exclude colnames from `sestats_all`. Default `'P.val|^hit '`
  excludes P-values and hit columns.

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
