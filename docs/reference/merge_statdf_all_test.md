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
  [`se_contrast_stats()`](se_contrast_stats.md) which used "all rows" of
  the source data, along with `define_hits=FALSE`.

- sestats_test:

  `list` with sestats data from
  [`se_contrast_stats()`](se_contrast_stats.md) which used a subset of
  rows, usually "detected rows" for the statistical test.

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
