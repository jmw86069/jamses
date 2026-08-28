# Combine SEStats objects by contrast,signal

Combine SEStats objects by contrast,signal

## Usage

``` r
combine_sestats(sestats_list, create_stats_df = TRUE, ...)
```

## Arguments

- sestats_list:

  `list` with two or more `list` objects in sestats format, as returned
  by
  [`se_contrast_stats()`](https://jmw86069.github.io/jamses/reference/se_contrast_stats.md).

- create_stats_df:

  `logical`, default TRUE, whether to create the 'stats_df' `data.frame`
  object by merging corresponding 'stats_dfs' `data.frame` objects
  together. This step assumes that columns of class "character" or
  "factor" will have identical values across all `data.frame` objects.
  Specifically, it assumes that at least one column has the same name,
  with row identifiers sufficient to merge each result. This assumption
  also implies that one should not combine sestats objects when the rows
  are not equivalent.

- ...:

  additional arguments are ignored.

## Value

`list` equivalent to output from
[`se_contrast_stats()`](https://jmw86069.github.io/jamses/reference/se_contrast_stats.md)
with minor exceptions:

- List elements 'idesign','icontrasts','normgroup' are each returned as
  a `list` with corresponding data from each sestats object. The reason
  is that these values cannot logically be combined into one new object,
  and therefore the next level of support is to retain the underlying
  data.

The 'hit_array' element will contain the unique set of values for each
dimension: Cutoffs, Contrasts, Signal.

## Details

This function concatenates two or more objects 'sestats' as returned by
[`se_contrast_stats()`](https://jmw86069.github.io/jamses/reference/se_contrast_stats.md).
It assumes that contrast:signal are unique across all entries, however
the same contrast can be present if associated with different signal
values.

The cutoffs are not considered when determining which results can be
combined, in other words only the 'contrast' and 'signal' values are
used to define a "unique key". Only unique contrast:signal combinations
are valid as input.

In other words, input will not accept multiple sestats objects which
contain the same contrast:signal and two different cutoffs.

## See also

Other jamses utilities:
[`choose_annotation_colnames()`](https://jmw86069.github.io/jamses/reference/choose_annotation_colnames.md),
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
[`shrink_matrix()`](https://jmw86069.github.io/jamses/reference/shrink_matrix.md),
[`sort_samples()`](https://jmw86069.github.io/jamses/reference/sort_samples.md),
[`strsplitOrdered()`](https://jmw86069.github.io/jamses/reference/strsplitOrdered.md),
[`sub_split_vector()`](https://jmw86069.github.io/jamses/reference/sub_split_vector.md),
[`update_function_params()`](https://jmw86069.github.io/jamses/reference/update_function_params.md),
[`update_list_elements()`](https://jmw86069.github.io/jamses/reference/update_list_elements.md)
