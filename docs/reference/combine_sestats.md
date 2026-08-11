# Combine SEStats objects by contrast,signal

Combine SEStats objects by contrast,signal

## Usage

``` r
combine_sestats(sestats_list, create_stats_df = TRUE, ...)
```

## Arguments

- sestats_list:

  `list` with two or more `list` objects in sestats format, as returned
  by [`se_contrast_stats()`](se_contrast_stats.md).

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
[`se_contrast_stats()`](se_contrast_stats.md) with minor exceptions:

- List elements 'idesign','icontrasts','normgroup' are each returned as
  a `list` with corresponding data from each sestats object. The reason
  is that these values cannot logically be combined into one new object,
  and therefore the next level of support is to retain the underlying
  data.

The 'hit_array' element will contain the unique set of values for each
dimension: Cutoffs, Contrasts, Signal.

## Details

This function concatenates two or more objects 'sestats' as returned by
[`se_contrast_stats()`](se_contrast_stats.md). It assumes that
contrast:signal are unique across all entries, however the same contrast
can be present if associated with different signal values.

The cutoffs are not considered when determining which results can be
combined, in other words only the 'contrast' and 'signal' values are
used to define a "unique key". Only unique contrast:signal combinations
are valid as input.

In other words, input will not accept multiple sestats objects which
contain the same contrast:signal and two different cutoffs.

## See also

Other jamses utilities:
[`choose_annotation_colnames()`](choose_annotation_colnames.md),
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
[`shrink_df()`](shrink_df.md), [`shrink_matrix()`](shrink_matrix.md),
[`sort_samples()`](sort_samples.md),
[`strsplitOrdered()`](strsplitOrdered.md),
[`sub_split_vector()`](sub_split_vector.md),
[`update_function_params()`](update_function_params.md),
[`update_list_elements()`](update_list_elements.md)
