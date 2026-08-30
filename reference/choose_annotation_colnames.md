# Choose interesting annotation colnames from a data.frame

Choose interesting annotation colnames from a data.frame

## Usage

``` r
choose_annotation_colnames(
  df,
  min_reps = 2,
  min_values = 2,
  max_values = Inf,
  keep_numeric = FALSE,
  simplify = TRUE,
  max_colnames = 20,
  ...
)
```

## Arguments

- df:

  `data.frame` with annotations that could be interesting to display at
  the top or side of a heatmap.

- min_reps:

  `numeric` minimum number of replicates required for a column to be
  considered interesting. For example, `min_reps=3` would require any
  value in a column to be repeated at least `3` times for that column to
  be interesting. This filter is intended to remove columns whose values
  are all unique, such as row identifiers.

- min_values:

  `numeric` minimum number of unique values required for a column to be
  considered interesting.

- max_values:

  `numeric` maximum number of unique values required for a column to be
  considered interesting. Too many values and the interest is lost.
  Also, too many values, and the color key becomes unbearable with too
  many labels.

- keep_numeric:

  `logical` indicating whether to keep columns with `numeric` values.
  When `keep_numeric == TRUE` it will override the rules above.

- simplify:

  `logical` indicating whether to filter out columns whose data already
  matches another column with 1:1 cardinality. This step requires
  `platjam::cardinality()` until that function is moved into the `jamba`
  package.

- max_colnames:

  `numeric` maximum number of colnames to return. Note that columns are
  not sorted for priority, so they will be returned in the order they
  appear in `df` after applying the relevant criteria.

- ...:

  additional arguments are ignored.

## Value

`character` vector of colnames in `df` that meet the criteria. If no
colnames meet the criteria, this function returns `NULL`.

## See also

Other jamses utilities:
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
[`shrink_df()`](https://jmw86069.github.io/jamses/reference/shrink_df.md),
[`shrink_matrix()`](https://jmw86069.github.io/jamses/reference/shrink_matrix.md),
[`sort_samples()`](https://jmw86069.github.io/jamses/reference/sort_samples.md),
[`strsplitOrdered()`](https://jmw86069.github.io/jamses/reference/strsplitOrdered.md),
[`sub_split_vector()`](https://jmw86069.github.io/jamses/reference/sub_split_vector.md),
[`update_function_params()`](https://jmw86069.github.io/jamses/reference/update_function_params.md),
[`update_list_elements()`](https://jmw86069.github.io/jamses/reference/update_list_elements.md)

## Examples

``` r
df <- data.frame(
   threereps=paste0("threereps_", letters[c(1,1,1,3,5,7,7)]),
   time=paste0("time_", letters[c(1:7)]),
   tworeps=paste0("tworeps_", letters[c(12,12,14,14,15,15,16)]),
   num=sample(1:7),
   class=paste0("class_", LETTERS[c(1,1,1,3,5,7,7)]),
   blah=rep("blah", 7),
   maxvalues=c("one", "two", "three", "four", "five", "six", "six"))
df
#>     threereps   time   tworeps num   class blah maxvalues
#> 1 threereps_a time_a tworeps_l   5 class_A blah       one
#> 2 threereps_a time_b tworeps_l   4 class_A blah       two
#> 3 threereps_a time_c tworeps_n   7 class_A blah     three
#> 4 threereps_c time_d tworeps_n   2 class_C blah      four
#> 5 threereps_e time_e tworeps_o   3 class_E blah      five
#> 6 threereps_g time_f tworeps_o   1 class_G blah       six
#> 7 threereps_g time_g tworeps_p   6 class_G blah       six

choose_annotation_colnames(df)
#> Error in loadNamespace(x): there is no package called ‘platjam’
df[,choose_annotation_colnames(df)]
#> Error in loadNamespace(x): there is no package called ‘platjam’

choose_annotation_colnames(df, max_values=5)
#> Error in loadNamespace(x): there is no package called ‘platjam’
df[,choose_annotation_colnames(df, max_values=5)]
#> Error in loadNamespace(x): there is no package called ‘platjam’

choose_annotation_colnames(df, simplify=FALSE)
#>   threereps     tworeps       class   maxvalues 
#> "threereps"   "tworeps"     "class" "maxvalues" 
df[,choose_annotation_colnames(df, simplify=FALSE)]
#>     threereps   tworeps   class maxvalues
#> 1 threereps_a tworeps_l class_A       one
#> 2 threereps_a tworeps_l class_A       two
#> 3 threereps_a tworeps_n class_A     three
#> 4 threereps_c tworeps_n class_C      four
#> 5 threereps_e tworeps_o class_E      five
#> 6 threereps_g tworeps_o class_G       six
#> 7 threereps_g tworeps_p class_G       six

choose_annotation_colnames(df, min_reps=3)
#> Error in loadNamespace(x): there is no package called ‘platjam’

choose_annotation_colnames(df, min_reps=1)
#> Error in loadNamespace(x): there is no package called ‘platjam’

choose_annotation_colnames(df, keep_numeric=TRUE)
#> Error in loadNamespace(x): there is no package called ‘platjam’

choose_annotation_colnames(df, min_reps=1)
#> Error in loadNamespace(x): there is no package called ‘platjam’

choose_annotation_colnames(df, min_reps=1, keep_numeric=TRUE)
#> Error in loadNamespace(x): there is no package called ‘platjam’
```
