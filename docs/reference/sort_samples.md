# Sort biological sample labels for experimental design

Sort biological sample labels for experimental design

## Usage

``` r
sort_samples(
  x,
  control_terms = c("WT|wildtype", "normal|healthy|healthycontrol|^hc$",
    "control|ctrl|ctl", "(^|[-_ ])(NT|NTC)($|[-_ ]|[0-9])", "none|empty|blank",
    "untreated|untrt|untreat", "Vehicle|veh", "ETOH|ethanol", "scramble|mock|sham",
    "ttx|PBS", "knockout", "mutant"),
  sortFunc = jamba::mixedSort,
  pre_control_terms = NULL,
  post_control_terms = NULL,
  ignore.case = TRUE,
  boundary = TRUE,
  perl = boundary,
  keep_factor_order = TRUE,
  ...
)
```

## Arguments

- x:

  character vector or factor

- control_terms:

  vector of regular expression patterns used to determine control terms,
  where the patterns are matched and returned in order.

- pre_control_terms:

  vector or NULL, optional control terms or regular expressions to use
  before the `control_terms` above. This argument is used as a
  convenient prefix to the default terms.

- post_control_terms:

  vector or NULL, optional control terms or regular expressions to use
  after the `control_terms` above. This argument is used as a convenient
  suffix to the default terms.

- ignore.case:

  logical passed to
  [`jamba::provigrep()`](https://jmw86069.github.io/jamba/reference/provigrep.html)
  indicating whether to ignore case-sensitive matching.

- boundary:

  logical indicating whether to require a word boundary at either the
  start or end of the control terms. When TRUE, it uses `perl=TRUE` by
  default, and allows either perl boundary or an underscore `"_"`.

- perl:

  logical indicating whether to use Perl regular expression pattern
  matching.

- keep_factor_order:

  logical indicating whether to maintain factor level order, if `x` is
  supplied as a factor. If `keep_factor_order==TRUE` then only `sort(x)`
  is returned.

- ...:

  additional arguments are ignored.

## Value

character vector ordered such that control terms are preferentially
first before non-control terms.

## Details

This function sorts a vector of sample labels using typical heuristics
that order typical control groups terms before test groups. For example,
`"Vehicle"` would be returned before `"Treatment"` since `"Vehicle"` is
a recognized control term.

It also employs
[`jamba::mixedSort()`](https://jmw86069.github.io/jamba/reference/mixedSort.html)
for proper alphanumeric sorting, for example so `"Time_5hr"` would be
sorted before `"Time_12hr"`.

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
[`shrink_df()`](shrink_df.md), [`shrink_matrix()`](shrink_matrix.md),
[`strsplitOrdered()`](strsplitOrdered.md),
[`sub_split_vector()`](sub_split_vector.md),
[`update_function_params()`](update_function_params.md),
[`update_list_elements()`](update_list_elements.md)

## Examples

``` r
# the defaults perform well for clear descriptors
sort_samples(c("Trt_12h", "Trt_9h", "Trt_1h", "Trt_9h", "Vehicle"));
#> [1] "Vehicle" "Trt_1h"  "Trt_9h"  "Trt_9h"  "Trt_12h"

# custom terms can be added before the usual control terms
sort_samples(c("Trt_12h", "Trt_9h", "Trt_1h", "Trt_9h", "Fixated", "Vehicle"),
   pre_control_terms="fixate");
#> [1] "Fixated" "Vehicle" "Trt_1h"  "Trt_9h"  "Trt_9h"  "Trt_12h"

# custom terms can be added after the usual control terms
sort_samples(c("Trt_12h", "Trt_9h", "Trt_1h", "Trt_9h", "Fixated", "Vehicle"),
   post_control_terms="fixate");
#> [1] "Vehicle" "Fixated" "Trt_1h"  "Trt_9h"  "Trt_9h"  "Trt_12h"
```
