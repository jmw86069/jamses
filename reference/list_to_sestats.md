# Prepare SEStats from a list of stat data.frame

Prepare SEStats from a list of stat data.frame

## Usage

``` r
list_to_sestats(
  stats_dfs,
  use_assay_name = "norm",
  p_cutoff = 1,
  adjp_cutoff = 0.05,
  fold_cutoff = 1.5,
  mgm_cutoff = 0,
  hit_pattern = "hit",
  fold_pattern = "fold",
  lfc_pattern = "^logFC|^log2F|^logfold|^log.fold|^lfc",
  p_pattern = "p.{0,1}val(ue|)",
  adjp_pattern = "adj.{0,1}p.{0,1}val(ue|)|p.{0,1}adj|fdr|q.{0,1}val(ue|)",
  mgm_pattern = "mgm|max.{0,1}group.{0,1}mean",
  contrast_pattern = "^.*[ ]([^ ]+-[^ ]+)$",
  verbose = FALSE,
  ...
)
```

## Arguments

- stats_dfs:

  `list` with one of two formats:

  1.  Containing a `list` named by `assay_names` (signal), which then
      contains a `list` of `data.frame` objects.

  2.  Containing `data.frame` objects.

  Note: The `list` of `data.frame` objects can be named by contrast, and
  these contrast names will be used if the contrast is not already
  encoded in the colnames the corresponding `data.frame`.

- use_assay_name:

  `character` string, default "norm", used when argument 'stats_dfs' is
  supplied as a `list` of `data.frames` (option 2 described for
  'stats_dfs'). In this case, the assay name will be 'use_assay_name'.

- p_cutoff, adjp_cutoff, fold_cutoff, mgm_cutoff:

  `numeric` values used when there is no matching hit column in each
  `data.frame`, matching the `hit_pattern`.

- hit_pattern, fold_pattern, lfc_pattern, p_pattern, adjp_pattern,
  mgm_pattern:

  `character` string with regular expression pattern used to match each
  corresponding column type: Hit column, fold change, log fold change,
  P-value (unadjusted), adjusted P-value (or FDR or Q-value), max group
  mean (mgm).

  - The pattern is assumed to be at the start of each column name
    string.

- contrast_pattern:

  `character` string used to match the contrast encoded in any of the
  cutoff columns for the '\_pattern' arguments.

  - The pattern is essentially any non-whitespace string that contains a
    hyphen '-' inside it. If it also contains a colon ':' the string is
    assumed to be a comp (see
    [`comp2contrast()`](https://jmw86069.github.io/jamses/reference/contrast2comp.md))
    otherwise it is considered a contrast name.

  - When present, the contrast is extracted from the column name.

  - When not present, the name from the `list` of `data.frame` objects
    is assumed to be the contrast name.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are ignored.

## Value

`list` consistent with
[`se_contrast_stats()`](https://jmw86069.github.io/jamses/reference/se_contrast_stats.md)
output, except limited to elements: 'stats_df', 'hit_array', 'hit_list'.
The content should be sufficient for use in
[`heatmap_se()`](https://jmw86069.github.io/jamses/reference/heatmap_se.md)
and other related functions.

## Details

This function converts a list of stat `data.frame` objects to a `list`
object as returned by
[`se_contrast_stats()`](https://jmw86069.github.io/jamses/reference/se_contrast_stats.md).
In future this function will return a S4 `SEStats` object.

### Limitations

This function only populates certain components:

- `stats_dfs` as a `list` by assay_name (signal), with `list` by
  contrast name with `data.frame` objects.

- `hit_array` as a 3-dimension array with 'Cutoffs', 'Contrasts', and
  'Signal' dimensions.

- `hit_list` as a `list` by assay name (signal), with `list` by contrast
  names, containing `numeric` values named by measurement.

This function does not return certain data which is not described in a
`list` of `data.frame` objects:

- `idesign`, `icontrasts`, `normgroup` cannot be defined or inferred.

Other limitations:

- `stats_df` with the combined stats across contrasts as one wide
  `data.frame`, is not created, partly because it does not assume the
  contrast name is included in each corresponding stat colname in each
  `data.frame`. Potential future enhancement.

### Rules Used

In each `data.frame` certain colnames are recognized:

- The hit columns is recognized by one or more colnames that begin with
  `"hit "`. The hit column should have values of `-1`, `0`, or `1` to
  indicate whether the row met the statistical criteria. Typically the
  statistical criteria are also included in the column header, followed
  by the contrast name.

  - If there is no hit column, the arguments '\_cutoff' are applied to
    each row, and a new hit column is created accordingly. If there are
    multiple cutoff values, multiple hit columns are created.

- Fold change column is recognized by a column beginning `"fold "` (with
  space), typically followed by the contrast name. For example
  `"fold Dex-Ctrl"`.

  - If no fold change column exists, it searches for log fold column,
    specifically assumed to be log2 fold change even when `"logFC"`. By
    convention, edgeR and limma use `"logFC"` for log2 fold change.

  - If neither fold change, nor log fold change, columns are found, this
    function stops.

  - The log2 fold change is converted to normal fold change for
    filtering.

- The raw P-value column should begin `"P-Value "` or any variation of
  "pval", "p.val", "pvalue", "p.value", case-insensitive.

- The adjusted P-value column should begin `"adj-P-Val "` or any
  variation of "adjp", "padj", "pvalue", "p.value", case-insensitive.

## See also

Other jamses utilities:
[`choose_annotation_colnames()`](https://jmw86069.github.io/jamses/reference/choose_annotation_colnames.md),
[`combine_sestats()`](https://jmw86069.github.io/jamses/reference/combine_sestats.md),
[`contrast2comp_dev()`](https://jmw86069.github.io/jamses/reference/contrast2comp_dev.md),
[`fold_to_log2fold()`](https://jmw86069.github.io/jamses/reference/fold_to_log2fold.md),
[`intercalate()`](https://jmw86069.github.io/jamses/reference/intercalate.md),
[`list2im_opt()`](https://jmw86069.github.io/jamses/reference/list2im_opt.md),
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
