# Mark statistical hits by threshold cutoffs

Mark statistical hits by threshold cutoffs

## Usage

``` r
mark_stat_hits(
  x,
  adjp_cutoff = NULL,
  p_cutoff = NULL,
  fold_cutoff = NULL,
  mgm_cutoff = NULL,
  ave_cutoff = NULL,
  adjp_colname = "adj.P.Val",
  p_colname = "P.Value",
  logfc_colname = "logFC",
  mgm_colname = "mgm",
  ave_colname = "AveExpr",
  assign_value = c("sign", "fold", "logfc"),
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  `data.frame` containing one or more statistical columns

- adjp_cutoff, p_cutoff, fold_cutoff, mgm_cutoff, ave_cutoff:

  `numeric` value for each cutoff to be enforced, or `NULL` or `NA` to
  ignore each threshold. Each argument must have only 1 value assigned
  to be enforced.

- adjp_colname, p_colname, logfc_colname, mgm_colname, ave_colname:

  `character` string for each colname in `x` to be used for the
  appropriate statistical threshold.

- assign_value:

  `character` string indicating the value assigned to hits: `"sign"`
  uses the sign of the log2 fold change; `"fold"` uses the normal space
  fold change, by [`log2fold_to_fold()`](log2fold_to_fold.md); `"logfc"`
  uses the log2 fold change value. If there is no matching colname for
  `logfc_colname` then all hits are assigned `1`. In all cases, entries
  which are not hits are assigned `0`.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are ignored.

## Value

`numeric` vector with length `nrow(x)` and values defined by argument
`assign_value`. The order is identical to the order of rows in `x`
input. The output vector will be named by `rownames(x)` if rownames
exist.

## Details

This function is lightweight method of applying one or more statistical
thresholds to define "statistical hits". The thresholds are based upon
three questions:

- Is it "detected"? (above minimum signal)

- Is it "changing"? (above minimum fold change)

- Is it "significant"? (below defined P-value threshold)

The reasoning is roughly described:

- If the signal for a measurement is not above a noise threshold, or
  above a defined level of signature required for adequate confirmation
  experiments, the other statistical measurements are not relevant.

- If the change between two experimental groups is not sufficient for
  follow-up experiments, or is below a biologically meaningful level of
  change, the other statistical measures are not relevant.

- If the signal is detected, and the change is potentially sufficient
  for follow-up confirmation experiments, and/or to induce biologically
  meaningful effects, it must also be statistically robust as defined by
  the relevant adjusted P-value.

These thresholds are dependent upon the experiment itself, and each
threshold, if used, must be well-defined and defensible.

That is, in order to define a signal threshold, one should evaluate the
level of noise below which a measured value is no longer sufficient for
follow-up experiments, or no longer reliable based upon the technology
being used.

In order to impost a minimum fold change threshold, one should have some
clear indication of any limitations in follow-up assay techniques, and
some indication of the magnitude of change expected for a biologically
meaningful response. In some cases, a biologically meaningful change may
be defined in other experiments, ideally showing small changes not
associated with biologically meaningful effects and changes which are
associated with biollogically meaningful effects.

A useful technique to review statistical thresholds is a volcano plot,
which depicts the relationship of log fold change versus adjusted
P-value. The plot can indicate the range of fold changes for which the
statistical model found significance. Some technologies or protocols
naturally compress the effective fold change, yielding a very narrow
volcano plot, while others with high variability may result in a
relatively short and wide volcano plot. The range of fold changes with
no significant P-value may indicate a reasonable expectation for
inherent variability, thus a fold change threshold may be defined above
that observed for the majority of non-significant entries.

Entries which meet the statistical criteria are marked:

- `-1` for entries that meet all criteria, with negative fold change

- `0` for entries that do not meet all thresholds

- `1` for entries that meet all criteria, with positive fold change

"Detected" is defined either by the "max group mean", representing the
highest group mean signal intensity, or by "average signal",
representing the average group mean signal intensity across all groups.
These columns should already be present in the input data `x`.

"Changing" is defined by the log2 fold change, which must meet the
criteria defined by `fold_cutoff` which is in normal space. For example,
`fold_cutoff=1.5` represents 1.5-fold change, and would be applied
`abs(log2fc) >= log2(fold_cutoff)`. Most statistical results are
reported using log2 fold change, but scientists usually define fold
changes in normal space.

"Significant" is define using the P-value, and/or adjusted P-value, and
requires entries to be at or below the threshold.

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
