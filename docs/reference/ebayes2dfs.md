# Convert limma eBayes fit to data.frame with annotated hits

Convert limma eBayes fit to data.frame with annotated hits

## Usage

``` r
ebayes2dfs(
  lmFit3 = NULL,
  lmFit1 = NULL,
  lmFit4 = NULL,
  define_hits = TRUE,
  adjp_cutoff = 0.05,
  p_cutoff = NULL,
  fold_cutoff = 1.5,
  int_adjp_cutoff = adjp_cutoff,
  int_p_cutoff = p_cutoff,
  int_fold_cutoff = fold_cutoff,
  mgm_cutoff = NULL,
  ave_cutoff = NULL,
  confint = FALSE,
  use_cutoff_colnames = TRUE,
  rename_headers = TRUE,
  return_fold = TRUE,
  merge_df = FALSE,
  include_ave_expr = FALSE,
  include_group_means = TRUE,
  transform_means = c("none", "exp2signed", "10^"),
  rowData_df = NULL,
  collapse_by_gene = FALSE,
  rename_contrasts = FALSE,
  sep = " ",
  int_grep = "[(].+-.+-.+[)]|-.+-",
  trim_colnames = c("t", "B", "F", "sca.t"),
  posthoc_test = c("none", "DEqMS"),
  verbose = FALSE,
  ...
)
```

## Arguments

- lmFit3:

  object returned by
  [`limma::eBayes()`](https://rdrr.io/pkg/limma/man/ebayes.html).

- lmFit1:

  object returned by
  [`limma::lmFit()`](https://rdrr.io/pkg/limma/man/lmFit.html),
  optional.

- lmFit4:

  object returned by `posthoc_test="DEqMS"` in
  [`run_limma_replicate()`](run_limma_replicate.md).

- define_hits:

  `logical` indicating whether to define hits using the statistical
  thresholds.

- adjp_cutoff, p_cutoff, fold_cutoff, mgm_cutoff, ave_cutoff:

  `numeric` values representing the appropriate statistical threshold,
  or `NULL` when a threshold should not be applied.

- int_adjp_cutoff, int_p_cutoff, int_fold_cutoff:

  `numeric` thresholds to apply only to interaction contrasts.

- confint:

  `logical` passed to
  [`limma::topTable()`](https://rdrr.io/pkg/limma/man/toptable.html),
  which defines whether to return confidence intervals for each log2
  fold change.

- use_cutoff_colnames:

  `logical` whether to include the statistical thresholds abbreviated in
  the `"hit"` colname, when `define_hits=TRUE`.

- rename_headers:

  `logical` indicating whether to rename statistical colnames returned
  by [`limma::topTable()`](https://rdrr.io/pkg/limma/man/toptable.html)
  to the colnames include the contrast name.

- return_fold:

  `logical` whether to return an additional column with the signed fold
  change, see [`log2fold_to_fold()`](log2fold_to_fold.md).

- merge_df:

  `logical` indicating whether to merge the final `data.frame` list into
  one `data.frame`.

- include_ave_expr:

  `logical` indicating whether to retain the column `"AveExpr"`. This
  column can be misleading, especially if the `mgm` (max group mean)
  threshold is used when determining statistical hits. This column is
  mainly useful in reviewing limma output, since it uses the `"AveExpr"`
  values to apply its moderated variance statistic.

- include_group_means:

  `logical` indicating whether to include each group mean along with the
  relevant contrast. These values are helpful, in that they should
  exactly represent the reported `logFC` value. Sometimes it is helpful
  and comforting to see the exact values used in that calculation.

- rowData_df:

  `data.frame` representing optional rowData annotation to be retained
  in the resulting stat `data.frame`. This argument is usually defined
  using `rowData_colnames` in
  [`se_contrast_stats()`](se_contrast_stats.md), which uses
  corresponding columns from `rowData(se)`.

- collapse_by_gene:

  `logical` indicating whether to apply `collapse_stats_by_gene` which
  chooses one "best" exemplar per gene when there are multiple rows that
  represent the same gene.

- rename_contrasts:

  `logical` (inactive) which will in future allow for automated renaming
  of contrasts.

- sep:

  `character` string used as a delimiter in certain output colnames.

- int_grep:

  `character` string used to recognize contrasts which are considered
  "interaction contrasts". The default pattern recognizes any contrasts
  that contain multiple fold changes, recognized by the presence of more
  than one hypen `"-"` in the contrast name.

- verbose:

  `logical` indicating whether to print verbose output.

## Value

`list` with one `data.frame` per contrast defined in the input `lmFit3`
object. When `define_hits=TRUE` there will be one column per statistical
threshold, named `"hit"` followed by an abbreviation of the statistical
thresholds which were applied. When `merge_df=TRUE` the returned data
will be one `data.frame` object.

## Details

This function is called by
[`run_limma_replicate()`](run_limma_replicate.md) as an extension to
[`limma::topTable()`](https://rdrr.io/pkg/limma/man/toptable.html), that
differs in that it is performed for each contrast in the input `lmFit3`
object.

By default the columns include the contrast, so that each `data.frame`
is self-described.

When `define_hits=TRUE`, then statistical thresholds are applied to
define a set of statistical hits. The thresholds available include:

1.  `adjp_cutoff` - applied to `"adj.P.Val"` for adjusted P-value.

2.  `p_cutoff` - applied to `"P.Value"` for raw, unadjusted P-value.

3.  `fold_cutoff` - normal space fold change, applied to `"logFC"` by
    using `log2(fold_cutoff)`.

4.  `mgm_cutoff` - max group mean, applied to the highest group mean
    value involved in each specific contrast.

5.  `ave_cutoff` - applied to `"AveExpr"` which represents the mean
    value across all sample groups.

Note that `mgm_cutoff` requires input `lmFit1` which stores the group
mean values used in the limma workflow.

Note also there are optional arguments specific to interaction
contrasts, which in this context is assumed to be a "fold change of fold
changes" style of contrast, for example:
`(groupA-groupB)-(groupC-groupD)`. The purpose is distinct interaction
thresholds is to enable reasonable data mining, sometimes with somewhat
more lenient thresholds for interaction contrasts. For example, one may
use `adjp_cutoff=0.01` and `int_adjp_cutoff=0.05`, or `fold_cutoff=2`
and `int_fold_cutoff=1.5`.

By default, `rename_headers=TRUE` causes colnames to include the
contrast, for example renaming colname `"logFC"` to `"logFC contrastA"`.
This change helps reinforce the source of the statistical results, and
allows the `data.frame` results to be merged together using
[`base::merge()`](https://rdrr.io/r/base/merge.html).

Indeed, `merge_df=TRUE` will cause all `data.frame` results to be merged
into one large `data.frame`, using
[`jamba::mergeAllXY()`](https://jmw86069.github.io/jamba/reference/mergeAllXY.html).

## See also

Other jamses stats: [`format_hits()`](format_hits.md),
[`handle_na_values()`](handle_na_values.md),
[`hit_array_to_list()`](hit_array_to_list.md),
[`process_sestats_to_hitim()`](process_sestats_to_hitim.md),
[`run_limma_replicate()`](run_limma_replicate.md),
[`save_sestats()`](save_sestats.md),
[`se_contrast_stats()`](se_contrast_stats.md),
[`sestats_to_df()`](sestats_to_df.md),
[`sestats_to_dfs()`](sestats_to_dfs.md), [`voom_jam()`](voom_jam.md)
