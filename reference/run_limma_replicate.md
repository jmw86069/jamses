# Run limma contrasts with optional probe replicates

Run limma contrasts with optional probe replicates

## Usage

``` r
run_limma_replicate(
  imatrix,
  idesign,
  icontrasts,
  weights = NULL,
  robust = FALSE,
  adjust.method = "BH",
  confint = FALSE,
  trim_colnames = c("t", "B", "F", "sca.t"),
  adjp_cutoff = 0.05,
  p_cutoff = NULL,
  fold_cutoff = 1.5,
  int_adjp_cutoff = adjp_cutoff,
  int_p_cutoff = p_cutoff,
  int_fold_cutoff = fold_cutoff,
  mgm_cutoff = NULL,
  ave_cutoff = NULL,
  block = NULL,
  rowData_df = NULL,
  collapse_by_gene = FALSE,
  correlation = NULL,
  define_hits = TRUE,
  posthoc_test = c("none", "DEqMS"),
  posthoc_args = list(DEqMS = list(PSM_counts = NULL, fit.method = "loess")),
  seed = 123,
  verbose = FALSE,
  ...
)
```

## Arguments

- confint:

  `logical` passed to
  [`limma::topTable()`](https://rdrr.io/pkg/limma/man/toptable.html),
  which defines whether to return confidence intervals for each log2
  fold change.

- adjp_cutoff, p_cutoff, fold_cutoff, mgm_cutoff, ave_cutoff:

  `numeric` values representing the appropriate statistical threshold,
  or `NULL` when a threshold should not be applied.

- int_adjp_cutoff, int_p_cutoff, int_fold_cutoff:

  `numeric` thresholds to apply only to interaction contrasts.

- rowData_df:

  `data.frame` representing optional rowData annotation to be retained
  in the resulting stat `data.frame`. This argument is usually defined
  using `rowData_colnames` in
  [`se_contrast_stats()`](https://jmw86069.github.io/jamses/reference/se_contrast_stats.md),
  which uses corresponding columns from `rowData(se)`.

- collapse_by_gene:

  `logical` indicating whether to apply `collapse_stats_by_gene` which
  chooses one "best" exemplar per gene when there are multiple rows that
  represent the same gene.

- correlation:

  `numeric` or `NULL` passed to
  [`limma::lmFit()`](https://rdrr.io/pkg/limma/man/lmFit.html). Note
  that when `block` is defined (and non-empty), and when
  `correlation=NULL`, the correlation will be calculated by calling
  [`limma::duplicateCorrelation()`](https://rdrr.io/pkg/limma/man/dupcor.html).

- define_hits:

  `logical` indicating whether to define hits using the statistical
  thresholds.

- seed:

  `numeric` value used to define
  [`set.seed()`](https://rdrr.io/r/base/Random.html) for
  reproducibility. To avoid setting seed, use `seed=NULL`.

- verbose:

  `logical` indicating whether to print verbose output.

## Value

`list` with the following entries:

- "stats_df": `data.frame` with all individual `data.frame` per
  contrast, merged together.

- "stats_df": `list` of individual `data.frame` per contrast, each
  result is the output from
  [`ebayes2dfs()`](https://jmw86069.github.io/jamses/reference/ebayes2dfs.md).

- "rep_fits": `list` of various intermediate model fits, dependent upon
  whether limma, limma-voom, or limma-DEqMS were used.

## Details

This function is called by
[`se_contrast_stats()`](https://jmw86069.github.io/jamses/reference/se_contrast_stats.md)
to perform the comparisons defined as contrasts. The
[`se_contrast_stats()`](https://jmw86069.github.io/jamses/reference/se_contrast_stats.md)
function operates on a `SummarizedExperiment` object, this function
operates on the `numeric` `matrix` values directly.

This function also calls
[`ebayes2dfs()`](https://jmw86069.github.io/jamses/reference/ebayes2dfs.md)
which extracts each contrast result as a `data.frame`, whose column
names are modified to include the contrast names.

This function optionally (not yet ported from previous implementation)
detects replicate probes, and performs the internal correlation
calculations recommended by `limma user guide` for replicate probes. In
that case, it detects each level of probe replication so that each can
be properly calculated. For example, Agilent human 4x44 arrays often
contain a large number of probes with 8 replicates; a subset of probes
with 4 replicates; then the remaining probes (the majority overall) have
only one replicate each. In that case, this function splits data into
8-replicate, 4-replicate, and 1-replicate subsets, calculates
correlations on the 8-replicate and 4-replicate subsets separately, then
runs limma calculations on the three subsets independently, then merges
the results into one large table. The end result is that the final table
contains one row per unique probe after adjusting for probe replication
properly in each scenario. As the Agilent microarray layout is markedly
less widely used that in past, the priority to port this methodology is
quite low.

## See also

Other jamses stats:
[`ebayes2dfs()`](https://jmw86069.github.io/jamses/reference/ebayes2dfs.md),
[`format_hits()`](https://jmw86069.github.io/jamses/reference/format_hits.md),
[`handle_na_values()`](https://jmw86069.github.io/jamses/reference/handle_na_values.md),
[`hit_array_to_list()`](https://jmw86069.github.io/jamses/reference/hit_array_to_list.md),
[`process_sestats_to_hitim()`](https://jmw86069.github.io/jamses/reference/process_sestats_to_hitim.md),
[`save_sestats()`](https://jmw86069.github.io/jamses/reference/save_sestats.md),
[`se_contrast_stats()`](https://jmw86069.github.io/jamses/reference/se_contrast_stats.md),
[`sestats_to_df()`](https://jmw86069.github.io/jamses/reference/sestats_to_df.md),
[`sestats_to_dfs()`](https://jmw86069.github.io/jamses/reference/sestats_to_dfs.md),
[`voom_jam()`](https://jmw86069.github.io/jamses/reference/voom_jam.md)
