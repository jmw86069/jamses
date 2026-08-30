# Limma-voom customized for Jam

Limma-voom customized for Jam

## Usage

``` r
voom_jam(
  counts,
  design = NULL,
  lib.size = NULL,
  normalize.method = "none",
  block = NULL,
  correlation = NULL,
  weights = NULL,
  span = 0.5,
  plot = FALSE,
  save.plot = TRUE,
  verbose = FALSE,
  ...
)
```

## Arguments

- counts:

  a numeric `matrix` containing raw counts, or an `ExpressionSet`
  containing raw counts, or a `DGEList` object. Counts must be
  non-negative and NAs are not permitted.

- design:

  design matrix with rows corresponding to samples and columns to
  coefficients to be estimated. Defaults to
  `model.matrix(~group, data=counts$samples)` if `counts` is a DGEList,
  otherwise defaults to the unit vector meaning that all samples are
  treated as replicates.

- lib.size:

  numeric vector containing the library sizes for each sample. Defaults
  to the columnwise count totals if `counts` is a matrix or to
  `normLibSizes(counts)` if `counts` is a `DGEList`.

- normalize.method:

  the microarray-style normalization method to be applied to the logCPM
  values. Choices are as for the `method` argument of
  `normalizeBetweenArrays` when the data is single-channel.

- block:

  vector or factor specifying a blocking variable on the samples. Has
  length equal to the number of samples (`ncol(counts)`).

- correlation:

  the intrablock correlation. Normally a single numeric value between -1
  and 1, but a vector of genewise correlations is also allowed.

- weights:

  prior weights. Can be a numeric matrix of individual weights of same
  dimensions as the `counts`, or a numeric vector of sample weights with
  length equal to `ncol(counts)`, or a numeric vector of gene weights
  with length equal to `nrow(counts)`.

- span:

  width of the smoothing window used for the lowess mean-variance trend.
  Expressed as a proportion between 0 and 1.

- plot:

  logical, should a plot of the mean-variance trend be displayed?

- save.plot:

  logical, should the coordinates and line of the plot be saved in the
  output?

## Details

This function is based directly upon
[`limma::voom()`](https://rdrr.io/pkg/limma/man/voom.html) with a small
adjustment to handle the presence of `NA` values, which otherwise causes
the [`stats::lowess()`](https://rdrr.io/r/stats/lowess.html) output to
be clearly incorrect. The correction removes `NA` values during this
step, producing a result as expected.

## See also

Other jamses stats:
[`ebayes2dfs()`](https://jmw86069.github.io/jamses/reference/ebayes2dfs.md),
[`format_hits()`](https://jmw86069.github.io/jamses/reference/format_hits.md),
[`handle_na_values()`](https://jmw86069.github.io/jamses/reference/handle_na_values.md),
[`hit_array_to_list()`](https://jmw86069.github.io/jamses/reference/hit_array_to_list.md),
[`process_sestats_to_hitim()`](https://jmw86069.github.io/jamses/reference/process_sestats_to_hitim.md),
[`run_limma_replicate()`](https://jmw86069.github.io/jamses/reference/run_limma_replicate.md),
[`save_sestats()`](https://jmw86069.github.io/jamses/reference/save_sestats.md),
[`se_contrast_stats()`](https://jmw86069.github.io/jamses/reference/se_contrast_stats.md),
[`sestats_to_df()`](https://jmw86069.github.io/jamses/reference/sestats_to_df.md)
