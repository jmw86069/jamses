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

## Details

This function is based directly upon
[`limma::voom()`](https://rdrr.io/pkg/limma/man/voom.html) with a small
adjustment to handle the presence of `NA` values, which otherwise causes
the [`stats::lowess()`](https://rdrr.io/r/stats/lowess.html) output to
be clearly incorrect. The correction removes `NA` values during this
step, producing a result as expected.

## See also

Other jamses stats: [`ebayes2dfs()`](ebayes2dfs.md),
[`format_hits()`](format_hits.md),
[`handle_na_values()`](handle_na_values.md),
[`hit_array_to_list()`](hit_array_to_list.md),
[`process_sestats_to_hitim()`](process_sestats_to_hitim.md),
[`run_limma_replicate()`](run_limma_replicate.md),
[`save_sestats()`](save_sestats.md),
[`se_contrast_stats()`](se_contrast_stats.md),
[`sestats_to_df()`](sestats_to_df.md),
[`sestats_to_dfs()`](sestats_to_dfs.md)
