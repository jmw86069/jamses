# Quick conversion of hit_array to hit_list

Quick conversion of hit_array to hit_list

## Usage

``` r
hit_array_to_list(
  hit_array,
  contrast_names = NULL,
  cutoff_names = NULL,
  assay_names = NULL,
  ...
)
```

## Arguments

- hit_array:

  `SEStats` or `list` with 'hit_array' element, or `array` with
  [`dimnames()`](https://rdrr.io/r/base/dimnames.html): 'Cutoffs',
  'Contrasts', 'Signals'.

- contrast_names:

  `character` vector of contrasts.

## Value

`list` named by `contrast_names`, that contains unique statistical hits
by combining entries across the `cutoff_names` and `assay_names` for
each contrast.

## Details

This function is mainly useful when there are multiple dimensions
unresolved in a hit_array, in which case this function will combine hits
across the different cutoffs and signals contained in the `hit_array` of
`sestats` output from [`se_contrast_stats()`](se_contrast_stats.md).

## See also

Other jamses stats: [`ebayes2dfs()`](ebayes2dfs.md),
[`format_hits()`](format_hits.md),
[`handle_na_values()`](handle_na_values.md),
[`process_sestats_to_hitim()`](process_sestats_to_hitim.md),
[`run_limma_replicate()`](run_limma_replicate.md),
[`save_sestats()`](save_sestats.md),
[`se_contrast_stats()`](se_contrast_stats.md),
[`sestats_to_df()`](sestats_to_df.md),
[`sestats_to_dfs()`](sestats_to_dfs.md), [`voom_jam()`](voom_jam.md)
