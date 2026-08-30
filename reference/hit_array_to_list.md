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
`sestats` output from
[`se_contrast_stats()`](https://jmw86069.github.io/jamses/reference/se_contrast_stats.md).

## See also

Other jamses stats:
[`ebayes2dfs()`](https://jmw86069.github.io/jamses/reference/ebayes2dfs.md),
[`format_hits()`](https://jmw86069.github.io/jamses/reference/format_hits.md),
[`handle_na_values()`](https://jmw86069.github.io/jamses/reference/handle_na_values.md),
[`process_sestats_to_hitim()`](https://jmw86069.github.io/jamses/reference/process_sestats_to_hitim.md),
[`run_limma_replicate()`](https://jmw86069.github.io/jamses/reference/run_limma_replicate.md),
[`save_sestats()`](https://jmw86069.github.io/jamses/reference/save_sestats.md),
[`se_contrast_stats()`](https://jmw86069.github.io/jamses/reference/se_contrast_stats.md),
[`sestats_to_df()`](https://jmw86069.github.io/jamses/reference/sestats_to_df.md),
[`voom_jam()`](https://jmw86069.github.io/jamses/reference/voom_jam.md)
