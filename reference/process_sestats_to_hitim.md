# Process sestats input into a hit incidence matrix

Process sestats input into a hit incidence matrix

## Usage

``` r
process_sestats_to_hitim(
  sestats,
  cutoff_names = NULL,
  contrast_names = NULL,
  assay_names = NULL,
  contrast_suffix = NULL,
  rename_contrasts = FALSE,
  rows = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- sestats:

  one of the following:

  - `list` object output from
    [`se_contrast_stats()`](https://jmw86069.github.io/jamses/reference/se_contrast_stats.md),
    containing `"hit_array"`

  - `SEStats` object output from
    [`se_contrast_stats()`](https://jmw86069.github.io/jamses/reference/se_contrast_stats.md),

  - `array` in format `"hit_array"` with dimnames
    `"Cutoffs","Contrasts","Signals"`.

  - `list` of `character` vectors representing `rownames(se)` for the
    parent
    [`heatmap_se()`](https://jmw86069.github.io/jamses/reference/heatmap_se.md)
    function.

  - `list` of `numeric` vectors named by `rownames(se)`.

- cutoff_names, contrast_names, assay_names:

  `character` or `numeric` passed to
  [`hit_array_to_list()`](https://jmw86069.github.io/jamses/reference/hit_array_to_list.md)
  when the input is `sestats` or `hit_array`.

- contrast_suffix:

  `character` optional suffix appended to the end of each contrast name.

- rename_contrasts:

  `logical` indicating whether to rename contrasts by calling
  [`contrast2comp()`](https://jmw86069.github.io/jamses/reference/contrast2comp.md)

- rows:

  `character` or `NULL` with optional fixed set of rownames expected in
  the output matrix. When `rows=NULL` all rows are returned using data
  from `sestats`. Otherwise, only rows defined by `rows` are returned.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are passed to
  [`contrast2comp()`](https://jmw86069.github.io/jamses/reference/contrast2comp.md)
  if relevant.

## See also

Other jamses stats:
[`ebayes2dfs()`](https://jmw86069.github.io/jamses/reference/ebayes2dfs.md),
[`format_hits()`](https://jmw86069.github.io/jamses/reference/format_hits.md),
[`handle_na_values()`](https://jmw86069.github.io/jamses/reference/handle_na_values.md),
[`hit_array_to_list()`](https://jmw86069.github.io/jamses/reference/hit_array_to_list.md),
[`run_limma_replicate()`](https://jmw86069.github.io/jamses/reference/run_limma_replicate.md),
[`save_sestats()`](https://jmw86069.github.io/jamses/reference/save_sestats.md),
[`se_contrast_stats()`](https://jmw86069.github.io/jamses/reference/se_contrast_stats.md),
[`sestats_to_df()`](https://jmw86069.github.io/jamses/reference/sestats_to_df.md),
[`sestats_to_dfs()`](https://jmw86069.github.io/jamses/reference/sestats_to_dfs.md),
[`voom_jam()`](https://jmw86069.github.io/jamses/reference/voom_jam.md)
