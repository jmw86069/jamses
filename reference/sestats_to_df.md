# Convert sestats to table summary

Convert sestats to table summary

## Usage

``` r
sestats_to_df(
  sestats,
  style = c("text", "integer"),
  dimname_order = c(3, 2, 1),
  rename_contrasts = FALSE,
  ...
)
```

## Arguments

- sestats:

  `SEStats` or `list` output from
  [`se_contrast_stats()`](https://jmw86069.github.io/jamses/reference/se_contrast_stats.md)

- style:

  `character` string indicating what values to use:

  - `"text"`: number of hits (number up, number down)

  - `"integer"`: only the integer number of hits

- rename_contrasts:

  `logical` indicating whether to rename contrasts using
  [`contrast2comp()`](https://jmw86069.github.io/jamses/reference/contrast2comp.md).
  The main benefit is reduction in length of the string that describes
  each contrast.

- ...:

  additional arguments are passed to
  [`contrast2comp()`](https://jmw86069.github.io/jamses/reference/contrast2comp.md),
  which is relevant when `rename_contrasts=TRUE`, relevant arguments
  include:

  - `contrast_delim="-"` to

  - `contrast_factor_delim="_"` to customize the delimiter between
    factors in the input group names, for example `"Treatment1_Time1"`
    uses "\_".

  - `comp_factor_delim=":"` to customize the delimiter between factors

  - `factor_order=NULL` to customize the order of factor comparisons

## Value

`data.frame` with a summary of statistical hits per contrast,
assay_name, and threshold.

## Details

Note: This function is intended to provide a simple data.frame summary
with the number of hits for each contrast, signal, cutoff. It is still
being tested, and updated for usability.

TODO: The order of `dimnames(hit_array)` should be user-customizable.
The series of dimnames in each lapply should use this order.

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
[`voom_jam()`](https://jmw86069.github.io/jamses/reference/voom_jam.md)

## Examples

``` r
if (FALSE) {
hitdf <- sestats_to_df(list(hit_array=hit_array));
hitdf_rowindex <- table(hitdf[[1]])[unique(hitdf[[1]])]
jamba::kable_coloring(
   hitdf[, -1, drop=FALSE],
   row.names=FALSE) %>%
   kableExtra::pack_rows(index=hitdf_rowindex)
}
```
