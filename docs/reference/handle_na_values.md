# Handle NA values in a numeric matrix

Handle NA values in a numeric matrix

## Usage

``` r
handle_na_values(
  x,
  idesign,
  handle_na = c("full1", "full", "partial", "none", "all"),
  na_value = 0,
  na_weight = 0,
  return_weights = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  `numeric` matrix

- idesign:

  `numeric` matrix with `rownames(idesign)` equal to `colnames(x)`,
  containing `0` or `1` to fill the design matrix.

- handle_na:

  `character` string to determine the method used to handle NA values in
  `x`.

- na_value:

  `numeric` or `NA` used to handle NA values.

- na_weight:

  `numeric` weight between `0` and `1` used for `NA` values in the
  weight matrix, used when `return_weights=TRUE`.

- return_weights:

  `logical` indicating whether to include a weight matrix as an
  attribute with name `"weights"`.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are ignored.

## Value

`numeric` matrix with equal dimensions as input `x`, where `NA` values
have been handled as defined by `handle_na`.

## Details

This function provides reasonable alternatives intended to manage the
presence of missing data encoded as `NA` values in a numeric matrix.

The alternatives are defined by argument `handle_na`:

- `"full"`: Retain `NA` values, except when an entire group is `NA` it
  is replaced with `na_value`. This option is intended for sparse data
  where non-NA values are accepted as real measurements for each group,
  and where a group with all NA values should be retained for
  statistical contrasts by assigning `na_value`. This method essentially
  keeps the data as-is, except when groups are otherwise entirely `NA`
  the value `na_value` is used in order to retain any relevant contrasts
  that involve this group.

- `"full1"`: Similar to `"full"`, retain `NA` values, except when an
  entire group is `NA`, then replace only one entry with `na_value`.
  This option is intended to keep data as-is, except to retain groups
  that are otherwise entirely `NA`. In these groups, only one `na_value`
  is used in order to prevent the group from contributing toward group
  variability or dispersion calculations.

- `"partial"`: Replace `NA` values with `na_value`, except when an
  entire group is `NA` the entire group is kept at `NA`.

- `"all"`: Replace all `NA` values with `na_value`.

- `"none"`: Perform no replacement of `NA` values.

## See also

Other jamses stats: [`ebayes2dfs()`](ebayes2dfs.md),
[`format_hits()`](format_hits.md),
[`hit_array_to_list()`](hit_array_to_list.md),
[`process_sestats_to_hitim()`](process_sestats_to_hitim.md),
[`run_limma_replicate()`](run_limma_replicate.md),
[`save_sestats()`](save_sestats.md),
[`se_contrast_stats()`](se_contrast_stats.md),
[`sestats_to_df()`](sestats_to_df.md),
[`sestats_to_dfs()`](sestats_to_dfs.md), [`voom_jam()`](voom_jam.md)
