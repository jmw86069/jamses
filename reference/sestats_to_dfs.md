# Extract stats as data.frame from SEStats results

Extract stat as data.frame from SEStats results

## Usage

``` r
sestats_to_dfs(
  sestats,
  assay_names = NULL,
  contrast_names = NULL,
  data_content = c("contrasts", "hits"),
  hits_use_lfc = FALSE,
  rename_contrasts = TRUE,
  se = NULL,
  rowData_colnames = NULL,
  row_type = "gene_name",
  verbose = FALSE,
  ...
)
```

## Arguments

- sestats:

  `SEStats` or `list` output from
  [`se_contrast_stats()`](https://jmw86069.github.io/jamses/reference/se_contrast_stats.md)

- assay_names:

  `character` string indicating which assay names to save, stored in
  `dimnames(sestats@hit_array)$Signal`. When `NULL` then all assay names
  are saved.

- contrast_names:

  `character` string indicating which contrasts to save, stored in
  `dimnames(sestats@hit_array)$Contrasts`. The default `NULL` will save
  all contrasts.

- data_content:

  `character` string describing the data content to include:

  - `"contrasts","hits"` - include worksheets per `contrast_names`, then
    assemble one `"hit sheet"` across all contrasts. One hit sheet is
    created for each value in `assay_names`.

  - `"contrasts"` - (default) include worksheets per `contrast_names`

  - `"hits"` - include only one `"hit sheet"` per value in
    `assay_names`.

- hits_use_lfc:

  `logical` default FALSE, indicating whether values in `"hits"` columns
  should use the log2 fold change.

  - `FALSE` (default) assigns `c(-1, 0, 1)` to indicate directionality
    after applying stat thresholds.

  - `TRUE` assigns the actual log2 fold change *only for hits* as
    defined by the stat thresholds.

- rename_contrasts:

  `logical` indicating whetheer to apply `contrasts2comp()` to shorten
  long contrast names. Currently this option only applies to the list
  element names, not the column headers.

- se:

  `SummarizedExperiment`, default NULL, used when `rowData_colnames` is
  defined.

- rowData_colnames:

  `character`, default NULL, with optional colnames used only when `se`
  is also provided. When defined, it provides additional annotations for
  each row as defined by `rowData(se)`.

- row_type:

  `character` with custom column name to use for the primary row
  identifier. The default `"probes"` is often not accurate, though this
  may not be problematic in practice. When defined, the first column is
  renamed to `row_type`.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are passed to
  [`save_sestats()`](https://jmw86069.github.io/jamses/reference/save_sestats.md).

## Value

`list` of `data.frame` based upon `data_content`:

- `data_content="contrasts"` will return `list` with one `data.frame`
  for each contrast, and each assay_name. Each list element is named by
  the contrast name.

- `data_content="hits"` will return `list` with one `data.frame` for
  each `assay_name`, including one column for each contrast. Each list
  element is named by "hits" plus the assay_name, for example
  `"hits quantile_counts"`.

- `data_content=c("contrasts", "hits")` will return `list` which
  includes both the options above.

## Details

The purpose is to output `data.frame` or `list` of `data.frame` from
SEStats results. This function is essentially a wrapper to
[`save_sestats()`](https://jmw86069.github.io/jamses/reference/save_sestats.md)
with `type="list"`.

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
[`sestats_to_df()`](https://jmw86069.github.io/jamses/reference/sestats_to_df.md),
[`voom_jam()`](https://jmw86069.github.io/jamses/reference/voom_jam.md)
