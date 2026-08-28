# Collapse SummarizedExperiment data by row

Collapse SummarizedExperiment data by row

## Usage

``` r
se_collapse_by_row(
  se,
  rows = rownames(se),
  row_groups,
  assay_names = NULL,
  group_func_name = c("sum", "mean", "weighted.mean", "geomean", "none"),
  rowStatsFunc = NULL,
  rowDataColnames = NULL,
  keepNULLlevels = FALSE,
  delim = "[ ]*[;,]+[ ]*",
  data_transform = c("none", "log2p+sqrt", "log2+sqrt", "log2p", "log2"),
  verbose = FALSE,
  ...
)
```

## Arguments

- se:

  `SummarizedExperiment`

- rows:

  `character` vector of `rows(se)` to use for analysis. When `rows=NULL`
  the default is to use all `rows(se)`.

- row_groups:

  `character` vector representing groups of rows to be combined.

- assay_names:

  `character` vector of `names(assays(se))` to use for the collapse
  operation. When `assay_names=NULL` the default is to use all
  `assays(se)`.

- group_func_name:

  `character` name of function used to aggregate measurement data within
  `row_groups`.

  - `sum` - takes the [`sum()`](https://rdrr.io/r/base/sum.html) of each
    value in the group. This option should be used together with
    `data_transform` when there has been any data transformation, so
    that the data is inverse-transformed prior to calculating the
    [`sum()`](https://rdrr.io/r/base/sum.html), after which data is
    re-transformed to its original state. This method is appropriate for
    log2p `log2(1 + x)` transformed abundance measurements for example.

  - `mean` - calculates the mean value per group. Note that in this case
    is it usually recommended not to define `data_transform` so that
    values are averaged in the appropriately transformed numeric space.

  - `weighted.mean` - calculates
    [`weighted.mean()`](https://rdrr.io/r/stats/weighted.mean.html)
    where weights `w` are defined by the values used. This method may be
    appropriate and effective with normal space abundance values derived
    from proteomics mass spec quantitation.

  - `geomean` - calculates geometric mean of values in each group.

  - `none` -

- rowStatsFunc:

  `function` optional function used instead of `group_func_name`.

- rowDataColnames:

  `character` subset of colnames in `rowData(se)` to be retained in the
  output data. Multiple values are combined usually by comma-delimited
  concatenation within `row_groups`, therefore it may be beneficial to
  include only relevant columns in that output.

- keepNULLlevels:

  `logical` indicating whether to drop unused factor levels in
  `row_groups`, this argument is passed to
  [`jamba::rowGroupMeans()`](https://jmw86069.github.io/jamba/reference/rowGroupMeans.html).

- delim:

  `character` string indicating a delimiter.

- data_transform:

  `character` string indicating which transformation was used when
  preparing the assay data. The assumption is that all assays were
  transformed by this method. During processing, data is
  inverse-transformed prior to applying the `group_func_name` or
  `rowStatsFunc` if supplied. After that function is applied, data is
  transformed using this function. The purpose is to enable taking the
  [`sum()`](https://rdrr.io/r/base/sum.html) in proper measured absolute
  units (in normal space for example) where relevant, after which is
  original numeric transformation is re-applied.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are passed to
  [`jamba::rowGroupMeans()`](https://jmw86069.github.io/jamba/reference/rowGroupMeans.html).

## Value

`SummarizedExperiment` object with these changes:

- rows will be collapsed by `row_groups`, for each `assays(se)`
  `numeric` matrix defined by `assay_names`. The collapse may optionally
  apply a data transformation defined in `data_transform` in order to
  apply an appropriate `numeric` summary calculation.

- `rowData(se)` will also be collapsed by
  [`shrinkDataFrame()`](https://jmw86069.github.io/jamses/reference/shrinkDataFrame.md)
  to combine unique values from each row annotation.

## Details

Purpose is to collapse rows of a `SummarizedExperiment` object, where
measurements for a given entity, usually a gene, are split across
multiple rows in the source data. The output of this function should be
measurements appropriately summarized to the gene level.

The key arguments are `group_func_name`, and `data_transform`. Note that
data is inverse-transformed based upon `data_transform`, prior to
calculating group summary values defined by `group_func_name`. The
reason is to enable using `group_func_name="sum"` on normal space
abundance values, when input data has already been transformed with
`log2(1 + x)` for example. In this case it is most appropriate to take
the `sum` of normal space abundance values, then to re-apply the
transformation afterwards.

However, when using `group_func_name="mean"` it is usually recommended
to use `data_transform="none"` so that data is maintained in
appropriately transformed state.

The driving use case is proteomics mass spectrometry data, where
measurements are described in terms of peptide sequences, with or
without optional post-translational modification (PTM), and the peptide
sequences are annotated to a source protein or gene. This function can
be used to:

- collapse peptide-PTM data to the peptide level

- collapse peptide data to the protein level

In future it may be used to collapse multiple microarray probe
measurements to the gene level, although that process is more likely to
be useful and recommended after performing probe-level statistical
analysis.

### Proteomics mass spectrometry analysis

For proteomics mass spectrometry data, proteins are inconsistently
fragmented into smaller peptides of varying sizes. The peptides are
usually separated on a chromatography column, from which aliquot
fractions are taken and measured by mass spectrometry. The total signal
derived from the original protein is therefore some combination of the
measured peptide parts.

In some upstream data processing tools, such as Proteomics Discoverer,
and PEAKS, the peptide data may be annotated with observed modification
events (PTM). In this scenario, peptide measurements are split across
multiple rows of data, where each row represents an observed combination
of peptide and PTMs.

### Collapse methods

It is fairly straightforward to observe peptide-PTM measurement data is
correlated with overall protein quantification, and that the specific
combination of peptide fragments may be inconsistent across samples.
That is, one may observe five peptides of protein A in one sample, and
may observe seven peptides of protein A in another sample. The
quantities of each peptide may be inconsistent, due to variability in
protein fragmentation across samples. However, the general sum of
peptide measurements is typically fairly stable across samples,
especially for proteins of moderate to high abundance which are known to
have stable abundance per cell.

Choice of method to collapse measurements is not trivial, and is
therefore configurable. In general, proteomics abundances are analyzed
after `log2( 1 + x )` transformation. However, measurements cannot be
summed in log2 form, which would be equivalent to multiplying
measurements in normal form. Measurements can be summed but only after
exponentiating the data, for example the reciprocal `( 2 ^ x ) - 1` is
sufficient.

## See also

Other jamses SE utilities:
[`geomx_to_se()`](https://jmw86069.github.io/jamses/reference/geomx_to_se.md),
[`make_se_test()`](https://jmw86069.github.io/jamses/reference/make_se_test.md),
[`se_collapse_by_column()`](https://jmw86069.github.io/jamses/reference/se_collapse_by_column.md),
[`se_detected_rows()`](https://jmw86069.github.io/jamses/reference/se_detected_rows.md),
[`se_normalize()`](https://jmw86069.github.io/jamses/reference/se_normalize.md),
[`se_rbind()`](https://jmw86069.github.io/jamses/reference/se_rbind.md),
[`se_to_assay_data()`](https://jmw86069.github.io/jamses/reference/se_to_assay_data.md),
[`se_to_assay_names()`](https://jmw86069.github.io/jamses/reference/se_to_assay_names.md),
[`se_to_rowcoldata()`](https://jmw86069.github.io/jamses/reference/se_to_rowcoldata.md)
