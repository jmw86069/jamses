# Convert NanoStringGeoMxSet to SummarizedExperiment

Convert NanoStringGeoMxSet to SummarizedExperiment, combining sData and
phenoData into colData

## Usage

``` r
geomx_to_se(x, assay_names = NULL, verbose = FALSE, ...)
```

## Arguments

- x:

  `NanoStringGeoMxSet` object

- assay_names:

  `character` vector with one or more assay names. The default NULL will
  use all available assay names recognized by
  [`se_to_assay_names()`](https://jmw86069.github.io/jamses/reference/se_to_assay_names.md).

- ...:

  additional arguments are ignored.

## Value

`SummarizedExperiment` with `rowData`, `colData`, and `assays`.

## Details

The main purpose of this function is to convert
`GeomxTools::NanoStringGeoMxSet-class` to `SummarizedExperiment` while
also handling the oddity that the GeoMx data stores some sample (column)
annotation in slot `sData` and others in slot `phenoData`. These two
tables are combined into one `colData` entry in the output object.

Each data matrix in `assayData` from GeoMx is stored into the
corresponding `assays` entry in the output object.

## See also

Other jamses SE utilities:
[`make_se_test()`](https://jmw86069.github.io/jamses/reference/make_se_test.md),
[`se_collapse_by_column()`](https://jmw86069.github.io/jamses/reference/se_collapse_by_column.md),
[`se_collapse_by_row()`](https://jmw86069.github.io/jamses/reference/se_collapse_by_row.md),
[`se_detected_rows()`](https://jmw86069.github.io/jamses/reference/se_detected_rows.md),
[`se_normalize()`](https://jmw86069.github.io/jamses/reference/se_normalize.md),
[`se_rbind()`](https://jmw86069.github.io/jamses/reference/se_rbind.md),
[`se_to_assay_data()`](https://jmw86069.github.io/jamses/reference/se_to_assay_data.md),
[`se_to_assay_names()`](https://jmw86069.github.io/jamses/reference/se_to_assay_names.md),
[`se_to_rowcoldata()`](https://jmw86069.github.io/jamses/reference/se_to_rowcoldata.md),
[`simulate_se_test()`](https://jmw86069.github.io/jamses/reference/simulate_se_test.md)
