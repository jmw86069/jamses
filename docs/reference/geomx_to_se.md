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
  [`se_to_assay_names()`](se_to_assay_names.md).

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

Other jamses SE utilities: [`make_se_test()`](make_se_test.md),
[`se_collapse_by_column()`](se_collapse_by_column.md),
[`se_collapse_by_row()`](se_collapse_by_row.md),
[`se_detected_rows()`](se_detected_rows.md),
[`se_normalize()`](se_normalize.md), [`se_rbind()`](se_rbind.md),
[`se_to_assay_data()`](se_to_assay_data.md),
[`se_to_assay_names()`](se_to_assay_names.md),
[`se_to_rowcoldata()`](se_to_rowcoldata.md)
