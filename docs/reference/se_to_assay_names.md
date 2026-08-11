# Get SE assay names

Get SE assay names consistently for `SummarizedExperiment`,
`SingleCellExperiment`, `Seurat` data, `NanoStringGeoMxSet`, and generic
Biobase `eSet` compatible objects that provide
[`assayData()`](https://rdrr.io/pkg/Biobase/man/assayData.html).

## Usage

``` r
se_to_assay_names(se, ...)
```

## Arguments

- se:

  `SummarizedExperiment` or other recognized data type inheriting from
  either `SummarizedExperiment`, `ExpressionSet`, `eSet`, or `Seurat`.

- ...:

  additional arguments are ignored.

## Value

`character` vector with corresponding assay names

## See also

Other jamses SE utilities: [`geomx_to_se()`](geomx_to_se.md),
[`make_se_test()`](make_se_test.md),
[`se_collapse_by_column()`](se_collapse_by_column.md),
[`se_collapse_by_row()`](se_collapse_by_row.md),
[`se_detected_rows()`](se_detected_rows.md),
[`se_normalize()`](se_normalize.md), [`se_rbind()`](se_rbind.md),
[`se_to_assay_data()`](se_to_assay_data.md),
[`se_to_rowcoldata()`](se_to_rowcoldata.md)
