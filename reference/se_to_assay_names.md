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

Other jamses SE utilities:
[`geomx_to_se()`](https://jmw86069.github.io/jamses/reference/geomx_to_se.md),
[`make_se_test()`](https://jmw86069.github.io/jamses/reference/make_se_test.md),
[`se_collapse_by_column()`](https://jmw86069.github.io/jamses/reference/se_collapse_by_column.md),
[`se_collapse_by_row()`](https://jmw86069.github.io/jamses/reference/se_collapse_by_row.md),
[`se_detected_rows()`](https://jmw86069.github.io/jamses/reference/se_detected_rows.md),
[`se_normalize()`](https://jmw86069.github.io/jamses/reference/se_normalize.md),
[`se_rbind()`](https://jmw86069.github.io/jamses/reference/se_rbind.md),
[`se_to_assay_data()`](https://jmw86069.github.io/jamses/reference/se_to_assay_data.md),
[`se_to_rowcoldata()`](https://jmw86069.github.io/jamses/reference/se_to_rowcoldata.md),
[`simulate_se_test()`](https://jmw86069.github.io/jamses/reference/simulate_se_test.md)
