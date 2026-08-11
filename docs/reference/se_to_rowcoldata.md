# Get SE colData and rowData

Get SE colData and rowData consistently for SummarizedExperiment,
SingleCellExperiment, Seurat data, and generic Biobase `eSet` compatible
objects that provide
[`featureData()`](https://rdrr.io/pkg/Biobase/man/featureData.html) and
[`phenoData()`](https://rdrr.io/pkg/Biobase/man/phenoData.html).

## Usage

``` r
se_to_rowcoldata(se, verbose = FALSE, ...)
```

## Arguments

- se:

  recognized object:

  - `SummarizedExperiment` or any object that inherits from this class

  - `SingleCellExperiment` because it inherits from
    `SummarizedExperiment`

  - `Seurat` which is converted to `SingleCellExperiment`

  - `ExpressionSet` and any object that inherits from this class, using
    [`Biobase::featureData()`](https://rdrr.io/pkg/Biobase/man/featureData.html)
    and
    [`Biobase::phenoData()`](https://rdrr.io/pkg/Biobase/man/phenoData.html).

  - `NanoStringGeoMxSet` which uses
    [`GeomxTools::sData()`](https://rdrr.io/pkg/NanoStringNCTools/man/NanoStringRccSet-class.html)
    or
    [`NanoStringNCTools::sData()`](https://rdrr.io/pkg/NanoStringNCTools/man/NanoStringRccSet-class.html)
    to define `colData_se`, otherwise is handled equivalent to
    `ExpressionSet`.

- ...:

  additional arguments are ignored

## Value

`list` with two components:

- `"colData_se"` as a `data.frame` with column metadata from `se`

- `"rowData_se"` as a `data.frame` with row metadata from `se`

## Details

This function provides a straightforward way to return the equivalent of
`data.frame(check.names=FALSE, rowData(se))` and
`data.frame(check.names=FALSE, colData(se))` for several types of object
types.

- It also defines [`rownames()`](https://rdrr.io/r/base/colnames.html)
  and [`colnames()`](https://rdrr.io/r/base/colnames.html) if either are
  missing.

- When `rowData` has no annotation columns, it defines one column
  `"rows"` using `rownames(se)`.

- If slot name `"rowRanges"` exists, and `"rowData"` either does not
  exist or has zero columns, it will use `rowRanges()`.

- When `colData` has no annotation columns, it defines one column
  `"columns"` using `colnames(se)`.

- When input `se` class is `"Seurat"`, it converts the object with
  [`Seurat::as.SingleCellExperiment()`](https://satijalab.org/seurat/reference/as.SingleCellExperiment.html)

- Any other object uses
  [`Biobase::featureData()`](https://rdrr.io/pkg/Biobase/man/featureData.html)
  for `"rowData_se"`, and
  [`Biobase::phenoData()`](https://rdrr.io/pkg/Biobase/man/phenoData.html)
  for `"colData_se"`.

To verify the logic used at each step, set `verbose=TRUE`.

For Class `"NanoStringGeoMxSet"` defined in `GeomxTools`, the
`colData_se` is defined using
[`GeomxTools::sData()`](https://rdrr.io/pkg/NanoStringNCTools/man/NanoStringRccSet-class.html)
in order to return the combined `data.frame` with
[`protocolData()`](https://rdrr.io/pkg/Biobase/man/protocolData.html)
and [`pData()`](https://rdrr.io/pkg/Biobase/man/phenoData.html)
together. Accordingly, it is possible to have duplicated colnames, which
becomes a problem for many downstream tools, so the second instance of
any duplicated colnames have `"_v1"` added by using
`jamba::makeNames(x, renameFirst=FALSE)`. The first instance of each
duplicated colname is not renamed. Use `verbose=TRUE` to confirm when
duplicated columns are detected and renamed.

## See also

Other jamses SE utilities: [`geomx_to_se()`](geomx_to_se.md),
[`make_se_test()`](make_se_test.md),
[`se_collapse_by_column()`](se_collapse_by_column.md),
[`se_collapse_by_row()`](se_collapse_by_row.md),
[`se_detected_rows()`](se_detected_rows.md),
[`se_normalize()`](se_normalize.md), [`se_rbind()`](se_rbind.md),
[`se_to_assay_data()`](se_to_assay_data.md),
[`se_to_assay_names()`](se_to_assay_names.md)
