# Collapse SummarizedExperiment data by column

Collapse SummarizedExperiment data by column

## Usage

``` r
se_collapse_by_column(
  se,
  columns = colnames(se),
  column_groups,
  assay_names = NULL,
  colDataColnames = colnames(SummarizedExperiment::colData(se)),
  keepNULLlevels = FALSE,
  groupFunc = jamba::rowGroupMeans,
  noise_floor = 0,
  noise_floor_value = 0,
  rmOutliers = FALSE,
  madFactor = 5,
  useMedian = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- se:

  `SummarizedExperiment` object

- columns:

  `character` vector of `colnames(se)` to include in the process.

- column_groups:

  `character` vector of column groupings, or `character` vector of
  `colnames(colData(se))` used to define the column groupings.

- assay_names:

  `character` vector with one or more `assayNames(se)` to apply the
  column grouping calculation defined in `groupFunc`. By default, all
  assay names in `assayNames(se)` are used.

- colDataColnames:

  `character` vector of `colData(se)` colnames to be included in the
  returned `SummarizedExperiment` after the column grouping. This
  argument is used to subset the columns, in cases where some columns do
  not need to be combined and returned in the output data.

- keepNULLlevels:

  `logical` indicating whether to return empty columns when there are
  not factor levels present in the data. This option is intended when
  `column_group` references a `factor` type, whose factor levels are not
  present in the current data, using `columns`. When
  `keepNULLlevels=TRUE` any missing levels will be present with `NA`
  values, which can be helpful for generating a consistent output.

- groupFunc:

  `function` used to perform row group calculations on a `numeric`
  matrix. The default is passed to
  [`jamba::rowGroupMeans()`](https://jmw86069.github.io/jamba/reference/rowGroupMeans.html),
  but can be substituted with another row-based function. It must accept
  arguments `x` and `groups`, but the other arguments are passed only if
  `groupFunc` permits these argument names, or `...`:

  - `x` as a `numeric` matrix (required),

  - `groups` as a `character` vector of column groups, in order of
    `colnames(x)` (required)

  - `rmOutliers` a `logical` indicating whether to apply outlier
    removal, though the function can ignore this value (optional).

  - `madFactor` a `numeric` value indicating the MAD threshold used when
    `rmOutliers=TRUE`; though again the function can ignore this value
    (optional).

  - `useMedian=FALSE` is `logical` and when `useMedian=FALSE` it
    disables calculating the
    [`median()`](https://rdrr.io/r/stats/median.html) value per group,
    and instead takes the group
    [`mean()`](https://rdrr.io/r/base/mean.html) value.

  - `...` additional arguments in `...` will be passed only if permitted
    by `groupFunc`.

- noise_floor:

  `numeric` value indicating the minimum numeric value permitted, *at or
  below* this value will be replaced with `noise_floor_value`. The
  default value `noise_floor=0` will therefore change all values at or
  below zero to `noise_floor_value=0` by default. Another alternative is
  to change abnormally low values such as zero `0` to `NA` so these
  values are not treated as actual measurements during the group summary
  calculation. This value and the replacement should be adjusted with
  caution. Use `noise_floor=NULL` or `noise_floor=-Inf` to disable this
  step.

- noise_floor_value:

  `numeric` or `NA` used as a replacement for `numeric` values *at or
  below* `noise_floor`, which occurs prior to calling the `groupFunc`
  summary calculation.

- rmOutliers, madFactor:

  `logical` and `numeric`, respectively, passed to `groupFunc` which by
  default is
  [`jamba::rowGroupMeans()`](https://jmw86069.github.io/jamba/reference/rowGroupMeans.html).

- useMedian:

  `logical` passed to argument `groupFunc()`, intended to be used by
  [`jamba::rowGroupMeans()`](https://jmw86069.github.io/jamba/reference/rowGroupMeans.html)
  to specify taking the mean and not the median value per row group.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are passed through `groupFunc`.

## Value

`SummarizedExperiment` object with these changes:

- columns will be collapsed by `column_groups`, for each `assays(se)`
  `numeric` matrix defined by `assay_names`.

- `colData(se)` will also be collapsed by
  [`shrinkDataFrame()`](https://jmw86069.github.io/jamses/reference/shrinkDataFrame.md)
  to combine unique values from each column annotation.

## Details

Purpose is to collapse columns of a `SummarizedExperiment` object, where
measurements for a given entity, usually a gene, are split across
multiple rows in the source data. The output of this function should be
measurements appropriately summarized to the gene level.

The driving use case is slightly different than with
[`se_collapse_by_row()`](https://jmw86069.github.io/jamses/reference/se_collapse_by_row.md),
in this case the function is mostly convenient method to calculate group
mean values in context of a `SummarizedExperiment` object, so it can be
used with
[`jamses::heatmap_se()`](https://jmw86069.github.io/jamses/reference/heatmap_se.md)
for example.

This function retains associated column annotations `colData(se)`, after
combining multiple values in an appropriate manner.

Optionally, this function will detect and remove individual outlier
values before calculating the group mean.

## See also

Other jamses SE utilities:
[`geomx_to_se()`](https://jmw86069.github.io/jamses/reference/geomx_to_se.md),
[`make_se_test()`](https://jmw86069.github.io/jamses/reference/make_se_test.md),
[`se_collapse_by_row()`](https://jmw86069.github.io/jamses/reference/se_collapse_by_row.md),
[`se_detected_rows()`](https://jmw86069.github.io/jamses/reference/se_detected_rows.md),
[`se_normalize()`](https://jmw86069.github.io/jamses/reference/se_normalize.md),
[`se_rbind()`](https://jmw86069.github.io/jamses/reference/se_rbind.md),
[`se_to_assay_data()`](https://jmw86069.github.io/jamses/reference/se_to_assay_data.md),
[`se_to_assay_names()`](https://jmw86069.github.io/jamses/reference/se_to_assay_names.md),
[`se_to_rowcoldata()`](https://jmw86069.github.io/jamses/reference/se_to_rowcoldata.md)
