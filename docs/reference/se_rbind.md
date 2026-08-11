# Combine SummarizedExperiment objects by row

Combine SummarizedExperiment objects by row, using
[`rbind()`](https://rdrr.io/r/base/cbind.html) logic.

## Usage

``` r
se_rbind(
  se_list,
  colnames_from = "_(n|p|neg|pos)_",
  colnames_to = "_X_",
  colnames_keep = NULL,
  colData_action = c("identical", "all"),
  colData_sep = ";",
  verbose = FALSE,
  ...
)
```

## Arguments

- se_list:

  `list` of `SummarizedExperiment` objects.

- colnames_from:

  `character` vector of patterns used with
  [`gsub()`](https://rdrr.io/r/base/grep.html) to convert
  [`colnames()`](https://rdrr.io/r/base/colnames.html) for each object
  in `se_list` to an identifier that will be shared across all objects
  in `se_list`.

- colnames_to:

  `character` vector of replacements used with
  [`gsub()`](https://rdrr.io/r/base/grep.html) alongside each entry in
  `colnames_from` to convert
  [`colnames()`](https://rdrr.io/r/base/colnames.html) for each object
  in `se_list` to an identifier that will be shared across all objects
  in `se_list`.

- colData_action:

  `character` string indicating the action used to combine
  [`colData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  across `se_list`:

  - `"identical"`: retain only those columns in
    [`colData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
    which are identical in all `se_list` objects.

  - `"all"`: retain all columns, but convert columns with mismatched
    values to store comma-delimited values.

- colData_sep:

  `character` string used as delimiter when `colData_action="all"` and
  when values in a column in
  [`colData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  differs across objects in `se_list`. Only values that differ are
  delimited, to minimize redundancy.

- ...:

  additional arguments are ignored.

## Value

`SummarizedExperiment` object whose
[`colData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
has been processed according to `colData_action` - either keeping only
columns with identical values, or keeping all values delimited as a
`character` string when values differ.

## Details

This function is intended to help the process of calling
`SummarizedExperiment::rbind()`.

The process:

1.  Convert [`colnames()`](https://rdrr.io/r/base/colnames.html) for
    each entry in `se_list` using `colnames_from` and `colnames_to`.
    This step is useful when each object in `se_list` may be using a
    different set of
    [`colnames()`](https://rdrr.io/r/base/colnames.html). For example
    `"sample_p_12"` and `"sample_n_12"` might be equivalent, so renaming
    them with `colnames_from=c("_[np]_")` and `colnames_to=c("_X_")`
    would convert both values to `"sample_X_12"`.

2.  Subset each object in `se_list` using only shared
    [`colnames()`](https://rdrr.io/r/base/colnames.html).

3.  Determine how to handle
    [`colData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
    columns that are not identical:

    - `colData_action="identical"`: will only keep columns whose values
      are identical across all objects in `se_list`.

    - `colData_action="all"`: will keep columns in
      [`colData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html),
      however non-identical columns will be converted to `character` and
      values will be comma-delimited.

4.  Perform [`rbind()`](https://rdrr.io/r/base/cbind.html).

### TODO:

- Write equivalent `se_cbind()` - it will wait until there is a driving
  use case.

- Consider retaining only shared
  [`assayNames()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
  across `se_list`.

- Consider optionally retaining user-defined
  [`assayNames()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html).
  (Alternatively, the user can subset the assayNames upfront, though it
  might be tedious). The recommended pattern in that case:

    se <- se_rbind(
       se_list=lapply(se_list, function(se){
          assays(se) <- assays(se)[assay_names];
          return(se)
       })
    )

## See also

Other jamses SE utilities: [`geomx_to_se()`](geomx_to_se.md),
[`make_se_test()`](make_se_test.md),
[`se_collapse_by_column()`](se_collapse_by_column.md),
[`se_collapse_by_row()`](se_collapse_by_row.md),
[`se_detected_rows()`](se_detected_rows.md),
[`se_normalize()`](se_normalize.md),
[`se_to_assay_data()`](se_to_assay_data.md),
[`se_to_assay_names()`](se_to_assay_names.md),
[`se_to_rowcoldata()`](se_to_rowcoldata.md)

## Examples

``` r
m1 <- matrix(rnorm(100), ncol=10);
colnames(m1) <- paste0("sample_p_", 1:10);
rownames(m1) <- paste0("row_", 1:10);
m2 <- matrix(rnorm(100), ncol=10);
colnames(m2) <- paste0("sample_n_", 1:10);
rownames(m2) <- paste0("row_", 11:20);
sample_id <- gsub("_[np]_", "_X_", colnames(m1));
m1
#>        sample_p_1 sample_p_2   sample_p_3 sample_p_4 sample_p_5 sample_p_6
#> row_1   0.8376299 -1.1661847 -0.003580308  1.3052615 -1.2621945 -0.1488164
#> row_2  -1.4434929  0.4079462 -1.495826814  0.8760961 -0.5518704  0.1124638
#> row_3  -0.2085702 -0.8630042 -0.768417027  0.4637961 -1.1827995  0.7246762
#> row_4  -0.4385635  0.3040420  0.408488505  0.4771142  0.6206636 -1.1874861
#> row_5  -0.2185938 -0.1464275  1.900136335 -0.4914053  0.4463130 -0.4996002
#> row_6   1.4599659 -1.4335622  0.110009123 -1.3193853  0.4218847 -1.0736430
#> row_7  -0.5820599 -0.7906078  1.140386825  1.2954258  0.4424648  1.0572402
#> row_8  -0.7830976  0.8851125  0.768081305 -1.4202195  0.5572457  1.2790726
#> row_9  -1.5196540  0.9030761 -1.168091622 -0.9388959  0.6393565  0.7876767
#> row_10 -0.8056981  2.0055733 -0.171112652  0.6289650 -1.9686616 -1.2224034
#>        sample_p_7  sample_p_8 sample_p_9  sample_p_10
#> row_1   0.4519521  0.57833494  0.5084848 -0.003988944
#> row_2   1.1504492  1.36467278 -0.1163584  0.847842769
#> row_3   0.1679410 -1.70157980  0.9255461 -0.100116526
#> row_4  -0.5661093 -0.28067628  0.6482298 -0.279629907
#> row_5  -1.0861182  0.06506802 -0.1502094  0.784438245
#> row_6  -0.6653028  0.57858929  1.0403770 -1.584616645
#> row_7   0.7148484 -1.16920662  0.2925587  0.478366148
#> row_8  -0.4316611  0.80618486  0.6687514  0.393566373
#> row_9   0.2276149  0.30739008 -0.5941776 -2.695329369
#> row_10  1.2949458  0.26380601  1.5804318  0.368377329
m2
#>          sample_n_1  sample_n_2 sample_n_3 sample_n_4 sample_n_5    sample_n_6
#> row_11 -2.168417747 -1.74102202  0.8687933 -0.0768659  1.1025652 -0.5870856158
#> row_12  0.659804377 -1.99258577  1.3693517  0.6873641 -0.5766189  0.0007641864
#> row_13 -0.453913733  0.55127421  0.7626511  0.1716315 -1.8516917  2.2144653193
#> row_14 -0.694936825 -0.03474206  0.4211472 -0.8301086 -0.1128632  0.9694343957
#> row_15 -0.006846303  1.85057170 -0.8682240 -0.2901591  1.3210693  0.7680077137
#> row_16  1.373052045  0.57367511  0.7295604 -1.3191257  0.6622543 -1.1083279118
#> row_17 -0.635323077  0.84969589  0.5002659 -0.9670319  0.4413832 -0.7862359200
#> row_18  0.558103294  1.33438359  0.6342503 -0.1446111  1.1837459  2.2841164803
#> row_19  0.341157868 -0.50071910  0.4236450 -1.7981326 -0.7715014 -1.0933007640
#> row_20 -1.179518629  0.51009793 -0.2018380 -1.6885425  0.7296892  0.2144793753
#>        sample_n_7 sample_n_8   sample_n_9 sample_n_10
#> row_11  0.8925711  0.3661144 -0.275890475  1.53242362
#> row_12  1.0187580 -0.8747814  0.682315245 -1.35799783
#> row_13  1.0891120  1.0244749 -0.117290715 -0.19961905
#> row_14 -0.1631290  0.9047589 -0.344675864  0.63152313
#> row_15 -0.8209867 -0.2382487  0.111620498  1.76202090
#> row_16 -0.3072572 -1.5578549 -0.283405315  0.42601436
#> row_17 -0.9020980  0.7613099 -0.591017164 -0.01375342
#> row_18  0.6270687  1.1291444 -0.315936931 -0.30755691
#> row_19  1.1203550 -0.2951078 -0.008152152  0.41430816
#> row_20  2.1272136  0.5362428  0.207495141  0.98905792
se1 <- SummarizedExperiment::SummarizedExperiment(
   assays=list(counts=m1),
   rowData=data.frame(measurement=rownames(m1)),
   colData=data.frame(sample=colnames(m1),
      sample_id=sample_id))
se2 <- SummarizedExperiment::SummarizedExperiment(
   assays=list(counts=m2),
   rowData=data.frame(measurement=rownames(m2)),
   colData=data.frame(sample=colnames(m2),
      sample_id=sample_id))
# this step fails because colnames are not shared
# do.call(SummarizedExperiment::rbind, list(se1, se2))

# keep only identical colData columns
se12 <- se_rbind(list(se1, se2))
SummarizedExperiment::colData(se12)
#> DataFrame with 10 rows and 1 column
#>               sample_id
#>             <character>
#> sample_X_1   sample_X_1
#> sample_X_2   sample_X_2
#> sample_X_3   sample_X_3
#> sample_X_4   sample_X_4
#> sample_X_5   sample_X_5
#> sample_X_6   sample_X_6
#> sample_X_7   sample_X_7
#> sample_X_8   sample_X_8
#> sample_X_9   sample_X_9
#> sample_X_10 sample_X_10

# keep all colData columns
se12all <- se_rbind(list(se1, se2),
   colData_action="all")
SummarizedExperiment::colData(se12all)
#> DataFrame with 10 rows and 2 columns
#>                             sample   sample_id
#>                        <character> <character>
#> sample_X_1   sample_p_1;sample_n_1  sample_X_1
#> sample_X_2   sample_p_2;sample_n_2  sample_X_2
#> sample_X_3   sample_p_3;sample_n_3  sample_X_3
#> sample_X_4   sample_p_4;sample_n_4  sample_X_4
#> sample_X_5   sample_p_5;sample_n_5  sample_X_5
#> sample_X_6   sample_p_6;sample_n_6  sample_X_6
#> sample_X_7   sample_p_7;sample_n_7  sample_X_7
#> sample_X_8   sample_p_8;sample_n_8  sample_X_8
#> sample_X_9   sample_p_9;sample_n_9  sample_X_9
#> sample_X_10 sample_p_10;sample_n.. sample_X_10
```
