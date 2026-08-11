# Normalize SummarizedExperiment data

Normalize SummarizedExperiment data

## Usage

``` r
se_normalize(
  se,
  method = c("quantile", "jammanorm", "limma_batch_adjust", "TMM", "TMMwsp", "RLE"),
  assay_names = NULL,
  output_method_prefix = NULL,
  output_assay_names = NULL,
  genes = NULL,
  samples = NULL,
  params = list(quantile = list(ties = TRUE), jammanorm = list(controlGenes = NULL,
    minimum_mean = 0, controlSamples = NULL, centerGroups = NULL, useMedian = FALSE,
    noise_floor = NULL, noise_floor_value = NULL), limma_batch_adjust = list(batch =
    NULL, group = NULL), TMM = list(refColumn = NULL, logratioTrim = 0.3, sumTrim = 0.05,
    doWeighting = TRUE, Acutoff = NULL), TMMwsp = list(refColumn = NULL, logratioTrim =
    0.3, sumTrim = 0.05, doWeighting = TRUE, Acutoff = NULL), RLE = list(refColumn =
    NULL, logratioTrim = 0.3, 
     sumTrim = 0.05, doWeighting = TRUE, Acutoff = NULL)),
  normgroup = NULL,
  floor = 0,
  enforce_norm_floor = TRUE,
  output_sep = "_",
  override = TRUE,
  populate_mcols = TRUE,
  verbose = FALSE,
  ...
)
```

## Arguments

- se:

  `SummarizedExperiment` object

- method:

  `character` vector indicating which normalization method(s) to apply.

  - `"quantile"`: quantile normalization via
    [`limma::normalizeQuantiles()`](https://rdrr.io/pkg/limma/man/normalizequantiles.html)

  - `"jammanorm"`: log-ratio normalization via
    [`jamma::jammanorm()`](https://rdrr.io/pkg/jamma/man/jammanorm.html)

  - `"limma_batch_adjust"`: batch adjustment via
    [`limma::removeBatchEffect()`](https://rdrr.io/pkg/limma/man/removeBatchEffect.html),
    recommended for data visualization, but not recommended for
    downstream statistical comparisons.

  - `"TMM"`: trimmed mean of M-values via
    [`edgeR::calcNormFactors()`](https://rdrr.io/pkg/edgeR/man/calcNormFactors.html)

  - `"TMMwsp"`: TMM with singleton pairing via
    [`edgeR::calcNormFactors()`](https://rdrr.io/pkg/edgeR/man/calcNormFactors.html)

  - `"RLE"`: relative log expression via
    [`edgeR::calcNormFactors()`](https://rdrr.io/pkg/edgeR/man/calcNormFactors.html)

- assay_names:

  `character` vector or one or more `names(assays(se))` that indicates
  which numeric matrix to use during normalization. When multiple values
  are provided, each matrix is normalized independently by each
  `method`.

- output_method_prefix:

  `character` vector (optional) with custom method prefix values to use
  when creating the new `assay_name` for each normalization. It must
  have length equal to `length(method)`, to be applied to each method in
  order. Note that `output_assay_names` takes priority, and when it is
  defined the `output_method_prefix` entries are ignored.

  Consider these arguments:

      assay_name="counts",
      method="limma_batch_adjust",
      output_method_prefix="lba"

  The assay_name created during normalization will be `"lba_counts"`.

- output_assay_names:

  `character` vector (optional) which overrides the default method for
  defining assay names for normalized data. This vector length must
  equal `length(method) * length(assay_names)`, and will be applied in
  the order data is normalized:

  1.  `assay_names` are iterated.

  2.  For each value in `assay_names`, each normalization in `method` is
      applied.

  Therefore the order of `output_assay_names` could follow this order:
  `method1_assay1`, `method1_assay2`, `method2_assay1`,
  `method2_assay2`.

- genes:

  `character` vector (optional) used to define a subset of gene rows in
  `se` to use for normalization. Values must match `rownames(se)`.

- samples:

  `character` vector (optional) used to define a subset of sample
  columns in `se` to use for normalization. Values must match
  `colnames(se)`.

- params:

  `list` (optional) parameters specific to each normalization method,
  passed to [`matrix_normalize()`](matrix_normalize.md). Any value which
  is not defined in the `params` provided will use the default value in
  [`matrix_normalize()`](matrix_normalize.md), for example
  `params=list(jammanorm=list(minimum_mean=2))` will use
  `minimum_mean=2` then use other default values relevant to the
  `jammanorm` normalization method.

- normgroup:

  `character` or equivalent vector that defines subgroups of `samples`
  to be normalized indendently of each normgroup, or `character` vector
  with one or more colnames in `colData(se)`.

  - When `normgroup` is `NULL` all data is normalized together by
    default.

  - When supplying colnames in `colData(se)`, values are combined using
    [`jamba::pasteByRow()`](https://jmw86069.github.io/jamba/reference/pasteByRow.html)
    which uses `_` underscore delimiter between non-empty fields in each
    column.

  - When supplying a vector of normgroup values, it should be equal to
    `length(samples)` or should be named using values in `samples` such
    that all `samples` are present in `names(normgroup)`.

  - Note that when data are normalized in independent normgroups, these
    normgroups are not normalized relative to each other. Data are only
    normalized within each independent normgroup.

- output_sep:

  `character` string used as a delimited between the `method` and the
  `assay_names` to define the output assay name, for example when
  `assay_name="counts"`, `method="quantile"`, and `output_sep="_"` the
  new assay name will be `"quantile_counts"`.

- override:

  `logical` indicating whether to override any pre-existing matrix
  values with the same output assay name. When `override=FALSE` and the
  output assay name already exists, the normalization will not be
  performed.

- populate_mcols:

  `logical`, default TRUE, whether to populate normalization details
  into `mcols(assays(se))`, including the normalization `method`, the
  source `assay_name` used during normalization, and values from
  `params`.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are passed to
  [`matrix_normalize()`](matrix_normalize.md).

## Value

`SummarizedExperiment` object where the normalized output is added to
`assays(se)` using the naming format `method_assayname`.

## Details

This function applies one or more data normalization methods to an input
`SummarizedExperiment` object. The normalization is applied to one or
more matrix data stored in `assays(se)`, each one is run independently.

Note that supplying `genes` and `samples` will apply normalization to
only those `genes` and `samples`, and this data will be stored in the
full `SummarizedExperiment` object `se` with `NA` values used to fill
any values not present in `genes` or `samples`.

For example if `assay_names` contains two assay names, and `method`
contains two methods, the output will include four normalizations, where
each assay name is normalized two ways. The output assay names will be
something like `"assay1_method1"`, `"assay1_method2"`,
`"assay2_method1"`, `"assay2_method2"`. It is not always necessary to
normalize data by multiple different methods, however when two methods
are similar and need to be compared, the `SummarizedExperiment` object
is a convenient place to store different normalization results for
downstream comparison. Further, the method
[`se_contrast_stats()`](se_contrast_stats.md) is able to apply
equivalent statistical contrasts to each normalization, and returns an
array of statistical hits which is convenient for direct comparison of
results.

This method calls [`matrix_normalize()`](matrix_normalize.md) to perform
each normalization step, see that function description for details on
each method.

## See also

Other jamses SE utilities: [`geomx_to_se()`](geomx_to_se.md),
[`make_se_test()`](make_se_test.md),
[`se_collapse_by_column()`](se_collapse_by_column.md),
[`se_collapse_by_row()`](se_collapse_by_row.md),
[`se_detected_rows()`](se_detected_rows.md),
[`se_rbind()`](se_rbind.md),
[`se_to_assay_data()`](se_to_assay_data.md),
[`se_to_assay_names()`](se_to_assay_names.md),
[`se_to_rowcoldata()`](se_to_rowcoldata.md)

## Examples

``` r
if (jamba::check_pkg_installed("farrisdata")) {

   # se_normalize
   # suppressPackageStartupMessages(library(SummarizedExperiment))
   GeneSE <- farrisdata::farrisGeneSE;
   samples <- colnames(GeneSE);
   genes <- rownames(GeneSE);

   GeneSE <- se_normalize(GeneSE,
      genes=genes,
      samples=samples,
      assay_names=c("raw_counts", "counts"),
      method="jammanorm",
      params=list(jammanorm=list(minimum_mean=5)))
   SummarizedExperiment::mcols(SummarizedExperiment::assays(GeneSE))
   names(SummarizedExperiment::assays(GeneSE))

   # review normalization factor values
   round(digits=3, attr(
      SummarizedExperiment::assays(GeneSE)$jammanorm_raw_counts, "nf"))

   # the data in "counts" was already normalized
   # so the normalization factors are very near 0 as expected
   round(digits=3,
      attr(SummarizedExperiment::assays(GeneSE)$jammanorm_counts, "nf"))


   # note that housekeeper genes are supplied in params
   # also this demonstrates output_method_prefix
   set.seed(123);
   hkgenes <- sample(rownames(GeneSE), 1000)
   GeneSE <- se_normalize(GeneSE,
      genes=genes,
      samples=samples,
      assay_names=c("raw_counts"),
      method="jammanorm",
      output_method_prefix="hkjammanorm",
      params=list(jammanorm=list(minimum_mean=5,
         controlGenes=hkgenes)))
   SummarizedExperiment::mcols(SummarizedExperiment::assays(GeneSE))

   # example showing quantile normalization
   GeneSE <- se_normalize(GeneSE,
      assay_names=c("raw_counts"),
      method="quantile")
   SummarizedExperiment::mcols(SummarizedExperiment::assays(GeneSE))

   # example showing quantile normalization with custom output_assay_names
   GeneSE <- se_normalize(GeneSE,
      assay_names=c("raw_counts"),
      method="quantile",
      output_assay_names="newquantile_raw_counts")
   SummarizedExperiment::mcols(SummarizedExperiment::assays(GeneSE))
}
#> DataFrame with 7 rows and 6 columns
#>                                    assay_name normalization_method
#>                                   <character>          <character>
#> counts                                 counts                   NA
#> raw_counts                         raw_counts                   NA
#> jammanorm_raw_counts     jammanorm_raw_counts            jammanorm
#> jammanorm_counts             jammanorm_counts            jammanorm
#> hkjammanorm_raw_counts hkjammanorm_raw_counts            jammanorm
#> quantile_raw_counts       quantile_raw_counts             quantile
#> newquantile_raw_counts newquantile_raw_counts             quantile
#>                        source_assay_name minimum_mean              controlGenes
#>                              <character>    <numeric>                    <list>
#> counts                                NA           NA                        NA
#> raw_counts                            NA           NA                        NA
#> jammanorm_raw_counts          raw_counts            5                        NA
#> jammanorm_counts                  counts            5                        NA
#> hkjammanorm_raw_counts        raw_counts            5 Aldh3b1,Gm8194,Gm7776,...
#> quantile_raw_counts           raw_counts           NA                        NA
#> newquantile_raw_counts        raw_counts           NA                        NA
#>                             ties
#>                        <logical>
#> counts                        NA
#> raw_counts                    NA
#> jammanorm_raw_counts          NA
#> jammanorm_counts              NA
#> hkjammanorm_raw_counts        NA
#> quantile_raw_counts         TRUE
#> newquantile_raw_counts      TRUE

```
