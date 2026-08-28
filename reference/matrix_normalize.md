# Normalize a numeric data matrix

Normalize a numeric data matrix

## Usage

``` r
matrix_normalize(
  x,
  method = c("quantile", "jammanorm", "limma_batch_adjust", "TMM", "TMMwsp", "RLE"),
  apply_log2 = c("ifneeded", "no", "always"),
  floor = 0,
  floor_value = floor,
  enforce_norm_floor = TRUE,
  params = list(quantile = list(ties = TRUE), jammanorm = list(controlGenes = NULL,
    minimum_mean = 0, controlSamples = NULL, centerGroups = NULL, useMedian = FALSE,
    noise_floor = NULL, noise_floor_value = NULL), limma_batch_adjust = list(batch =
    NULL, group = NULL), TMM = list(refColumn = NULL, logratioTrim = 0.3, sumTrim = 0.05,
    doWeighting = TRUE, Acutoff = NULL), TMMwsp = list(refColumn = NULL, logratioTrim =
    0.3, sumTrim = 0.05, doWeighting = TRUE, Acutoff = NULL)),
  normgroup = NULL,
  normgroup_rows = NULL,
  subset_columns = NULL,
  debug = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  `numeric` matrix with sample columns, and typically gene rows, but any
  measured assay row will meet the assumptions of the method.

- method:

  `character` string indicating which normalization method to apply.

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

- apply_log2:

  `character` string indicating whether to apply log2 transformation:
  `"ifneeded"` will apply log2 transform when any absolute value is
  greater than 40; `"no"` will not apply log2 transformation; `"always"`
  will apply log2 transform. Note the log2 transform is applied with
  `jamba::log2signed(x, offset=1)` which is equivalent to `log(1 + x)`
  except that negative values are also transformed using the absolute
  value, then multiplied by their original sign.

- floor:

  `numeric` value indicating the lowest accepted numeric value, below
  which values are assigned to this floor. The default `floor=0`
  requires all values are `0`, and any values below `0` are assigned
  `0`. Note that the `floor` is applied after log2 transform, when the
  log2 transform is performed.

- floor_value:

  `numeric` or `NA` used to replace values in `x` when they are at or
  below `floor`, default `floor_value=floor`, however it can be useful
  to assign `NA` to replace zero in circumstances when that is
  preferable.

- enforce_norm_floor:

  `logical` indicating whether to enforce the `floor` for the normalized
  results, default is `TRUE`. For example, when `floor=0` any values at
  or below `0` are set to `0` before normalization. After normalization
  some of these values will be above or below `0`. When
  `enforce_norm_floor=TRUE` these values will again be set to `0`
  because they are considered to be below the noise threshold of the
  protocol, and adjustments are not relevant; also any normalized values
  below the `floor` will also be set to `floor`.

- params:

  `list` of parameters relevant to the `method` of normalization. The
  `params` should be a `list` named by the `method`, whose values are a
  list named by the relevant method parameter. See examples.

- normgroup:

  `character` or equivalent vector that defines subgroups of `samples`
  to be normalized indendently of each normgroup.

  - When `normgroup` is `NULL` all data is normalized together as
    default.

  - The `normgroup` vector is expected to be in the order of
    `colnames(x)`.

  - Note that when data are normalized in independent normgroups, these
    normgroups are not normalized relative to each other. Data are only
    normalized within each independent normgroup.

- normgroup_rows:

  `list` with `character` values which must match `rownames(x)`, and
  where `names(normgroup_rows)` must also match the values in
  `normgroup`.

  - The purpose is to permit using independent rows for each normgroup,
    for example if the detected genes (rows) in each matrix differ for
    each normgroup, the normalization should only use each appropriate
    set of detected genes (rows).

  - When `normgroup_rows` are applied, data for rows not present in each
    particular normgroup will be set to `NA`, since those values were
    not included in the normalization for that normgroup.

- subset_columns:

  `integer` intended for internal use when `normgroup` is provided. This
  argument is used to instruct each normalization method to use an
  appropriate subset of `params` based upon the subset of columns being
  analyzed.

- verbose:

  `logical` indicating whether to print verbose output.

## Value

`numeric` matrix with the same dimensions as the input matrix `x`.
Additional information may be returned as `attributes(x)`:

- `"norm_method"`: a `character` string with the method used.

- `"nf"`: a `numeric` vector with normalization factors, returned only
  by `"jammanorm"`, `"TMM"`, and `"TMMwsp"`.

- `"hk"`: a `character` vector of `rownames(x)` used as housekeeper
  `controlGenes` by `"jammanorm"`.

## Details

This function is a wrapper for several relevant normalization methods
that operate on a numeric matrix.

## Normalization Methods Implemented:

### method='quantile'

Quantile-normalization performed by
[`limma::normalizeQuantiles()`](https://rdrr.io/pkg/limma/man/normalizequantiles.html).
This method has one parameter `"ties"` passed to
[`limma::normalizeQuantiles()`](https://rdrr.io/pkg/limma/man/normalizequantiles.html),
the default here `ties=TRUE` which handles tied numeric expression
values in a robust way to avoid unpredictability otherwise. This option
is especially relevant with expression count data, where integer counts
cause a large number of values to be represented multiple times.

### method='jammanorm'

Median-normalization performed by
[`jamma::jammanorm()`](https://rdrr.io/pkg/jamma/man/jammanorm.html).
This method shifts expression data as shown on MA-plots, so the median
expression is zero across all samples, using only the rows that meet the
relevant criteria.

Some relevant criteria to define rows used for normalization:

- `controlGenes` defines specific genes to use for normalization, such
  as housekeeper genes. It may also be useful to use detected genes
  here, so the normalization only considers those genes defined as
  detected by the protocol.

- `minimum_mean` sets a numeric threshold and requires the mean
  expression (shown on the x-axis of the MA-plot) to be at least this
  value.

Note that when both `controlGenes` and `minimum_mean` are defined, both
criteria are enforced. So the `controlGenes` are also required to have
expression of at least `minimum_mean`.

Also note that all rows of data are normalized by this method, only the
subset of rows defined by `controlGenes` and `minimum_mean` are used to
compute the normalization factor.

### method='limma_batch_adjust'

Batch adjustment performed by
[`limma::removeBatchEffect()`](https://rdrr.io/pkg/limma/man/removeBatchEffect.html)
which is intended to apply batch-adjustment as a form of normalization,
but which does not represent full normalization itself. There are two
relevant parameters: `"batch"` which is a vector of batch values in
order of `colnames(x)`, and `"group"` which is a vector of sample groups
in order of `colnames(x)`.

## Normalization groups via `normgroup`

The `normgroup` argument is intended as a convenient method to apply a
normalization method to each independent `normgroup`. This situation is
especially useful when a study contains multiple tissue types, or
multiple data types, that may not be appropriate to normalize directly
relative to one another.

For example, one could normalize total RNA-seq and nascent 4sU-seq data
independently, without expectation that the two would ever have a common
frame of reference to normalize one relative to another. However, both
may be amenable to `"quantile"` or `"jammanorm"` median normalization.

Similarly, one could normalize each tissue type independently, which may
be appropriate when analyzing data that contains very different
mammalian tissue organ samples, such as muscle and brain.

It would generally not be appropriate to use quantile normalization
across muscle and brain samples, since the overall pattern and
distribution of expression values is not expected to be similar.
Quantile normalize assumes (and imposes) a common distribution, by
adjusting mean expression signal at each quantile to a common mean
expression across all samples.

For a rough approximation of cross-tissue normalization, one could apply
`"quantile"` normalization within each `normgroup` defined by tissue
type, then apply `"jammanorm"` median normalization to apply a linear
adjustment of signal across tissue types. The median normalization does
not affect distribution, thus will not affect intra-tissue contrasts,
except by adjusting its overall signal which may change downstream
assumptions regarding signal thresholds.

It is recommended not to compare directly across tissue types. In some
cases a two-way contrast may be appropriate, where fold change within
one tissue type is compared to the fold change within another tissue
type. However, even in that case the two tissue types do not need to be
normalized relative to each other upfront - the within-tissue fold
change serves as one method of normalizing the observations across
tissue types.

## Other useful parameters

Note the `floor` and `enforce_norm_floor` have recommended default
values `floor=0` and `enforce_norm_floor=TRUE`.

- `floor` is applied prior to normalization, typically to minimize
  effects of low, noisy signal on the normalization process itself.
  Specifically, this floor is used to remove negative values, which may
  be by-products of upstream signal processing. "A measured signal at or
  below the noise floor of a platform technology is effectively the same
  as a signal at the noise floor."

- `enforce_norm_floor` is applied after normalization, typically as a
  convenience, also to prevent low, noisy signal from contributing to
  downstream analysis steps.

These defaults will set any assay value at or below `0` to `0`, and
after normalization any values whose input values were at or below `0`
will also be set to `0` to prevent normalizing a value of `0` to
non-zero. Any normalized value at or below `0` will also be set to `0`
to prevent results from containing negative normalized values.

The assumption for this default is that a value of zero is not a
measurement but represents the lack of a measurement. Similarly, the
intent of `floor` is a numeric threshold at or below there is no
confidence in the reported measurement, therefore values at or below
this threshold are treated as equivalent to the threshold for the
purpose of downstream analyses.

Some platforms like QPCR for example, have substantially lower
confidence at high CT values, where expression values using the equation
`2^(40-CT)` might impose a noise threshold at expression 32 or lower.
This noise threshold for QPCR means any expression measurement of 32 or
lower is as likely to be `32` as it is to be `2`, and therefore any
differences between reported expression of `32` and `2` should not be
considered relevant. Applying `floor=32` in this case accomplishes this
goal by setting all values at or below `32` to `32`. Of course when
using this method `matrix_normalize()` the data should be log2
transformed, which means the `floor` should also be log2 transformed,
e.g. `floor=log2(32)` which is `floor=5`.

One alternative might be to set values at or below zero to `NA` prior to
normalization, and before calling `matrix_normalize()`. In this case,
only non-NA values will be used during normalization according to the
`method` being used.

## See also

Other jamses utilities:
[`choose_annotation_colnames()`](https://jmw86069.github.io/jamses/reference/choose_annotation_colnames.md),
[`combine_sestats()`](https://jmw86069.github.io/jamses/reference/combine_sestats.md),
[`contrast2comp_dev()`](https://jmw86069.github.io/jamses/reference/contrast2comp_dev.md),
[`fold_to_log2fold()`](https://jmw86069.github.io/jamses/reference/fold_to_log2fold.md),
[`intercalate()`](https://jmw86069.github.io/jamses/reference/intercalate.md),
[`list2im_opt()`](https://jmw86069.github.io/jamses/reference/list2im_opt.md),
[`list_to_sestats()`](https://jmw86069.github.io/jamses/reference/list_to_sestats.md),
[`log2fold_to_fold()`](https://jmw86069.github.io/jamses/reference/log2fold_to_fold.md),
[`make_block_arrow_polygon()`](https://jmw86069.github.io/jamses/reference/make_block_arrow_polygon.md),
[`mark_stat_hits()`](https://jmw86069.github.io/jamses/reference/mark_stat_hits.md),
[`merge_statdf_all_test()`](https://jmw86069.github.io/jamses/reference/merge_statdf_all_test.md),
[`point_handedness()`](https://jmw86069.github.io/jamses/reference/point_handedness.md),
[`point_slope_intercept()`](https://jmw86069.github.io/jamses/reference/point_slope_intercept.md),
[`shortest_unique_abbreviation()`](https://jmw86069.github.io/jamses/reference/shortest_unique_abbreviation.md),
[`shrinkDataFrame()`](https://jmw86069.github.io/jamses/reference/shrinkDataFrame.md),
[`shrink_df()`](https://jmw86069.github.io/jamses/reference/shrink_df.md),
[`shrink_matrix()`](https://jmw86069.github.io/jamses/reference/shrink_matrix.md),
[`sort_samples()`](https://jmw86069.github.io/jamses/reference/sort_samples.md),
[`strsplitOrdered()`](https://jmw86069.github.io/jamses/reference/strsplitOrdered.md),
[`sub_split_vector()`](https://jmw86069.github.io/jamses/reference/sub_split_vector.md),
[`update_function_params()`](https://jmw86069.github.io/jamses/reference/update_function_params.md),
[`update_list_elements()`](https://jmw86069.github.io/jamses/reference/update_list_elements.md)

## Examples

``` r
# use farrisdata real world data if available
if (requireNamespace("farrisdata", quiet=TRUE)) {

   suppressPackageStartupMessages(suppressWarnings(
      library(SummarizedExperiment)))

   # test matrix_normalize()
   GeneSE <- farrisdata::farrisGeneSE;
   imatrix <- assays(GeneSE)$raw_counts;
   genes <- rownames(imatrix);
   samples <- colnames(imatrix);
   head(imatrix);

   # matrix_normalize()
   # normalize the numeric matrix directly
   imatrix_norm <- matrix_normalize(imatrix,
      genes=genes,
      samples=samples,
      method="jammanorm",
      params=list(minimum_mean=5))
   names(attributes(imatrix_norm))

   # review normalization factors
   round(digits=3, attr(imatrix_norm, "nf"));

   # example for quantile normalization
   imatrix_quant <- matrix_normalize(imatrix,
      genes=genes,
      samples=samples,
      method="quantile")
   names(attributes(imatrix_quant))
}
#> Loading required namespace: farrisdata


# simulate reasonably common expression matrix
set.seed(123);
x <- matrix(rnorm(9000)/4, ncol=9);
colnames(x) <- paste0("sample", LETTERS[1:9]);
rownames(x) <- paste0("gene", jamba::padInteger(seq_len(nrow(x))))
rowmeans <- rbeta(nrow(x), shape1=2, shape2=5)*14+2;
x <- x + rowmeans;
for (i in 1:9) {
   x[,i] <- x[,i] + rnorm(1);
}

# display MA-plot with jamma::jammaplot()
jamma::jammaplot(x)


# normalize by jammanorm
xnorm <- matrix_normalize(x, method="jammanorm")
jamma::jammaplot(xnorm, maintitle="method='jammanorm'")


# normalize by jammanorm with housekeeper genes
hk_genes <- sample(rownames(x), 10);
xnormhk <- matrix_normalize(x,
   method="jammanorm",
   params=list(jammanorm=list(controlGenes=hk_genes)))

jamma::jammaplot(xnormhk,
   maintitle="method='jammanorm' with housekeeper genes",
   highlightPoints=list(housekeepers=hk_genes),
   highlightColor="green");


xnormhk6 <- matrix_normalize(x,
   method="jammanorm",
   params=list(jammanorm=list(
      controlGenes=hk_genes,
      minimum_mean=6)))
hk_used <- attr(xnormhk6, "hk")[[1]];
jamma::jammaplot(xnormhk6,
   maintitle="method='jammanorm' with housekeeper genes, minimum_mean=6",
   highlightPoints=list(housekeepers=hk_genes,
      hk_used=hk_used),
   highlightColor=c("red", "green"));


# normalize by quantile
xquant <- matrix_normalize(x, method="quantile")
jamma::jammaplot(xquant,
   maintitle="method='quantile'")


# simulate higher noise for lower signal
rownoise <- rnorm(prod(dim(x))) * (3 / ((rowmeans*1.5) - 1.5));
xnoise <- x;
xnoise <- xnoise + rownoise;
jamma::jammaplot(xnoise,
   maintitle="simulated higher noise at lower signal");


# simulate non-linearity across signal
# sin(seq(from=pi*4/10, to=pi*7/10, length.out=100))-0.8
rowadjust <- (sin(pi * jamba::normScale(rowmeans, from=3.5/10, to=5.5/10)) -0.9) * 20;
xwarp <- xnoise;
xwarp[,3] <- xnoise[,3] + rowadjust;
jamma::jammaplot(xwarp,
   maintitle="signal-dependent noise and non-linear effects");


# quantile-normalization is indicated for this scenario
xwarpnorm <- matrix_normalize(xwarp,
   method="quantile");
jp <- jamma::jammaplot(xwarpnorm,
   maintitle="quantile-normalized: signal-dependent noise and non-linear effects");

```
