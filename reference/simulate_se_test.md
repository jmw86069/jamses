# Simulate SummarizedExperiment tests

Simulate SummarizedExperiment tests, with adjustable batch effect.

## Usage

``` r
simulate_se_test(
  ngroups = 4,
  nreps = 5,
  nrow = 250,
  multiplier = 1.5,
  offset = 7,
  batch_ex = 0.15,
  hit_exp = 2,
  hit_fraction = 1/2,
  hit_max = 2.8,
  noise_factor = 1,
  seed = 123,
  assay_name = "counts",
  sparsity = 0,
  verbose = FALSE,
  ...
)
```

## Arguments

- ngroups:

  `integer` number of experimental groups

- nreps:

  `integer` number of replicates per group, can be used to provide the
  number of replicates for each group in order.

- nrow:

  `integer` number of rows (measurements)

- multiplier:

  `numeric` value multiplied by
  [`rnorm()`](https://rdrr.io/r/stats/Normal.html) to adjust the
  magnitude of values produced, default 1.

- offset:

  `numeric` value added to the output of
  [`rnorm()`](https://rdrr.io/r/stats/Normal.html), default 7.

- batch_ex:

  `numeric` adjustment for batch effect, default `0.15` is a moderate
  batch effect, it would require adjustment or blocking or covariate for
  effective analysis.

- hit_exp:

  `numeric` default 2, adjust the "sharpness" of differential changes,
  higher values are sharper, with fewer strong hits.

- hit_fraction:

  `numeric` value between 0 and 1 indicating the fraction of rows to
  simulate as having a fold change, default 1/2.

- hit_max:

  `numeric` maximum value for a simulated fold change, default 2.8,
  intended to be interpreted as log2 of a 7-fold change.

- noise_factor:

  `numeric` multiplied by
  [`rnorm()`](https://rdrr.io/r/stats/Normal.html) to add additional
  noise.

- seed:

  `numeric` passed to [`set.seed()`](https://rdrr.io/r/base/Random.html)
  when provided, for reproducible random output.

- assay_name:

  `character` name to use for the assay name in the output.

- sparsity:

  `numeric` value from 0 to 1, default 0, indicating the fraction of
  values that are converted to `NA` to simulate sparse data
  measurements. It can be provided as a vector and applied to each group
  in order. In some proteomics datasets, control samples may be
  substantially more sparse than other groups, for example if the
  control is non-targeted IP, or negative control. In this case, data
  can be simulated by using `sparsity=c(0.6, 0, 0)` for three groups.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are ignored.

## Value

`SummarizedExperiment` object

## Details

This function extends
[`make_se_test()`](https://jmw86069.github.io/jamses/reference/make_se_test.md)
in two ways:

1.  For four-group scenario, it will generate two-factor changes.

2.  Batch effects are generated, adjustable with `batch_ex`. Use
    `batch_ex=0` for no batch effect. See examples for demonstration of
    centering within batch as visual cue that batch effect is
    well-defined.

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
[`se_to_assay_names()`](https://jmw86069.github.io/jamses/reference/se_to_assay_names.md),
[`se_to_rowcoldata()`](https://jmw86069.github.io/jamses/reference/se_to_rowcoldata.md)

## Examples

``` r
# batch effects
seb <- simulate_se_test();
#> Warning: first element used of 'each' argument
#> Warning: first element used of 'each' argument
#> Warning: data length [1249] is not a sub-multiple or multiple of the number of rows [250]
#> Warning: first element used of 'length.out' argument
#> Warning: first element used of 'length.out' argument
#> Warning: first element used of 'length.out' argument
#> Warning: first element used of 'length.out' argument
hmb <- heatmap_se(seb,
   controlSamples=colnames(seb)[c(1:5, 11:15)],
   column_title="With batch effect\nglobal-centered")
#> Error in loadNamespace(x): there is no package called ‘platjam’
hmb
#> Error: object 'hmb' not found

# when centering versus a sparse control group, some values can be lost:
hmbg <- heatmap_se(seb,
   controlSamples=colnames(seb)[c(1:5, 11:15)],
   centerby_colnames="batch",
   column_title="With batch effect\ncentered by batch")
#> Error in loadNamespace(x): there is no package called ‘platjam’
hmb + hmbg
#> Error: object 'hmb' not found

hmbc <- heatmap_se(seb,
   correlation=TRUE,
   controlSamples=colnames(seb)[c(1:5, 11:15)],
   column_title="With batch effect\nglobal-centered")
#> Error in loadNamespace(x): there is no package called ‘platjam’
hmbc
#> Error: object 'hmbc' not found

hmbgc <- heatmap_se(seb,
   correlation=TRUE,
   controlSamples=colnames(seb)[c(1:5, 11:15)],
   centerby_colnames="batch",
   column_title="With batch effect\ncentered by batch")
#> Error in loadNamespace(x): there is no package called ‘platjam’
hmbgc
#> Error: object 'hmbgc' not found
hmbc + hmbgc
#> Error: object 'hmbc' not found
```
