# Make SummarizedExperiment test data

Make SummarizedExperiment test data

## Usage

``` r
make_se_test(
  ngroups = 2,
  nreps = 3,
  nrow = 50,
  multiplier = 1,
  offset = 7,
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

- nrow:

  `integer` number of rows (measurements)

- multiplier:

  `numeric` value multiplied by
  [`rnorm()`](https://rdrr.io/r/stats/Normal.html) to adjust the
  magnitude of values produced, default 1.

- offset:

  `numeric` value added to the output of
  [`rnorm()`](https://rdrr.io/r/stats/Normal.html), default 7.

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

- mreps:

  `integer` number of replicates per group, can be used to provide the
  number of replicates for each group in order.

## Value

`SummarizedExperiment` object

## See also

Other jamses SE utilities:
[`geomx_to_se()`](https://jmw86069.github.io/jamses/reference/geomx_to_se.md),
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
se <- make_se_test();
heatmap_se(se)

heatmap_se(se, controlSamples=1:3)


se <- make_se_test(ngroups=3, assay_name="expression");
heatmap_se(se, apply_hm_column_title=TRUE)

heatmap_se(se, controlSamples=1:3, control_label="versus A", apply_hm_column_title=TRUE)


se <- make_se_test(ngroups=3, assay_name="expression", nreps=c(3, 5, 4));
heatmap_se(se, controlSamples=1:3, control_label="versus A", apply_hm_column_title=TRUE)


se <- make_se_test(ngroups=3, assay_name="expression", nreps=c(3, 5, 4), hit_max=4);
heatmap_se(se, controlSamples=1:3, control_label="versus A", apply_hm_column_title=TRUE)


sedesign <- groups_to_sedesign(SummarizedExperiment::colData(se)[, "group", drop=FALSE])
plot_sedesign(sedesign)

sestats <- se_contrast_stats(se=se, sedesign=sedesign, assay_names="expression")
plot_sedesign(sedesign, sestats=sestats,
   contrast_style="none", sestats_style="label")

heatmap_se(se, controlSamples=1:3, control_label="versus A",
   column_split="group",
   apply_hm_column_title=TRUE, sestats=sestats, color_max=5)


if (jamba::check_pkg_installed("venndir")) {
   venndir::venndir(hit_array_to_list(sestats), overlap_type="each")
   venndir::venndir(hit_array_to_list(sestats), overlap_type="each", proportional=TRUE)
   venndir::venndir(hit_array_to_list(sestats), overlap_type="each", show_labels="ncs")
}




# demonstrate sparsity
se2 <- make_se_test(sparsity=c(0.5, 0));
hm2a <- heatmap_se(se2,
   column_title="global centered\nall values shown")
# when centering versus a sparse control group, some values can be lost:
hm2b <- heatmap_se(se2, controlSamples=1:3,
   column_title="centered vs A\nsome values become NA")
hm2a + hm2b
#> Warning: Heatmap/annotation names are duplicated: centered expression

# use naControlFloor
hm2c <- heatmap_se(se2,
   column_title="centered vs A\nnaControlFloor=7",
   naControlFloor=7, naControlAction="floor", controlSamples=1:3)
hm2a + hm2b + hm2c
#> Warning: Heatmap/annotation names are duplicated: centered expression
#> Warning: Heatmap/annotation names are duplicated: centered expression, centered
#> expression

```
