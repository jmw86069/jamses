# jamses
Jam SummarizedExperiment Stats (jamses)

This package is under active development, as these functions
make up a core set of methods used across multiple
Omics analysis projects. A summary of goals and relevant
features are described below.

## Goal

The core goal is to make data analysis of
`SummarizedExperiment` objects straightforward for
common scenarios:

* Define a design matrix based upon sample groups

   * by default using `~ 0 + group` syntax.

* Define statistical contrasts for factor comparisons

   * one or more design factors are recognized
   * goal to enable one-way and two-way contrasts*
   * factor comparisons are applied in proper order
   
* Automate analysis of a series of contrasts

   * allow multiple data matrices stored as `assays`
   * intended to help compare data normalization/processing steps


> \*Note: [The Limma User's Guide](https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/)
describes alternate approaches for one-way and two-way contrasts.
The approach used in `jamses` uses the `~0 + x` style of grouping,
which defines each experiment group with independent replicates.
In this context, a two-way contrast is defined as testing the
fold change of one-way fold changes. The `jamses` methods ensure
that two-way contrasts compare the same factors in proper order,
for example `(A_treated - A_control)` can be compared to
`(B_treated - B_control)`, but will not be compared to
`(B_treated - B_knockout)`. Similarly, one-way contrasts will
only compare one factor change at a time, and would not
generate a contrast `(A_treated - B_control)`.


## Define experiment design and contrasts

* `groups_to_sedesign()` takes by default a `data.frame` where each column
represents an experiment factor, and performs the following steps:

   * defines appropriate experiment `design` matrix
   * defines a `contrast` matrix limited to changes in single
   experiment factors at a time, combining equivalent contrasts across
   secondary factors where appropriate, up to `max_depth` factor depth.
   * defines an appropriate `samples` vector typically defined by `rownames()`
   of the input `data.frame`.

* Output is `SEDesign` as described:

   * `SE@samples` - contains the vector of samples.
   * `SE@design` contains the design matrix.
   * `SE@contrasts` contains the contrasts matrix.

* The input `data.frame` columns represent design factors,
and column values represent factor levels.
If columns are `character` they are coerced to `factor`,
otherwise the order of existing `factor` levels is maintained.
* Contrasts are defined such that the first level is the control group
in each comparison.


## New object class: `SEDesign`

`SEDesign` is an S4 object that contains the following slots:

* `"samples"`: a `character` vector of sample names, derived from
`colnames()` of the `SummarizedExperiment` object.
* `"design"`: a `matrix` representing the design matrix, with specific
constraints:

   * `rownames()` are equal to `colnames()` of the `SummarizedExperiment`
   object. This requirement helps ensure that the design and assay data
   are integrated.
   * `colnames()` are defined as sample groups, typically equivalent to
   using `model.matrix( ~ 0 + groups )`, based upon the 
   [The Limma User's Guide](https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/)
   from the `limma` Bioconductor package.

* `"contrasts"`: a `matrix` representing the contrast matrix used for
statistical comparisons, with constraints:

   * `rownames(contrasts)` must be equal to `colnames(contrasts)` in order
   to ensure all data is properly synchronized.

### Use cases for `SEDesign`:

* `SEDesign` can be subset by `samples` using the form: `SE[samples, ]`

   * Subsetting `samples` may be useful in order to remove outlier samples,
   or to focus analysis on a subset of samples within best statistical
   practices (outside scope of this package).
   * Subsetting `SEDesign` by `samples` will force appropriate subsetting
   of the `design` matrix. This subset will cascade to `contrasts`.
   * Argument `min_reps=1` defines the minimum samples required in a
   `design` group in order for the group to be retained. It may be beneficial
   to require at least `n=3` replicates for analysis for example.

* `SEDesign` can be subset by `groups` using the form: `SE[, groups]`

   * The `design` matrix will be subset accordingly.
   * The `samples` vector will be subset to remove samples that are not
   present in any remaining design groups.
   * The `contrasts` matrix will be subset to remove all contrasts where one
   or more groups are no longer present.

* `SEDesign` can be subset by `contrasts`.

   * The `contrasts` matrix will be subset accordingly.
   * The `design` matrix will be subset to include only those groups
   required by the contrasts.

### Accessor functions to `SEDesign`:

* `samples()` will list the samples, equivalent to `rownames(design)`.
* `samples()<-` will set values in slot `samples`, and update `rownames(design)`.
* `design()` will return the `design` matrix which defines experiment groups.
* `design()<-` will set the `design` matrix, also ensure that:

   * `colnames(design)` are in the same order as `rownames(contrasts)`, and
   * `rownames(design)` are in the same order as `samples`.
   
* `contrasts()` will return the `contrasts` matrix.
* `contrasts()<-` will set the contrast matrix, and verify the
`rownames(contrasts)` are in the same order as `colnames(design)`.
* `groups()` will export `colnames(design)` representing the design groups.
* `groups()<-` will update group names with the corresponding vector.


## Data normalization

`SummarizedExperiment` objects are normalized by:

* `se_normalize()` - lightweight wrapper to one of several normalization
functions:

   * `jamma::jammanorm()` which normalizes the median log fold change to zero.
   This method is also used in `DESeq2::estimateSizeFactors()` for example.
   When using `jamma::jammaplot()` for MA-plots, the `jammanorm()` method
   is conceptually equivalent to shifting the y-axis values to zero, so
   the median log fold change is zero. (Assumptions apply.)
   * `limma::normalizeQuantiles()` for quantile normalization.
   * `limma::removeBatchEffect()` for adjustment of batch effects,
   recommended for visualization and clustering, not typically recommended
   prior to statistical comparisons. Instead, we suggest using a blocking
   factor, which can be passed to `se_contrast_stats()`.

* `matrix_normalize()` - is the core function for `se_normalize()` and operates
on individual numeric data matrices.


## Heatmaps

`heatmap_se()` is a convenient wrapper for `ComplexHeatmap::Heatmap()`.
Some useful arguments and features are described below:

* `se`: the `SummarizedExperiment` data for the heatmap
* `assay_name`: define a specific data matrix to display
* `rows`: choose a specific subset of rows
* `sestats`: output from `se_contrast_stats()` to display a hit matrix
alongside the heatmap. Also `rows` are automatically subset to
show statistical hits.
* `top_colnames`: display column annotations across the top. When
absent, it will auto-detect columns which may have interesting
annotations, based upon the cardinality of each column annotation.
* `rowData_colnames`: display row/gene annotations on the left
* `sample_color_list`: supply pre-defined set of colors for top and
left annotations
* `centerby_colnames`: optional data centering sub-groups, where data
is centered within each centerby group.
* `normgroup_colnames`: optional normalization groups, also used
for independent data centering. Columns are also "split" by this value,
visually separating columns by each normgroup.

Other options can be passed to `ComplexHeatmap::Heatmap()`:

* `row_split`: split by `rowData(se)` column annotations, or by number
of dendrogram sub-tree clusters.


## Statistical comparisons

### `se_contrast_stats()` is the central function

* Applies design and contrasts defined by `SEDesign`, or given
specific `design` and `contrasts` matrix objects.
* Initially calls `limma` functions on a series of `assays`,
for one or more data matrices in the `SummarizedExperiment` object.

   * Optionally applies `limma-voom` methodology for count data, for
   example RNA-seq or Nanostring, to apply the `voom` methodology of
   estimating dispersion and applying weights to contrast stats.
   * Optionally applies `limma-DEqMS` methodology for proteomics mass
   spec data, to apply error model based upon PSM counts per row.

* Applies statistical thresholds to define statistical hits for follow-up:

   * `adjp_cutoff`: adjusted P-value threshold
   * `fold_cutoff`: normal space fold change threshold
   * `mgm_cutoff`: "mgm" is an abbreviation for "max group mean";
   which requires at least one experiment group involved in the
   contrast to have group mean at or above this threshold. The cutoff
   is useful to require at least one group to have signal above
   a noise threshold.

* Additional threshold options:

   * `p_cutoff`: unadjusted P-value threshold
   * `ave_cutoff`: average expression threshold, using the equivalent of
   `limma` column `AveExpr` with the mean group expression. Note this
   `AveExpr` value is calculated mean per row across all groups,
   and is subject to skewing.
   * `int_adjp_cutoff`, `int_fold_cutoff`, `int_mgm_cutoff` are
   optionally used for interaction contrasts, for example sometimes a two-way
   contrast threshold is more lenient during data mining experiments,
   eg. `int_adjp_cutoff=0.1` for example.

* Options:

   * `use_voom=TRUE`: enable limmavoom workflow for count data
   * `posthoc_test="DEqMS"`: enable DEqMS post-hoc adjustment to the limma
   empirical Bayes model fit.
   * `floor_min`, `floor_value`: logic to handle numeric values below a
   pre-defined noise threshold, values below this threshold are assigned
   `floor_value`. This filtering may allow converting values from `0` zero
   to `NA`, and as such the values do not contribute to replicate
   measurements. For data with substantial missing measurements, this
   option may be beneficial.
   * `block`: optional blocking factor, which enables the
   `limma::duplicateCorrelation()` calculation, which also applied and
   used during model fit per the
   [The Limma User's Guide](https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/).

* Returns a `list` referred to as `sestats`:

   * `"hit_array"`: `array` of named numeric vectors indicating direction of
   statistical hits after applying the relevant threshold cutoffs:
   
      * `1` up-regulated
      * `-1` down-regulated
      
   * `"stats_dfs"`: `list` of `data.frame` for each contrast,
   where column headers also include the statistical contrast to
   help confirm the identity of each analysis result.
   * `"stats_df"`: one super-wide `data.frame` with results across
   all contrasts, assay data matrices, and hit thresholds. This
   matrix is intended to help filter for hits across the various
   different contrasts and thresholds.

* Save to excel with `save_sestats()`

   * This function helps automate saving every statistical contrast
   to its own Excel worksheet.
   * Each worksheet is styled consistently, making it easy to flip
   through each worksheet.


* Future work:

   * Port the function `plot_comparisons()` which creates a visual
   plot of each one-way and two-way contrast, with axes defined by the
   experiment factors.
   * Enable equivalent analysis workflow steps using `DESeq2` methodology.

