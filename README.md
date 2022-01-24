# jamses
Jam SummarizedExperiment Stats (jamses)

Note this package is in active development, functions have
been in use for quite some time in offline context, and are
being ported here for convenient re-use.

The core goal is to help make the process of analyzing data
stored in `SummarizedExperiment` objects straightforward for
a few specific scenarios:

* Defining sample groups to use in statistical comparisons
* Defining statistical contrasts with one or more design factors,
specifically with the goal to enable one-way and two-way contrasts.
* Automate the process of running a series of pairwise comparisons
for several contrasts, sometimes using multiple data matrices stored
in the `SummarizedExperiment` slot `assays`. This step is intended
to analyze data normalized/processed differently so the results can
be easily and rapidly compared.

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
   using `model.matrix( ~ 0 + groups )`, based upon the Limma User Guide
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

## Data normalization

`SummarizedExperiment` objects are normalized by:

* `se_normalize()` - lightweight wrapper to one of several normalization
functions:

   * `jamma::jammanorm()` which normalizes the median log fold change to zero.
   * `limma::normalizeQuantiles()` for quantile normalization.
   * `limma::removeBatchEffect()` for adjustment of batch effects, recommended
   mainly for visualization and clustering, and not ideally prior to
   statistical comparisons.

* `matrix_normalize()` - is the core function for `se_normalize()` and operates
on individual numeric data matrices.
