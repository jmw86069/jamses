# jamses 0.0.8.900

## changes to existing functions

* `groups_to_sedesign()` new behavior

   * two-way interaction contrasts now retain the order of first-order
   contrasts, preferring to alter the reverse `rev(factor_order)`
   in order to build interaction contrasts on top of existing contrasts.
   * some verbose output now requires `(verbose >= 2)`
   * added verbose output for two-way contrasts removed because they
   are equivalent, but written with different order of comparisons.


# jamses 0.0.7.900

## changes to existing functions

* `groups_to_sedesign()` new argument
`default_sort=c("asis", "sort_samples", "mixedSort")`:

   * Note in all cases, a `factor` column is sorted by its factor levels.
   * `"asis"`: will convert any non-`factor` column to `factor`, with levels
   defined in the order they appear in the input data. Any existing `factor`
   column will remain as-is.
   * `"sort_samples"`: will call `sort_samples()` in order to recognize
   control terms, and preferentially place them first, for example
   "Control", or "Wildtype". All other terms are sorted by `jamba::mixedSort()`
   for alphanumeric sorting.
   * `"mixedSort"`: will use `jamba::mixedSortDF()` to sort group columns:
   
      * `factor` columns are sorted according to factor levels
      * `character` columns are alphanumeric sorted using `jamba::mixedSort()`
   
   * To influence the sort order, `data.frame` columns should have
   `factor` values with defined levels in order.
   Otherwise by default, the first entries in the `data.frame` provided
   will become the first factor levels.

* `validate_sedesign()` default argument `verbose=FALSE`, also silencing
other internal messages which were previously always turned on.
* `matrix_normalize()` new argument `normgroup` normalizes subsets of
matrix columns independently, then re-assembles the original matrix.

   * Purpose is convenience of maintaining one data matrix, while also
   applying normalization independently on each sub-matrix. The
   attributes returned from the normalization method are also combined,
   for example the normalization factors `"nf"` reported by
   `jamma::jammanorm()` are combined in proper order.
   * This feature is distinct from batch adjustment, in that the
   data in each sub-matrix are never normalized relative to each
   other.
   
      * For example, one could normalize total RNA-seq and nascent 4sU-seq
      data independently, without expectation that the two would ever
      have a common frame of reference to normalize one relative to another.
      * Similarly, one could normalize each tissue type independently,
      which may be appropriate when analyzing data that contains very
      different mammalian tissue organ samples, such as muscle and brain.
      It would generally not be appropriate to use quantile normalization
      across muscle and brain samples, since the overall pattern and
      distribution of expression values is not expected to be similar.
      Quantile normalize assumes (and imposes) a common distribution,
      by adjusting mean expression signal at each quantile to a common
      mean expression across all samples.
      * For a rough approximation of cross-tissue normalization, one
      could apply `"quantile"` normalization within each `normgroup` defined
      by tissue type, then apply `"jammanorm"` median normalization to
      apply a linear adjustment of signal across tissue types. The median
      normalization does not affect distribution, thus will not affect
      intra-tissue contrasts, except by adjusting its overall signal
      which may change downstream assumptions regarding signal thresholds.
      It is still not advise to compare directly across tissue types,
      except in some cases a two-way contrast may be appropriate
      (comparing intra-tissue fold change to intra-tissue fold change.)

* `se_normalize()` also gained new argument `normgroup`, although
this argument is merely passed through to `matrix_normalize()`.

## new functions

* `sestats_to_df()` - simple method to convert `sestats` output from
`se_contrast_stats()` to a `data.frame` summary with number of hits
for each comparison. It should be suitable for `kable()` output,
with row groups by cutoff.


# jamses 0.0.6.900

## changes to existing functions

* `se_contrast_stats()`

   * new argument `sedesign` allows providing `SEDesign` output from
   `groups_to_design()`. When `isamples` is also provided, it will be
   used to subset the `sedesign` as appropriate, using
   `validate_design()`.

* `validate_sedesign()` was updated to prevent the annoying warning
when comparing vectors of two different sizes.

## new functions

* `save_sestats()` is a wrapper function intended to help export stat
summary tables into Excel using `jamba::writeOpenxlsx()`. This function
defines useful column style functions for statistical data types, such
as P-values, fold changes, log fold changes, etc. It also saves
all statistical contrasts, with one contrast in each worksheet.
The worksheet names are limited to 31 characters, the limitation
of MS Excel.

# jamses 0.0.5.900

Bumped version, to sync the transition of slicejam functions into jamses.

# jamses 0.0.4.900

## new functions

* `se_collapse_by_row()` is a new function intended to collapse
rows, for example peptides into protein abundance.


# jamses 0.0.3.900

# updates to existing functions

* `se_contrast_stats()` new arguments `block`, and `correlation` intended
to allow blocking factors in experiment design.
* `run_limma_replicate()` new arguments `block`, and `correlation` recognized
from `se_contrast_stats()`.


# jamses 0.0.2.900

# updates to existing functions

* `ebayes2dfs()` was updated to handle data that does not have a column
with specific probe or gene annotations. In this situation `rownames()`
will be used from output of `limma::topTable()`.


# jamses 0.0.1.900

## initial package

Methods are being ported from another internal R package.

## New functions

* `matrix_normalize()` - wrapper function to normalize numeric matrix data
* `se_normalize()` - function to normalize SummarizedExperiment assay data

* `groups_to_sedesign()` - function to convert experiment group
into `SEDesign` with proper design matrix, and contrast matrix.
This function was ported from `splicejam::groups2contrasts`.

* `se_contrast_stats()` - wrapper function to perform numerous statistical
contrasts, in this case calling `limma` methods with or without limma-voom.
* `voom_jam()` - lightly extended variation of `limma::voom()` solely
to enable presence of `NA` values without distortion of the model fit.

* `handle_na_values()` - function to manage presence of NA values using
one of several appropriate strategies.

* `run_limma_replicate()` - wrapper function to perform `limma` moderated
t-test steps, and to produce consistent, fully-annotated summary `data.frame`.
* `ebayes2dfs()` - internal wrapper function to convert `limma::eBayes()`
output to a list of `data.frame` tables.

* `mark_stat_hits()` - wrapper function used to apply one or more
statistical thresholds to mark statistical hits among rows of
statistical contrast results.

* `update_function_params()` - internal function used to apply user-specific
parameters to a pre-defined set of parameters stored in function formal
arguments.
* `update_list_elements()` - internal function called by
`update_function_params()` as a generic way to update individual elements
in a `list` with user-defined values.
* `log2fold_to_fold()` and `fold_to_log2fold()`


## New S4 object class

`SEDesign`: SummarizedExperiment design object

* `design`: numeric matrix indicating the experiment design
* `contrasts`: numeric matrix indicating statistical contrasts
* `samples`: `character` vector of sample replicates, also equivalent
to `rownames(design)`

Some constraints:

* `rownames(design)` contains sample replicate identifiers,
intended to correspond to `colnames(SummarizedExperiment)`
* `colnames(design)` contains the name of each sample group
* `rownames(contrasts)` contains the name of each sample group
* Enforced constraint: `colnames(design) == rownames(contrasts)`

Related accessor and assignment methods:

* `[` subset operations use `[samples, groups]` to enable re-ordering,
or subsetting data by sample and/or groups. For example, one may want
to remove a technical outlier sample from an analysis, or include
only a subset of groups in downstream statistical contrasts.
* `validate_sedesign()` contains the full logic for constraints
and subset operations:

   * the order of `samples` must match `rownames(design)`
   * the order of `colnames(design)` must match `rownames(contrasts)`
   * when a subset of `samples` is defined, `design` is updated
   to remove groups which are no longer represented by `min_reps`
   replicates. By default, a group with zero replicates is no
   longer a valid group
   * When a subset of groups is defined, perhaps by using a subset of
   samples as described above, the `contrasts` entries are evaluated
   to see which entries require the groups to be removed. For each contrast
   that requires a group being removed, the contrast is also removed.
   This step prevents retaining an incomplete contrast.

* `samples()` will return the `object@samples`
* `samples()<-` to assign new values to `object@samples`, which will
also update `rownames(design)` accordingly. Note these values are
expected to match `colnames(SummarizedExperiment)`.
* `groups()` will return `colnames(object@design)`, representing all
sample groups defined.
* `groups()<-` will assign values to `colnames(object@design)`,
which will also update `rownames(object@contrasts)` if there are
statistical contrasts already defined in the `object`.
* `contrastnames()` will return `colnames(object@contrasts)` representing
the name of each contrast.
* `contrastnames()<-` will assign values to `colnames(object@contrasts)`,
which may be useful to define a sensible label.
* `design()` will return the design matrix `object@design`.
* `design()<-` will assign a new design matrix to `object@design`,
then will call `validate_sedesign()` to enforce the relevant
constraints. For example `colnames(object@design)` must match
`rownames(object@contrasts)` if contrasts are defined in the object.
* `contrasts()` will return the contrast matrix `object@contrasts`.
* `contrasts()<-` will assign a new contrast matrix to `object@contrasts`,
then will call `validate_sedesign()` to enforce the relevant
constraints. For example `colnames(object@design)` must match
`rownames(object@contrasts)`.

