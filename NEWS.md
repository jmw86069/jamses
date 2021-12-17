# slicejam 0.0.3.900

# updates to existing functions

* `se_contrast_stats()` new arguments `block`, and `correlation` intended
to allow blocking factors in experiment design.
* `run_limma_replicate()` new arguments `block`, and `correlation` recognized
from `se_contrast_stats()`.


# slicejam 0.0.2.900

# updates to existing functions

* `ebayes2dfs()` was updated to handle data that does not have a column
with specific probe or gene annotations. In this situation `rownames()`
will be used from output of `limma::topTable()`.


# slicejam 0.0.1.900

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

