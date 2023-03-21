# jamses 0.0.36.900

## bug fixes

* `heatmap_se()`

   * `row_split` passed as integer value when `cluster_rows=FALSE`
   caused an error:
   "Error: Length or nrow of `row_split` should be same as nrow of `matrix`."
   The error was caused because rows cannot be split into subgroups
   in the absence of a row clustering.
   * The same error for `column_split` was also present, and was fixed.
   * Resolution is that `row_split` or `column_split` is converted to `NULL`
   when there is no respective row or column dendrogram or clustering
   function.

## changes to existing functions

* `save_sestats()`

   * Default format for mgm and mean is "num" instead of "int",
   which assumes log2-scaled values that are generally at or less than 20.

* `heatmap_se()`

   * the default gap between row annotations is defined by the
   ComplexHeatmap option `"ROW_ANNO_PADDING"`, so the gap between heatmap
   and annotation is the same as the gap between two annotations.
   In future the gap width might become a separate option but not just yet.


# jamses 0.0.35.900

## changes to existing functions

Dependency on `jamma` bumped to version `0.0.28.900` for new arguments
to `naControlAction`, see below.

* `heatmap_se()`

   * New arguments for data centering: `controlFloor`,
   `naControlAction`, `naControlFloor`
   * `controlFloor` is a `numeric` value for optional noise floor used
   for control groups during centering. If the control group mean/median
   value is below `controlFloor`, it is set to `controlFloor` to define
   a minimum baseline for centering. It is mostly useful when the effective
   noise floor for a platform is above zero.
   * `naControlAction` and `naControlFloor` handle the specific case where all
   `controlSamples` have `NA` values for a particular row.
   By default, all centered values would become `NA` thereby making
   platform measurements disappear from the heatmap.
   Sometimes it is preferable to define another reference
   so the non-NA values can be indicated in the heatmap.
   * See `jamma::centerGeneData()` for more information.
   * The most common alternatives:
   
      * `naControlAction="row"`: uses the remaining non-NA values from that row
      * `naControlAction="floor"`: assigns the `numeric` value
      `naControlFloor`, which is useful when the reference value
      should be a known baseline limit of detection.

* `shrink_df()`

   * Argument `numShrinkFunc` now includes `na.rm=TRUE` by default.

# jamses 0.0.34.900

## changes to existing functions

* `heatmap_se()`

   * new argument `use_raster=TRUE` which allows turning it off.
   When R package `magick` is not installed, on some systems it causes
   an error when trying to render ComplexHeatmap with `use_raster=TRUE`.
   Using `use_raster=FALSE` avoids the error.

* `groups_to_sedesign()`

   * small update to enforce `mixedSortDF(..., honorFactor=TRUE)`
   to ensure that `data.frame` columns that contain `factor` values
   will use the `levels` to define the proper sort order.

# jamses 0.0.33.900

## package dependencies

* `ComplexHeatmap` now requires version `2.13.4` or higher, and points
to the Github repository instead of Bioconductor. This change is required
for R-3.6.1 whose matching Bioconductor version 3.10 only offers
ComplexHeatmap version 2.2.0, and throws an error with `heatmap_se()`
when rendering row annotations.
* `jamba` dependency now requires version `0.0.88.900`, to pick up the subtle
but important change to `mixedSort()`
* `colorjam` was bumped to version `0.0.23.900` for consistency.

## changes to existing functions

* `groups_to_sedesign()`

   * 

# jamses 0.0.32.900

## new functions

* `se_collapse_by_column()` which is the equivalent of `se_collapse_by_row()`
except that the default math is much simpler, to perform mean calculation
with optional outlier removal prior to calculating the mean.
* new experimental function `shrink_df()` as potential drop-in replacement
for `shrinkDataFrame()`. This new function is intended to have simpler
processing, with less complicated history than `shrinkDataFrame()`.

## changes to existing functions

* minor updates to help text for `se_collapse_by_row()`, added `@family`
to this and other functions.
* `se_contrast_stats()` was changed from using `matrixStats::rowMaxs()`
instead to use `apply(x, 1, max, na.rm=TRUE)` - due to persistent
Segfaults caused by `matrixStats::rowMaxs()` on multiple architectures
(OSX, linux) and with recent package versions, at least for R-3.6.1.

# jamses 0.0.31.900

## new functions

* `names_contrast2comp()` and `names_comp2contrast()` are analogous to
`contrast2comp()` and `comp2contrast()` except they operate on the
object names, making it convenient to wrapper around `list` objects
where the list names are the contrasts or comps.

## changes to existing functions

* `se_contrast_stats()`

   * argument `normgroup` now accepts `colnames(colData(se))` to define
   normgroups by column annotations.

## bug fixes

* fixed apparent regression in `groups_to_sedesign()` for interaction contrasts

   * version 0.0.30.900 apparently introduced regression bug in interaction
   contrasts that returns only the interactions in `list` form instead
   of all contrasts in `sedesign` form. Short-lived, fortunately.
   * Some debug potential was added with `verbose=TRUE`, that prints
   iterative contrast name generation, with indent to help show when the
   function is iteratively calling itself.

* `heatmap_se()`

   * fixed bug with user-supplied `row_names_gp` or `column_names_gp`
   that does not also define `fontsize`.
   * new arguments `cluster_column_slices=FALSE`, `cluster_row_slices=FALSE`
   to allow this behavior to be overridden when relevant

## changes to existing functions

* `heatmap_se()`

   * moved substantial logic out, and into `process_sestats_to_hitim()`,
   which handled the hit incidence matrix, and deciding which rows
   are included in the heatmap with default `rows=NULL`.
   * quality of life change: when there is only one `assayNames(se)`, then
   `assay_name` defaults to the only possible value, instead of acting
   like it couldn't possibly know what I mean. Haha.
   * `legend_at` is smarter when `centerby_colnames=FALSE` and there is
   no data centering, the legend only includes negative values when
   data contains negative values.
   * `legend_labels` are smarter when `color_max <= 1`, it uses the numerical
   value without transformation.
   * vastly expanded the examples and the param help text.

## new functions

* `process_sestats_to_hitim()`

   * previously internal logic to `heatmap_se()`, it surpassed the amount
   of logic hidden inside a function, especially when it's used twice.

# jamses 0.0.30.900

## changes to existing functions

* `se_contrast_stats()`

   * removed duplicated `merge()` called across stats `data.frame` results.
   * implemented `normgroups` for independent model fit among subsets
   of samples.


# jamses 0.0.29.900

## bug fixes

* `heatmap_se()` did not properly check `cluster_columns` for `function`
input, so it did not allow using a custom clustering function.

## changes to existing functions

* `heatmap_se()`

   * new arguments: `row_names_gp` and `column_names_gp` allow custom
   settings for row and column labels, respectively.


# jamses 0.0.28.900

## bug fixes

* `se_contrast_stats()` with `posthoc_test="DEqMS"` and only one contrast
threw an error
`"Error in fit$sca.t[results.table$gene, coef_col] : subscript out of bounds"`.
The error actually occurs in `ebayes2dfs()` the function that converts
statistical model fit into a list of `data.frame` objects.

   * The `DEqMS` package produces `numeric` matrix results that it also
   summarizes into a `data.frame` via `DEqMS::outputResult()`. When there
   is only one contrast, many of the `numeric` matrix objects also have
   one column, and R apparently dropped the `colnames()` for the `sca.t`
   matrix, causing `DEqMS::outputResult()` to fail when argument `coef_col`
   is given the `character` name of the contrast instead of the `integer`
   column number.
   * The contrast name is preserved in the `coefficients` matrix,
   so the contrast is converted to column number, then passed as `coef_col`,
   and resolves the error.
   * Calls to `DEqMS::outputResult()` use `integer` values for argument
   `coef_col`.
   * Some internal step in `DEqMS` does not properly handle data with only
   one contrast, probably forgot to include `drop=FALSE` in matrix subsetting.

## changes to existing functions

* `groups_to_sedesign()` updates

   * input `SummarizedExperiment` is now recognized, using `group_colnames`
   to define the experimental design `data.frame`. Data can also be subset
   using `isamples` for example.
   * new argument `contrast_names` to provide specific contrasts upfront,
   therefore skipping all internal logic of defining contrasts.
   * internal values `max_depth` and `factor_order` are validated using
   the design `data.frame` such that it will not pursue `max_depth` nor
   `factor_order` higher than the number of design columns available.
   * When there are no valid contrasts for the given design groups,
   the function returns empty contrasts instead of throwing an error.
   One could argue that an error is actually preferred, since it forces
   the input data to be valid, and forces the user to review the data.

* `contrast2comp()` and `comp2contrast()` were reworked to handle alternate
delimiters, which mainly involved wrapping the delimiter in proper
regular expression patterns. During the process, obscure testing showed
certain delimiters cause problems with regular expressions, so those
problems were pre-emptively fixed. The help docs were expanded somewhat.
* `sestats_to_df()` was refactored to be cleaner, correcting the issue
of colname mislabling when `dimname_order` was different than expected.
The function now works smoothly with any `dimname_order` values, and accepts
`character` or `integer` input to specify the dimname ordering.

## new functions

* `format_hits()` is an internal function, but it might be a useful
convenience function so it is exported. It is called by `sestats_to_df()`
to provide a summary `data.frame` of statistical results.
It converts a hit vector into either:

   * a summary string with hits/up/down,
   * one integer count of hits, or
   * a vector of integers for c("hits", "up", "down")`.


# jamses 0.0.27.900

## bug fixes

* `groups_to_sedesign()` threw an error in some circumstances when the
input `ifactors` only contained one column, but the column contained
levels with delimiters recognized by `factor_sep`, and some levels without
the delimiter. The error was: `"Error in combn(iMatch, 2) : n < m"`

   * The workaround was to supply levels with consistent
   delimiters, or a design that used more than one column initially.
   * However, the R code was adjusted to handle these situations properly
   without workarounds.

# jamses 0.0.26.900

## updates to existing functions

* `heatmap_se()` updates:

   * still in progress: several new arguments intended to customize exact
   font size and grid size for color legends, annotations, around the heatmap.
   * new arguments: `show_heatmap_legend`, `show_top_legend`, `show_left_legend`
   intended to allow hiding various color legends, mostly to save plot space.
   * `all_sample_colors` sometimes does not match left or top annotation
   names, and now will create categorical colors in more scenarios.
   * new argument `mark_rows` to enable optional `ComplexHeatmap::anno_mark()`
   row labels for a subset of heatmap rows. `mark_labels_gp` to customize the
   displayed font.
   * new argument `right_annotation` for custom annotation
   * new arguments `show_top_annotation_name`, `show_left_annotation_name`,
   `left_annotation_name_rot` to customize these label positions.
   * new arguments `mark_rows` and `mark_labels_gp` used to define
   `anno_mark()` call-out labels for a subset of rows in the heatmap. This
   step is difficult to define upfront when the data in `se` is subset
   by some other aspect of `heatmap_se()`, since `anno_mark()` requires
   numeric index positions for each label, and the order of rows may not
   be known upfront.


## minor bug fixes

* `hit_array_to_list()` was returning `NA` values when the `sestats$hit_array`
also included `NA` values. The function now removes any `NA` values before
returning each vector.
* `heatmap_se()` argument `alt_sestats` was throwing an error when supplied
with a `numeric` matrix instead of `sestats` object, this input type
was been corrected, and the code slightly refactored for clarity.
* `heatmap_se()` was passing multiple values to the left_annotation argument
`show_legend`, which only accepts one `logical` value. This bug was corrected.


## new dev function

* `contrasts2comp_dev()` is an experiment in converting the code to be
vectorized, though the exceptions to handle unbalanced contrasts may
still prevent this function from being as fast as intended.


# jamses 0.0.25.900

## updates to existing functions

* `heatmap_se()` updates:

   * new argument `data_type` used in color legend title, and heatmap title.
   * new arguments `legend_at`, `legend_labels` to customize heatmap color
   gradient legend positions and labels.
   * data centering can be turned off with `centerby_colnames=FALSE`.
   * `sestats` and `alt_sestats` can be supplied as incidence matrices directly,
   without requiring them to be in `sestats` or `hit_array` formats.
   * `top_colnames` can be hidden with `top_colnames=FALSE`; also when
   `colData(se)` is empty, `top_annotation` is hidden.
   * `assay_name` may contain multiple values to define gene hits in `sestats`;
   however only the first matching `assay_name` in `names(assays(se))` is
   used for the heatmap data.
   


# jamses 0.0.24.900

## updates to existing functions

* `heatmap_se()` new argument `row_subcluster`:

   * intended to help drill-down into a subset of heatmap rows by
   using elements of `row_split` as returned by
   `jamba::heatmap_row_order()`.
   * The rows are split by `row_split` using the full heatmap data,
   and can be split by...
   
      * cluster number using `row_split=5` integer values,
      * annotation columns in `rowData(se)`,
      * `character` vector, or
      * `data.frame`
   
   * In any case, `heatmap_row_order()` returns a `list` of `rownames(se)`,
   which is used to create the heatmap. In fact, rows are displayed in
   the same order, but without any associated dendrogram even if there
   is a dendrogram associated with the subset of rows.
   * `row_subcluster` can be an `integer` vector referring to the row
   list element in order, or `character` vector of one or more `names()`
   associated with the row list elements.

* `heatmap_se()` arguments `cutoff_name` and `alt_cutoff_name` default
values changed to `NULL`, so hits in `sestats` and `alt_sestats` will
use all cutoffs by default, unless specified otherwise.

# jamses 0.0.23.900

## updates to existing functions

* `heatmap_se()` was updated

   * `isamples` is applied to subset the columns displayed in the heatmap,
   however it no longer subsets data during the data centering step.
   This change allows centering by patient with time 0 as the control point,
   which produces a stripe of `0` values for each patient at time zero,
   all other times become difference from time zero. The heatmap can then
   use `isamples` with only samples not at time zero.
   * new argument `subset_legend_colors=TRUE` will now filter colors in
   the `top_annotation` and `left_annotation` color legends to remove
   categorical colors which are not represented in the data being displayed.

      * This update prevents showing a huge list of categorical colors that
      are not relevant to the heatmap. Note that even when two heatmaps
      are added together (with ComplexHeatmap `+` or `%v%`) the amazing
      coding in ComplexHeatmap will merge together color legends by
      annotation color name, showing only the unique set of colors.
      * It seems that ComplexHeatmap retains all colors in the heatmap object,
      but hides colors that are unused. The `heatmap_se()` code was actually
      encoding the `at` color labels for all colors regardless of which
      colors were present.
      So it is unclear the benefit of this option, but the function
      itself was improved.
      
   * color legends are displayed either in order of factor levels,
   or sorted using `jamba::mixedSort()`
   * new argument `row_anno_fontsize` to control the specific fontsize
   for row annotation labels. These labels appear for row annotations,
   and appears beside the column labels on the heatmap.

## bug fixes

* `heatmap_se()` when `sestats` is supplied but `se` does not
contain all rows present in `sestats` hit array, it produced an error.
* `heatmap_se()` when annotations with factor columns were subset to remove
a color from the legend, the color legend was still encoding `at` labels
for the full set of colors, causing `ComplexHeatmap::HeatmapAnnotation()`
to print a warning message.

## new functions

* `choose_annotation_colnames()` is a wrapper function to hold the simple
logic for choosing annotation columns to include for `heatmap_se()`.
By default it chooses columns with at least one repeated value, at
least two unique values, and only includes the first column with
matching cardinality. The last criterion ensures it won't display
six columns which show the same pattern of information for all samples.

## removed functions

* `call_fn_ellipsis()` was removed since it was migrated to `jamba`.
(Previous change.)


# jamses 0.0.22.900

## bug fixes

* `hit_array_to_list()` was returning `matrix` for single contrasts,
due to `apply()` tendency to simplify output to vector. The method
was changed to use a different approach guaranteed to return `list`.
* `heatmap_se()` had a typo in an edge condition when using `alt_sestats`,
which has been corrected.

# jamses 0.0.21.900

## bug fixes

* `heatmap_se()` was not properly auto-detecting `show_row_names`,
which only worked when `rows` was supplied specifically. Now it
works when `rows` is not supplied and/or only `sestats` is also supplied.

# jamses 0.0.20.900

## changes to existing functions

* `sestats_to_df()` was updated:

   * now reports `""` whenever the `hit_array` value is NA,
   because the NA in this context means
   the cutoff was not tested for statistical hits. This situation
   happens when using `int_adjp_cutoff`,`int_fold_cutoff`,`int_mgm_cutoff`
   that differ from the main contrast cutoffs. Therefore, reporting
   `0` implies the threshold was tested and no hits were found, which
   is incorrect.
   * also whenever there are zero hits it reports `"0 hits"` and not
   `"0 hits (0 up, 0 down)"`. This change saves unnecessary text, also
   makes it easier to see which entries have hits.

# jamses 0.0.19.900

## changes to existing functions

* `heatmap_se()` was updated:

   * The logic of obtaining a hit list for a contrast must account for
   multiple cutoffs, particularly when interaction contrast cutoffs
   differ from pairwise cutoffs. This logic was pushed into new function
   `hit_array_to_list()`.
   * The default color legend for `left_annotation` and `top_annotation`
   now includes `color_bar="discrete"` for numeric annotation columns.
   The legend will display discrete color breaks instead of a continuous
   color bar. The breaks are defined in the color function, in
   most cases the recommended method is calling `platjam::design2colors()`,
   which typically defines reasonable color breaks.

## new functions

* `hit_array_to_list()`

   * Input is `hit_array` from `sestats`, as output from `se_contrast_stats()`.
   * It optionally takes a subset of the following, which match the `dimnames()`
   in `hit_array`:
   
      * `cutoff_names`
      * `contrast_names`
      * `assay_names`
   
   * Then it iterates each `contrast_names`:
   
      * It returns one named vector of hits for that contrast, named
      by the measurement, with values `-1` or `1` indicating the direction.
      * When multiple `assay_names` or `cutoff_names` are involved, it
      combines results into one set representing the union of hits.
   
   * The intent is to retrieve a `list` of hits across `contrast_names`,
   regardless if there is one or more `cutoff_names` or `assay_names`.

# jamses 0.0.18.900

## changes to existing functions

* `sestats_to_df()` was updated:

   * when `se_contrast_stats()` uses different interaction statistical
   thresholds from the pairwise thresholds, the output `hit_array`
   contains NA values.
   
      * previously any non-empty vector simply counted the length
      * new behavior is to count non-NA entries
      * future "fix" will be for the cell to have length=0, with no `NA`.
      This situation only happens when `*_cutoff` and `int_*_cutoff`
      define different thresholds.


# jamses 0.0.17.900

## functions removed

* `call_fn_ellipse()` was removed and placed in the `jamba` package,
package dependency was bumped up to `jamba(>= 0.0.81.900)`

## bug fixes

* `heatmap_se()` error when not supplied `rows` nor `sestats`, now
defaults to show all `rownames(se)`.
* `heatmap_se()` error when `colorjam` package not previously loaded,
added missing package prefix.

## changes to existing functions

* `heatmap_se()` argument `normgroup_colname=NULL` and
`centerby_colnames=NULL`, instead of previous defaults which
should not be used unless relevant to the design.

## changes to existing functions

* `save_sestats()` was updated to prioritize sorting `row_type`
column headers first in the output Excel table, previously the
method looked for known patterns, which obviously failed when
`row_type` differed from those known patterns. Now the `row_type`
value itself is used in addition to the expected patterns.

# jamses 0.0.16.900

## functions removed

* `call_fn_ellipse()`

## changes to existing functions

* `heatmap_se()` new argument `correlation`

   * when `correlation=TRUE` it calculates a sample correlation
   matrix, and will display a correlation heatmap.
   * The option was added to this function in order to share much
   of the remaining logic involving data centering, and heatmap
   column annotations.
   * Note that it is possible to create a correlation heatmap using
   a subset of statistical hits, which although typically
   creates a nice-looking correlation structure, may not be that
   informative since it only includes statistical hits.

## bug fixes

* `heatmap_se()` error when `row_split` was an integer number of clusters,
and `cluster_rows` is a `function`, which is not permitted by
`ComplexHeatmap::Heatmap()`. Instead, when `cluster_rows` or `cluster_columns`
is a function, the function is evaluated upfront, so the resulting
dendrogram or hclust is passed to `Heatmap()`, which does permit
integer row and column split.

   * Note that when `row_split` or `column_split` are not single numeric
   values, the `function` must not be evaluated, otherwise that process
   fails. Haha. Makes sense. Now both scenarios are handled, hopefully
   this covers all scenarios.


# jamses 0.0.15.900

## changes to existing functions

* `heatmap_se()`

   * Argument `row_split` now also accepts named vector where names
   correspond to rownames in the heatmap.
   * Added help text for all function arguments.

* save_sestats()`

   * new argument `use_assay_suffix=TRUE` will try to add abbreviated
   `assay_names` to the end of each Excel sheet name, when there is
   more than one unique `assay_names` value.
   * The Excel sheet names are required to be unique, and have no
   more than `max_nchar_sheetname` characters in each string.
   The sheet names are progressively shortened as necessary, then
   are made unique when required by calling `jamba::makeNames()`.

# jamses 0.0.14.900

## changes to existing functions

This update attempts to extend `se_contrast_stats()` with methods
from the Bioconductor package DEqMS. This package models error as
a function of PSM_count, which is the number of peptide score matrix
entries that were combined to formulate an abundance estimate for
each protein. The approach is analogous to using `voom` as an extension
of `limma`. In fact, DEqMS is another extension of the `limma`
error model.

* `se_contrast_stats()`

   * argument `posthoc_test=c("none", "DEqMS")` to allow post-hoc tests
   such as `DEqMS` to be applied after the core `limma` steps.
   This argument is pushed into `run_limma_replicates()`.
   * argument `posthoc_args` is a list named by the `posthoc_test` value,
   containing a list of named elements used as arguments for the relevant
   post-hoc test functions.

* `run_limma_replicate()`

   * argument `posthoc_test=c("none", "DEqMS")` to allow post-hoc tests
   such as `DEqMS` to be applied after the core `limma` steps.
   This argument is also passed to `ebayes2dfs()`.
   * argument `posthoc_args` is a list named by the `posthoc_test` value,
   containing a list of named elements used as arguments for the relevant
   post-hoc test functions.
   * argument `trim_colnames` added `"sca.t"` as additional column to
   omit by default.

* `ebayes2dfs()`

   * argument `lmFit4` is optional and intended to contain additional
   output from post-hoc methods such as `DEqMS::spectraCounteBayes()`.
   * argument `trim_colnames` added `"sca.t"` as additional column to
   omit by default.
   * The order of colnames in the `data.frame` returned was updated to
   be more consistent for different formats. Specifically, when
   `posthoc_test="DEqMS"` the P-value columns sort with the `"sca.P.Value"`
   and `"sca.adj.pval"` columns before the usual limma columns.
   * when calling `mark_stat_hits()` it updates the expected colnames
   when `posthoc_test="DEqMS"`, using `adjp_colname="sca.adj.pval"` and
   `p_colname="sca.P.Value"`. In this way, hits are filtered by the
   DEqMS adjusted P-value and P-value, respectively. All other statistics
   are carried over from `lmFit3` and work as before.

* `heatmap_se()`

   * The incidence matrix data for `alt_sestats` was updated to handle
   proper hit matrix, in haste to make the `left_annotation` extensible
   the sestats hit matrix was used twice. This is why the package is
   private, for now.

# jamses 0.0.13.900

## changes to existing functions

* `heatmap_se()` was modified to handle updated `sample_color_list`
that contains `function` output from `circlize::colorRamp2()`.

# jamses 0.0.12.900

## changes to existing functions

* `heatmap_se()` was modified to accomodate other object types,
such as `"MethyLumiSet"` which does not have `rowData()`,`colData()`,`assays()`,
accessor functions, instead uses a common alternative for
Bioconductor objects: `featureData()`, `phenoData()` and `assayData()`.
The `heatmap_se()` changes will try to use those accessor functions
for any object whose class does not contain `"SummarizedExperiment"`.

# jamses 0.0.11.900

## new functions

* `se_detected_rows()` is a new function to assist the process of
defining "detected" rows which are suitable for downstream analysis.
The rules are essentially simple heuristics based upon the total
number of replicates, group replicates, fraction of group replicates,
and number of sample groups where a row has "valid" measurements.
In this case "valid" is defined as meeting a minimum abundance
threshold, and is not `NA`.

# jamses 0.0.10.900

## new functions

* `contrast2comp()`, `comp2contrast()` are two reciprocal functions
intended to convert long contrast names to a short comparison form
called a `"comp"`. The `"comp"` can be converted back to the full
contrast name. In cases where a two-way contrast can be written
in two equivalent forms, argument `factor_order` can be used to
define the desired order. The statistical results are the same,
but sometimes it helps to have specific ordering. See help docs
and examples for details, and worked examples.
Main driving motivation:

   * `ComplexHeatmap::Heatmap()` labels are too long to fit onscreen
   without making other non-ideal adjustments to plot dimensions.
   * `venndir::venndir()` labels are also too long.
   * Excel worksheets are limited to 31 characters, a typical two-way
   contrast with 3-characters for each factor level uses 35 characters,
   the `"comp"` for uses only 15.

      * contrast: `"(Aaa_Bbb-Ccc_Bbb)-(Aaa_Ddd-Ccc_Ddd)"` (35 characters)
      * comp: `"Aaa-Cca:Bba-Dda"` (15 characters)

## changes to existing functions

* argument `rename_contrasts` was added to:

   * `heatmap_se()`
   * `sestats_to_df()`
   * `save_sestats()`

# jamses 0.0.9.900

## changes to existing functions

* `save_sestats()` was updated to use `jamba-version-0.0.79-900` new
`writeOpenxlsx()`, which uses a substantially faster method of saving
worksheets, and allows building a multi-sheet Workbook prior to saving
to a file. This change avoids the time-consuming, repeated steps:
loading, saving, updating, saving, loading, updating, saving, etc.
Previously, each successive worksheet added to an Excel file took
longer to load and save. For a 15,000 row, 8 column `data.frame`,
saving once took 3 seconds, but saving 12 worksheets took 15 minutes!
The new approach took about 25 seconds.

## new functions

* `heatmap_se()` is a work in progress, but attempts to automate
display of expression data as a heatmap using `ComplexHeatmap::Heatmap()`.
By default it uses statistical hits defined in `sestats` output
from `se_contrast_stats()`.
* `call_fn_ellipsis()` is copied from `multienrichjam::call_fn_ellipsis()`.
This function allows passing ellipsis argument `...` to another
function, where that function does not allow `...` and therefore
all arguments must be strictly defined. Any arguments in `...`
which are not defined in the other function are removed before
calling that function.

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

