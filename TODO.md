
# TODO for jamses

## 08aug2025

* `se_contrast_stats()`

   * Figure out how to test a subset of rows, but retain the total
   rows in the output.
   Purpose is to mimic DESeq2 output, where it tests only a subset of rows,
   but returns summary information about all rows anyway.
   
      * Perhaps option to return all rows?
      * Should the group mean and logFC values be reported without P-values?
      Probably yes.

   * DONE. Update `se_contrast_stats()` to handle `detected_genes` as list
   alongside normgroups.
   Purpose is to apply independent detected genes for corresponding normgroups.

## 07aug2025



## 11apr2025

* `heatmap_se()`

   * Consider "trimmed down" method to indicate controls and centerGroups:
   
      * controlSamples: Thin black bar at the top or bottom.
      * centerGroups: Thin colored/greyscale rectangles?
   
   * Consider adding something like `column_sort`, `row_sort`.
   Either `character` vector matching `colnames(colData(se))`, or
   vector to use for sorting directly, or
   named vector to use where names match `colnames(se)`.

## 02apr2025

* Add `print()`/`show()` for `SEDesign` class:

   * group names, contrast names, number of samples

* S4 object `SEStats`

   * Slots:
   
      * hit_array
      * stat_dfs
      * design
      * contrasts
      * normgroup (or within "metadata" slot?)

* `plot_sedesign()` enhancements

   * Consider option to "highlight" one or more contrasts:

      * Emphasized: group boxes, contrast border, contrast color,
      contrast font.
      * Un-emphasized: greyed group boxes, greyed border/fill, hide label.
      * Hidden: optionally hide the "un-emphasized" contrasts, and factors.

   * Consider subsetting by factor levels.
   Probably not necessary if the option to emphasize/hide is available.
   * Consider easier method to assign group colors
   
      * Fix `colorset`, assign colors by factor or contrast, something like
      `sample_color_list` (list named by factor, with vector of colors named
      by factor level), or `colorSub`/`color_sub` just vector of colors
      named by factor levels.
      * Optional bonus points: when group colors are supplied, optionally
      create group colors by blending the factor level colors.
      (It might become too muddy, but would be fun to test.)

* `heatmap_se()`

   * Make it work with `DESeqDataSet` - optional log2(1 + x) transform counts.
   * Consider "easy option" to add row mean (or centered mean) row annotation.


* `groups_to_sedesign()` enhancements

   * Consider option not to use `limma::eBayes()`, passed through to
   `run_limma_replicate()`
   * Make it work for more input types:
   `SummarizedExperiment`, `SingleCellExperiment`, optionallly `Seurat`,
   `DESeq2::DESeqDataSet`, `Biobase::ExpressionSet`.
   * Is `groups_to_sedesign()` the best name? What about `make_sedesign()`?
   * Consider recognizing `numeric` columns to apply as a covariate
   in the `design`.

* `se_contrast_stats()` enhancements

   * Accept more input data types:
   `SummarizedExperiment`, `SingleCellExperiment`, optionallly `Seurat`,
   `DESeq2::DESeqDataSet`, `Biobase::ExpressionSet`.
   * Consider adding DESeq2 methods equivalent to limma:
   ```
   dds <- makeExampleDESeqDataSet(n=100, m=18)
   dds$genotype <- factor(rep(rep(c("I","II","III"),each=3),2))
   design(dds) <- ~ genotype + condition + genotype:condition
   dds <- DESeq(dds)
   resultsNames(dds)
   ```
   * Confirm that design and contrast matrices can be supplied.
   * Use `DESeq2::results()` as equivalent to `limma::topTable()`


## 18mar2025

* Vignette: "How to Create Insightful Heatmaps"

   * Centering, grouping, advanced centering
   * Subsetting by cluster, etc.

## 24jan2025

* `heatmap_se()`

   * Consider some way to indicate control samples used for centering.
   Asterisk beside label? Asterisk above/below heatmap?

* For `se_normalize()` consider accepting `batch` and `group` params
for `limma_batch_adjust` using colnames of `colData(se)`.
* Consider some "unifying framework" by which:

   * Groups are recognized.
   * Contrasts are defined using factor level changes.
   * Groups are assigned colors, based upon factor design, and factor levels.
   * Contrasts are assigned colors, using group colors. (Some rational way.)
   * Venn diagram comparisons are defined, using contrast factor changes.
   Bonus points for including twoway contrasts with component oneway contrasts.

* Extend `SEDesign` for `block` and `normgroup`

   * Optional slots which could be empty?
   * In practice, `normgroup` subdivides analyses.
   * In practice, `block` would be used during `limma` analysis.

* `SEStats` S4 object for output from `se_contrast_stats()`
   
   * Main drivers:
   
      * Need to combine multiple SEStats objects together (for convenience).
      * Need to extract component pieces easily.
      * Need safer print summary of object itself.
      * Improve default behaviors, via generic functions.

   * Proposed slots:
   
      * `"hit_array"`
      * `"stats_dfs"` - each contrast in `data.frame` format
      * `"stats_df"` - overall merged `data.frame`
      * `"sedesign"` - (with `"block"`, `"normgroup"`) - for reproducibility
      * `"metadata"` - effectively "miscellaneous" supporting data
   
   * Methods:
   
      * `contrast_names()`, `assayNames()`, `hit_names()` - core dimnames
      * `sestats_to_list()` - return hits as `list`, e.g. `hit_array_to_list()`
      * `sestats_to_im()` - above but returns signed incidence matrix
      * `sestats_to_sedesign()`, `sedesign()` - returns `SEDesign`
      * `sestats_to_statlist()` to convert to older `list` format for
      backward compatibility
      * `statlist_to_sestats()` to convert older `list` format to `SEStats`
      for forward compatibility
   
   * Concatenation by `c()` or `rbind()`, `cbind()`:
   
      * Simplest approach: Do not maintain `sedesign` or else convert to
      `list()` and store as metadata.
      * Combine `hit_array`
      * Combine `sedesign`? Or not, in case two are not compatible.
   
   * `hit_array()` - access to the array of statistical hits by dimensions:
   
      * `cutoff_name`
      * `contrast_name`
      * `assay_name`
      * `method_name` - Add this dimension to enable alternative methods
   
   * `hit_im()` - incidence `matrix` for specific dimensions in `hit_array`
   * `hit_list()` - `list` of stat hit direction, named by entity
   * `sestats_to_df()` - `data.frame` suitable for RMarkdown and `kable()`
   * accessors: `assay_names()`, `contrast_names()`, `cutoff_names()`,
   `method_names()`

## 05dec2024

* `plot_sedesign()`
   
   * Debug arguments `group_fill`, `group_border`, which cause warning/error.
   * Fix bug/warning/errors, some caused by comparing multiple value `if()`:
   `"Warning in max(subset(contrast_group_df, angle %in% c(0, 180))$bump):"`
   `"no non-missing arguments to max; returning -Inf"`

## 13nov2024

* Fix the README.Rmd.
* `heatmap_se()` - Consider recognizing `table` input.

   * 2-dimensional table could be converted to `matrix`, then assigning
   `dimnames()` to `column_title` and `row_title` if not already defined.
   * 3-dimensional table could be converted to 2-dimensional `matrix`,
   stacking each matrix slice with `rbind()`, then using `row_split`
   to represent the third dimension?

* `heatmap_se()`

   * Consider convenient way to define `controlSamples`

      * Current: Must provide `colnames(se)`, then define `control_label`
      * Goal: Define groups, control group, and consistent `control_label`
      * Suggested patterns:
      
         * First group per centering group (`centerby_colnames`).
         * Pattern-match group name within each centering group.
      
      * Define column(s) with grouping
      * Default takes first ordered group
      * Optionally specify group(s) to use as controls
      * Benefit is that `control_label` can be defined also,
      `versus [group_name]`

   * Consider some way to indicate control samples used for centering.
   
      * Some ideas:
      
         * Bold text
         * `control_suffix="*"` - add an asterisk `*` to control samples
         `"Sample B"` would become `"* Sample B"`
   
   * Consider new argument: `column_label_colname`, to match `row_label_colname`


## 28oct2024

* DONE. Prepare `pkgdown` package docs.

## 07oct2024

* Update `README.Rmd`

   * DONE. Include visuals for `plot_sedesign()`, and `sestats()`
   * Include visuals for `heatmap_se()`
   * Include examples for `heatmap_column_group_labels()`
   * Include description for `shrinkDataFrame()`
   * Consider examples using `colorjam::col_div_xf()` within a customized
   `ComplexHeatmap::Legend()`.

* `save_sestats()`

   * DONE. Consider option to populate the "hits" worksheet using `logFC` values
   instead of only using `c(1, 0, -1)`.

* Consider new helper functions

   * `contrast_list_by_factor()`: subdivides a set of contrasts by which
   factor is being compared
   * `sestats_to_dfs()`: wrapper for `save_sestats()` which returns
   the list of `data.frames`.

* Consider adding conversion of `Seurat` to `SingleCellExperiment`

   * Key addition: maintain all assay matrix objects using their original
   names: `layername_assayname`
   * Option to apply `log2(1 + x)` if range exceeds a threshold (e.g. 50)
   * Option to convert sparse `Matrix` objects to vanilla `matrix`.

* Add tests

   * `heatmap_se()` basic workflows
   * Conversion of `Seurat` to `SingleCellExperiment`

## 18sep2024

* Adapt `heatmap_column_group_labels()` for `Seurat` and `SingleCellExperiment`
* `plot_sedesign()`, `groups_to_sedesign()`

   * Consider option to return/filter "useful information":
   
      * Filter for minimum number of replicates in a group to be permitted
      in a contrast.
      * Return number of replicates per group in a contrast.

## 16sep2024

* DONE. Adapt `heatmap_se()` for `SingleCellExperiment` objects.

   * DONE. By proxy, it could also work for `Seurat` objects.

## 10sep2024

* DONE. Fix error with `save_sestats()` when `rowData_colnames` contains
columns already present, causing them to be added twice.
The error is caused with `writeOpenxlsx()` by duplicate colnames.

## 14aug2024

* `plot_sedesign()` - Improve the method for bumping arrows.

   * Currently all contrasts that share the same y-position or x-position
   are bumped relative to each other. The goal is to bump contrasts
   only when they overlap another contrast.
   * When B-A:X and D-C:X are on the same y-position, but one does not cross
   the other, they should not bump each other.
   * The new algorithm may need to account for x,y,intercept and bump number,
   so that a contrast is only bumped when it matches all three values.
   This way, when one contrast is bumped, it should free some space for other
   contrasts. It may also need to be adjusted iteratively.
   * Interesting to consider whether to begin with longest or shortest
   contrasts (by Eucliden distance) in order to control the symmetry.

* Consider method to combine two `SEStats` objects

   * Typically expected to append contrasts across two objects
   with the same hit thresholds, assay_names. Not required, however.

* `heatmap_se()`: Consider option to add item count to `row_title`.

   * New argument? `tabulate_row_split=TRUE`
   * Get the row split information, tabulate number of rows per split,
   define `row_title` to include the number in parentheses.
   * Label should be "(1,234 `row_type`)".
   * Bonus points: Trim trailing "s" if there is only one row.


* DONE. `save_sestats()`

   * DONE. Debug misalignment of sheet name with contrasts. Sigh.
   * DONE. Consider optional methods to shorten the sheet names
   
      * custom function to edit the sheet name before
      truncating to `max_nchar_sheetname` number of characters.
      For example to edit things like `"male", "female"` to `"m", "f"`
      to save space.
      * option to abbreviate each term to use only the first letter
      or first N unique letters?

## 08aug2024

* DONE. `matrix_normalize()` and `se_normalize()`

   * DONE. Consider option to define `reference_samples` during normalization.
   CORRECTION: Test existing option `controlSamples` to confirm it works
   as intended for this purpose. CONFIRMED.

      * The goal is to normalize only the `controlSamples` to themselves,
      then normalize all other samples relative to that subset.

## 24jun2024

* Consider function to take `SEStats` and reverse/flip contrasts
using a preferred set of contrasts.

   * Use case is when receiving data with (groupA-groupB) but the preferred
   order is (groupB-groupA).
   * Function iterates `SEStats` and finds matching contrasts to reverse.
   * All hit signs and fold changes are multiplied by `-1` to flip the sign.
   * All matching contrast names are reversed.
   * All colnames which include the contrast also have the contrast flipped.

* Consider new function to take `list` of `data.frame` with statistical
results and create `SEStats` object.

   * It should parse each `data.frame` and create `"hit "` column if needed.
   * It should define `hit_list` for each `data.frame`

* Consider function to manipulate dimnames in `SEStats`: contrast_names,
assay_names, cutoff_names.

   * Basic example is to convert `contrast` to `comp` or vice versa.
   * Could be through accessor `contrast_names(SEStats)` and
   `contrast_names(SEStats)<-`
   * Rename `assay_names` as needed.
   * Rename `cutoff_names`?? Perhaps not rename.
   

* Consider function to "validate" `SEStats` object, suggested rules:

   * `hit_array` contrast_names must match `hit_list` and `stats_dfs`.
   * `hit_array` assay_names must match `hit_list` and `stats_dfs`.
   * `hit_array` cutoff_names must match `hit_list` and `stats_dfs`.

* Consider function to `c()` multiple `SEStats` objects together.

   * It would combine `hit_list`, throw error whenever the input
   and new data both have: `assay_name`, `contrast_name`, and `cutoff_name`.

* Consider function to subset `SEStats` by dimensions:
`assay_name`, `contrast_name`, `cutoff_name`

   * For example: `SEStats[1, 4:8, 2]`

## 31may2024

* Add testing for all variations of `se_contrast_stats()`

   * create test case with `NA` values, to test `handle_na`

## 30may2024

* DONE. Debug `se_contrast_stats()` discrepancies when `isamples` is provided
in different orders. All outputs should be identical regardless of input
order. Bug was `handle_na_values()` and was fixed, unclear when it was
introduced.

## 15may2024

* DONE. `save_sestats()`

   * include a summary table with the `"hit"` column from each contrast

## 23apr2024

* Design idea for `plot_sedesign()`

   * Allow and document how to create a multi-panel plot, with this
   plot as one panel in the output?
   
      * The driving use case is to display contrasts alongside a Venn diagram
      or heatmap.

   * Option to plot one contrast, showing the arrow (or arrows) and
   displaying at the top the full contrast name, optionally
   the abbreviated "comp".
   * Potential to create a small "logo" for a single contrast, that
   could be displayed as a thumbnail in another slide.
   * Potential to create a small "logo" for a few contrasts,
   for example alongside a Venn diagram comparing two or three contrasts.
   The thumbnail could represent which contrasts are being displayed.

## 16apr2024

* Implement `testthis` unit tests for `se_contrast_stats()`

   * Highest priority: test each option of `handle_na`
   * Test `use_voom`.
   * Test `block`.
   * Test `normgroup`.
   * Test `block` and `use_voom` together, `voom_block_twostep=TRUE`,
   which should call correlation twice.
   * Test `block` and `normgroup` notably when one `normgroup` contains 2+
   unique `block` values, and another `normgroup` has only 1 `block`.
   * Bonus points: Test using a subset of `isamples` that no longer contains
   a valid contrast. It should fail - though in future it could potentially
   return the group mean values, then stop short of performing contrasts.

* Extend `SEDesign` for `block` and `normgroup`:

   * Pressing need to maintain `normgroup` and `block` with `sedesign`.
   Basic goal: when `sedesign` is subset, `normgroup` and `block` are also
   subset consistently.
   * New slot names:
   
      * `"normgroup"` - `character` vector to match `samples()`, and
      `rownames()` of the design matrix.
      Or could it be `data.frame` to permit multiple column values?
      Upon use each row would be concatenated to make one value per sample.
      When `sedesign` is subset, `"normgroup"` is also subset consistently.
      * `"block"` - `data.frame` with one or more columns, indicating.
      Its primary purpose is to maintain values per `samples()`, so they
      can be subset and maintained consistently.
      It will be pushed to `se_contrast_stats()` how to deal with the
      actual values:
      
         * `character` values are considered `block` covariates for `limma`.
         For now, only `character` will be permitted.
         * `numeric` values can be encoded as a scalar covariate, and
         would be appended as a new column in `design` - and therefore must
         be added as a new row in `contrasts` with empty `0` values.
         * `integer`, `factor` values are encoded as an ordinal covariate,
         converted to rank integer values.
   
   * `normgroup()`, `normgroup()<-` - set/get functions for `sedesign@normgroup`
   * `block()`, `block()<-` - set/get functions for `sedesign@block`

* `groups_to_sedesign()` to handle `normgroup`, `block`

   * When `normgroup` is supplied, contrasts should be limited to those
   within each `normgroup`.
   It may be accomplished by including `normgroup` as factor columns,
   but not included with `factor_order` so that comparisons cannot
   involve multiple `normgroup` values.
   * Store `normgroup` in the `SEDesign` object, see above.

## 12mar2024

* New function idea: `heatmap_to_xlsx()` or `heatmap_to_df()`
to convert `heatmap_se()` output to `data.frame` or save with
`jamba::writeOpenxlsx()`.

   * Idea is to extract data from the heatmap as displayed, in order to
   save it to a file (e.g. Excel or tab-delimited) for later review.
   * Use techniques to extract the `left_annotation` and `top_annotation` data.
   * Annotation colors can be extracted and used to colorize corresponding
   cells in Excel.
   * It should mainly target features made available by `heatmap_se()`,
   in other words any custom `HeatmapAnnotation()` functions would likely
   not be supported. Things like annotation bar charts, line plots, etc. would
   not be easily supported for export.
   * Main features:
   
      * Rows and columns match the order as they appear in the heatmap.
      * Left annotations: `rowData_colnames`, `sestats`, `sestats_alt`
      * Top annotations: `top_colnames`
      
         * As below, it is unclear how best to handle multiple header rows,
         it makes the output file much more difficult to use,
         whether it is saved as tab-delimited text, or Excel.
      
      * Optional: Extract heatmap color function, apply to each cell.
      (This step sounds interesting, but is likely to be a bad idea.
      Applying a color per cell in Excel would be very time-consuming
      and inefficient, and would produce a much larger file.)
      * Optional `column_labels`, `row_labels`
      * Row split encoded as a row annotation, using values from `row_title`.
      * Column split - unclear the best approach to use:
      
         * Adding column annotations requires saving multiple header rows,
         which makes the exported file substantially more difficult to use.
         * One option is to use column split as a prefix, e.g.
         "Group 1 - colname1", "Group 1 - colname2", "Group 2 - colname3", etc.
      
      * Option to supply the `se` object, and append additional data
      via `rowData_colnames`.

   * Required: Ability to save `HeatmapList` objects, already drawn using
   `ComplexHeatmap::draw()`.
   
      * Bonus points: recognize when multiple heatmaps are present, e.g. from
      `ComplexHeatmap::draw(hm1 + hm2 + hm3)`.
      * It needs to iterate each heatmap in the `HeatmapList`, potentially also
      obtaining `left_annotation` and `top_annotation` from each `Heatmap`.
      Ugh. Each heatmap could have different `top_annotations` which means
      the `top_annotation` would need to have a column with row names for
      each heatmap's unique `top_annotation`.

## 07mar2024

* `groups_to_sedesign()`

   * Consider adding argument `normgroup` to restrict contrasts within
   each `normgroup`.

* `heatmap_se()`: Consider some way to indicate `controlSamples`.

   * Note that sometimes `controlSamples` are not displayed on the heatmap,
   this might be an exception where it is not shown.
   The rule might be: if not all `controlSamples` are in `isamples` then do
   not indicate `controlSamples`.
   * One option is to render an asterisk just above each column? `"*"`

## 06mar2024

* `se_normalize()`

   * DONE. Consider allowing user-specified `assay_name` to store the
   resulting data? Some mechanism to allow customizing the output `assay_name`.
   New arguments:
   
      * DONE. `output_method_prefix` which by default uses each value in `method`
      to formulate each output `assay_name`
      * DONE. `output_assay_names` 
      
         * specify completely user-defined output `assay_name` values,
         overriding `method` and `output_method_prefix` when defining
         the output `assay_name`.
         * `character` vector, must be equal to the number of normalizations.
         Note that each normalization is applied to each `assay_names` in order:
         `method1_assay1`, `method1_assay2`, `method2_assay1`, `method2_assay2`.
   
   * DONE. Consider using `mcols(assays(se))` to store annotation regarding
   the normalization:
   
      * DONE. store values: `"normalization_method"`, `"source_assay_name"`,
      `"params"` (params may not be easy to use, since params is a `list`
      whose elements may include `character` vectors, `numeric` values with
      many decimal places, etc.
      * DENIED. Consider comparing new `mcols()` entries to existing `mcols()`
      entries to decide whether to overwrite the existing entry, or to
      create a new `assays()` entry with versioned name and different params.
      Decision: Do **not** make this change, push this back onto users
      to define `output_method_prefix` or `output_assay_names`, with the
      idea that the user should specifically request creating a new
      `assay_name`.
   
   * DENIED. If `method="jammanorm"` is called with different params, should it
   create a new `assay_name` entry associated with those params?
   **Decision: No.**

* `heatmap_column_group_labels()`

   * Consider permitting a subset of `se` columns being passed, to
   allow labeling only a subset of columns per operation.
   Ignore columns which are not provided in the `se` object.
   * DONE. Consider option to hide the group line, when the group label
   is also empty or whitespace.

* `heatmap_row_group_labels()`

   * New function completely analogous to `heatmap_column_group_labels()`.
   * Default `rot` text rotation is 180 degrees (sideways).
   * Consider adding something like `hm_title_buffer` to adjust the
   buffer to the left of `left_annotation`, `HeatmapBody`?

## 26feb2024

* `se_normalize()`

   * DONE. new methods from `edgeR::calcNormFactors()`:
   
      * "TMM": trimmed mean of M-values
      * "TMMwsP": modified TMM with singleton pairing, removes zeros
      * "RLE": relative log expression
      * "upperquartile": uses upper 75% (by default) or user-defined threshold

   * DENIED. Consider composite `assay_name` with param values?
   
      * for example running the same normalization with two sets of `params`

* DENIED. consider renaming `"limma_batch_effect"`

   * something shorter, and without underscores
   * `"lba"`, `"limmabatchadjust"`
   * (Instead user can customize the resulting assay_name.)

* DELAYED. Consider something like `list_normalization_methods()`

   * purpose is to list available normalizations as a programmatic reference
   * returns `data.frame`: abbreviation, full_name, description
   * the `abbreviation` ultimately becomes the new `assay_name`
   * Currently `jamba::jargs(se_normalize, "method")` will print available
   options, except without descriptive information. `? se_normalize`

## 23feb2024

* `heatmap_se()`

   * When an entry in `sample_color_list` contains a color function, the
   breaks are taken from `attr(x, "breaks")`, and it uses the same
   values as labels.
   
      * Consider recognizing optional attribute with "labels" to use:
      `attr(x, "labels")`. It must be the same length as `"breaks"`.

## 25jan2024

* `se_contrast_stats()`

   * when using blocking factor (argument `block`) the group mean values
   slightly differ from what could be calculated manually, even when
   applying limma batch adjustment.
   * Excellent overview of experiment design and contrast models in R:
   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7873980/
   doi: 10.12688/f1000research.27893.1


## 10jan2024

* `heatmap_se()`

   * DONE. Consider using `"legend_title"` instead of `"name"` as the heatmap
   legend title, so that `"name"` can be defined with a specific value.
   The driving need is when adding two heatmaps together, if they both
   have the same `"name"` the name appears twice when calling
   `ComplexHeatmap::list_components()`. Instead, we should allow `name`
   to be defined uniquely for each heatmap, even while both heatmaps
   share the same legend title.
   * Consider convenient method and default to add separation between
   heatmap and annotations, e.g. `ht_opt(ROW_ANNO_PADDING=grid::unit(4, "mm"))`
   and `ht_opt(COLUMN_ANNO_PADDING=grid::unit(4, "mm"))`.
   It should set the value, then revert to previous value afterward.

## 08jan2024

* DONE: debug `se_contrast_stats()` with `handle_na="full1"`

## 11dec2023

* `contrasts_to_venn_setlists()`

   * DONE. Assign list names based upon the contents of each Venn setlist,
   instead of the default names that are not user-friendly.
   * Consider "two-way" Venn subsets that include:
   two-way contrast,
   corresponding one-way contrasts.

* `plot_sedesign()`

   * When `contrast_depths=2` and supplying `sestats` it labels the one-way
   and two-way contrasts, but should only label one-way contrasts.
   The workaround is to alter `sestats$hit_array` to include only the
   contrasts being displayed.

* `heatmap_se()`

   * consider some method to "group" the `sestats` contrasts, for example
   sub-grouping them by which factor(s) are being compared.
   Use `contrasts_to_factors()` to generate a table, then split by which
   factors have comparisons (with delimiter `"-"`).

## 04dec2023

* `contrast2comp()` and `comp2contrast()`

   * consider method that can convert contrasts inside stat colnames:
   
      * `"logFC factorA_factorB-factorC-factorB"` to
      `"logFC factorA-factorD:factorB"`
      * `"hit mgm5 adjp0.05 fc1 PNU_Control_Nano147car1-Veh_Control_Nano147car1"`
      to `"hit mgm5 adjp0.05 fc1 PNU-Veh:Control:Nano147car1"`
      * It therefore needs to detect the location of the contrast or comp,
      convert that substring, then replace with the new value.

* migrate `platjam::design2colors()`

   * Consider adding SE-specific method
   
      * assign colors to `rowData()` and `colData()` in one step
      * bonus points for sharing color assignments for shared terms 

## 09nov2023

* add documentation for generic methods, `contrasts()`, `contrast_names()`,
`groups()`, `samples()`, etc.
* consider DESeq2 support

   * Research the equivalent design model, contrasts model for equivalent
   comparisons in DESeq2 as used with limma-voom. Often people seem to
   include only groups relevant to each contrast in the DESeq2 workflow,
   unclear if this is recommended guidance.
   
      * Guidance is similar to limma-voom, include each sample group as
      relevant to the overall comparisons, using `normgroup` to separate
      subsets that are expected to differ substantially by sample type,
      or where each subset may represent different variabilities.
      * Therefore, this option would initially follow the limma-voom default
      with user-defined, optional use of `normgroup` which can define
      independent subsets of sample groups.
      * There should probably be an option to perform each contrast in
      its own independent subset, which would not be equivalent to `normgroup`,
      for example "untreated", "treated1", "treated2" would imply three
      one-way contrasts, and each contrast would be performed in its own
      unique group-to-group analysis.
   
   * Include options specific to DESeq2:
   
      * lfcShrink method
      * normalization method (or none)
      * optional outlier removal (e.g. Cook's distance method)

   * Decide how to handle stats column headers which differ from limma,
   so they are easily used in downstream methods: `"P.Value"`, `"adj.P.Val"`,
   `"logFC"`, etc.
   * Consider returning dispersion data sufficient to produce the dispersion
   QC plot to review. (Consider doing similar for `voom()` internals.)
   * Consider including "method" as a dimension in `hit_array`,
   so that limma-voom, DESeq2 could both be used and directly compared.
   Similarly the `posthoc_method="DEqMS"` may be included so that
   its effects can be compared to limma without the custom posthoc adjustment.


## 18oct2023

* `plot_sedesign()`

   * Option to filter by `max_depth` (to show only oneway contrasts)
   or `contrast_depths` (to show only oneway, or only twoway contrasts).
   * When `group_buffer` is adjusted, the replicate `"n=3"` labels are
   not shifted, so they can appear outside the square.
   * Consider option to adjust aspect ratio.
   * Consider option to adjust the offset between contrasts.

* `filter_contrast_names()`

   * Consider option to enable sequential comparisons, for example:
   Time1-Time0, Time2-Time1, Time3-Time2. This option may be more useful
   when enabled only for certain factors. Not as useful for something
   like "Treatment" where each treatment usually does not have a sequential
   order.

## 06oct2023

* `heatmap_se()`

   * Add optional argument to apply heatmap title to `column_title`,
   so that it appears above the heatmap automatically.
   **It will override showing `column_split` titles.** So it may be best
   practice to use `heatmap_column_group_labels()`
   * Consider argument to add whitespace padding to the heatmap title,
   to make it more convenient to use `heatmap_column_group_labels()`.

## 04oct2023

* Done. fix bug in `contrasts_to_venn_setlists()` resulting in some sets with
>4 entries.
* Consider integrating `sestats` with `contrasts_to_venn_setlists()` so that
the output is actually a `list` suitable for `venndir::venndir()`.

## 02oct2023

* `contrast_names()<-`

   * validate the supplied `value` and print error message when incompatible.

* `plot_sedesign()`

   * consider easy filter to hide two-way contrasts
   * consider scaling the arrow width (decreased) when there is
   a large number of contrasts to display, especially when large number
   are "bumped" per factor level.
   * consider making `arrow_ex`,`head_ex` arguments to `make_block_arrow()`
   visible in `plot_sedesign()`. Perhaps change `twoway_lwd` to `twoway_cex`,
   then applying some global `cex` to all contrast arrows and connectors.
   * Done. Fix bug with custom contrasts, two-way contrast color `NA` throws error.

* `validate_sedesign()`

   * argument `contrasts` is not properly subsetting contrasts by name,
   however `contrast_names()<-` appears to be working fine. Same mechanism.

* Done. `filter_contrast_names()`

   * new function: take long list of "all versus all" contrasts and subset
   for specific control factor levels.

## 20sep2023

* `heatmap_se()`

   * consider indicating `controlSamples` when a subset of samples are used,
   as opposed to relying upon `control_label`. Especially useful when
   a subset of potential samples are used, or when there are multiple
   normgroups displayed, each with their own `controlSamples`.
   * consider option to apply `hm_title` to `column_title` which may
   help when adding two heatmaps together, they would have the
   column_title defined within each heatmap itself.
   Downside: no column titles defined by `column_split`.

* `SEDesign`: consider associating blocking factor to one or more contrasts

   * Unclear how: blocking factor is per-sample annotation
   * `se_contrast_stats()` does not indicate blocking factor - this is the
   real need, so output represents the comparison.

* `SEDesign`: consider associating factor name with each contrast

   * Background: Each contrast compares one or more factors.
   It would be useful to "know" the factors being compared during
   automated analysis.
   
      * Contrasts could be subset by factor(s) being compared
      * Center data using other factors in the design:
      when showing Treatment contrasts, could center by Time or Genotype
      * Venn diagrams make the most sense when comparing within factor:
      Treatment contrasts to compare across time points or genotypes.
      Gene knockout comparisons across treatments (Vehicle, DEX, Etoposide)

   * Proposed changes to `SEDesign`
   
      * new slot: `"contrast_factors"`, `list` named by `colnames(contrasts)`
      values contain zero or more `colnames(factors)`
      * new slot: `"factors"`, `data.frame`
      
         * `colnames(factors)` are experimental factors (Treatment, Time)
         * values are factor levels inferred during `groups_to_sedesign()`
         * rows are samples, `rownames(factors) == rownames(design)`

      * subsetting SEDesign

         * by sample would also subset `factors`
         * by group would also subset `contrasts` and `contrast_factors`

* `SEDesign`: consider supporting `~Treatment + Time` style design matrix

   * Current design and contrast matrices use `~ 0 + Treatment + Time`,
   so that contrasts indicate distinct experiment groups.
   * `~Treatment + Time` format changes:
   
      * `design` includes `"(Intercept)"` then factor levels
      * instead of groupA-groupB contrasts, it encodes coefficients
      using factor levels.

   * `plot_sedesign()` needs to handle this format differently,
   visualizing factor comparisons differently as well.
   * Two-factor designs might compare one factor, which implies using
   the other factor as a covariate term, equivalent to modeling
   one factor with the second factor as a blocking factor.
   * Unclear how to represent a specific interaction term if encoded as
   `~Treatment + Time + Treatment:Time`

## 31aug2023

* Method to `rbind()` or `cbind()` SummarizedExperiment objects.

   * Main goal is to automate the process of aligning samples
   * Keep or remove metadata columns that do not match
   * When metadata does not match:
   
      1. Remove the metadata column altogether (easiest option), or
      2. Keep metadata, but "combine" values by comma-delimiting (easiest),
      or by some numerical operation (complex to design R function interface).

## 25aug2023

* Create `SEStats` object as more formal S4 object output `se_contrast_stats`()

   * slotNames

      * `stats_dfs`: each `data.frame` from each `contrast_name`, `assay_name`
      * `hit_array`: N-dimensional array of stat hits with direction:

         * `assay_name`
         * `contrast_name`
         * `cutoff_name`
         * `method_name`? (to compare limmavoom, limma, DESeq2, edgeR?) adding
         this dimension could be fairly invisible

      * `stats_df`: **omit** - separate function to merge `data.frame`
      * `metadata`: method, parameters, etc.
   
   * methods
   
      * `hit_list()`:       calls `hit_array_to_list()`
      * `to_df()`:          calls `sestats_to_df()` to create table summary of counts
      * `hits()`:           converts `hit_list()` into incidence matrix?
      * `contrast_names()`: extracts `contrast_name` vector
      * `assay_names()`:    extracts `assay_name` vector
      * `reapply_cutoffs()`: (new) re-calculate `hit_array`,
      by iterating `stats_dfs` and applying stat cutoffs.
      * `rbind_sestats()`:  combines multiple `SEStats` objects
   
   * Other related todo:

      * consider backward compatibility?

         * `SEStats_to_list()` to convert `SEStats` to previous `list` format
         * `list_to_SEStats()` to convert previous format to new `SEStats`
      
      * `heatmap_se()` should accept `SEStats` input
      * `hit_array_to_list()` should accept `SEStats` input
      * `sestats_to_df()` should accept `SEStats` input

* `heatmap_se()`

   * use new `SEStats` object input
   * consider expanding `sestats` to show separate `cutoff_name` and
   `method_name` entries as distinct stripes in the incidence matrix,
   rather than including hits across `cutoff_name` entries together.
   The labels could become quite long.

* Expand `se_contrast_stats()` to call corresponding DESeq2 methods.

## 22aug2023

* `se_contrast_stats()`

   * For very large data volume, the method seems to take more memory than
   absolutely necessary and could be trimmed:
   
      * Remove `rownames()` for `stats_dfs` and `stats_df`.
      Row identifiers are roughly 50% the size of each `data.frame`,
      so they should not be stored in a column and as rownames. The
      rownames are less "safe" to R manipulation, so values will be
      retained in a specific column (first column).
      Other functions that utilize `data.frame` objects must use column
      values and not rely upon rownames, in theory should already be true.
      * Option to omit `stats_df`, roughly 25% overall object size.
      Matter of fact, this object could be replaced by a function that
      converted `stats_dfs` into `stats_df` dynamically.
      * option to omit `stats_dfs`, roughly 75% overall object size,
      but used to create volcano plots, and to review specific results.
   
   * When using `block` the process becomes substantially slower with large
   data.
      
      * Description of the scenario, and supporting evidence:
      
         * The slow step occurs during `limma::lmFit()` when supplied with `block`
         and when `correlation=NULL`, so it is calculated inside `lmFit()`.
         * The actual rate-limiting step is `limma::duplicateCorrelation()`
         which is run to determine `correlation`.
         * In one test case with 220k rows, 30 columns, this step took 7 minutes.
         * A subset with a random 10k rows took 15 seconds, and estimated
         the same summary correlation value.
         * When provided the pre-calculated `correlation` the `lmFit()` using
         220k rows took 1.5 seconds.

      * Potential workarounds:
      
         1. When `block` is defined, and `correlation` is not supplied,
         calculate `correlation` using a subset of up to N rows (e.g. 10000).
         The correlation could even be calculated 10 times using random
         subsets of 1000 rows, then take the average.
         The max rows could be a new argument `max_correlation_rows=10000`
         to make this process explicit, and customizable.
         2. The process above could be "pushed back" to any functions calling
         `se_contrast_stats()` to avoid having an invisible difference
         between using this function, and using `limma::lmFit()`
   
   * Need to store arguments such as `correlation`, and `block` alongside
   the returned `sestats` object.

* Debug why loading `jamses` causes the warning to the effect:
`"design() has already been defined"`.

## 14aug2023

* `plot_sedesign()`

   * Option to size block arrows by the relative number of hits?
   * Option to define block arrow sizes as a vector, one per arrow,
   which probably means one per contrast name.

* `heatmap_se()`

   * Consider option to cluster rows by stats hit matrix data. Unclear how,
   but would be useful to sub-group hits.

* `se_normalize()`: use `mcols(assays(se))`:

   * Driving use case:
   
      * somewhere to store parameters used during analysis.
      * `metadata(assays(se)[[1]])` cannot be used, since metadata on
      the matrix itself is easily lost upon any sort of manipulation,
      subsetting of the matrix.
      * secondary use case: ability to filter for normalized data
      
   * upon creating a new entry in `assays(se)` it should also populate
   `mcols(assays(se))` with columns of annotation. These are optional
   annotations, but could be convenient for including things like
   arguments to the methods used.
   * Consider updating `platjam` data import methods similarly.
   * Suggested annotation names:
   
      * `"method"`: name of the normalization method applied
      * `"params"`: `list` of parameters for the given method,
      obtained from the argument with `params[[method]]`
      * `"parent_assay_name"`: may not be practical. The `assay_name`
      can be edited, which would invalidate the relationship.
      * `preferred`: optional flag to indicate the preferred `assay_name`
      for downstream analysis?

* `se_normalize()`: change default `output_sep="_"` to `output_sep="."`?

   * benefit is that method names (containing underscores) are more
   easily distinguished when concatenated:
   
      * `"jammanorm.limma_batch_adjust.totalIntensity"` versus
      * `"jammanorm_limma_batch_adjust_totalIntensity"`

## 27jul2023

* Remaining TODO for `plot_sedesign()`:

   * Fix error when any of axis1,axis2,axis3,axis4 are assigned
   empty values.
   * Options for rendering two-way contrasts:
   
      * PARTIAL. Consider option to "loop" two-way contrasts around the one-way
      contrast, similar to what happens when two contrasts are too
      close to each other.
   
         * Loops are implemented but did not meet the favor of artistic review.
         * For contrasts with very small gap between them, loop seems the best.
         * For contrasts with room for an S-shaped swoop, that looks nicer.
   
      * Option to adjust the midpoint of S-shaped "swoop",
      so the middle of the swoop would not be the diagonal
      midpoint between the end of contrast 1 and start of contrast 2.
      It  could be shifted closer to contrast 1 or contrast 2.
      It could help when trying to avoid label overlaps.

   * Call `sedesign_to_factors()` instead of calculating internally.
   * Consider refactoring the drawing order so that contrasts are drawn,
   then all labels are (potentially) drawn atop the graphics so they are
   more consistently visible. (Could be done if porting to grid graphics.)
   * Accept contrast_names as input, instead of requiring an `sedesign` object.
   * Consider using grid graphics with
   with the `units="snpc"` pattern to maintain fixed aspect ratio:
    ```R
    pushViewport(
    viewport(
       x=0.5, y=0.5,
       width=unit(min(1, diff(xlim)/diff(ylim)), "snpc"),
       height=unit(min(1, diff(ylim)/diff(xlim)), "snpc"),
       xscale=xlim,
       yscale=ylim))
    ```
   * Consider drawing boxes which reflect the `groupedAxis()` grouped
   regions along each axis, which may involve slightly nested boxes.
   It may be useful to shade the boxes light grey.

## 12jul2023

* enhance `sedesign`

   * Consider adding `factor_labels` to have design labels for each
   factor in the group label.

* DONE (initial implementation): new function `plot_sedesign()`

   * Simplified version of `plotComparisonTable()` from past work.
   * Factor levels should be "known".
   * Layout should follow logic similar to (or directly extending)
   `vcd::mosaic()` which defines factors on x-axis, y-axis, then
   subdivides each axis in order with sub-factors.
   * Layout should define every observed combination of factor levels.
   * Then block arrows (or regular arrows) can be used to indicate each
   pairwise contrast.
   * Block arrows are "bumped" when they would overlap another block arrow.
   * Bonus points: Block arrow labels can include the contrast name,
   or comp name (shortened version of the contrast), and optionally
   the number of statistical hits.
   * Bonus points: Indicate the number of replicates per group.
   * Two-way contrasts would indicate both one-way contrasts, with
   a "ribbon" connecting the end of the first contrast to the
   beginning of the next contrast. Something like a
   beta-sheet in protein structure graphics.
   * Most contrasts should only be vertical or horizontal, therefore
   non-standard contrasts would be angled, indicating that they
   are comparing more than one factor at a time.

* Completed TODO for `plot_sedesign()`:

   * DONE. Option to indicate "n=8" number of replicates per group.
   * DONE. Provide or add custom label; option to hide labels.
   * DONE. Show subset of contrasts based upon contrast name.
   * DONE. Consider option to "flip" two-way contrasts:
   (A-B)-(C-D) is equivalent to (A-C)-(B-D) and can be shown either way.
   * DONE. Decide how to assign colors to contrasts. For now, colors
   are assigned to one-way contrasts, then two-way contrasts inherit
   two colors as a gradient.
   * DONE. Adjust axis labels to be adjacent to the square boxes used per group.
   This step used argument `pos` which places the axis at a fixed coordinate
   inside the plot.
   * DONE. Allow custom position for contrast labels, scaled from 0 to 1
   indicating the position along each contrast (from-to) to place the
   text label. This step also applies text justification, which helps
   prevent text from spilling beyond the relevant end of the arrow.
   * DONE. Add argument `sestats` to print the hits per contrast.

      * DONE. Define how to handle multiple `assay_names`, `cutoff_names`.
      It passes arguments to `hit_array_to_list()`.
      Anything more complicated requires the user to pass argument
      with custom labels.

   * DONE. Sort contrasts:
   
      * by axis display order; currently axes are ordered in the order
      that factor levels appear in the contrasts.
      * then by "distance" along each axis,
      so longer contrasts are "bumped" consistently.

   * DONE. Consider making hit direction optional when `sestats` is used.

## 27jun2023

* `heatmap_se()`

   * Consider automating color assignment when `top_colnames` or
   `rowData_colnames` have no colors assigned in `sample_color_list`.
   It may involve calling `platjam::design2colors()` with empty
   argument for `group_colnames`, which may involve moving that
   function into `colorjam`.

* DONE: `heatmap_se()`

   * DONE: when `correlation=TRUE` the label with number of rows should
   indicate number of columns? See below, both dimensions are indicated.
   Or number of rows used to calculate correlation of N columns?
   * DONE: Consider new argument `column_type="samples"` and default heatmap
   title also indicates the number of columns (samples).

## 05jun2023

* DONE: `se_collapse_by_column()`

   * DONE: consider changing default `noise_floor_value=NA` to
   `noise_floor_value=0`.
   * DONE: consider new argument `useMedian=FALSE` so that the default
   behavior is to take the mean and not the median value per group.
   * DONE: change to use `jamba::call_fn_ellipse()` so that calls to row
   group functions will only pass arguments accepted by that function.

## 19may2023

* DONE: add `heatmap_column_group_labels()`

   * Custom function to augment `heatmap_se()` for a specific scenario,
   but the scenario is used often enough to warrant making the
   function available here.
   * The function draws group labels above a rendered heatmap, with
   underline drawn by default indicating the heatmap columns under each
   label. It draws multiple levels of group labels when defined.
   * It uses `grid` coordinates *after rendering* to determine where to
   position labels, which requires drawing the heatmap first.
   * Also the group labels cannot be included in the `Heatmap` object,
   which is not ideal, but is also a known limitation.

## 05may2023

* `se_contrast_stats()`

   * DONE: add optional `rowData()` colnames to stat `data.frame` output,
   for example adding `rowData_colnames=c("SYMBOL", "GENENAME")` would
   keep gene symbol and gene name alongside microarray probe IDs.

   * in future, optionally run DESeq2 equivalent steps to limma/limma-voom,
   by replacing the `run_limma_replicate()` step with optional function
   to wrapper DESeq2 steps.

* tests for `se_contrast_stats()` with various `handle_na` values,
and with/without `rowData_colnames`.

## 05apr2023

* DONE: use `testthat` unit testing
* `sedesign` object

   * DONE: add method `contrastNames()` (or `contrast_names()`)
   * consider adding `comps()` as shortcut for `contrast2comp()`

* `contrast2comp()`

   * optimize performance, it is surprisingly slow, but functional
   * consider embedding the factor order into the output for two-way
   contrasts, to ensure the output exactly matches input even when
   alternate outputs are mathematically equivalent.

* add `plot_sedesign()`, `plot_contrasts()`

   * display design with similar layout orientation to `vcd::mosaic()`
   * indicate contrasts with block arrows drawn between groups
   * optionally label groups by number of replicates
   * optionally label contrasts by number of statistical hits

* DONE: add `contrast_names_to_sedesign()`

   * convenience function to produce `sedesign` from only `contrastNames`
   * mostly helpful for demo purposes, development of `plot_sedesign()`


## 28mar2023

* `heatmap_se()`

   * Add optional title above `sestats` incidence matrix display,
   which could address need to display which `assay_name`, and `cutoff_name`
   was used.
   * Might be useful to have additional option to provide a label to display
   atop the `sestats` incidence matrix.
   * Probably needs other options such as font sizing, label rotation, etc. Oi.

## 12mar2023

* `heatmap_se()` improvements

   * DONE: Some method to hide `column_title` labels. Use `column_title=" "`.
   * Consider `hmgrouplabel` function `heatmap_column_group_labels()`
   which enables cleaner column labels when the heatmap is split by one
   or more variables in `colData(se)`.

## 22nov2022

* Add concept of `"normgroup"` to `sestats_to_df()`

   * Implied: `"normgroup"` should also be added to `sestats`.

* Add per-gene logic and functions developed for the DM-JDM manuscript.
* Bonus points: optionally print commands as they are being run, for example:

   * general idea is that `se_contrasts()` is a wrapper around limma,
   so it should be able to print the equivalent calls to `limma::lmFit()`
   for example, so a user can observe the progression of analysis steps.
   * command to convert SE to expression matrix / and ExpressionSet used by limma
   * command to run optional `voom_jam()`; then the vanilla `limma::voom()`
   * command to run `run_limma_replicate()`
   * command to run `limma::lmFit()`, `limma::eBayes()`

* COMPLETE: previous update to `groups_to_sedesign()` introduced a regression (error)
for input that generates two-way contrasts.

## 20oct2022

* `heatmap_se()`

   * when `centerby_colnames=FALSE` no data centering is performed,
   the color legend should match the range of data, and hide negative
   values if there are no negative values.
   * allow custom `col` color function in form of `circlize::colorRamp2()`
   * deprecate `sample_color_list` into `color_list`; or consider adding
   `row_color_list` so row colors can contain the same colnames
   with different color assignments compared with `sample_color_list`. Hmmm.

## 13oct2022

* `se_contrast_stats()`

   * COMPLETE: argument `normgroup` to enforce independent statistical analyses
   within each unique normgroup.

## 21sep2022

* `se_contrast_stats()` to include additional gene annotation data

   * similarly, the first column currently hardcoded `"probes"`, should
   use the appropriate column name, or allow it to be defined by argument.
   * rowname is currently the only identifier included in results.
   * implementation ideas: it could either use input `rowData(se)`, in
   the limma model fit, or annotation can be added while each `stats_df`
   stat `data.frame` is created.

* `groups_to_sedesign()`

   * (COMPLETE) Mechanism to supply specific contrast names to be used
   in place of auto-generated contrasts.
   * (COMPLETE) Allow `SummarizedExperiment` input, with `group_colnames` used
   to define sample groups.

* `sestats` object output from `se_contrast_stats()`:

   * proper slot names: stats_df, stats_dfs, hit_array, hit_list, sedesign
   * proper print method that calls `sestats_to_df()` by default.
   * functions:
   
      * `hit_array(sestats)` with arguments assay_name, cutoff, contrasts

* `sestats_to_df()` bug:

   * apparently the colnames do not match the dimnames, `"cutoff"` is displayed
   for `"assay_name"` values. Probably something was reversed during recursive
   `list` navigation.

## 12sep2022

* `heatmap_se()` does not have option to customize arguments
`row_names_gp` nor `column_names_gp`, which could be used to
colorize, highlight, boldface, individual labels in the heatmap.

   * Currently `row_names_gp` is defined internally, in order to
   define `fontsize` based upon the number of rows and columns in
   the heatmap. The `fontsize` could be applied to user-defined
   `row_names_gp` or `column_names_gp`, however sometimes the
   rows and columns are defined dynamically - making it difficult
   to sync a vector of `grid::gpar(col=c("red", "black"))` values
   to the exact rownames.
   * To accommodate dynamic rows, it might need another method
   such as named vectors for known `grid::gpar()` attributes: `col`,
   `fontsize`, `fill`, `alpha`, `fontfamily`, `fontface`, `cex`, `font`
   (`font` is an alias for `fontface`, should be passed as `fontface`).

## 03aug2022

* `groups_to_sedesign()` implement normalization groups with these rules:

   * No one-way contrast, direct comparison, should be permitted which
   compares to different normalization groups.
   * Two-way contrasts are permitted across normalization groups, only
   those that involve one-way comparisons within a normalization group.

## 18jul2022

* Some form of power calculation that leverages the same methods used
by `se_contrast_stats()`. It could leverage the `ssizeRNA` package
for RNA-seq data.
* `contrast2comp()` is somehow very slow for even 10 to 20 contrast names.
* `contrast2comp()` adjustments to tolerate having label prefix.

   * for example `"fold (A_c-B_c)(A_d-B_d)"` could be recognized as
   `"label contrast"` and return `"fold A-B:c-d"`.
   * main assumption is no spaces in contrast name
   * may be useful to specify whether to look for label prefix or suffix,
   in the event a label contains a hyphen or dash `"-"` character.

## 11jul2022

* `heatmap_se()` font size customizations sufficient for manuscript prep.

## 27jun2022

* COMPLETE: `heatmap_se()` when not supplied `rows` nor `sestats` results in an error.
* COMPLETE: `heatmap_se()` should have option not to center data.
* COMPLETE: `heatmap_se()` should allow custom incidence matrix through
`sestats` and `alt_sestats`, to avoid having to provide `sestats` or
`hit_array`.

## 22jun2022

* COMPLETE: `heatmap_se()` - ability to "drill-down" into row clusters

## 21jun2022

* COMPLETE: Bug in `heatmap_se()` when `sestats` is supplied but `se` does not
contain all rows present in `sestats` hit array.
* COMPLETE: `heatmap_se()` needs control over the annotation name fontsize.

## 16jun2022

* `heatmap_se()`

   * COMPLETE: argument `isamples` should be useful to define samples
   to display in the heatmap.
   However, it would be useful to perform centering before
   subsetting by samples, in order to produce more useful graphs with
   paired data. For example, center data by each patient at time zero,
   then display the other timepoints (since the patient time zero would
   always be exactly zero, it only contributes a blank stripe to the
   heatmap).
   * COMPLETE: `top_annotation` and `left_annotation` shows color key
   for all colors, not just the those which are displayed in the heatmap
   annotation.
   
      * The potential downside is that it might affect results when two
      heatmaps are added together using `+` or `%v%` - in those cases
      equivalent color keys are sometimes merged together. Unsure how
      it would be handled when they are only partially identical.
      * UPDATE: It actually merges color keys by name, and displays
      unique colors for each name. That's so awesome.
   
   * COMPLETE: create new function to choose interesting annotation colnames,
   with logic that removes columns with 1:1 cardinality compared to
   other chosen columns. Columns that only repeat the information
   are no longer interesting.

## 24may2022

* `sestats_to_df()` consider making a wider output format intended
for `kable`, with relevant columns grouped by `assay_name`.
This way the output includes one row per contrast.
* COMPLETE: `sestats_to_df()` should report blank cell whenever NA values are
present, to convey that no cutoff was applied for that scenario,
rather than implying the cutoff was applied and there were no hits.

## 26apr2022

* `se_detected_rows()` was added, however

   * in future would be nice to apply constraints based upon contrasts.
   (In hindsight, I'm not sure what use case I had in mind!)

* COMPLETE: Migrate the `volcano_plot()` from `slicejam::volcano_plot()`

   * currently also draws block arrows in the plot margins
   * Migrated into `jamma` package, alongside MA-plots.


## 18apr2022

* New class `sestats` to replace output from `se_contrast_stats`

   * slotNames:
   
      * hit_array
      * hit_list
      * stat_dfs
      * stat_df
   
   * subset: [signal, contrast, cutoff]
   
      * intended to help filter for specific results
      * when supplying `sestats[matrix(ncol=3)]` it will subset for
      each element in the `hit_array`
   
   * accessors:
   
      * `summary()`, `print()` prints the `data.frame` summary of hit counts
      * `hits()` will return `list` of `list`
      * `hit_array()` will return the full `sestats@hit_array`
      * `hit_im()` will return an incidence matrix of hits
      
   * converters
   
      * `as.list()` to convert to previous `list` format
      * `list2sestats()` to convert from previous `list` format to `sestats`


* COMPLETE: Some method to rename two-way contrasts to save character space.

   * Notes:
   
      * Contrasts should be easily distinguished, since `se_contrast_stats()`
      encodes the contrast into column headers, for example
      
         * `"logFC CellA_Treated-CellA_Control"`
         * `"adj.P.Value CellA_Treated-CellA_Control"`
         * `"P.value CellA_Treated-CellA_Control"`
         * `"hit mgm5 adjP0.01 fc1.5 CellA_Treated-CellA_Control"`
         * Currently "easy" to find the contrast by taking the last
         non-whitespace string from the column header.
         * Contrasts should have no whitespace, use `":"` delimiter?
         
      * Renamed contrasts ideally do not use parentheses, since the goal
      is to reduce characters.
      
         * Technically `(CellA_Treated-CellA_Control)-(CellB_Treated-CellB_Control)`
         * Equivalent: `CellA_Treated-CellA_Control:CellB_Treated-CellB_Control`
         * Also equivalent: `Treated-Control:CellA-CellB`
   
   * Two-way contrast:
      * `(CellA_Treated-CellA_Control)-(CellB_Treated-CellB_Control)`
      (59 characters, 27 per contrast)

   * Alternative two-way syntax:
      * `Treated-Control x CellA-CellB`
      * `Treated-Control:CellA-CellB`

   * Two-way with one unchanging factor:
      * `(CellA_Treated_WT-CellA_Control_WT)-(CellB_Treated_WT-CellB_Control_WT)`
      (71 characters, 33 per contrast)

   * Alternative two-way with extra factor:
      * `CellA-CellB:Treated-Control:WT` (32 characters)
      * `Treated-Control CellA-CellB WT`

   * One-way contrast with extra factor:
      * `CellA_Treated_WT-CellA_Control_WT`
      (33 characters)
   
   * Alternative one-way with extra factor:
      * `Treated-Control:CellA:WT` (24 characters)
   
   
## 14apr2022

Simpler methods for common visualizations

* COMPLETE: expression heatmap - using `ComplexHeatmap::Heatmap()`
* Venn diagram of hits - using `venndir::venndir()`

   * design idea: use second order contrasts to determine which
   contrasts are "compatible" to be used in the same Venn diagrams.
   * contrast in cell type A; compared to same contrast in cell type B
   * one treatment-control in cell type A; another treatment-control cell type A
   * implied: it would not generally choose treatmentA-control cellA; treatmentB-control cellB
   * it could compare two-way contrasts using third-order contrast logic


* COMPLETE: volcano plots - migrate function from slicejam into `jamma` package.


## 13apr2022

* COMPLETE: `save_sestats()` is super slow, for 6 worksheets,
~15 columns, 25k rows, took about 5 minutes. Should be much faster.

   * Likely imposes changes to `jamba::writeOpenxlsx()`
   * Is conditional formatting is the slow step? If so, skip it.
   * Is categorical formatting is the slow step? If so, skip it.
   * Changes were made in `jamba::writeOpenxlsx()` to operate on
   an open Workbook without saving, passing the Workbook to each
   internal step also without saving. Each worksheet is added to
   the Workbook, and it is only saved at the end. Saved about 20x time
   especially for large multi-sheet Workbooks.

## 11apr2022

* COMPLETE: New function: `save_sestats()` or something similar.

   * Saves statistical results to Excel `.xlsx` file.
   * One contrast per worksheet.
   * Uses `jamba::writeOpenxlsx()` and defines each column type.
   * Also calls `jamba::set_xlsx_colwidths()` to set proper column widths.
   * TODO: Optionally saves the superwide `data.frame` stat table output.

* `se_contrast_stats()` enhancements

   * Consider object type `sestats` so it can have proper `print.sestats()`
   and `summary.sestats()` generic functions.
   * When using blocking factor with voom, follow guidance to
   apply voom/contrasts using the two-step process.
   * Some method to enforce `normgroup`
   
      * each group should be analyzed independently
      * results should be collated together into one `hit_array`

* `groups_to_sedesign()` needs some method to define first-order contrasts
distinct from second-order contrasts.

   * The second-order contrasts appear to use the wrong `factor_order`, e.g.
   
      * first order contrasts: `"A_Treat-A_Veh"` and `"B_Treat-B_Veh"`
      * second order contrast observed:`"(B_Treat-A_Treat)-(B_Veh-A_Veh)"`
      * second order contrast intended:`"(B_Treat-B_Veh)-(A_Treat-A_Veh)"`

* `se_contrast_stats()`, `se_normalize()`, `matrix_normalize()` need optional argument `normgroup`

   * COMPLETE: `se_normalize()` and `matrix_normalize()`
   * each subset of samples is processed independently within each `normgroup`.
   * intended to allow batch processing while also keeping distinct subsets
   of samples independent of one another.
   * basic workflow:
   
      * Apply method to entire data, optional, to create the full output matrix.
      * Iterate each `normgroup` subset, apply method, populate each result
      into the full output matrix.
      * In the case of batch adjustment, if `normgroup` subset does not work
      for batch adjustment (only one batch is represented) then copy data
      as-is into the output matrix.

* Consider some mechanism to hold blocking factor within `sedesign` object.


## contrast design plot

* x,y axes represent design factors, with factor levels at each position
* draw block arrow for each pairwise contrast, one group to another group
* Bonus points: use hit_array to plot hit counts inside each block arrow

* each block in the grid represents a group
* probably use fixed aspect=1 to ensure blocks are square

* if there is one factor, use it as x-axis, and use y-axis with no label
* if there are more than two factors, use them as a pseudo-z-axis with slight offset
* when two block arrows are co0linear, fan them out from the center

   * co-linear: identical slope (m) and intercept (b)
   * distance between parallel lines: abs(b1 - b2) / sqrt(1 + m^2)


## stats result object type?

* Goal would be to store output from `se_contrast_stats()` in
a proper S4 object:

   * stats_df
   * stats_dfs
   * hit_array
   * SEDesign (samples, design, contrasts used)

* Method to access stats data.frame by contrast
* Method to access statistical hits
* Method to review SEDesign used for the analysis
* Method to convert/aggregate probe-level hits to gene level, if relevant
* Method to use alongside SE data to visualize data
* Method to re-apply stat thresholds using stored results


## Add capability to `se_contrast_stats()` to call DESeq2

* wrap contrasts,design into DESeq2-compatible workflow
* store stat table output from DESeq2 similar to current method
* store relevant DESeq2-specific output


## Add capability to `se_contrast_stats()` to test isoforms/exons

* wraps `limma::diffSplice` or `DEXSeq` differential effects per gene



## Miscellaneous

* Vignette idea: Practical guide to creating gene expression
heatmaps using ComplexHeatmap.

