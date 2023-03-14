
# TODO for jamses

## 12mar2023

* `heatmap_se()` improvements

   * Some method to hide `column_title` labels.
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

