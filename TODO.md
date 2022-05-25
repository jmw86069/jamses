
# TODO for jamses

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

