
# TODO for jamses

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
   * Some method to enforce `normgroups`
   
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

