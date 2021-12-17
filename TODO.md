
# TODO for jamses

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

