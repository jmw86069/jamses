
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Jam SummarizedExperiment Stats (jamses)

This package is under active development. These functions are in
frequent use in active Omics analysis projects. A summary of goals and
relevant features are described below.

[Jamses Online Documentation](https://jmw86069.github.io/jamses)

## Goals of jamses

The core goal is to make data analysis and visualization of
`SummarizedExperiment` objects straightforward for common scenarios. It
also accepts `SingleCellExperiment` and `Seurat` objects.

- **[Make Effective
  Heatmaps](https://jmw86069.github.io/jamses/articles/how_to_heatmap.html)**

- **Apply Normalization / Adjustment**

- `SEDesign`: **Create Design and Contrast matrices**

  - Use `~ 0 + group` syntax, see below.
  - Integrate Samples, Groups, and Contrasts.
  - Visualize with `plot_sedesign()`.

- `SEStats`: **Analyze Multiple Contrasts**

  - Use limma / limma-voom / limma-DEqMS.

- **Convenient contrast labels**

       Contrast:
       `(Knockout_treated-Knockout_control)-(Wildtype_treated-Wildtype_control)`

       Comp:
       `Knockout-Wildtype:treated-control`

- **Integrate with other tools**.

  - `SummarizedExperiment`, `SingleCellExperiment`, `Seurat`,
    `ExpressionSet`, `DESeqDataset`, `DGElist`
  - `venndir::venndir()` - see Github `"jmw86069/venndir"` to create
    directional Venn diagrams.

## Make Effective Heatmaps

[How to Make a
Heatmap](https://jmw86069.github.io/jamses/articles/how_to_heatmap.html)

Jamses `heatmap_se()` reinforces heatmap principles, and applies some
opinions.

- Data are centered.
- Always use divergent colors for centered data.
- Data are **not scaled**.
- Row and Column annotations are “easy”.

<details>

<summary>

## Approach for Statistical Contrasts

</summary>

[The Limma User’s
Guide](https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/)
(LUG) is an amazing resource which describes numerous approaches for
one-way and two-way contrasts **which are mathematically equivalent**.
For a more thorough discussion please review these approaches to confirm
that `~0 + x` is mathematically identical to `~ x`, and only differs in
how estimates are reported.

- Jamses uses the `~0 + x` strategy.

- Each experiment group is defined using independent replicates.<br>

- This approach *does not imply* that there is “no intercept” during the
  model fit, see LUG for details.

- One-way contrasts compare only one factor per contrast.

  - Valid one-way contrast:

  <!-- -->

      (A_treated - A_control) # valid one-way contrast

  - Invalid one-way contrast:

  <!-- -->

      (A_treated - B_control) # not a valid one-way contrast

- Two-way contrasts in Jamses compare the fold change of two
  ***compatible*** one-way fold changes.

  - Two compatible one-way contrasts:

  <!-- -->

      (A_treated - A_control)  # one-way contrast
                               # "treated-control" for A
      (B_treated - B_control)  # compatible one-way contrast
                               # "treated-control" for B

  - Corresponding two-way contrast:

  <!-- -->

      (B_treated - B_knockout) # incompatible one-way contrast

</details>

## SEDesign: Design and Contrasts

- `groups_to_sedesign()` takes by default a `data.frame` where each
  column represents an experiment factor, and creates the following:
- Output is `SEDesign` as an S4 object with slot names:
  - `design(sedesign)` - the design matrix.
  - `contrasts(sedesign)` - the contrasts matrix.
  - `samples(sedesign)` - vector of samples.

### Example SEDesign object

The example below uses a `character` vector of group names per sample,
with two factors separated by underscore `"_"`. The same data can be
provided as a `data.frame` with two columns.

``` r
library(jamses)

igroups <- jamba::nameVector(paste(rep(c("WT", "KO"), each=6),
   rep(c("Control", "Treated"), each=3),
   sep="_"),
   suffix="_rep");
igroups <- factor(igroups, levels=unique(igroups));

jamba::kable_coloring(color_cells=FALSE,
   caption="Sample to group association",
   data.frame(groups=igroups))
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<caption>

Sample to group association
</caption>

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:left;">

groups
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

WT_Control_rep1
</td>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

WT_Control
</td>

</tr>

<tr>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

WT_Control_rep2
</td>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

WT_Control
</td>

</tr>

<tr>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

WT_Control_rep3
</td>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

WT_Control
</td>

</tr>

<tr>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

WT_Treated_rep1
</td>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

WT_Treated
</td>

</tr>

<tr>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

WT_Treated_rep2
</td>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

WT_Treated
</td>

</tr>

<tr>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

WT_Treated_rep3
</td>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

WT_Treated
</td>

</tr>

<tr>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

KO_Control_rep1
</td>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

KO_Control
</td>

</tr>

<tr>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

KO_Control_rep2
</td>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

KO_Control
</td>

</tr>

<tr>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

KO_Control_rep3
</td>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

KO_Control
</td>

</tr>

<tr>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

KO_Treated_rep1
</td>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

KO_Treated
</td>

</tr>

<tr>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

KO_Treated_rep2
</td>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

KO_Treated
</td>

</tr>

<tr>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

KO_Treated_rep3
</td>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

KO_Treated
</td>

</tr>

</tbody>

</table>

The resulting design and contrasts matrices are shown below:

``` r
sedesign <- groups_to_sedesign(igroups);
jamba::kable_coloring(
   colorSub=c(`-1`="dodgerblue", `1`="firebrick"),
   caption="Design matrix output from design(sedesign).",
   data.frame(check.names=FALSE, design(sedesign)));
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<caption>

Design matrix output from design(sedesign).
</caption>

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:right;">

WT_Control
</th>

<th style="text-align:right;">

WT_Treated
</th>

<th style="text-align:right;">

KO_Control
</th>

<th style="text-align:right;">

KO_Treated
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

WT_Control_rep1
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(255, 255, 255, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: firebrick !important;">1</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

</tr>

<tr>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

WT_Control_rep2
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(255, 255, 255, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: firebrick !important;">1</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

</tr>

<tr>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

WT_Control_rep3
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(255, 255, 255, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: firebrick !important;">1</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

</tr>

<tr>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

WT_Treated_rep1
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(255, 255, 255, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: firebrick !important;">1</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

</tr>

<tr>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

WT_Treated_rep2
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(255, 255, 255, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: firebrick !important;">1</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

</tr>

<tr>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

WT_Treated_rep3
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(255, 255, 255, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: firebrick !important;">1</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

</tr>

<tr>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

KO_Control_rep1
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(255, 255, 255, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: firebrick !important;">1</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

</tr>

<tr>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

KO_Control_rep2
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(255, 255, 255, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: firebrick !important;">1</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

</tr>

<tr>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

KO_Control_rep3
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(255, 255, 255, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: firebrick !important;">1</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

</tr>

<tr>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

KO_Treated_rep1
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(255, 255, 255, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: firebrick !important;">1</span>
</td>

</tr>

<tr>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

KO_Treated_rep2
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(255, 255, 255, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: firebrick !important;">1</span>
</td>

</tr>

<tr>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

KO_Treated_rep3
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(255, 255, 255, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: firebrick !important;">1</span>
</td>

</tr>

</tbody>

</table>

``` r

jamba::kable_coloring(
   colorSub=c(`-1`="dodgerblue", `1`="firebrick"),
   caption="Contrast matrix output from contrasts(sedesign).",
   data.frame(check.names=FALSE, contrasts(sedesign)));
```

<table class="table" style="margin-left: auto; margin-right: auto;">

<caption>

Contrast matrix output from contrasts(sedesign).
</caption>

<thead>

<tr>

<th style="text-align:left;">

</th>

<th style="text-align:right;">

KO_Control-WT_Control
</th>

<th style="text-align:right;">

KO_Treated-WT_Treated
</th>

<th style="text-align:right;">

WT_Treated-WT_Control
</th>

<th style="text-align:right;">

KO_Treated-KO_Control
</th>

<th style="text-align:right;">

(KO_Treated-WT_Treated)-(KO_Control-WT_Control)
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

WT_Control
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(255, 255, 255, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: dodgerblue !important;">-1</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(255, 255, 255, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: dodgerblue !important;">-1</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(255, 255, 255, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: firebrick !important;">1</span>
</td>

</tr>

<tr>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

WT_Treated
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(255, 255, 255, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: dodgerblue !important;">-1</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(255, 255, 255, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: firebrick !important;">1</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(255, 255, 255, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: dodgerblue !important;">-1</span>
</td>

</tr>

<tr>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

KO_Control
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(255, 255, 255, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: firebrick !important;">1</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(255, 255, 255, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: dodgerblue !important;">-1</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(255, 255, 255, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: dodgerblue !important;">-1</span>
</td>

</tr>

<tr>

<td style="text-align:left;border-left:1px solid #DDDDDD;white-space: nowrap;">

KO_Treated
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(255, 255, 255, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: firebrick !important;">1</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(0, 0, 0, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: transparent !important;">0</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(255, 255, 255, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: firebrick !important;">1</span>
</td>

<td style="text-align:right;border-left:1px solid #DDDDDD;white-space: nowrap;">

<span style="     color: rgba(255, 255, 255, 255) !important;border-radius: 4px; padding-right: 4px; padding-left: 4px; background-color: firebrick !important;">1</span>
</td>

</tr>

</tbody>

</table>

For convenience, SEDesign can be visualized using `plot_sedesign()`:

``` r
# plot the design and contrasts
plot_sedesign(sedesign);
title(main="plot_sedesign(sedesign):")
```

![](man/figures/README-plot_sedesign-1.png)<!-- -->

- One-way contrasts are shown with a wide block arrow.

- Two-way contrasts are shown by connecting two block arrows with a
  “squiggly curved line”.

  - It connects the end of one contrast\
    to the beginning of the next contrast.
  - The order indicates that the first contrast is subtracted by the
    second contrast.
