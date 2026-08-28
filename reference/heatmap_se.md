# Heatmap for SummarizedExperiment data

Heatmap for SummarizedExperiment data

## Usage

``` r
heatmap_se(
  se,
  sestats = NULL,
  hm_name = NULL,
  hm_title = NULL,
  rows = NULL,
  row_type = "rows",
  column_type = "samples",
  data_type = "expression",
  correlation = FALSE,
  assay_name = NULL,
  contrast_names = NULL,
  contrast_suffix = "",
  cutoff_name = NULL,
  alt_sestats = NULL,
  alt_assay_name = assay_name,
  alt_contrast_names = NULL,
  alt_contrast_suffix = "",
  alt_cutoff_name = NULL,
  isamples = colnames(se),
  normgroup_colname = NULL,
  centerby_colnames = NULL,
  controlSamples = NULL,
  control_label = "",
  noise_floor = NA,
  noise_floor_value = noise_floor,
  transformation = NULL,
  controlFloor = NA,
  naControlAction = c("na", "row", "floor", "min"),
  naControlFloor = 0,
  top_colnames = NULL,
  top_annotation = NULL,
  top_annotation_name_gp = grid::gpar(),
  rowData_colnames = NULL,
  rowData_side = c("left", "right"),
  left_annotation = NULL,
  left_annotation_name_gp = grid::gpar(),
  left_annotation_name_rot = 90,
  right_annotation = NULL,
  simple_anno_size = grid::unit(8, "mm"),
  legend_title_gp = grid::gpar(fontsize = 10),
  legend_labels_gp = grid::gpar(fontsize = 10),
  legend_grid_cex = 1,
  row_names_gp = NULL,
  row_split = NULL,
  row_subcluster = NULL,
  row_title_rot = 0,
  sample_color_list = NULL,
  legend_at = NULL,
  legend_labels = NULL,
  subset_legend_colors = TRUE,
  row_cex = 0.8,
  column_cex = 1,
  row_anno_fontsize = 11,
  useMedian = FALSE,
  show_row_names = NULL,
  show_row_dend = length(rows) < 2000,
  mark_rows = NULL,
  mark_labels_gp = grid::gpar(),
  column_title = character(0),
  apply_hm_column_title = FALSE,
  hm_title_buffer = 0,
  show_heatmap_legend = TRUE,
  show_top_legend = TRUE,
  show_left_legend = TRUE,
  legend_border_color = "black",
  show_top_annotation_name = c("right", "left", "none"),
  show_left_annotation_name = c("bottom", "top", "none"),
  row_label_colname = NULL,
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  cluster_rows = function(x, ...) {
     amap::hcluster(jamba::rmNA(naValue = 0, x), ...,
    method = "euclidean", link = "ward")
 },
  cluster_row_slices = FALSE,
  column_names_gp = NULL,
  column_split = NULL,
  column_split_sep = ",",
  color_max = 3,
  color_floor = 0,
  lens = 2,
  rename_contrasts = TRUE,
  rename_alt_contrasts = TRUE,
  use_raster = TRUE,
  verbose = FALSE,
  debug = FALSE,
  ...
)
```

## Arguments

- se:

  `SummarizedExperiment` by default, or one of the following:

  - `SummarizedExperiment` with accessor functions `rowData()`,
    `colData()`, and `assays()`. It will use `values(rowRanges())` if no
    slot `rowData` exists.

  - `SingleCellExperiment` with accessor functions `rowData()`,
    `colData()`, and `assays()`. It will use `values(rowRanges())` if no
    slot `rowData` exists.

  - `Seurat` object, which is coerced to `SingleCellExperiment` and
    handled accordingly

  - `ExpressionSet` or compatible object with accessor functions
    `featureData()`, `phenoData()`, and `assayData()`.

- sestats:

  one of the following types of data:

  - `list` output from
    [`se_contrast_stats()`](https://jmw86069.github.io/jamses/reference/se_contrast_stats.md),
    which specifically contains `hit_array` as a 3-dimensional array of
    hits with dimensions "Cutoffs", "Contrasts", "Signal".

  - `numeric` matrix intended to represent an incidence matrix, where a
    value `0` indicates absence, and non-zero indicates presence. This
    format is useful for supplying any incidence matrix, such as
    gene-by-pathway (for example Github package
    "jmw86069/multienrichjam" provides `mem$memIM` with gene-by-pathway
    matrix), or gene-by-class (see Github package "jmw86069/pajam" for
    examples using ProteinAtlas protein classification, including
    membrane-bound, secreted, transcription factors, etc.), or any
    incidence matrix defined by Github "jmw86069/venndir" function
    `list2im_value()` or `list2im()` which converts input to a Venn
    diagram into an incidence matrix.

  - When `sestats` is supplied, data is converted to incidence matrix,
    then columns are matched with `contrast_names`. All rows with
    non-zero entry in those columns are included in the heatmap. When
    `rows` is also supplied, then the intersection of incidence matrix
    rows and `rows` is displayed in the heatmap.

  - Note that `alt_sestats` does not subset rows displayed in the
    heatmap.

- hm_name:

  `character` string, or `NULL` (default) which uses the `data_type`
  value. Note that the legend title uses the `data_type`, and is also
  used for `hm_name` when `hm_name=NULL`. The `hm_name` is most useful
  to customize because this string is used as the prefix for grid
  graphical components, for example seen with
  [`ComplexHeatmap::list_components()`](https://rdrr.io/pkg/ComplexHeatmap/man/list_components.html).
  When two heatmaps or a `HeatmapList` is drawn, the names can be used
  to define specific grid regions of each heatmap. If the heatmaps share
  the same `hm_name` then the regions will also have identical name and
  cannot be addressed distinctly.

- hm_title:

  `character` string, or `NULL` (default) which generates a heatmap
  title using the dimensions, `assay_name`, `data_type`, and a string
  which describes the data centering. When provided as a `character`
  string, it is used as-is. (In future this value may accept variable
  names.)

- rows:

  `character` vector of `rownames(se)` to define a specific set of rows
  to display. When `sestats` is supplied, then the intersection of
  `rows` with genes defined by `sestats` is displayed. Note that rows
  are required to be in `rownames(se)`, all other rows are dropped.

- row_type:

  `character` string used in the title of the heatmap which indicates
  how many rows are displayed. For example
  `"1,234 genes detected above background"` or
  `"1,234 DEGs by limma-voom"`. When `row_type=""` or `row_type=NULL`
  this information is not included in the heatmap title.

- column_type:

  `character` string used in the title of the heatmap which indicates
  how many column are displayed. For example `"12 samples"` or
  `"12 biological replicates"`. When `column_type=""` or
  `column_type=NULL` this information is not included in the heatmap
  title.

- data_type:

  `character` string used as title of the heatmap color gradient legend,
  for example `"expression"` indicates the data contains gene expression
  measurements. Notes:

  - The prefix `"centered"` is automatically appended whenever the data
    is also centered for the heatmap. Set `centerby_colnames=FALSE` to
    display data that is not centered.

  - The prefix `"correlation of"` is automatically appended when
    `correlation=TRUE` which displays correlation of whatever data is
    included in the heatmap.

- correlation:

  `logical` indicating whether to calculate sample correlation, and plot
  a sample-by-sample correlation heatmap. This option is included here
  since many of the same arguments are required for data centering, and
  sample annotations. Note that `color_max` is forced to a maximum value
  of `1.0`, representing the maximum correlation value.

- assay_name:

  `character` string indicating the name in `assays(se)` to use for data
  to be displayed in the heatmap.

  - When multiple `assay_name` values are supplied, the first assay_name
    that matches `names(assays(se))` will be used in the heatmap. In
    this way, multiple `assay_names` can be supplied to define
    statistical hits in `sestats`, which calls
    [`hit_array_to_list()`](https://jmw86069.github.io/jamses/reference/hit_array_to_list.md)
    to combine hits across `assay_name` entries; but only the first
    `assay_name` found in `se` is used for the heatmap values.

  - When there is only one value for `assayNames(se)`, then `assay_name`
    will default to this value, instead of acting like it couldn't
    possibly know what was intended. Haha.

  - Lastly, `assay_name` can be a `numeric` index, helpful in case
    `assays(se)` contains no names - not recommended but it can happen.

- contrast_names:

  `character` vector of contrasts in `sestats$hit_array` to use for the
  heatmap. When `contrast_names=NULL` then all contrasts are displayed,
  which is the default.

- contrast_suffix:

  `character` string with optional suffix to append to the end of each
  contrast name label for `sestats` hit incidence matrix beside the
  heatmap. This suffix may be useful when comparing two methods for the
  same set of contrast names, with `sestats` and `alt_sestats`.

- cutoff_name:

  `character` or `integer` index used to define the specific statistical
  cutoffs to use from `sestats$hit_array`. This argument is passed to
  [`hit_array_to_list()`](https://jmw86069.github.io/jamses/reference/hit_array_to_list.md)
  as `cutoff_names`.

- alt_sestats, alt_assay_name, alt_contrast_names, alt_contrast_suffix,
  alt_cutoff_name:

  arguments analogous to those described above for `sestats` which are
  used when `alt_sestats` is supplied.

- isamples:

  `character` vector of `colnames(se)` used to visualize a subset of
  samples used for the data centering step. Note that data centering
  uses all columns supplied in `se`, and after centering, the subset of
  columns defined in `isamples` is displayed in the heatmap. This
  distinction makes it possible to center data by some control group,
  then optionally not display the control group data.

- normgroup_colname:

  `character` vector of colnames in `colData(se)` used during data
  centering. When supplied, samples are centered independently within
  each normgroup grouping. These values are equivalent to using
  `centerby_colnames`.

- centerby_colnames:

  either:

  - `character` vector of colnames in `colData(se)` used during data
    centering. When supplied, samples are centered independently within
    each centerby grouping. It is typically used for things like cell
    lines, to center each cell line by a time point control, or
    untreated control.

  - `NULL` to perform centering across all columns in `se`.

  - `FALSE` to disable centering.

- controlSamples:

  `character` optional vector of samples to use as the reference during
  data centering. Note that samples are still centered within each
  normgroup and centerby grouping, and within that grouping samples are
  centered to the `controlSamples` which are present in that grouping.
  Any center group for which no samples are defined in `controlSamples`
  will use all samples in that center group. Typically, `controlSamples`
  is used to define a specific group as the reference for centering, so
  changes are displayed relative to that group. Make sure to define
  `control_name` to include an appropriate label in the heatmap title.

- control_label:

  `character` string used in heatmap title to describe the control used
  during data centering, relevant when `controlSamples` is also
  supplied. Recommended format: `"versus Wildtype"` or `"vs. Wildtype"`.
  The heatmap title will include data centering and `control_label` in
  this format: `"centered within {centerby_colnames}, {control_label}"`,
  for example `"centered within Genotype/Time, versus Vehicle"`.

- noise_floor, noise_floor_value:

  `numeric` values used to adjust data *before* data centering. It is
  most useful when changing `0` to `NA`, which is appropriate when zero
  `0` represents 'not measured'.

  - `noise_floor`: threshold at or below which all `numeric` values are
    set to `noise_floor_value`.

    - For example, `noise_floor=0` and `noise_floor_value=NA` will set
      any value at or below `0` to `NA`. This strategy is appropriate
      when values of zero `0` are considered 'not measured', and
      therefore should not be considered when centering data. Values of
      `NA` will be colored using the argument 'na_col' whose default is
      `na_col="grey"` in
      [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html).

    - Another option is to set `noise_floor=0` or some technical noise
      floor pertinent to the platform, then
      `noise_floor_value=noise_floor`. This approach sets any value at
      or below the noise floor to the noise floor. The approach
      therefore assumes any measurement at or below the effective limit
      of detection is at the limit of detection, since no measurement
      below that value is trustworthy. Some platforms (e.g. some tandem
      mass spectrometry proteomics) produces extremely high intensities
      when a peptide is detected, and zero otherwise. However, the
      effective noise floor is also high, even in log2 space it could be
      log2 of 20. Any measurement below log2 of 20 is considered
      equivalent to 20, since the workflow is not capable of measuring
      15 for example.

  - `noise_floor_value`: single value used as a replacement when any
    value meets the criteria, at or below `noise_floor`.

    - The default is to replace with `noise_floor`.

    - Another option is `NA` so the values do not affect data centering.

- transformation:

  `function` default NULL, optional `numeric` transformation applied
  *after* applying `noise_floor`, and *before* data is centered.

  - An example might be `transformation = jamba::log2signed()` which
    applies `log2(1 + x)` transformation, while also preserving negative
    values if they exist.

- controlFloor, naControlAction, naControlFloor:

  passed to
  [`jamma::centerGeneData()`](https://rdrr.io/pkg/jamma/man/centerGeneData.html)
  to customize control group values during data centering:

  - `controlFloor` imposes an optional noise floor to control group
    mean/median values, so the summary value during centering is at
    least `controlFloor`. Useful for defining an effective noise floor
    for a platform technology. Use `NA` to apply no control floor.

  - `naControlAction` defines the action taken only when values for all
    control samples are `NA`:

    - "na": all values become NA, since the control group is NA.

    - "row": non-NA values are centered to the row non-NA values.

    - "floor": non-NA values are centered versus `naControlFloor`

    - "min": non-NA values are centered to non-NA min across all rows

  - `naControlFloor` is a `numeric` value used when
    `naControlAction="floor"`, which causes the group reference value to
    use the value provided in `naControlFloor`.

- top_colnames:

  one of the following types:

  - `character` vector of colnames to use from `colData(se)` as
    annotations to display in `top_annotation` above the heatmap.

  - `NULL`, will call
    [`choose_annotation_colnames()`](https://jmw86069.github.io/jamses/reference/choose_annotation_colnames.md)
    to detect reasonable colnames: columns with more than one unique
    value; columns with at least one duplicated value.

  - `FALSE` will hide the `top_colnames`, which also occurs when
    `colData(se)` is empty.

- top_annotation:

  specific heatmap annotation as defined by
  [`ComplexHeatmap::HeatmapAnnotation()`](https://rdrr.io/pkg/ComplexHeatmap/man/HeatmapAnnotation.html).
  When supplied, the `top_colnames` described above is not used.

- top_annotation_name_gp:

  [`grid::gpar`](https://rdrr.io/r/grid/gpar.html) object to customize
  the annotation name displayed beside the top annotation.

- rowData_colnames:

  `character` vector of colnames in `rowData(se)` to use for heatmap
  annotations displayed on the left side of the heatmap. Specific colors
  can be included in `sample_color_list` as a named `list` of color
  vectors or color functions. The names of this list must match colnames
  to be displayed, otherwise
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
  will define its own color function.

- rowData_side:

  `character` string with the side to place the annotations defined by
  `rowData_colnames`. Default 'left', although 'right' is available as
  of 0.0.74.900.

  - When `rowData_side='right'` the annotation data are prepared as
    usual, and assigned to `right_annotation`.

  - When `rowData_side='right'` and either `right_annotation` or
    `mark_rows` are enabled, the right annotation is therefore added to
    combine them, resulting in a `HeatmapList` output, thereby also
    causing the output of this function to be `HeatmapList`. It can be
    drawn as usual with `ComplexHeatmap::draw())`, however it may be
    difficult to combine heatmaps vertically using `%v%`.

- left_annotation:

  specific heatmap annotation as defined by
  [`ComplexHeatmap::rowAnnotation()`](https://rdrr.io/pkg/ComplexHeatmap/man/rowAnnotation.html).
  When supplied, the `rowData_colnames` and `sestats` row annotations
  are not displayed. In order to supply custom row annotations and not
  lose `left_annotation` defined above, supply the row annotations as
  `right_annotation`.

- left_annotation_name_gp:

  [`grid::gpar`](https://rdrr.io/r/grid/gpar.html) object to customize
  the annotation name displayed beside the left annotation.

- left_annotation_name_rot:

  `numeric` rotation of left annotation label, in degrees, where `0`
  indicates normal text, and `90` is rotated vertically.

- right_annotation:

  specific heatmap annotation as defined by
  [`ComplexHeatmap::HeatmapAnnotation()`](https://rdrr.io/pkg/ComplexHeatmap/man/HeatmapAnnotation.html).

  - This element is created automatically when `mark_rows` is supplied.

  - Note: When `right_annotation` and `mark_rows` are both defined, the
    annotations are added, which produces `HeatmapList` instead of
    `HeatmapAnnotation`, and as a result, it is also added to the
    heatmap for correct results, causing the output to become
    `HeatmapList` instead of `Heatmap`. It can still be drawn, but
    cannot easily be combined vertically using `%v%` for example.

- simple_anno_size:

  [`grid::unit`](https://rdrr.io/r/grid/unit.html) size used to define
  heatmap annotation sizes (height or width of each line) for any simple
  annotations.

- legend_title_gp:

  [`grid::gpar`](https://rdrr.io/r/grid/gpar.html) to customize the
  legend title fonts, applied to each legend: top annotation, left
  annotation, main heatmap.

- legend_labels_gp:

  [`grid::gpar`](https://rdrr.io/r/grid/gpar.html) to customize the
  legend label fonts, applied to each legend: top annotation, left
  annotation, main heatmap.

- legend_grid_cex:

  `numeric` multiplied to adjust the relative size of each legend grid
  unit, applied to each relevant metric.

- row_names_gp:

  `gpar` to define custom column name settings. When `"fontsize"` is not
  defined, the automatic font size calculation is added to the
  `row_names_gp` supplied.

- row_split:

  is used to define heatmap split by row, ultimately passed to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
  argument `row_split`. However, the input type can vary:

  - `integer` number of row splits based upon row clustering. If
    `row_split` is greater than the number of rows, it will be set to
    the number of rows.

  - `character` value or values in colnames of `rowData(se)` to split
    using row annotation in `se`.

  - `data.frame` whose
    [`rownames()`](https://rdrr.io/r/base/colnames.html) must contain
    all rows to be displayed in the heatmap. This argument is passed
    directly to
    [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
    to apply the split appropriately.

  - `character` or `factor` vector named by `rownames(se)` with another
    custom row split, passed directly to
    [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
    argument `row_split`, with proper order for rows being displayed

- row_subcluster:

  `integer` or `character` vector representing one or more elements
  returned by `row_split` to use as a drill-down sub-cluster heatmap.
  This argument is experimental, and is intended to make it easy to
  "drill down" into specific row clusters.

  - The process internally creates a full heatmap using all arguments as
    defined, then extracts the
    [`jamba::heatmap_row_order()`](https://jmw86069.github.io/jamba/reference/heatmap_row_order.html)
    which contains row split data in a `list` of rownames vectors. The
    `list` elements that match `row_subcluster` are extracted and used
    again for a subsequent heatmap, and are displayed in the same order
    in which they appear in the original full heatmap - which means
    `cluster_rows=FALSE` is defined at this point. However `row_split`
    is retained for this subset of rows, to indicate the original row
    split annotation.

  - Note that `row_subcluster` must match the
    [`names()`](https://rdrr.io/r/base/names.html) returned by
    [`jamba::heatmap_row_order()`](https://jmw86069.github.io/jamba/reference/heatmap_row_order.html)
    for the full heatmap, or should include a `numeric` index for the
    `list` element or elements to use.

  - In principle this process would be run in two stages: First, view a
    heatmap with `row_split=6`, then re-run the same heatmap with
    `row_subcluster=4` to see cluster number 4 from the full heatmap.

- row_title_rot:

  `numeric` value indicating text rotation in degrees to use for row
  titles.

- sample_color_list:

  named `list` of color vectors or color functions, where names
  correspond to colnames in either `colData(se)` or `rowData(se)`, and
  which are passed to corresponding left or top annotation functions.
  When colors are not defined,
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
  will define colors using its own internal function.

- legend_at, legend_labels:

  `numeric` and `character`, respectively, to define custom values for
  the heatmap color gradient legend.

  - When `legend_at` is supplied, it is used as provided.

  - When `legend_labels` is supplied, it is used only when its length
    equals `length(legend_at)`, in which case it is used as provided.

  - When `centerby_colnames=FALSE` and the matrix data does not contain
    negative values, `legend_at` uses integers from `0` to `color_max`,
    to avoid presenting a color legend with unnecessary negative values.
    However, when `color_max <= 1` it uses `pretty(c(0, color_max))`,
    removing extraneous values, then ensuring the maximum value is
    `color_max`. For example when `color_max=0.85`, the `legend_at` is
    likely to be `c(0, 0.2, 0.4, 0.6, 0.8, 0.85)`.

  - When `centerby_colnames` is not `FALSE`, and/or data contains
    negative values, the `legend_at` is symmetric above and below zero.
    When `color_max <= 1` the label is created using
    `pretty(c(-color_max, color_max))`, as described above, so
    `color_max` is used as the minimum and maximum value. When
    `color_max > 1` the `legend_at` uses integer steps.

  - When `color_max <= 1` the `legend_labels` are presented as-is with
    no transformation.

  - When `color_max > 1` the `legend_labels` are transformed with
    `exp2signed(x)` which is the inverse of `log2(1 + x)`. This inverse
    tranform displays normal space values, in the case of centered data,
    the values represent normal space fold changes. For example the
    `legend_at=c(-2, -1, 0, 1, 2)` would result in
    `legend_labels=c("-4", "-2", "1", "2", "4")`.

  - When `correlation=TRUE` the `legend_labels` by default use
    `legend_at`, following rules for `color_max <= 1` above. Otherwise,
    `legend_labels` values inverse transformed from `log2(1 + x)` in
    order to display normal space fold change values,

  - To override any of this behavior, supply both `legend_at` and
    corresponding `legend_labels`.

- subset_legend_colors:

  `logical` indicating whether to subset colors shown in the color key
  defined by `sample_color_list`, which is useful when the heatmap only
  represents a subset of categorical color values.

  - When `subset_legend_colors == TRUE`, the color key will only include
    colors shown in the `top_annotation`.

  - When `subset_legend_colors == FALSE` all colors defined in
    `sample_color_list` will be included for each relevant column.

- row_cex, column_cex:

  `numeric` values used to adjust the row and column name font size,
  relative to the automatic adjustment that is already done based upon
  the number of rows and columns being displayed.

- row_anno_fontsize:

  `numeric` base font size for row annotation labels. This value is only
  used when `left_annotation_name_gp` is not supplied. Note these labels
  appears underneath row annotations, alongside column labels, and
  therefore they are also adjusted by multiplying `column_cex` so these
  labels are adjusted together.

- useMedian:

  `logical` passed to
  [`jamma::centerGeneData()`](https://rdrr.io/pkg/jamma/man/centerGeneData.html)
  during data centering.

- show_row_names, show_row_dend:

  `logical` indicating whether to display row names, and row dendrogram,
  respectively. With more than 2,000 rows this step can become somewhat
  slow.

- mark_rows:

  `character` vector of values in `rownames(se)` that should be labeled
  using
  [`ComplexHeatmap::anno_mark()`](https://rdrr.io/pkg/ComplexHeatmap/man/anno_mark.html)
  in call-out style.

  - Usually this argument is used when `show_row_names=FALSE`, hiding
    the row labels, but is not required.

  - Values in `mark_rows` are intersected with rows included in the
    heatmap, therefore only matching entries will be labeled.

  - Note when `mark_rows` and another of `right_annotation` or
    `rowData_side='right'` is active, the mark labels are always placed
    nearest the heatmap, consistent with ComplexHeatmap convention.
    Also, in this case the output of this function becomes `HeatmapList`
    instead of `Heatmap`, due to how ComplexHeatmap adds multiple
    annotations together.

- mark_labels_gp:

  [`grid::gpar`](https://rdrr.io/r/grid/gpar.html) to customize the font
  used by labels when `mark_rows` is supplied.

- column_title:

  `character` optional title to include at the top of the heatmap. It
  can include a single value, or multiple values representing each
  `column_split` in the order they appear.

  - Note: This argument is ignored when `apply_hm_column_title=TRUE`.

  - When `column_title=character(0)` (default) or `column_title=""`, the
    [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
    uses its usual default behavior, which is to assign `column_title`
    using `column_split` values when they are being used.

- apply_hm_column_title:

  `logical` (default FALSE) whether to apply the heatmap title to
  `column_title`. This option makes it convenient to display the title
  atop the heatmap without additional effort, however it hides any other
  `column_title` created by using `column_split`. When using both
  `column_split` and `apply_hm_column_title=TRUE` it may be useful to
  call
  [`heatmap_column_group_labels()`](https://jmw86069.github.io/jamses/reference/heatmap_column_group_labels.md).

- hm_title_buffer:

  `numeric` number of whitespace lines to add to the heatmap title
  `(attr(hm, "hm_title")` between the title and the heatmap below it.
  This whitespace can be useful when also calling
  [`heatmap_column_group_labels()`](https://jmw86069.github.io/jamses/reference/heatmap_column_group_labels.md),
  to provide enough space to draw the additional annotations.

- show_heatmap_legend, show_left_legend, show_top_legend:

  `logical` indicating whether each legend should be displayed.
  Sometimes there are too many annotations, and the color legends can
  overwhelm the figure. Note that `show_left_legend` is applied in a
  specific order, with these rules:

  - `show_left_legend` is extended to at least length 2, then values are
    used in order for: `sestats`, `rowData_colnames`, in order, using
    whichever is defined.

  - If `sestats` is defined, the first value in `show_left_legend` is
    used for this annotation, then the remaining values are used for
    `rowData_colnames`. Setting the first `show_left_legend` value to
    `FALSE` will ensure the legend for `sestats` is not displayed.

  - If `rowData_colnames` is defined, then the remaining values in
    `show_left_legend` are recycled for all columns in
    `rowData_colnames`, and applied in order. In this way, individual
    columns can have the legend displayed or hidden.

  - If `alt_sestats` is defined, the legend is always hidden, in favor
    of showing only the legend for `sestats` without duplicating this
    legend.

- legend_border_color:

  `character` color used as border color tofor be used as a border color
  for the various legend colors. Note this argument recognizes only the
  first color provided, and does not recycle different colors across the
  various legend borders.

- show_top_annotation_name, show_left_annotation_name:

  `character`, default 'right' for top annotation, 'bottom' for left
  annotation. Indicating which side to place the label, or 'none' to
  hide the label.

  - For compatibility with previous work, and ComplexHeatmap in general,
    `logical` is accepted, where TRUE uses the default side, and FALSE
    hides the label.

- row_label_colname:

  `character` string used as a row label, where this value is a colname
  in `rowData(se)`. It is useful when rownames are some identifier that
  is not user-friendly, and where another column in the data may provide
  a more helpful label, for example `"SYMBOL"` to display gene symbol
  instead of accession number.

- cluster_columns, cluster_rows:

  `logical`, `function`, `hclust` or `dendrogram`. Anything except
  `FALSE` will cluster columns or rows, respectively. Example input by
  type:

  - `TRUE` uses 'Jam' default, which replaces NA with zero `0`, then
    applies: `amap::hcluster(x, method="euclidean", link="ward")`.

  - `FALSE` performs no clustering.

  - `function` must take a `numeric` matrix input, and return either
    `hclust` or `dendrogram`.

  - `hclust` or `dendrogram` input must exactly match the columns or
    rows, as appropriate. It does not support use of `column_km` or
    `row_km` alongside clustering, however these arguments can be used
    alongside `function`.

- cluster_column_slices, cluster_row_slices:

  `logical` default FALSE, whether to cluster columns and rows,
  respectively, when also using `column_split` or `row_split`.

  Note that `TRUE` implies the respective `cluster_columns` or
  `cluster_rows` is TRUE or is a `function` that returns `hclust` or
  `dendrogram` output.

- column_names_gp:

  `gpar` to define custom column name settings. When `"fontsize"` is not
  defined, the automatic font size calculation is added to the
  `column_names_gp` supplied.

- column_split:

  `character` or `integer` vector used to define heatmap column split.

- column_split_sep:

  `character` string used as delimited when `column_split` defines
  multiple split levels.

- color_max:

  `numeric` value passed to
  [`colorjam::col_div_xf()`](https://jmw86069.github.io/colorjam/reference/col_div_xf.html)
  which defines the upper limit of color gradient used in the heatmap.

- color_floor:

  `numeric` value passed to
  [`colorjam::col_div_xf()`](https://jmw86069.github.io/colorjam/reference/col_div_xf.html)
  argument `floor` which defines the minimum non-zero numeric value for
  a color to be applied. This option is available to prevent coloring
  values below the `color_floor` which can be useful in some
  circumstances.

- lens:

  `numeric` value passed to
  [`colorjam::col_div_xf()`](https://jmw86069.github.io/colorjam/reference/col_div_xf.html)
  to control the intensity of color gradient applied to the numeric
  range.

- rename_contrasts, rename_alt_contrasts:

  `logical` (default TRUE) whether to rename long contrast names in
  `sestats` and `alt_sestats` using
  [`contrast2comp()`](https://jmw86069.github.io/jamses/reference/contrast2comp.md).

- use_raster:

  `logical` passed to
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
  to determine whether heatmaps should be converted to raster images,
  which effectively turns each heatmap panel into a single graphical
  object. Recommend `use_raster=TRUE` and also installing R package
  `magick` which greatly enhances speed and quality of rasterized
  heatmap output. When `magick` is not available, it may be best to use
  `use_raster=FALSE`. When `use_raster=FALSE` each pixel square of a
  heatmap is its own graphical object. For heatmaps with very large
  dimensions, having each pixel as an object can make the heatmap
  extremely large in memory, and sometimes pixels can overlap others
  because the minimum pixel size of the output graphics device does not
  reflect the actual size of each pixel.

- verbose:

  `logical` indicating whether to print verbose output.

- debug:

  `logical` indicating debug mode, data is returned in a `list`:

  - `hm` object
    [`ComplexHeatmap::Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)

  - `top_annotation` object
    [`ComplexHeatmap::HeatmapAnnotation`](https://rdrr.io/pkg/ComplexHeatmap/man/HeatmapAnnotation.html)
    for columns

  - `left_annotation` object
    [`ComplexHeatmap::HeatmapAnnotation`](https://rdrr.io/pkg/ComplexHeatmap/man/HeatmapAnnotation.html)
    for rows

  - `hm_title` object `character` string with the heatmap title.

- ...:

  additional arguments are passed to supporting functions.

## Value

`Heatmap` object, capable of being drawn using
[`ComplexHeatmap::draw()`](https://rdrr.io/pkg/ComplexHeatmap/man/draw-dispatch.html).

- In special scenarios when `right_annotation` is provided together with
  `mark_rows` or `rowData_side='right'`, the output is `HeatmapList`.

- In either case, the row and column order can be returned using
  [`jamba::heatmap_row_order()`](https://jmw86069.github.io/jamba/reference/heatmap_row_order.html)
  and
  [`jamba::heatmap_column_order()`](https://jmw86069.github.io/jamba/reference/heatmap_column_order.html),
  respectively.

- When `row_label_colname` is provided, the output of
  [`jamba::heatmap_row_order()`](https://jmw86069.github.io/jamba/reference/heatmap_row_order.html)
  will include named `character` vectors. The names will include the
  labels as displayed, while the values will match the data being
  plotted. Typically, data are available using `hm@matrix` for `Heatmap`
  objects, and `hm@ht_list[[1]]@matrix` for `HeatmapList` objects, where
  the first heatmap in the list is typically the main heatmap, with
  additional heatmaps on the right indicating the right_annotation data
  where applicable.

## Details

Note: Still a work in progress. This function is the basis for the
majority of heatmaps created for Omics data.

This function is a bold attempt to simplify the intricate task of
creating an expression heatmap, using
[`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html),
given a `SummarizedExperiment` object.

It attempts to enable:

- selection of `assays(se)` to use in the heatmap

- use of `rowData(se)` or `colData(se)` to produce row and column
  annotations, respectively.

- re-use of defined colors for annotations, see
  `platjam::design2colors()`

- define and adjust heatmap color gradient and scale

- data centering by row: versus all columns, or specific controls,
  optionally within independent centering groups

- filtering rows to show only the statistical hits

- display annotation of statistical hits beside the heatmap

- split rows or columns using `rowData(se)` and `colData(se)`,
  respectively

- heatmap title to display key options used, for easy reference

### Additional Features

- data centering can be disabled with `centerby_colnames=FALSE`.

- alternative hits can be displayed using `alt_sestats`. It does not
  subset heatmap rows, it inherits rows from `sestats`.

- display a subset of columns after row centering, useful to hide the
  control group for certain figures.

- option to display correlation heatmap, using the same data centering,
  then calculates Pearson correlation across sample columns.

- labels and legend grids can be customized to exact sizes with
  [`grid::gpar()`](https://rdrr.io/r/grid/gpar.html) and
  [`grid::unit()`](https://rdrr.io/r/grid/unit.html) definitions, for
  manuscript figures.

- mark annotations option to label a subset of rows

- row subclusters can be visualized using `row_subcluster` to drill down
  into specific subclusters from hierarchical clustering, k-means
  clustering, or any `row_split`.

### Data Centering

The intent is to display expression values from `assays(se)`, centered
across all columns, or with customization defined by `centerby_colnames`
and `normgroup_colnames`. The resulting centered data can be subsetted
by argument `isamples`, which occurs after centering in order to
decouple the centering step from the display of resulting data. To
subset samples involved in centering itself, either subset the input
`se` data, or supply `controlSamples` to define a subset of samples used
as the baseline in centering. See
[`jamma::centerGeneData()`](https://rdrr.io/pkg/jamma/man/centerGeneData.html)
for more details.

Paired data, also called repeated measures data, can be visualized by
including the pairing as `centerby_colnames` so that centering is
calculated within each pairing subgroup. In this case if also using
`controlSamples` to define a "time zero" or "baseline", then all
baseline samples will have exactly zero, if there is only one replicate
per pairing group at the baseline. In this case, it may be useful to
create the full heatmap once to confirm the centering is performing as
intended, then create a second heatmap using `isamples` to show only the
non-baseline samples - thus removing the large chunk of values with 0.

Note: data centering can be disabled with `centerby_colnames=FALSE`.

### Heatmap Title

A heatmap title is returned as an attribute `attr(hm, "hm_title")`,
which describes:

- total rows displayed, with `row_type` indicating the measured entity
  (gene, probe, DEGs, etc.)

- total columns displayed, with `column_type` indicating the sampled
  entity (samples, total replicates, etc.)

- the `assay_name` for the data being displayed

- relevant options for data centering, for example `"global-centered"`
  (by default) or `"Centered within Cell Line, versus Wildtype"`

To include the heatmap title:

`ComplexHeatmap::draw(hm, column_title=attr(hm, "hm_title))`

### Top and Left Annotations

The top heatmap annotations use `colData(se)` with user-supplied
`top_colnames` or by auto-detecting those colnames that apply to
multiple `colnames(se)`. Colors can be supplied using argument
`sample_color_list`, as described below.

The an incidence matrix of statistical hits can be displayed on the left
of the heatmap, using arguments `sestats` and `alt_sestats`. These
arguments can accept either the output of
[`se_contrast_stats()`](https://jmw86069.github.io/jamses/reference/se_contrast_stats.md),
or they can be a `numeric` matrix with values `c(-1, 0, 1)`, indicating
statistical hits down, no change, and up, respectively. The contrasts
can optionally be subset with `contrast_names`, which corresponds to
columns in the matrix if supplied in that format.

When `sestats` is supplied, it will subset all heatmap rows to include
only rows with at least one non-zero value in the incidence matrix. If
argument `rows` is supplied, then all `rownames(se)` matching `rows` are
displayed, regardless of statistical hits.

For comparison across other `sestats` results, argument `alt_sestats` is
treated similar to `sestats` except that the heatmap is not subset based
upon these values. That means the heatmap will be subset to match hits
defined by `sestats` but not `alt_sestats`. The `alt_sestats` incidence
matrix is displayed to the far left of the `sestats` incidence matrix.
For clarity, it can be useful to add `alt_sestats_suffix` to add a
suffix to each contrast label, for example if `sestats` represents limma
hits, use `sestats_suffix=" limma"`, and if `alt_sestats` represents
limma-voom hits, use `alt_sestats_suffix=" limmavoom"`.

Argument `rowData_colnames` can be supplied, which enables display of
`rowData(se)` annotations in the `left_annotation` of the heatmap.
Colors can be supplied using argument `sample_color_list`.

Argument `sample_color_list` is a `list` named by each annotation column
to be displayed as top or left annotation. Each list element is either:

- a `character` vector of R colors named by `character` value, or

- a `function` defined by
  [`circlize::colorRamp2()`](https://rdrr.io/pkg/circlize/man/colorRamp2.html)
  to be applied for `numeric` column values. In this case the `breaks`
  used to define the color function are used to define the color legend.

The function `platjam::design2colors()` can be used to create
`sample_color_list` starting with a `data.frame` of annotations, and
will soon be moved into this package.

A custom `left_annotation` can be supplied, but this method currently
prevents the other annotations described above from being displayed. To
display automated annotations with `rowData_colnames` and custom row
annotations, supply custom annotations with `right_annotation`. Note
that annotations must be supplied in exact row order, which is usually
easiest when supplying `rows` with specific set of rows.

### Compatible Input Formats

Data provided in `se` is expected to be `SummarizedExperiment`, however
other Bioconductor data types are accepted that provide accessor
functions: `featureData()`, `phenoData()`, and `assayData()`, including
for example the `"MethyLumiSet"` class.

Note that `matrix` input is currently not supported, however it can be
converted to `SummarizedExperiment` like this:

    se <- SummarizedExperiment::SummarizedExperiment(
       assays=list(data=matrix),
       rowData=data.frame(Gene=rownames(matrix)),
       colData=data.frame(Sample=colnames(matrix)))

## See also

Other jamses heatmaps:
[`detect_heatmap_components()`](https://jmw86069.github.io/jamses/reference/detect_heatmap_components.md),
[`heatmap_column_group_labels()`](https://jmw86069.github.io/jamses/reference/heatmap_column_group_labels.md),
[`heatmap_profile_plot()`](https://jmw86069.github.io/jamses/reference/heatmap_profile_plot.md)

## Examples

``` r
se <- make_se_test(nrow=1000, ngroups=4, nreps=8)

# optionally define factor levels to force the order of labels
SummarizedExperiment::rowData(se)$Class <- factor(
   sample(head(LETTERS, 5), size=nrow(se), replace=TRUE))

# basic heatmap
hm <- heatmap_se(se, rowData_colnames="Class")
#> 'magick' package is suggested to install to give better rasterization.
#> 
#> Set `ht_opt$message = FALSE` to turn off this message.

# draw by printing hm, or call draw() to add useful options
ComplexHeatmap::draw(hm,
   column_title=attr(hm, "hm_title"),
   merge_legends=TRUE)


# define specific colors
sample_color_list <- list(
   group=colorjam::group2colors(
      unique(SummarizedExperiment::colData(se)$group)),
   Class=colorjam::group2colors(
      unique(SummarizedExperiment::rowData(se)$Class)))

heatmap_se(se,
   rowData_colnames="Class",
   show_left_annotation_name="top",
   left_annotation_name_rot=150,
   sample_color_list=sample_color_list)
#> 'magick' package is suggested to install to give better rasterization.
#> 
#> Set `ht_opt$message = FALSE` to turn off this message.


# split rows by "Class"
heatmap_se(se,
   rowData_colnames="Class",
   row_split="Class",
   sample_color_list=sample_color_list)
#> 'magick' package is suggested to install to give better rasterization.
#> 
#> Set `ht_opt$message = FALSE` to turn off this message.


# let's have some fun now
hm2 <- heatmap_se(se,
   column_split=c("group"),
   column_title_rot=90,
   row_split=c("Class"),
   rowData_colnames=c("Class"),
   cluster_row_slices=FALSE,
   sample_color_list=sample_color_list)
#> 'magick' package is suggested to install to give better rasterization.
#> 
#> Set `ht_opt$message = FALSE` to turn off this message.
hm2drawn <- ComplexHeatmap::draw(hm2, merge_legends=TRUE)


# as an example, extract the row order
# technically you should use hm2drawn, but usually hm2 is enough
hro <- jamba::heatmap_row_order(hm2drawn);
jamba::sdim(hro)
#>   rows     class
#> A  196 character
#> B  210 character
#> C  199 character
#> D  215 character
#> E  180 character
lapply(hro, head, 7)
#> $A
#>   row_0640   row_0405   row_0417   row_0730   row_0847   row_0582   row_0858 
#> "row_0640" "row_0405" "row_0417" "row_0730" "row_0847" "row_0582" "row_0858" 
#> 
#> $B
#>   row_0087   row_0075   row_0421   row_0624   row_0238   row_0878   row_0090 
#> "row_0087" "row_0075" "row_0421" "row_0624" "row_0238" "row_0878" "row_0090" 
#> 
#> $C
#>   row_0834   row_0593   row_0670   row_0823   row_0066   row_0902   row_0422 
#> "row_0834" "row_0593" "row_0670" "row_0823" "row_0066" "row_0902" "row_0422" 
#> 
#> $D
#>   row_0349   row_0755   row_0633   row_0784   row_0208   row_0265   row_0382 
#> "row_0349" "row_0755" "row_0633" "row_0784" "row_0208" "row_0265" "row_0382" 
#> 
#> $E
#>   row_0948   row_0115   row_0874   row_0165   row_0227   row_0598   row_0781 
#> "row_0948" "row_0115" "row_0874" "row_0165" "row_0227" "row_0598" "row_0781" 
#> 
# (the names will differ from values when `row_labels` are customized)

# center by WildType samples
# - controlSamples
# - control_label
hm2 <- heatmap_se(se,
   controlSamples=rownames(subset(
      SummarizedExperiment::colData(se), group %in% "groupA")),
   control_label="vs groupA",
   column_split=c("group"),
   column_title_rot=90,
   row_split=c("Class"),
   rowData_colnames=c("Class"),
   cluster_row_slices=FALSE,
   sample_color_list=sample_color_list)
#> 'magick' package is suggested to install to give better rasterization.
#> 
#> Set `ht_opt$message = FALSE` to turn off this message.
hm2drawn <- ComplexHeatmap::draw(hm2,
   column_title=attr(hm2, "hm_title"),
   merge_legends=TRUE)


# add "callout" labels for a subset of rows
mark_rows <- c(sample(jamba::heatmap_row_order(hm2drawn)[[1]], size=5),
   sample(jamba::heatmap_row_order(hm2drawn)[[1]], size=3));

# turn off ComplexHeatmap warning when using RStudio
ComplexHeatmap::ht_opt(message=FALSE)

hm3 <- heatmap_se(se,
   mark_rows=mark_rows,
   controlSamples=rownames(
      subset(SummarizedExperiment::colData(se), group %in% "groupA")),
   control_label="vs groupA",
   column_split=c("group"),
   column_title_rot=90,
   row_split=c("Class"),
   rowData_colnames=c("Class"),
   cluster_row_slices=FALSE,
   sample_color_list=sample_color_list)
ComplexHeatmap::draw(hm3,
   column_title=attr(hm3, "hm_title"),
   merge_legends=TRUE)


# sestats can accept list, incidence matrix, hit_array, or sestats
# this example defines random set of hits
sestats_list <- list(
   contrast1=setNames(sample(c(1, -1), replace=TRUE, size=50),
      sample(rownames(se), size=50)),
   contrast2=setNames(sample(c(1, -1), replace=TRUE, size=50),
      sample(rownames(se), size=50)))
hm4 <- heatmap_se(se,
   controlSamples=rownames(
      subset(SummarizedExperiment::colData(se), group %in% "groupA")),
   control_label="vs groupA",
   sestats=sestats_list,
   column_split=c("group"),
   row_split=c("Class"),
   rowData_colnames=c("Class"),
   cluster_row_slices=FALSE,
   sample_color_list=sample_color_list)
ComplexHeatmap::draw(hm4,
   column_title=attr(hm4, "hm_title"),
   merge_legends=TRUE)


# it doesn't take much effort to run stats really quick
sedesign <- groups_to_sedesign(SummarizedExperiment::colData(se)[, "group", drop=FALSE])
contrast_names(sedesign) <- jamba::vigrep("-groupA", contrast_names(sedesign))
sestats <- se_contrast_stats(se=se,
   fold_cutoff=4,
   sedesign=sedesign, assay_name="counts")
hm4s <- heatmap_se(se,
   controlSamples=rownames(
      subset(SummarizedExperiment::colData(se), group %in% "groupA")),
   control_label="vs groupA",
   sestats=sestats,
   column_split=c("group"),
   row_split=6,
   rowData_colnames=c("Class"),
   cluster_row_slices=FALSE,
   sample_color_list=sample_color_list)
ComplexHeatmap::draw(hm4s,
   column_title=attr(hm4s, "hm_title"),
   merge_legends=TRUE)


# for fun, "drill down" into cluster 5
hm4s_4 <- heatmap_se(se,
   controlSamples=rownames(
      subset(SummarizedExperiment::colData(se), group %in% "groupA")),
   control_label="vs groupA",
   sestats=sestats,
   column_split=c("group"),
   row_split=6,
   row_subcluster=4,
   rowData_colnames=c("Class"),
   cluster_row_slices=FALSE,
   sample_color_list=sample_color_list)
#> Warning: The heatmap has not been initialized. You might have different results
#> if you repeatedly execute this function, e.g. when row_km/column_km was
#> set. It is more suggested to do as `ht = draw(ht); row_order(ht)`.
ComplexHeatmap::draw(hm4s_4,
   column_title=attr(hm4s_4, "hm_title"),
   merge_legends=TRUE)



# sestats can be provided as an incidence matrix
if (jamba::check_pkg_installed("venndir")) {
# convert sestats to list
sestats_hitlist <- hit_array_to_list(sestats)
# convert sestats hitlist to incidence matrix
# - for fun, use only the first two contrasts
sestats_hitim <- venndir::list2im_value(sestats_hitlist[1:2])
print(head(sestats_hitim));

# convert sestats_list to signed incidence matrix
sestats_im <- venndir::list2im_value(sestats_list)
print(head(sestats_im, 10));
# if the list has items (no direction) use venndir::list2im_opt()

hm5 <- heatmap_se(se,
   controlSamples=rownames(
      subset(SummarizedExperiment::colData(se), group %in% "groupA")),
   control_label="vs groupA",
   sestats=sestats_hitim,
   column_split=c("group"),
   rowData_colnames=c("Class"),
   cluster_row_slices=FALSE,
   sample_color_list=sample_color_list)
ComplexHeatmap::draw(hm5,
   column_title=attr(hm5, "hm_title"),
   merge_legends=TRUE)
}
#>          groupB-groupA groupC-groupA
#> row_0022            -1            -1
#> row_0030            -1             0
#> row_0066             1             0
#> row_0075             1             0
#> row_0080            -1             0
#> row_0087             1             0
#>          contrast1 contrast2
#> row_0230         1         0
#> row_0428         1         0
#> row_0080         1         1
#> row_0839        -1         0
#> row_0342        -1         0
#> row_0370        -1         0
#> row_0562         1         0
#> row_0783        -1         0
#> row_0333         1         0
#> row_0687        -1         0



# customize column label fonts using column_names_gp
column_bold <- ifelse(
   SummarizedExperiment::colData(se)$group %in% "groupA",
   2, 1);
hm6 <- heatmap_se(se,
   controlSamples=rownames(
      subset(SummarizedExperiment::colData(se), group %in% "groupA")),
   control_label="vs WildType",
   column_names_gp=grid::gpar(col=sample_color_list$group[
      as.character(SummarizedExperiment::colData(se)$group)],
      font=column_bold),
   column_split=c("group"),
   row_split=c("Class"),
   rowData_colnames=c("Class"),
   cluster_row_slices=FALSE,
   sample_color_list=sample_color_list)
ComplexHeatmap::draw(hm6,
   column_title=attr(hm6, "hm_title"),
   merge_legends=TRUE)


# correlation=TRUE, any heatmap becomes a sample correlation heatmap
hm6corr <- heatmap_se(se,
   correlation=TRUE,
   apply_hm_column_title=TRUE,
   controlSamples=rownames(
      subset(SummarizedExperiment::colData(se), group %in% "groupA")),
   control_label="vs groupA",
   column_names_gp=grid::gpar(col=sample_color_list$group[
      as.character(SummarizedExperiment::colData(se)$group)],
      font=rep(c(1, 2, 1), c(3, 5, 24))),
   column_split=c("Group"),
   sample_color_list=sample_color_list)
ComplexHeatmap::draw(hm6corr,
   merge_legends=TRUE)


## Final heatmap:
# 1. Applies heatmap title automatically.
# 2. Hides the top_colnames
# 3. Adds fancy grouped labels above the heatmap.
#
# apply_hm_column_title=TRUE
#    convenient way to define a title,
#    but it does not also display column_split labels
#
# hm_title_buffer=4
#    convenient way to insert some whitespace lines
#
# heatmap_column_group_labels()
#    adds to a drawn heatmap - it must already be drawn
#
SummarizedExperiment::colData(se)$Genotype <- rep(c("WT", "KO"), each=16);
SummarizedExperiment::colData(se)$Treatment <- rep(c("Control", "Dex"), each=8);
hm7 <- heatmap_se(se,
   apply_hm_column_title=TRUE,
   hm_title_buffer=3,
   controlSamples=rownames(
      subset(SummarizedExperiment::colData(se), group %in% "groupA")),
   control_label="vs groupA",
   sestats=sestats_list,
   top_colnames=FALSE,
   column_split=c("group"),
   row_split=c("Class"),
   rowData_colnames=c("Class"),
   cluster_row_slices=FALSE,
   sample_color_list=sample_color_list)
hm7_drawn <- ComplexHeatmap::draw(hm7,
   merge_legends=TRUE)

# now add fancy labels
heatmap_column_group_labels(
   hm_group_list=c("Treatment", "Genotype"),
   se=se,
   hm_drawn=hm7_drawn)

# Note: this step does not work consistently inside RStudio plot pane,
# in that case call dev.new() then run the step above to create hm7_drawn,
# then repeat the step below
#
# adjust the height of labels with argument y_offset_lines
# with positive values (upward), or negative values (downward).
```
