
#' Expression heatmap of sestats hits
#'
#' Expression heatmap of sestats hits
#'
#' Note: Still a work in progress.
#'
#' This function is a bold attempt to automate the intricate task
#' of creating an expression heatmap, using `ComplexHeatmap::Heatmap()`.
#' This function attempts to apply reasonable defaults to the many
#' possible varieties of heatmaps that can be produced, including:
#'
#' * Optional top and left annotations based upon data in `se`
#' from `colData(se)` and `rowData(se)`, respectively. Categorical
#' colors are re-used from `sample_color_list` supplied, otherwise
#' sensible default values are defined.
#' * Data can use any `assays(se)`.
#' * Data is centered by default, unless `centerby_colnames=FALSE`.
#' * Optional `sestats` can be supplied to subset the heatmap to
#' statistical hits defined in `sestats$hit_array`. A hit incidence
#' matrix is displayed as left annotation, indicating up- and
#' down-regulation.
#' * Another optional `alt_sestats` can be supplied to
#' compare statistical hits from `sestats` and `alt_sestats`, although
#' the heatmap rows are only subset based upon `sestats`.
#' * Optionally display a subset of sample columns after centering,
#' so data can be centered by control group, then display only the
#' non-control group samples.
#' * Optional sample correlation heatmap, which re-uses the same data
#' centering, then calculates Pearson correlation across sample columns.
#' * Rows and columns can be split by annotations, using annotations in
#' `rowData()` and `colData()`, respectively, or by `character` or `factor`,
#' or by `integer` which will cluster then split data into that many
#' subclusters.
#' * Various labels and legend grids can be customized to exact sizes
#' based upon `grid::gpar()` and `grid::unit()` definitions, useful
#' for publication figures.
#' * Row mark annotations can be used to label only a subset of rows,
#' useful when the heatmap includes far too many labels to read them all.
#' * Specific row split subclusters can be visualized using
#' `row_subcluster` to define the specific `row_split` entry, useful
#' for drilling into a specific subcluster from hierarchical clustering
#' without the manual effort to extract the subset of rows in that cluster
#' then re-running a new `heatmap_se()`.
#'
#'
#' The intent is to display expression values from `assays(se)`,
#' centered across all columns, or with customization defined by
#' `centerby_colnames` and `normgroup_colnames`. The resulting centered
#' data can be subsetted by argument `isamples`, which occurs after
#' centering in order to decouple the centering step from the display
#' of resulting data. To subset samples involved in centering itself,
#' either subset the input `se` data, or supply `controlSamples` to
#' define a subset of samples used as the baseline in centering.
#' See `jamma::centerGeneData()` for more details.
#'
#' Note: data centering can be disabled with `centerby_colnames=FALSE`.
#'
#' The top heatmap annotations use `colData(se)` with user-supplied
#' `top_colnames` or by auto-detecting those colnames that apply
#' to multiple `colnames(se)`.
#'
#' The `hit_array` data is used to define an incidence matrix of up/down
#' hits, which is displayed to the left of the heatmap. The contrasts
#' can optionally be subset with `contrast_names`.
#'
#' The heatmap title is returned as an `attr(hm, "hm_title")` that
#' describes the assay data used, data centering, and total rows
#' displayed. Draw the resulting heatmap like this:
#'
#' `ComplexHeatmap::draw(hm, column_title=attr(hm, "hm_title))`
#'
#' For comparison across other `sestats` results, argument `alt_sestats`
#' allows supplying an alternative hit array. These hit arrays are placed
#' as `left_annotation`, alongside optional data defined by `rowData_colnames`.
#'
#' When `rowData_colnames` is supplied, data in the corresponding colnames
#' of `rowData(se)` are also displayed in `left_annotation`. Colors can
#' be defined in `sample_color_list`.
#'
#' Argument `sample_color_list` is a `list` named by each annotation column
#' to be displayed as top or left annotation. Each list element is a vector
#' of R colors named by `character` value, or for `numeric` columns is a
#' color `function` as produced by `circlize::colorRamp2()`.
#' The function `platjam::design2colors()` is intended to create
#' `sample_color_list`, and will soon be moved into this package.
#'
#' A custom `left_annotation` can be supplied, but this method currently
#' prevents the other annotations described above from being displayed.
#' Currently the best way to supply custom row annotations in addition
#' to those described above, supply `right_annotation` to be displayed
#' on the right side of the heatmap.
#'
#' Data provided in `se` is expected to be `SummarizedExperiment`, however
#' it also accepts other Bioconductor data types that provide
#' accessor functions `featureData()`, `phenoData()`, and `assayData()`,
#' including for example `"MethyLumiSet"` class.
#'
#' @param se `SummarizedExperiment` object with accessor functions:
#'    `rowData()`, `colData()`, and `assays()`;
#'    or another suitable Bioconductor object with accessor functions:
#'    `featureData()`, `phenoData()`, and `assayData()`.
#' @param sestats one of the following types of data:
#'    * `list` output from `se_contrast_stats()`, which
#'    specifically contains `hit_array` as a 3-dimensional array of hits
#'    with dimensions "Cutoffs", "Contrasts", "Signal".
#'    * `numeric` matrix intended to represent an incidence matrix,
#'    where a value `0` indicates absence, and non-zero indicates presence.
#'    This format is useful for supplying any incidence matrix, such as
#'    gene-by-pathway (for example Github package "jmw86069/multienrichjam"
#'    provides `mem$memIM` with gene-by-pathway matrix),
#'    or gene-by-class (see Github package "jmw86069/pajam"
#'    for examples using ProteinAtlas protein classification, including
#'    membrane-bound, secreted, transcription factors, etc.), or any
#'    incidence matrix defined by Github "jmw86069/venndir" function
#'    `list2im_value()` or `list2im()` which converts input to a Venn diagram
#'    into an incidence matrix.
#'    * When `sestats` is supplied, data is converted to incidence matrix,
#'    then columns are matched with `contrast_names`. All rows with non-zero
#'    entry in those columns are included in the heatmap.
#'    When `rows` is also supplied, then the intersection of incidence
#'    matrix rows and `rows` is displayed in the heatmap.
#'    * Note that `alt_sestats` does not subset rows displayed in the
#'    heatmap.
#' @param rows `character` vector of `rownames(se)` to define a specific
#'    set of rows to display. When `sestats` is supplied, then the
#'    intersection of `rows` with genes defined by `sestats` is displayed.
#' @param row_type `character` string used in the title of the heatmap
#'    which indicates how many rows are displayed. For example
#'    `"1,234 genes detected above background"` or
#'    `"1,234 DEGs by limma-voom"`.
#' @param data_type `character` string used as title of the heatmap
#'    color gradient legend, for example `"expression"` indicates
#'    the data contains gene expression measurements. Notes:
#'    * The prefix `"centered"` is automatically appended whenever
#'    the data is also centered for the heatmap. Set `centerby_colnames=FALSE`
#'    to display data that is not centered.
#'    * The prefix `"correlation of"` is automatically appended when
#'    `correlation=TRUE` which displays correlation of whatever data
#'    is included in the heatmap.
#' @param correlation `logical` indicating whether to calculate sample
#'    correlation, and plot a sample-by-sample correlation heatmap.
#'    This option is included here since many of the same arguments
#'    are required for data centering, and sample annotations.
#'    Note that `color_max` is forced to a maximum value of `1.0`,
#'    representing the maximum correlation value.
#' @param assay_name `character` string indicating the name in
#'    `assays(se)` to use for data to be displayed in the heatmap.
#'    When multiple `assay_name` values are supplied, the first
#'    assay_name that matches `names(assays(se))` will be used in the
#'    heatmap. In this way, multiple `assay_names` can be supplied to
#'    define statistical hits in `sestats`, which calls `hit_array_to_list()`
#'    to combine hits across `assay_name` entries; but only the first
#'    `assay_name` found in `se` is used for the heatmap values.
#' @param contrast_names `character` vector of contrasts in
#'    `sestats$hit_array` to use for the heatmap. When `contrast_names=NULL`
#'    then all contrasts are displayed, which is the default.
#' @param contrast_suffix `character` string with optional suffix to append
#'    to the end of each contrast name label for `sestats` hit incidence
#'    matrix beside the heatmap. This suffix may be useful when comparing
#'    two methods for the same set of contrast names, with `sestats` and
#'    `alt_sestats`.
#' @param cutoff_name `character` or `integer` index used to define the
#'    specific statistical cutoffs to use from `sestats$hit_array`. This
#'    argument is passed to `hit_array_to_list()` as `cutoff_names`.
#' @param alt_sestats,alt_assay_name,alt_contrast_names,alt_contrast_suffix
#'    arguments analogous to those described above for `sestats` which
#'    are used when `alt_sestats` is supplied.
#' @param isamples `character` vector of `colnames(se)` used to visualize a
#'    subset of samples used for the data centering step. Note that
#'    data centering uses all columns supplied in `se`, and after centering,
#'    the subset of columns defined in `isamples` is displayed in the heatmap.
#'    This distinction makes it possible to center data by some control group,
#'    then optionally not display the control group data.
#' @param normgroup_colname `character` vector of colnames in `colData(se)`
#'    used during data centering. When supplied, samples are centered
#'    independently within each normgroup grouping. These values are
#'    equivalent to using `centerby_colnames`.
#' @param centerby_colnames either:
#'    * `character` vector of colnames in `colData(se)`
#'    used during data centering. When supplied, samples are centered
#'    independently within each centerby grouping. It is typically used
#'    for things like cell lines, to center each cell line by a time
#'    point control, or untreated control.
#'    * `NULL` to perform centering across all columns in `se`.
#'    * `FALSE` to disable centering.
#' @param controlSamples `character` optional vector of samples to use as the
#'    reference during data centering. Note that samples are still
#'    centered within each normgroup and centerby grouping, and within
#'    that grouping samples are centered to the `controlSamples`
#'    which are present in that grouping. Any center group for which no
#'    samples are defined in `controlSamples` will use all samples in that
#'    center group. Typically, `controlSamples` is used to define a
#'    specific group as the reference for centering, so changes are displayed
#'    relative to that group. Make sure to define `control_name` to include
#'    an appropriate label in the heatmap title.
#' @param control_name `character` string used in heatmap title
#'    to describe the control used during data centering, relevant when
#'    `controlSamples` is also supplied.
#' @param top_colnames one of the following types:
#'    * `character` vector of colnames to use from
#'    `colData(se)` as annotations to display in `top_annotation` above
#'    the heatmap.
#'    * `NULL`, will call `choose_annotation_colnames()` to detect
#'    reasonable colnames: columns with more than one unique value;
#'    columns with at least one duplicated value.
#'    * `FALSE` will hide the `top_colnames`, which also occurs when
#'    `colData(se)` is empty.
#' @param top_annotation specific heatmap annotation as defined by
#'    `ComplexHeatmap::HeatmapAnnotation()`. When supplied, the `top_colnames`
#'    described above is not used.
#' @param top_annotation_name_gp `grid::gpar` object to customize the
#'    annotation name displayed beside the top annotation.
#' @param rowData_colnames `character` vector of colnames in `rowData(se)`
#'    to use for heatmap annotations displayed on the left side of
#'    the heatmap. Specific colors can be included in `sample_color_list`
#'    as a named `list` of color vectors or color functions. The names
#'    of this list must match colnames to be displayed, otherwise
#'    `ComplexHeatmap::Heatmap()` will define its own color function.
#' @param left_annotation specific heatmap annotation as defined by
#'    `ComplexHeatmap::rowAnnotation()`. When supplied, the `rowData_colnames`
#'    and `sestats` row annotations are not displayed. In order to supply
#'    custom row annotations and not lose `left_annotation` defined above,
#'    supply the row annotations as `right_annotation`.
#' @param left_annotation_name_gp `grid::gpar` object to customize the
#'    annotation name displayed beside the left annotation.
#' @param left_annotation_name_rot `numeric` rotation of left annotation
#'    label, in degrees, where `0` indicates normal text, and `90` is
#'    rotated vertically.
#' @param right_annotation specific heatmap annotation as defined by
#'    `ComplexHeatmap::HeatmapAnnotation()`. This element is created
#'    automatically when `mark_rows` is supplied.
#' @param simple_anno_size `grid::unit` size used to define heatmap
#'    annotation sizes (height or width of each line) for any simple
#'    annotations.
#' @param legend_title_gp `grid::gpar` to customize the legend title
#'    fonts, applied to each legend: top annotation, left annotation,
#'    main heatmap.
#' @param legend_labels_gp `grid::gpar` to customize the legend label
#'    fonts, applied to each legend: top annotation, left annotation,
#'    main heatmap.
#' @param legend_grid_cex `numeric` multiplied to adjust the relative
#'    size of each legend grid unit, applied to each relevant metric.
#' @param row_names_gp `gpar` to define custom column name settings.
#'    When `"fontsize"` is not defined, the automatic font size calculation
#'    is added to the `row_names_gp` supplied.
#' @param row_split is used to define heatmap split by row, ultimately
#'    passed to `ComplexHeatmap::Heatmap()` argument `row_split`. However,
#'    the input type can vary:
#'    * `integer` number of row splits based upon row clustering
#'    * `character` value or values in colnames of `rowData(se)` to split
#'    using row annotation in `se`.
#'    * `character` or `factor` vector named by `rownames(se)` with another
#'    custom row split, passed directly to `ComplexHeatmap::Heatmap()`
#'    argument `row_split`, with proper order for rows being displayed
#' @param row_subcluster `integer` or `character` vector representing one
#'    or more elements returned by `row_split` to use as a drill-down
#'    sub-cluster heatmap. This argument is experimental, and is intended
#'    to make it easy to "drill down" into specific row clusters.
#'    * The process internally creates a full heatmap using all arguments
#'    as defined, then extracts the `jamba::heatmap_row_order()` which
#'    contains row split data in a `list` of rownames vectors. The `list`
#'    elements that match `row_subcluster` are extracted and used again
#'    for a subsequent heatmap, and are displayed in the same order
#'    in which they appear in the original full heatmap - which means
#'    `cluster_rows=FALSE` is defined at this point. However `row_split`
#'    is retained for this subset of rows, to indicate the original
#'    row split annotation.
#'    * Note that `row_subcluster` must match the `names()` returned
#'    by `jamba::heatmap_row_order()` for the full heatmap, or should
#'    include a `numeric` index for the `list` element or elements to
#'    use.
#'    * In principle this process would be run in two stages: First,
#'    view a heatmap with `row_split=6`, then re-run the same heatmap
#'    with `row_subcluster=4` to see cluster number 4 from the full
#'    heatmap.
#' @param row_title_rot `numeric` value indicating text rotation in degrees
#'    to use for row titles.
#' @param sample_color_list named `list` of color vectors or color functions,
#'    where names correspond to colnames in either `colData(se)` or
#'    `rowData(se)`, and which are passed to corresponding left or top
#'    annotation functions. When colors are not defined,
#'    `ComplexHeatmap::Heatmap()` will define colors using its own internal
#'    function.
#' @param legend_at,legend_labels `numeric` and `character`, respectively,
#'    to define custom values for the heatmap color gradient legend.
#'    When `legend_at` is NULL, it uses `x=ceiling(color_max)` and defines
#'    breaks at each integer value from `-x` to `+x`.
#'    When `correlation=TRUE` the `legend_labels` by default use `legend_at`.
#'    Otherwise, `legend_labels` values inverse transformed from `log2(1 + x)`
#'    in order to display normal space fold change values, for example
#'    the `legend_at=c(-2, -1, 0, 1, 2)` would result in:
#'    `legend_labels=c("-4", "-2", "1", "2", "4")`.
#'    To override this behavior, supply both the custom `legend_at` values
#'    and corresponding `legend_labels`.
#' @param subset_legend_colors `logical` indicating whether to subset colors
#'    shown in the color key defined by `sample_color_list`, which is useful
#'    when the heatmap only represents a subset of color values.
#'    When `subset_legend_colors == TRUE`, the color key will only
#'    include colors shown in the `top_annotation`, when
#'    `subset_legend_colors == FALSE` all colors defined in
#'    `sample_color_list` will be included for each relevant column.
#' @param row_cex,column_cex `numeric` values used to adjust the row and
#'    column name font size, relative to the automatic adjustment that
#'    is already done based upon the number of rows and columns being
#'    displayed.
#' @param row_anno_fontsize `numeric` base font size for row
#'    annotation labels. This value is only used when `left_annotation_name_gp`
#'    is not supplied. Note these labels appears underneath row annotations,
#'    alongside column labels, and therefore they are also adjusted
#'    by multiplying `column_cex` so these labels are adjusted together.
#' @param useMedian `logical` passed to `jamma::centerGeneData()` during
#'    data centering.
#' @param show_row_names,show_row_dend `logical` indicating whether to
#'    display row names, and row dendrogram, respectively. With more than
#'    2,000 rows this step can become somewhat slow.
#' @param mark_rows `character` vector of values in `rownames(se)` that
#'    should be labeled using `ComplexHeatmap::anno_mark()` in call-out
#'    style. Usually this argument is used when `show_row_names=FALSE`,
#'    hiding the row labels, but is not required. Values in `mark_rows`
#'    are intersected with rows displayed in the heatmap, therefore only
#'    matching entries will be labeled.
#' @param mark_labels_gp `grid::gpar` to customize the font used by labels
#'    when `mark_rows` is supplied.
#' @param show_heatmap_legend,show_left_legend,show_top_legend `logical`
#'    indicating whether each legend should be displayed. Sometimes there
#'    are too many annotations, and the color legends can overwhelm the
#'    figure. Note that `show_left_legend` is applied in a specific order,
#'    with these rules:
#'    * `show_left_legend` is extended to at least length 2, then values
#'    are used in order for: `sestats`, `rowData_colnames`, in order,
#'    using whichever is defined.
#'    * If `sestats` is defined, the first value in `show_left_legend`
#'    is used for this annotation, then the remaining values are used
#'    for `rowData_colnames`. Setting the first `show_left_legend` value
#'    to `FALSE` will ensure the legend for `sestats` is not displayed.
#'    * If `rowData_colnames` is defined, then the remaining values in
#'    `show_left_legend` are recycled for all columns in
#'    `rowData_colnames`, and applied in order.
#'    In this way, individual columns can have the legend displayed or hidden.
#'    * If `alt_sestats` is defined, the legend is always hidden, in favor
#'    of showing only the legend for `sestats` without duplicating this legend.
#' @param show_top_annotation_name,show_left_annotation_name `logical`
#'    indicating whether to display the annotation name beside the top and
#'    left annotations, respectively.
#' @param row_label_colname `character` string used as a row label, where
#'    this value is a colname in `rowData(se)`. It is useful when rownames
#'    are some identifier that is not user-friendly, and where another column
#'    in the data may provide a more helpful label, for example `"SYMBOL"`
#'    to display gene symbol instead of accession number.
#' @param cluster_columns,cluster_rows `logical` indicating whether
#'    to cluster columns by hierarchical clustering; or `function` with
#'    a specific function that produces `hclust` or `dendrogram` output,
#'    given a `numeric` matrix. Note that `cluster_rows` default will replace
#'    `NA` values with zero `0` to avoid errors with missing data, and
#'    uses `amap::hcluster()` by default which is a one-step compiled
#'    process to perform distance calculation and hierarchical clustering.
#' @param column_names_gp `gpar` to define custom column name settings.
#'    When `"fontsize"` is not defined, the automatic font size calculation
#'    is added to the `column_names_gp` supplied.
#' @param column_split `character` or `integer` vector used to define
#'    heatmap column split.
#' @param column_split_sep `character` string used as delimited when
#'    `column_split` defines multiple split levels.
#' @param color_max `numeric` value passed to `colorjam::col_div_xf()`
#'    which defines the upper limit of color gradient used in the heatmap.
#' @param lens `numeric` value passed to `colorjam::col_div_xf()` to control
#'    the intensity of color gradient applied to the numeric range.
#' @param rename_contrasts,rename_alt_contrasts `logical` indicating
#'    whether to rename long contrast names in `sestats` and `alt_sestats`
#'    using `contrast2comp()`.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param debug `logical` indicating debug mode, data is returned in a `list`:
#'    * `hm` object `ComplexHeatmap::Heatmap`
#'    * `top_annotation` object `ComplexHeatmap::HeatmapAnnotation` for columns
#'    * `left_annotation` object `ComplexHeatmap::HeatmapAnnotation` for rows
#'    * `hm_title` object `character` string with the heatmap title.
#' @param ... additional arguments are passed to supporting functions.
#'
#' @export
heatmap_se <- function
(se,
   sestats=NULL,
   rows=NULL,
   row_type="rows",
   data_type="expression",
   correlation=FALSE,
   assay_name=NULL,
   contrast_names=NULL,
   contrast_suffix="",
   cutoff_name=NULL,
   alt_sestats=NULL,
   alt_assay_name=assay_name,
   alt_contrast_names=NULL,
   alt_contrast_suffix="",
   alt_cutoff_name=NULL,
   isamples=colnames(se),
   normgroup_colname=NULL,
   centerby_colnames=NULL,
   controlSamples=NULL,
   control_label="",
   top_colnames=NULL,
   top_annotation=NULL,
   top_annotation_name_gp=grid::gpar(),
   rowData_colnames=NULL,
   left_annotation=NULL,
   left_annotation_name_gp=grid::gpar(),
   left_annotation_name_rot=90,
   right_annotation=NULL,
   simple_anno_size=grid::unit(8, "mm"),
   legend_title_gp=grid::gpar(fontsize=10),
   legend_labels_gp=grid::gpar(fontsize=10),
   legend_grid_cex=1,
   row_names_gp=NULL,
   row_split=NULL,
   row_subcluster=NULL,
   row_title_rot=0,
   sample_color_list=NULL,
   legend_at=NULL,
   legend_labels=NULL,
   subset_legend_colors=TRUE,
   row_cex=0.8,
   column_cex=1,
   row_anno_fontsize=11,
   useMedian=FALSE,
   show_row_names=NULL,
   show_row_dend=length(rows) < 2000,
   mark_rows=NULL,
   mark_labels_gp=grid::gpar(),
   show_heatmap_legend=TRUE,
   show_top_legend=TRUE,
   show_left_legend=TRUE,
   show_top_annotation_name=TRUE,
   show_left_annotation_name=TRUE,
   row_label_colname=NULL,
   cluster_columns=FALSE,
   cluster_rows=function(x, ...){
      amap::hcluster(rmNA(naValue=0, x),
         ...,
         method="euclidean",
         link="ward")},
   column_names_gp=NULL,
   column_split=NULL,
   column_split_sep=",",
   color_max=3,
   lens=2,
   rename_contrasts=TRUE,
   rename_alt_contrasts=TRUE,
   verbose=FALSE,
   debug=FALSE,
   ...)
{
   #
   if (!jamba::check_pkg_installed("ComplexHeatmap")) {
      stop("This function requires Bioconductor package ComplexHeatmap.");
   }
   if (!jamba::check_pkg_installed("venndir")) {
      stop("This function requires Github package venndir from 'jmw86069/venndir'");
   }
   if (!suppressPackageStartupMessages(require(SummarizedExperiment))) {
      stop("This function requires Bioconductor package SummarizedExperiment.");
   }
   if (length(correlation) == 0) {
      correlation <- FALSE;
   }

   # row_subcluster
   if (length(row_subcluster) > 0) {
      if (verbose) {
         jamba::printDebug("heatmap_se(): ",
            "Preparing sub-cluster heatmap.");
      }
      hm_total <- heatmap_se(
         se=se,
         sestats=sestats,
         rows=rows,
         row_type=row_type,
         data_type=data_type,
         correlation=correlation,
         assay_name=assay_name,
         contrast_names=contrast_names,
         contrast_suffix=contrast_suffix,
         cutoff_name=cutoff_name,
         alt_sestats=alt_sestats,
         alt_assay_name=alt_assay_name,
         alt_contrast_names=alt_contrast_names,
         alt_contrast_suffix=alt_contrast_suffix,
         alt_cutoff_name=alt_cutoff_name,
         isamples=isamples,
         normgroup_colname=normgroup_colname,
         centerby_colnames=centerby_colnames,
         controlSamples=controlSamples,
         control_label=control_label,
         top_colnames=top_colnames,
         top_annotation=top_annotation,
         top_annotation_name_gp=top_annotation_name_gp,
         rowData_colnames=rowData_colnames,
         left_annotation=left_annotation,
         left_annotation_name_gp=left_annotation_name_gp,
         left_annotation_name_rot=left_annotation_name_rot,
         right_annotation=right_annotation,
         simple_anno_size=simple_anno_size,
         legend_title_gp=legend_title_gp,
         legend_labels_gp=legend_labels_gp,
         legend_grid_cex=legend_grid_cex,
         row_split=row_split,
         row_subcluster=NULL,
         row_title_rot=row_title_rot,
         sample_color_list=sample_color_list,
         legend_at=legend_at,
         legend_labels=legend_labels,
         subset_legend_colors=subset_legend_colors,
         row_cex=row_cex,
         column_cex=column_cex,
         row_anno_fontsize=row_anno_fontsize,
         useMedian=useMedian,
         show_row_names=show_row_names,
         show_row_dend=show_row_dend,
         mark_rows=mark_rows,
         mark_labels_gp=mark_labels_gp,
         show_heatmap_legend=show_heatmap_legend,
         show_top_legend=show_top_legend,
         show_left_legend=show_left_legend,
         show_top_annotation_name=show_top_annotation_name,
         show_left_annotation_name=show_left_annotation_name,
         row_label_colname=row_label_colname,
         cluster_columns=cluster_columns,
         cluster_rows=cluster_rows,
         column_split=column_split,
         column_split_sep=column_split_sep,
         color_max=color_max,
         lens=lens,
         rename_contrasts=rename_contrasts,
         rename_alt_contrasts=rename_alt_contrasts,
         verbose=verbose,
         debug=debug,
         ...)
      row_order <- jamba::heatmap_row_order(hm_total);
      if (length(names(row_order)) == 0) {
         names(row_order) <- as.character(seq_along(row_order))
      }

      if (is.numeric(row_subcluster)) {
         row_order_use <- row_order[seq_along(row_order) %in% row_subcluster];
      } else {
         row_order_use <- row_order[names(row_order) %in% as.character(row_subcluster)];
      }
      if (verbose) {
         jamba::printDebug("heatmap_se(): ",
            "row_subcluster: ", row_subcluster);
         jamba::printDebug("heatmap_se(): ",
            "selected: '",
            sep="', '",
            names(row_order_use),
            "'")
      }
      # set new argument values for the drill-down heatmap
      rows <- unlist(unname(row_order_use));
      cluster_rows <- FALSE;
      row_split <- nameVector(rep(names(row_order_use),
         lengths(row_order_use)),
         rows);
      # subset the se data?
      se <- se[rows,];

      if (debug) {
         return(invisible(row_order));
      }
   }

   # sestats - define rows to use
   gene_hitlist <- NULL;
   alt_gene_hitlist <- NULL;
   gene_hits_im <- NULL;
   gene_hits <- NULL;
   alt_gene_hits_im <- NULL;
   if (length(sestats) > 0) {
      if ("list" %in% class(sestats) && "hit_array" %in% names(sestats)) {
         hit_array <- sestats$hit_array;
      } else if ("matrix" %in% class(sestats)) {
         gene_hits_im <- sestats;
         hit_array <- NULL;
      } else {
         hit_array <- sestats;
      }
      if (length(hit_array) == 0) {
         if (verbose) {
            jamba::printDebug("heatmap_se(): ",
               "sestats is using a custom incidence matrix.");
         }
         if (length(contrast_names) > 0 &&
               any(contrast_names %in% colnames(gene_hits_im))) {
            contrast_names <- intersect(contrast_names,
               colnames(gene_hits_im));
            gene_hits_im <- gene_hits_im[, contrast_names, drop=FALSE];
         }
         gene_hits <- rownames(gene_hits_im);
      } else {
         if (verbose) {
            jamba::printDebug("heatmap_se(): ",
               "sestats is generating an incidence matrix.");
         }
         if (length(contrast_names) == 0) {
            contrast_names <- dimnames(hit_array)[[2]];
         }
         gene_hitlist <- hit_array_to_list(hit_array,
            cutoff_names=cutoff_name,
            contrast_names=contrast_names,
            assay_names=assay_name);
         gene_hits <- names(jamba::tcount(names(unlist(unname(
            gene_hitlist)))));
         # confirm all gene_hits are present in the data provided
         if (!all(gene_hits %in% rownames(se))) {
            gene_hits <- intersect(gene_hits, rownames(se));
         }
         gene_hits_im <- venndir::list2im_value(gene_hitlist,
            do_sparse=FALSE)[gene_hits,,drop=FALSE];
      }

      # optionally rename contrasts
      if (rename_contrasts) {
         colnames(gene_hits_im) <- tryCatch({
            # paste0(".          ", gsub(":", ":",
               contrast2comp(colnames(gene_hits_im))
            # ));
         }, error=function(e){
            colnames(gene_hits_im)
         });
      }

      if (length(contrast_suffix) > 0 && any(nchar(contrast_suffix)) > 0) {
         colnames(gene_hits_im) <- paste0(colnames(gene_hits_im),
            contrast_suffix);
      }
   }

   # rows as user-defined vector
   rows <- intersect(rows, rownames(se));
   if (length(rows) > 0) {
      if (length(sestats) > 0) {
         rows_im <- (gene_hits_im * 0)[rep(1, length(rows)), , drop=FALSE];
         rownames(rows_im) <- rows;
         gene_hits_rows <- intersect(rows, gene_hits);
         if (length(gene_hits_rows) > 0) {
            rows_im[match(gene_hits_rows, rows),] <-
               gene_hits_im[match(gene_hits_rows, rownames(gene_hits_im)),,drop=FALSE];
         }
         gene_hits_im <- rows_im;
      } else {
         gene_hits_im <- NULL;
      }
      gene_hits <- rows;
      rows <- gene_hits;
   } else {
      rows <- rownames(se);
      if (length(gene_hits) > 0) {
         rows <- intersect(gene_hits, rows);
      }
      gene_hits <- rows;
   }

   # alt_sestats only for rows and gene_hits defined from sestats
   if (length(sestats) > 0 && length(alt_sestats) > 0) {
      if ("list" %in% class(alt_sestats) && "hit_array" %in% names(alt_sestats)) {
         alt_hit_array <- alt_sestats$hit_array;
      } else if ("matrix" %in% class(alt_sestats)) {
         alt_gene_hits_im1 <- alt_sestats;
         alt_hit_array <- NULL;
      } else {
         alt_hit_array <- alt_sestats;
      }
      if (length(alt_hit_array) == 0) {
         if (verbose) {
            jamba::printDebug("heatmap_se(): ",
               "alt_sestats is using a custom incidence matrix.");
         }
         if (length(alt_contrast_names) > 0 &&
               any(alt_contrast_names %in% colnames(alt_gene_hits_im1))) {
            alt_contrast_names <- intersect(alt_contrast_names,
               colnames(alt_gene_hits_im1));
         } else {
            alt_contrast_names <- colnames(alt_gene_hits_im1);
         }
         alt_gene_hits_im1 <- alt_gene_hits_im1[, alt_contrast_names, drop=FALSE];
      } else {
         if (length(alt_contrast_names) == 0) {
            alt_contrast_names <- dimnames(alt_hit_array)[[2]];
         }
         alt_gene_hitlist <- hit_array_to_list(alt_hit_array,
            cutoff_names=alt_cutoff_name,
            contrast_names=alt_contrast_names,
            assay_names=alt_assay_name);
         alt_gene_hits_im1 <- venndir::list2im_value(alt_gene_hitlist,
            do_sparse=FALSE);
      }
      # start with empty gene_hits_im, expand to ncol required here
      alt_gene_hits_im <- (gene_hits_im * 0)[,rep(1, ncol(alt_gene_hits_im1)), drop=FALSE];
      colnames(alt_gene_hits_im) <- colnames(alt_gene_hits_im1);
      # determine shared genes to populate with values
      genes_shared <- intersect(rownames(alt_gene_hits_im),
         rownames(alt_gene_hits_im1))
      if (length(genes_shared) > 0) {
         alt_gene_hits_im[genes_shared,] <- alt_gene_hits_im1[genes_shared,];
      }

      # optionally rename contrasts
      if (rename_contrasts) {
         colnames(alt_gene_hits_im) <- tryCatch({
            contrast2comp(colnames(alt_gene_hits_im));
         }, error=function(e){
            colnames(alt_gene_hits_im)
         })
      }

      if (length(alt_contrast_suffix) > 0 && any(nchar(alt_contrast_suffix)) > 0) {
         colnames(alt_gene_hits_im) <- paste0(colnames(alt_gene_hits_im),
            alt_contrast_suffix);
      }
   }

   # validate sample_color_list
   # remove NA
   # convert color name to hex
   if (length(sample_color_list) > 0) {
      sample_color_list <- lapply(sample_color_list, function(i){
         if (is.function(i)) {
            i
         } else {
            i <- i[!is.na(i) & !is.na(names(i))];
            i_is_hex <- grepl("^#", i);
            if (any(!i_is_hex)) {
               i[!i_is_hex] <- jamba::rgb2col(alpha=FALSE,
                  col2rgb(i[!i_is_hex]))
            }
            i;
         }
      })
   }

   # pull colData and rowData as data.frame
   # to be tolerant of other data types
   # Note: This process does not subset by `rows` or `isamples` yet
   if (grepl("SummarizedExperiment", ignore.case=TRUE, class(se))) {
      rowData_se <- data.frame(check.names=FALSE,
         rowData(se));
      colData_se <- data.frame(check.names=FALSE,
         colData(se))
   } else {
      if (verbose) {
         jamba::printDebug("heatmap_se(): ",
            "using accessor functions: ",
            c("featureData()", "phenoData()"))
      }
      rowData_se <- as(Biobase::featureData(se), "data.frame");
      rownames(rowData_se) <- rownames(se);
      colData_se <- as(Biobase::phenoData(se), "data.frame")
      rownames(colData_se) <- colnames(se);
   }

   # normgroup for column split
   normgroup_colname <- intersect(normgroup_colname,
      colnames(colData_se));
   if (length(column_split) == 0) {
      if (length(normgroup_colname) > 0 &&
            nrow(unique(colData_se[isamples, normgroup_colname, drop=FALSE])) > 1) {
         column_split <- jamba::pasteByRow(
            colData_se[isamples, normgroup_colname, drop=FALSE],
            sep=column_split_sep);
         names(column_split) <- isamples;
      } else {
         column_split <- NULL;
      }
   } else {
      if (any(c("factor", "character") %in% class(column_split))) {
         if (!any(duplicated(column_split)) &&
               all(column_split %in% colnames(colData_se))) {
            column_split <- jamba::pasteByRowOrdered(
               data.frame(check.names=FALSE,
                  colData_se[isamples, column_split, drop=FALSE]),
               keepOrder=TRUE,
               sep=column_split_sep);
            names(column_split) <- isamples;
         } else if (all(names(column_split) %in% isamples)) {
            column_split <- column_split[isamples];
         } else if (length(column_split) == length(isamples)) {
            # leave as-is but add isamples as names
            names(column_split) <- isamples;
         } else {
            column_split <- NULL;
         }
      } else if (length(column_split) == 1 && is.numeric(column_split)) {
         # leave as-is
      } else {
         column_split <- NULL;
      }
   }

   # column font size
   column_fontsize <- jamba::noiseFloor(
      column_cex * 60/(length(isamples))^(1/2),
      ceiling=20,
      minimum=2);

   # row font size
   if (correlation) {
      row_fontsize <- jamba::noiseFloor(
         row_cex * (60*(14 / 10))/(length(isamples))^(1/2),
         minimum=1,
         ceiling=20);
   } else {
      row_fontsize <- jamba::noiseFloor(
         row_cex * (60*(14 / 10))/(length(gene_hits))^(1/2),
         minimum=1,
         ceiling=20);
   }

   # choose interesting top_annotation colnames when none are supplied
   if (length(top_colnames) == 0) {
      top_colnames <- choose_annotation_colnames(colData_se,
         ...);
      if (length(top_colnames) > 0 && verbose) {
         jamba::printDebug("heatmap_se(): ",
            "derived top_colnames: ",
            top_colnames);
      }
   }
   top_colnames <- intersect(top_colnames,
      colnames(colData_se));
   if (length(top_annotation) == 0 &&
         length(top_colnames) > 0 &&
         !any(top_colnames %in% FALSE)) {
      # subset color key by data shown in the heatmap
      top_color_list <- NULL;
      # subset any factor columns to limit colors shown in the legend
      top_df <- data.frame(check.names=FALSE,
         colData_se[isamples, top_colnames, drop=FALSE]);
      for (top_colname in top_colnames) {
         if (is.factor(top_df[[top_colname]])) {
            top_df[[top_colname]] <- factor(top_df[[top_colname]]);
         }
      }
      if (any(top_colnames %in% names(sample_color_list))) {
         if (subset_legend_colors) {
            top_color_list <- lapply(jamba::nameVector(top_colnames), function(top_colname){
               sample_colors <- sample_color_list[[top_colname]];
               if (!is.function(sample_colors)) {
                  if (is.factor(top_df[isamples, top_colname])) {
                     uniq_values <- levels(top_df[isamples, top_colname]);
                  } else {
                     uniq_values <- jamba::mixedSort(unique(
                        as.character(
                           top_df[isamples, top_colname])));
                  }
                  sample_colors <- sample_colors[uniq_values];
                  if (length(sample_colors) == 0) {
                     sample_colors <- rep(NA, length.out=length(uniq_values));
                  }
                  names(sample_colors) <- uniq_values;
                  if (any(is.na(sample_colors))) {
                     # fallback plan for missing values is to assign
                     # generic rainbow categorical colors
                     sample_colors[is.na(sample_colors)] <- colorjam::rainbowJam(
                        n=sum(is.na(sample_colors)),
                        ...);
                  }
               }
               sample_colors;
            })
         } else {
            top_color_list <- sample_color_list[intersect(top_colnames, names(sample_color_list))];
         }
      }
      top_param_list <- c(
         lapply(jamba::nameVector(top_colnames), function(iname){
            if (iname %in% names(top_color_list)) {
               if (is.function(top_color_list[[iname]])) {
                  if ("breaks" %in% names(attributes(top_color_list[[iname]]))) {
                     list(border=TRUE,
                        title_gp=legend_title_gp,
                        labels_gp=legend_labels_gp,
                        grid_height=grid::unit(4 * legend_grid_cex, "mm"),
                        grid_width=grid::unit(4 * legend_grid_cex, "mm"),
                        color_bar="discrete",
                        at=attr(top_color_list[[iname]], "breaks"))
                  } else {
                     list(border=TRUE,
                        title_gp=legend_title_gp,
                        labels_gp=legend_labels_gp,
                        grid_height=grid::unit(4 * legend_grid_cex, "mm"),
                        grid_width=grid::unit(4 * legend_grid_cex, "mm"),
                        color_bar="discrete")
                  }
               } else {
                  if (subset_legend_colors) {
                     list(border=TRUE,
                        title_gp=legend_title_gp,
                        labels_gp=legend_labels_gp,
                        grid_height=grid::unit(4 * legend_grid_cex, "mm"),
                        grid_width=grid::unit(4 * legend_grid_cex, "mm"),
                        color_bar="discrete",
                        at=jamba::rmNA(names(top_color_list[[iname]])))
                  } else {
                     list(border=TRUE,
                        grid_height=grid::unit(4 * legend_grid_cex, "mm"),
                        grid_width=grid::unit(4 * legend_grid_cex, "mm"),
                        title_gp=legend_title_gp,
                        labels_gp=legend_labels_gp)
                  }
               }
            } else {
               list(border=TRUE,
                  grid_height=grid::unit(4 * legend_grid_cex, "mm"),
                  grid_width=grid::unit(4 * legend_grid_cex, "mm"),
                  title_gp=legend_title_gp,
                  labels_gp=legend_labels_gp)
            }
         }));
      top_annotation <- ComplexHeatmap::HeatmapAnnotation(
         border=TRUE,
         df=top_df,
         annotation_name_gp=top_annotation_name_gp,
         annotation_legend_param=top_param_list,
         simple_anno_size=simple_anno_size,
         show_legend=show_top_legend,
         show_annotation_name=show_top_annotation_name,
         col=top_color_list);
   }

   # left_annotation
   if (length(left_annotation) == 0 && !correlation) {
      row_anno_fontsize <- jamba::noiseFloor(
         column_cex * row_anno_fontsize,
         ceiling=24,
         minimum=2);
      left_anno_list <- list();
      left_color_list <- list();
      left_param_list <- list();
      if (length(show_left_legend) == 0) {
         show_left_legend <- TRUE;
      }
      if (length(show_left_legend) < 2) {
         show_left_legend <- rep(show_left_legend,
            length.out=2);
      }
      show_left_legend_v <- logical(0);
      # sestats annotations
      if (length(sestats) > 0) {
         if (verbose) {
            jamba::printDebug("heatmap_se(): ",
               "Preparing sestats incidence matrix left_annotation.");
         }
         show_left_legend_v <- show_left_legend[1];
         show_left_legend <- tail(show_left_legend, -1);
         left_anno_list <- c(list(
            hits=gene_hits_im[gene_hits, , drop=FALSE]),
            left_anno_list);
         left_color_list <- c(list(
            hits=colorjam::col_div_xf(1.5)),
            left_color_list);
         left_param_list <- c(list(
            hits=list(
               at=c(-1, 0, 1),
               color_bar="discrete",
               title_gp=legend_title_gp,
               labels_gp=legend_labels_gp,
               grid_height=grid::unit(4 * legend_grid_cex, "mm"),
               grid_width=grid::unit(4 * legend_grid_cex, "mm"),
               border=TRUE,
               labels=c("down", "no change", "up"))),
            left_param_list);
      }
      # alt_sestats annotations
      if (length(alt_sestats) > 0) {
         if (verbose) {
            jamba::printDebug("heatmap_se(): ",
               "Preparing alt_sestats incidence matrix left_annotation.");
         }
         show_left_legend_v <- c(FALSE,
            show_left_legend_v);
         left_anno_list <- c(list(
            hits_alt=alt_gene_hits_im[gene_hits, , drop=FALSE]),
            left_anno_list);
         left_color_list <- c(list(
            hits_alt=colorjam::col_div_xf(1.5)),
            left_color_list);
         left_param_list <- c(list(
            hits_alt=list(
               at=c(-1, 0, 1),
               color_bar="discrete",
               title_gp=legend_title_gp,
               labels_gp=legend_labels_gp,
               grid_height=grid::unit(4 * legend_grid_cex, "mm"),
               grid_width=grid::unit(4 * legend_grid_cex, "mm"),
               border=TRUE,
               labels=c("down", "no change", "up"))),
            left_param_list);
      }
      # rowData annotations
      if (length(rowData_colnames) > 0) {
         if (verbose) {
            jamba::printDebug("heatmap_se(): ",
               "preparing left_annotation for rowData_colnames: ",
               rowData_colnames);
         }
         show_left_legend_v <- c(show_left_legend_v,
            rep(show_left_legend, length.out=length(rowData_colnames)))
         # subset any factor columns to limit colors shown in the legend
         left_df <- data.frame(check.names=FALSE,
            rowData_se[gene_hits, rowData_colnames, drop=FALSE]);
         for (rowData_colname in rowData_colnames) {
            if (is.factor(left_df[[rowData_colname]])) {
               left_df[[rowData_colname]] <- factor(left_df[[rowData_colname]]);
            }
         }
         # put it together
         left_anno_list <- c(list(
            df=left_df),
            left_anno_list);
         use_color_list_names <- intersect(rowData_colnames,
            names(sample_color_list));
         leftanno_color_list <- NULL;
         if (length(use_color_list_names) > 0) {
            if (subset_legend_colors) {
               leftanno_color_list <- lapply(jamba::nameVector(use_color_list_names), function(use_color_list_name){
                  sample_colors <- sample_color_list[[use_color_list_name]];
                  if (!is.function(sample_colors)) {
                     if (is.factor(left_df[gene_hits, use_color_list_name])) {
                        # for factors we honor the order of factor levels
                        uniq_values <- levels(left_df[gene_hits, use_color_list_name]);
                     } else {
                        # note we add mixedSort() so they are alphanumerical
                        uniq_values <- jamba::mixedSort(unique(
                           as.character(
                              left_df[gene_hits, use_color_list_name])));
                     }
                     sample_colors <- sample_colors[uniq_values];
                     names(sample_colors) <- uniq_values;
                     if (any(is.na(sample_colors))) {
                        # fallback plan for missing values is to assign
                        # generic rainbow categorical colors
                        sample_colors[is.na(sample_colors)] <- colorjam::rainbowJam(
                           n=sum(is.na(sample_colors)),
                           ...);
                     }
                  }
                  sample_colors;
               })
            } else {
               leftanno_color_list <- sample_color_list[use_color_list_names];
            }
            left_color_list <- c(
               jamba::rmNULL(
                  leftanno_color_list),
               left_color_list);
            if (debug) {
               jamba::printDebug("heatmap_se(): ",
                  "left_color_list:");
               print(left_color_list[!jamba::sclass(left_color_list) %in% "function"]);
            }
         }
         left_param_list <- c(
            lapply(jamba::nameVector(rowData_colnames), function(iname){
               if (iname %in% names(left_color_list)) {
                  if (is.function(left_color_list[[iname]])) {
                     if ("breaks" %in% names(attributes(left_color_list[[iname]]))) {
                        list(border=TRUE,
                           color_bar="discrete",
                           title_gp=legend_title_gp,
                           labels_gp=legend_labels_gp,
                           grid_height=grid::unit(4 * legend_grid_cex, "mm"),
                           grid_width=grid::unit(4 * legend_grid_cex, "mm"),
                           at=attr(left_color_list[[iname]], "breaks"))
                     } else {
                        list(border=TRUE,
                           title_gp=legend_title_gp,
                           labels_gp=legend_labels_gp,
                           grid_height=grid::unit(4 * legend_grid_cex, "mm"),
                           grid_width=grid::unit(4 * legend_grid_cex, "mm"),
                           color_bar="discrete")
                     }
                  } else {
                     if (subset_legend_colors) {
                        list(border=TRUE,
                           title_gp=legend_title_gp,
                           labels_gp=legend_labels_gp,
                           grid_height=grid::unit(4 * legend_grid_cex, "mm"),
                           grid_width=grid::unit(4 * legend_grid_cex, "mm"),
                           color_bar="discrete",
                           at=jamba::rmNA(names(left_color_list[[iname]])))
                     } else {
                        list(border=TRUE,
                           grid_height=grid::unit(4 * legend_grid_cex, "mm"),
                           grid_width=grid::unit(4 * legend_grid_cex, "mm"),
                           title_gp=legend_title_gp,
                           labels_gp=legend_labels_gp)
                     }
                  }
               } else {
                  list(border=TRUE,
                     grid_height=grid::unit(4 * legend_grid_cex, "mm"),
                     grid_width=grid::unit(4 * legend_grid_cex, "mm"),
                     title_gp=legend_title_gp,
                     labels_gp=legend_labels_gp)
               }
            }),
            left_param_list);
      }

      # put it all together
      if (length(left_anno_list) > 0) {
         if (length(left_annotation_name_gp) == 0) {
            left_annotation_name_gp <- grid::gpar(fontsize=row_anno_fontsize);
         }
         left_alist <- alist(
            simple_anno_size=simple_anno_size,
            col=left_color_list,
            annotation_legend_param=left_param_list,
            show_legend=show_left_legend_v,
            show_annotation_name=show_left_annotation_name,
            annotation_name_rot=left_annotation_name_rot,
            annotation_name_gp=left_annotation_name_gp,
            border=TRUE);
         if (debug > 1) {
            jamba::printDebug("heatmap_se(): ",
               "left_alist:");
            print(sdim(left_alist));
            print(left_alist);
            jamba::printDebug("heatmap_se(): ",
               "left_anno_list:");
            print(sdim(left_anno_list));
            print(left_anno_list);
         }
         left_arglist <- c(
            left_alist,
            left_anno_list);
         left_annotation <- do.call(ComplexHeatmap::rowAnnotation,
            left_arglist);
         if (debug > 1) {
            print(sdim(left_annotation@anno_list))
            for (i in left_annotation@anno_list){
               print(i@show_legend)
            }
         }
      }
   }

   # optional row_split
   if (length(row_split) > 0) {
      if (correlation) {
         # correlation uses colData for split
         if (any(c("factor", "character") %in% class(row_split))) {
            if (!any(duplicated(row_split)) &&
                  all(row_split %in% colnames(colData_se))) {
               row_split <- data.frame(check.names=FALSE,
                  colData_se[isamples, row_split, drop=FALSE]);
            } else if (all(names(row_split) %in% isamples)) {
               row_split <- row_split[isamples];
            } else if (length(row_split) == length(isamples)) {
               names(row_split) <- isamples;
            } else {
               row_split <- NULL;
            }
         } else if (length(row_split) == 1 && is.numeric(row_split)) {
            # leave as-is
         } else {
            row_split <- NULL;
         }
      } else {
         # non-correlation uses rowData for split
         if (any(c("factor", "character") %in% class(row_split))) {
            if (!any(duplicated(row_split)) &&
                  all(row_split %in% colnames(rowData_se))) {
               row_split <- data.frame(check.names=FALSE,
                  rowData_se[gene_hits, row_split, drop=FALSE]);
            } else if (all(names(row_split) %in% gene_hits)) {
               row_split <- row_split[gene_hits];
            } else {
               row_split <- NULL;
            }
         } else if (length(row_split) == 1 && is.numeric(row_split)) {
            # leave as-is
         } else {
            row_split <- NULL;
         }
      }
   } else if (correlation) {
      row_split <- column_split;
   }

   assay_name <- head(intersect(assay_name,
      names(assays(se))), 1);
   if (length(assay_name) == 0) {
      stop("assay_name must be supplied.")
   }
   if (verbose) {
      jamba::printDebug("heatmap_se(): ",
         "using assay_name '", assay_name, "'");
   }
   norm_label <- paste0(assay_name, " ", data_type);

   if (any(centerby_colnames %in% FALSE)) {
      centerby_colnames <- NULL;
      centerGroups <- FALSE;
      centerby_label <- "";
   } else {
      centerby_colnames <- intersect(centerby_colnames,
         colnames(colData_se));
      if (length(centerby_colnames) > 0) {
         centerby_label <- paste0("centered by ",
            jamba::cPaste(centerby_colnames,
               sep="/"));
         centerGroups <- jamba::pasteByRow(
            colData_se[,centerby_colnames]);
      } else {
         centerGroups <- NULL;
         centerby_label <- "global-centered";
      }
      if (length(control_label) > 0 && any(nchar(control_label)) > 0) {
         centerby_label <- paste(centerby_label,
            control_label);
      }
   }

   # row_labels
   if (correlation) {
      # for correlation, samples are shown on rows
      # so it must use colData_se
      if (length(row_label_colname) == 0 ||
            !all(row_label_colname %in% colnames(colData_se))) {
         row_labels <- isamples;
      } else if (length(row_label_colname) > 1) {
         row_labels <- jamba::pasteByRow(
            colData_se[isamples, row_label_colname, drop=FALSE],
            sep=",");
      } else {
         row_labels <- colData_se[isamples, , drop=FALSE][[row_label_colname]];
      }
   } else {
      # for expression data, rowData_se must be used
      if (length(row_label_colname) == 0) {
         row_labels <- gene_hits;
      } else if (length(row_label_colname) > 1) {
         row_labels <- jamba::pasteByRow(
            rowData_se[gene_hits, row_label_colname, drop=FALSE],
            sep=",");
      } else {
         row_labels <- rowData_se[gene_hits, , drop=FALSE][[row_label_colname]];
      }
   }
   if (length(show_row_names) == 0) {
      show_row_names <- (length(gene_hits) <= 500);
   }

   # heatmap legend labels
   if (length(legend_at) == 0) {
      if (correlation) {
         if (abs(color_max) > 1) {
            color_max <- 1;
         }
         legend_at <- seq(-ceiling(color_max),
            to=ceiling(color_max),
            by=0.25);
      } else {
         legend_at <- seq(-ceiling(color_max),
            to=ceiling(color_max));
      }
   }
   if (length(legend_labels) != length(legend_at)) {
      if (correlation) {
         legend_labels <- legend_at;
      } else {
         legend_labels <- round(jamba::exp2signed(legend_at+0.001, offset=0))
         if (any(duplicated(legend_labels))) {
            legend_labels <- round(10 * jamba::exp2signed(legend_at+0.00001, offset=0)) / 10;
         }
      }
   }

   # pull assay data separately so we can tolerate other object types
   # Note columns are not subset here so they can be used during centering.
   # After centering, isamples is used to subset columns as needed.
   if (grepl("SummarizedExperiment", ignore.case=TRUE, class(se))) {
      se_matrix <- assays(se[gene_hits, ])[[assay_name]];
   } else {
      se_matrix <- assayData(se[gene_hits, ])[[assay_name]];
   }

   # cluster_columns
   if (!is.function(cluster_columns) && cluster_columns %in% TRUE) {
      cluster_columns <- function(x, ...) {
         amap::hcluster(rmNA(naValue=0, x),
            ...,
            method="euclidean",
            link="ward")}
   }

   # define heatmap matrix
   # After centering, columns are subset using isamples.
   if (length(centerGroups) > 0 && any(centerGroups %in% FALSE)) {
      if (verbose) {
         jamba::printDebug("heatmap_se(): ",
            "Not centering data.");
      }
      se_matrix <- se_matrix[, isamples, drop=FALSE];
      hm_name <- data_type;
   } else {
      if (verbose) {
         jamba::printDebug("heatmap_se(): ",
            "Centering data.");
      }
      se_matrix <- jamma::centerGeneData(
         useMedian=useMedian,
         centerGroups=centerGroups,
         x=se_matrix,
         controlSamples=controlSamples,
         ...)[, isamples, drop=FALSE];
      hm_name <- paste0("centered\n", data_type);
   }
   if (correlation) {
      # call correlation function cor()
      se_matrix <- jamba::call_fn_ellipsis(cor,
         x=se_matrix,
         use="pairwise.complete.obs",
         ...);
      cluster_rows <- cluster_columns;
      if (length(centerGroups) > 0 && any(centerGroups %in% FALSE)) {
         hm_name <- paste0("correlation of\n", data_type);
      } else {
         hm_name <- paste0("correlation of\ncentered\n", data_type);
      }
   }

   # optional mark_rows
   mark_rows <- intersect(mark_rows, gene_hits);
   if (length(mark_rows) > 0) {
      mark_at <- match(mark_rows, rownames(se_matrix));
      right_annotation_mark <- ComplexHeatmap::rowAnnotation(
         foo=ComplexHeatmap::anno_mark(
            at=mark_at,
            labels=row_labels[mark_at],
            labels_gp=mark_labels_gp))
      if (length(right_annotation) == 0) {
         right_annotation <- right_annotation_mark;
      } else {
         right_annotation <- right_annotation + right_annotation_mark;
      }
   }

   # pre-calculate row clusters
   # This step is required to enable row_split as integer number of clusters,
   # which is not accepted when supplying a function.
   # This step does not work with character or data.frame row_split
   if (length(row_split) == 1 &&
         is.numeric(row_split) &&
         is.function(cluster_rows)) {
      cluster_rows <- cluster_rows(se_matrix);
   }
   if (length(column_split) == 1 &&
         is.numeric(column_split) &&
         is.function(cluster_columns)) {
      cluster_columns <- cluster_columns(se_matrix);
   }

   # optional customization of row and column names gp
   if (length(row_names_gp) == 0) {
      row_names_gp <- grid::gpar(fontsize=row_fontsize)
   } else {
      if (!"fontsize" %in% names(row_names_gp)) {
         row_names_gp$fontsize <- grid::gpar(fontsize=row_fontsize)$fontsize;
      }
   }
   if (length(column_names_gp) == 0) {
      column_names_gp <- grid::gpar(fontsize=column_fontsize)
   } else {
      if (!"fontsize" %in% names(column_names_gp)) {
         column_names_gp$fontsize <- grid::gpar(fontsize=column_fontsize)$fontsize;
      }
   }


   # define heatmap
   hm_hits <- jamba::call_fn_ellipsis(ComplexHeatmap::Heatmap,
      matrix=se_matrix,
      use_raster=TRUE,
      top_annotation=top_annotation,
      left_annotation=left_annotation,
      right_annotation=right_annotation,
      heatmap_legend_param=list(
         border=TRUE,
         color_bar="discrete",
         at=legend_at,
         labels=legend_labels,
         grid_height=grid::unit(4 * legend_grid_cex, "mm"),
         grid_width=grid::unit(4 * legend_grid_cex, "mm"),
         title_gp=legend_title_gp,
         labels_gp=legend_labels_gp
      ),
      clustering_method_rows="ward.D",
      column_split=column_split,
      row_split=row_split,
      row_title_rot=row_title_rot,
      cluster_column_slices=FALSE,
      border=TRUE,
      name=hm_name,
      show_row_names=show_row_names,
      show_row_dend=show_row_dend,
      show_heatmap_legend=show_heatmap_legend,
      row_labels=row_labels,
      row_names_gp=row_names_gp,
      column_names_gp=column_names_gp,
      col=colorjam::col_div_xf(color_max,
         lens=lens,
         ...),
      cluster_columns=cluster_columns,
      cluster_rows=cluster_rows,
      ...)
   hm_title <- paste0(
      formatInt(length(gene_hits)),
      " ", row_type,
      "\n", norm_label,
      ifelse(any(nchar(centerby_label) > 0),
         paste0(",\n", centerby_label),
         ""))
   attr(hm_hits, "hm_title") <- hm_title;
   if (debug) {
      ret_list <- list(
         hm=hm_hits,
         top_annotation=top_annotation,
         left_annotation=left_annotation,
         hm_title=hm_title
      );
      ret_list$gene_hits_im <- gene_hits_im;
      ret_list$alt_gene_hits_im <- alt_gene_hits_im;
      ret_list$gene_hitlist <- gene_hitlist;
      ret_list$alt_gene_hitlist <- alt_gene_hitlist;
      return(ret_list)
   }
   hm_hits
}

