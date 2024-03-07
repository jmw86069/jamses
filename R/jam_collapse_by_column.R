
#' Collapse SummarizedExperiment data by column
#'
#' Collapse SummarizedExperiment data by column
#'
#' Purpose is to collapse columns of a `SummarizedExperiment` object,
#' where measurements for a given entity, usually a gene, are split
#' across multiple rows in the source data. The output of this function
#' should be measurements appropriately summarized to the gene level.
#'
#' The driving use case is slightly different than with `se_collapse_by_row()`,
#' in this case the function is mostly convenient method to calculate
#' group mean values in context of a `SummarizedExperiment` object,
#' so it can be used with `jamses::heatmap_se()` for example.
#'
#' This function retains associated column annotations `colData(se)`,
#' after combining multiple values in an appropriate manner.
#'
#' Optionally, this function will detect and remove individual outlier
#' values before calculating the group mean.
#'
#' @return `SummarizedExperiment` object with these changes:
#'    * columns will be collapsed by `column_groups`, for each `assays(se)`
#'    `numeric` matrix defined by `assay_names`.
#'    * `colData(se)` will also be collapsed by `shrinkDataFrame()` to
#'    combine unique values from each column annotation.
#'
#' @family jamses utilities
#'
#' @param se `SummarizedExperiment` object
#' @param columns `character` vector of `colnames(se)` to include in the
#'    process.
#' @param column_groups `character` vector of column groupings, or
#'    `character` vector of `colnames(colData(se))` used to define the
#'    column groupings.
#' @param assay_names `character` vector with one or more `assayNames(se)`
#'    to apply the column grouping calculation defined in `groupFunc`.
#'    By default, all assay names in `assayNames(se)` are used.
#' @param colDataColnames `character` vector of `colData(se)` colnames to
#'    be included in the returned `SummarizedExperiment` after the column
#'    grouping. This argument is used to subset the columns, in cases where
#'    some columns do not need to be combined and returned in the output
#'    data.
#' @param keepNULLlevels `logical` indicating whether to return empty columns
#'    when there are not factor levels present in the data. This option is
#'    intended when `column_group` references a `factor` type, whose factor
#'    levels are not present in the current data, using `columns`. When
#'    `keepNULLlevels=TRUE` any missing levels will be present with `NA`
#'    values, which can be helpful for generating a consistent output.
#' @param groupFunc `function` used to perform row group calculations on a
#'    `numeric` matrix. The default is passed to `jamba::rowGroupMeans()`,
#'    but can be substituted with another row-based function.
#'    It must accept arguments `x` and `groups`, but the other arguments
#'    are passed only if `groupFunc` permits these argument names, or `...`:
#'    * `x` as a `numeric` matrix (required),
#'    * `groups` as a `character` vector of column groups, in order of
#'    `colnames(x)` (required)
#'    * `rmOutliers` a `logical` indicating whether to apply outlier removal,
#'    though the function can ignore this value (optional).
#'    * `madFactor` a `numeric` value indicating the MAD threshold used when
#'    `rmOutliers=TRUE`; though again the function can ignore this value
#'    (optional).
#'    * `useMedian=FALSE` is `logical` and when `useMedian=FALSE` it disables
#'    calculating the `median()` value per group, and instead takes the
#'    group `mean()` value.
#'    * `...` additional arguments in `...` will be passed only if permitted
#'    by `groupFunc`.
#' @param noise_floor `numeric` value indicating the minimum numeric value
#'    permitted, *at or below* this value will be replaced with
#'    `noise_floor_value`.
#'    The default value `noise_floor=0` will therefore change all values
#'    at or below zero to `noise_floor_value=0` by default.
#'    Another alternative is to change abnormally low
#'    values such as zero `0` to `NA` so these values are not
#'    treated as actual measurements during the group summary calculation.
#'    This value and the replacement should be adjusted
#'    with caution. Use `noise_floor=NULL` or `noise_floor=-Inf` to disable
#'    this step.
#' @param noise_floor_value `numeric` or `NA` used as a replacement for
#'    `numeric` values *at or below* `noise_floor`, which occurs prior to
#'    calling the `groupFunc` summary calculation.
#' @param rmOutliers,madFactor `logical` and `numeric`, respectively, passed
#'    to `groupFunc` which by default is `jamba::rowGroupMeans()`.
#' @param useMedian `logical` passed to argument `groupFunc()`, intended
#'    to be used by `jamba::rowGroupMeans()` to specify taking the mean
#'    and not the median value per row group.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed through `groupFunc`.
#'
#' @export
se_collapse_by_column <- function
(se,
 columns=colnames(se),
 column_groups,
 assay_names=NULL,
 colDataColnames=colnames(SummarizedExperiment::colData(se)),
 keepNULLlevels=FALSE,
 groupFunc=jamba::rowGroupMeans,
 noise_floor=0,
 noise_floor_value=0,
 rmOutliers=FALSE,
 madFactor=5,
 useMedian=FALSE,
 verbose=FALSE,
 ...)
{
   # Purpose is to collapse columns of a SummarizedExperiment object, intended
   # for technical replicates where an individual sample measurements are
   # represented in multiple columns.
   #
   # if (length(column_groups) == 1) {
   #    # ?
   # }
   ## Combine technical replicates, using MAD factor cutoff of 5 for outliers

   # validate assay_names
   if (length(assay_names) == 0) {
      assay_names <- SummarizedExperiment::assayNames(se);
   } else {
      assay_names <- intersect(assay_names,
         SummarizedExperiment::assayNames(se));
      if (length(assay_names) == 0) {
         stop("'assay_names' did not match any assayNames(se), please remedy or supply assay_names=NULL.")
      }
   }

   # confirm colnames(se) are defined
   if (length(colnames(se)) == 0) {
      colnames(se) <- paste("column_", jamba::padInteger(seq_len(ncol(se))))
   }

   # validate columns
   if (!all(columns %in% colnames(se))) {
      stop("'columns' must all be present in colnames(se).")
   }

   # validate column_groups
   if (length(column_groups) == 0) {
      stop("'column_groups' is required.")
   }
   if (all(column_groups %in% colnames(SummarizedExperiment::colData(se)))) {
      column_groups <- jamba::pasteByRowOrdered(
         data.frame(check.names=FALSE,
            SummarizedExperiment::colData(
               se[, columns, drop=FALSE])[, column_groups, drop=FALSE]),
         sep="_");
      names(column_groups) <- columns;
   } else if (length(column_groups) != length(columns)) {
      stop("'column_groups' must equal length(columns) or contain colnames(colData(se)).")
   }

   #################################################################
   # Iterate each assay_name and perform the collapse
   # Any assay_names not included will be dropped from the resulting object
   assays_grouped_list <- lapply(jamba::nameVector(assay_names), function(assay_name){
      if (verbose) {
         jamba::printDebug("se_collapse_by_column(): ",
            "Collapsing assay_name:",
            assay_name);
      }
      iMatrix <- SummarizedExperiment::assays(
         se[,columns, drop=FALSE])[[assay_name]];
      if (length(noise_floor) > 0) {
         noise_match <- (!is.na(iMatrix) & iMatrix <= noise_floor);
         if (any(noise_match)) {
            iMatrix[noise_match] <- noise_floor_value;
         }
      }
      iMatrixColGrp <- jamba::call_fn_ellipsis(groupFunc,
         iMatrix,
         rmOutliers=rmOutliers,
         madFactor=madFactor,
         useMedian=useMedian,
         groups=column_groups);

      ## Revert zero back to NA?
      # 29nov2022 - not necessary in favor of using noise_floor
      # iMatrixColGrp[iMatrixColGrp == 0] <- NA;

      iMatrixColGrp;
   });

   #################################################################
   ## Now try to be clever and create a new colData() which indicates
   ## the column groupings
   if (verbose) {
      jamba::printDebug("se_collapse_by_column(): ",
         "Collapsing colData(se)");
   }
   colDataShrunk <- shrinkDataFrame(
      x=SummarizedExperiment::colData(se[,columns])[, colDataColnames, drop=FALSE],
      groupBy=column_groups,
      includeNumReps=TRUE,
      verbose=FALSE,
      ...);

   # re-order rows to match levels(column_groups)
   colDataShrunk <- colDataShrunk[colnames(assays_grouped_list[[1]]), ,
      drop=FALSE]

   if (verbose) {
      jamba::printDebug("se_collapse_by_column(): ",
         "head(colDataShrunk):");
      print(head(colDataShrunk));
      jamba::printDebug("se_collapse_by_column(): ",
         "sdim(assays_grouped_list):");
      print(jamba::sdim(assays_grouped_list));
      jamba::printDebug("se_collapse_by_column(): ",
         "head(assays_grouped_list[[1]]):");
      print(head(assays_grouped_list[[1]]));
   }

   se_shrunk <- SummarizedExperiment::SummarizedExperiment(
      assays=assays_grouped_list,
      rowData=SummarizedExperiment::rowData(se),
      colData=colDataShrunk,
      metadata=list(
         genes=rownames(se),
         samples=rownames(colDataShrunk)));
   #colnames(SEshrunk) <- rownames(colDataShrunk);
   se_shrunk;
}
