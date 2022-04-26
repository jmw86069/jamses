
#' SummarizedExperiment heuristics to define detected rows
#'
#' SummarizedExperiment heuristics to define detected rows
#'
#' This function is intended to help apply common logical rules
#' to define valid, "detected" rows for downstream analysis.
#'
#' The rules:
#'
#' * minimum value at or above which a measurement is "valid"
#' * minimum total replicates with "valid" measurement, across all sample columns
#' * minimum replicates with "valid" measurement required in any sample group
#' * minimum percent replicates with "valid" measurement required in any sample group
#' * minimum sample groups with "valid" criteria above required
#'
#' ## Example
#'
#' Consider an experiment with 7 groups, and n=3 replicates, which
#' contains 21 total samples.
#'
#' Assume one row of data that contains 6 "valid" measurements.
#' * If these 6 "valid" measurements are found in only 2 groups,
#' both groups contain n=3 "valid" measurements. This row may have
#' sufficient data for a statistical comparison across these two groups.
#' * However, if the 6 "valid" measurements are also found across
#' 6 different groups, it may not be suitable for statistical testing.
#'
#' ## Use of normgroups
#'
#' Detection can be carried out within `"normgroups"`, which are
#' independent subsets of sample columns. In most cases this method
#' is not necessary, but is intended when the detected rows should
#' be independently calculated for two or more subsets of sample
#' columns.
#'
#' A specific example might be an experiment that measures treatment
#' effects in two very different tissue types, like lung and muscle.
#' The detected genes in lung may well not be the same as detected
#' genes in lung. And in fact, statistical comparisons may not be intended
#' to compare muscle and lung directly. (That judgement is left
#' to the analyst.) One may define a column in `colData(se)` that represents
#' tissue type, with values `"muscle"`, and `"lung"`, then define
#' this column with argument `normgroup_colname`. The detection will
#' be done within each independent normgroup, returned as a `list`
#' named `"detected_normgroup"`. The detected rows are also combined
#' into `"detected_rows"` which returns rows detected across
#' all normgroups.
#'
#'
#' @return `list` with the following elements:
#' * `detected_rows` is a `character` vector of detected `rownames(se)`
#' * `detected_normgroup` is a `list` of `logical` vectors for each normgroup,
#' where the vectors encode whether a row is detected within each normgroup.
#' * `detected_df` is a `data.frame` with summary information for each
#' normgroup.
#'
#' @param se `SummarizedExperiment` object
#' @param assay_name `character` or `integer` index, referring to the
#'    entry in `assays(se)` to use when determining valid measurements.
#' @param group_colnames `character` vector of colnames in `colData(se)` which
#'    defines sample grouping.
#' @param normgroup_colname `character` string with optional colname
#'    in `colData(se)` to use for normgroups.
#' @param detect_mincounts `numeric` value at or above which a measurement
#'    is considered "valid".
#' @param detect_totalreps `numeric` minimum total number of replicates
#'    that must contain "valid" measurements.
#' @param detect_minreps `numeric` minimum replicates which must contain
#'    "valid" measurements in any given sample group.
#' @param detect_minpct `numeric` minimum fraction of available replicates
#'    in a sample group that must contain "valid" measurements.
#' @param detect_mingroups `numeric` minimum number of sample groups that
#'    are considered "valid" based upon other criteria.
#' @param isamples `character` optional vector of `colnames(se)` to use
#'    during this analysis. This vector is useful for example, when excluding
#'    outlier samples that were defined by other methods.
#' @param ... additional arguments are ignored.
#'
#'
#' @export
se_detected_rows <- function
(se,
 assay_name=1,
 group_colnames,
 normgroup_colname=NULL,
 detect_mincounts=0,
 detect_totalreps=1,
 detect_minreps=2,
 detect_minpct=0.65,
 detect_mingroups=1,
 isamples=colnames(se),
 ...)
{
   #
   if (length(normgroup_colname) > 0) {
      isamples_normgroup <- split(isamples,
         SummarizedExperiment::colData(se[,isamples])[[normgroup_colname]]);
   } else {
      isamples_normgroup <- list(overall=isamples);
   }

   #
   detect_normgroup <- lapply(isamples_normgroup, function(isamples_i){
      # incidence matrix of detected values
      detect_im <- (
         !is.na(SummarizedExperiment::assays(se[,isamples_i])[[assay_name]]) &
            SummarizedExperiment::assays(se[,isamples_i])[[assay_name]] > detect_mincounts) * 1;

      # define groups for this normgroup
      igroups_i <- jamba::pasteByRow(sep="_",
         SummarizedExperiment::colData(se[,isamples_i])[,group_colnames]);

      # sum detected points by group
      se_Num_Samples_Detected <- rowSums(detect_im, na.rm=TRUE);
      se_Num_Group <- jamba::rowGroupMeans(detect_im,
         groups=igroups_i,
         rowStatsFunc=rowSums);

      # percent detected per group
      reps_per_group <- tcount(igroups_i);
      # first determine replicates per group
      se_Pct_Group <- se_Num_Group / reps_per_group[colnames(se_Num_Group)];
      # se_Pct_Group <- se_Num_Group /
      #    rep(apply(head(se_Num_Group, 10000), 2, function(i){
      #          max(i, na.rm=TRUE)
      #       }),
      #       each=nrow(detect_im));

      # apply all criteria to each row
      # at least mingroups groups with at least minreps valid reps
      # at least mingroups groups with at least minpct valid reps
      im_criteria <- (
         se_Num_Group >= detect_minreps &
            se_Pct_Group >= detect_minpct) * 1;
      valid_rows <- (rowSums(im_criteria, na.rm=TRUE) >= detect_mingroups &
            se_Num_Samples_Detected >= detect_totalreps);
      jamba::nameVector(valid_rows,
         rownames(detect_im));
   });

   detect_t <- jamba::renameColumn(
      data.frame(check.names=FALSE,
         t(sapply(detect_normgroup, table))),
      from=c("FALSE", "TRUE"),
      to=c("not detected", "detected"));
   detect_df <- data.frame(check.names=FALSE,
      normgroup=rownames(detect_t),
      detect_t)

   # detected_genes
   detected_rows <- names(which(Reduce("|", detect_normgroup)));
   return(list(
      detected_rows=detected_rows,
      detected_normgroup=detect_normgroup,
      detect_df=detect_df))
}
