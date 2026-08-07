
#' Merge stats data.frame from two sestats results
#'
#' Merge stats data.frame from two sestats results, specifically
#' when one object has "all rows" and is intended only to produce
#' the non-P-value statistics such as group mean, max group mean (mgm),
#' fold, and logFC.
#'
#' * `sestats_all` should have all rows, not just those tested.
#' It will specifically exclude stats columns with `'P.val'` or
#' `'^hit '` in the name, using argument `exclude_from_all`.
#' * `sestats_test` should have only the subset of rows tested.
#'
#' This function may be used internally within jamses to automate
#' the job of testing a subset of rows, but reporting summary
#' values for all rows - including those rows not formally tested.
#'
#' @return `list` object in the same format as `sestats_test`, with only
#' the elements 'stats_df' and 'stats_dfs' replaced with the expanded
#' complete set of data.
#'
#' @family jamses utilities
#'
#' @param sestats_all `list` with sestats data from `se_contrast_stats()`
#'    which used "all rows" of the source data, along with `define_hits=FALSE`.
#' @param sestats_test `list` with sestats data from `se_contrast_stats()`
#'    which used a subset of rows, usually "detected rows" for the statistical
#'    test.
#' @param exclude_from_all `character` pattern used with `jamba::unvigrep()`
#'    to exclude colnames from `sestats_all`. Default `'P.val|^hit '`
#'    excludes P-values and hit columns.
#' @param ... additional arguments are ignored.
#'
#' @export
merge_statdf_all_test <- function(
   sestats_all,
   sestats_test,
   exclude_from_all = "P.val|^hit ",
   ...
) {
   # iterate stats_dfs
   # each assay_name
   # validate SEStats
   if (!inherits(sestats_all, "SEStats")) {
      if (is.list(sestats_all)) {
         sestats_all <- SEStats(stats_dfs = sestats_all[["stats_dfs"]])
      } else {
         cli::cli_abort(paste0(
            "{.var sestats_all} must be {.cls SEStats} or ",
            "a {.cls list} with {.var stats_dfs}."
         ))
      }
   }
   # validate SEStats
   if (!inherits(sestats_test, "SEStats")) {
      if (is.list(sestats_test)) {
         if ("sedesign" %in% names(sestats_test)) {
            sedesign_test <- sestats_test[["sedesign"]]
         } else {
            sedesign_test <- SEDesign(
               design=sestats_test[["design"]],
               contrasts=sestats_test[["contrasts"]])
         }
         sestats_test <- SEStats(
            sedesign = sedesign_test,
            hit_array = sestats_test[["hit_array"]],
            stats_dfs = sestats_test[["stats_dfs"]],
            metadata = sestats_test[["metadata"]]
         )
      } else {
         cli::cli_abort(paste0(
            "{.var sestats_test} must be {.cls SEStats} or ",
            "a {.cls list} with {.var stats_dfs}."
         ))
      }
   }

   for (iassay_name in names(sestats_all@stats_dfs)) {
      # each contrast_name
      for (icontrast_name in names(sestats_all@stats_dfs[[iassay_name]])) {
         stats_df_all <- sestats_all@stats_dfs[[iassay_name]][[icontrast_name]]
         stats_df_test <- sestats_test@stats_dfs[[iassay_name]][[
            icontrast_name
         ]]
         add_from_all <- setdiff(stats_df_all[[1]], stats_df_test[[1]])
         colnames_from_all <- colnames(stats_df_all)
         if (length(exclude_from_all) == 1) {
            colnames_from_all <- jamba::unvigrep(
               exclude_from_all,
               colnames(stats_df_all)
            )
         }
         stats_df_new <- jamba::mergeAllXY(
            stats_df_test,
            subset(
               stats_df_all[, colnames_from_all, drop = FALSE],
               stats_df_all[[1]] %in% add_from_all
            )
         )[, colnames(stats_df_test), drop = FALSE]
         rownames(stats_df_new) <- stats_df_new[[1]]
         # add back into sestats_test
         sestats_test@stats_dfs[[iassay_name]][[icontrast_name]] <- stats_df_new
      }
   }

   # iterate stats_df
   # each assay_name
   for (iassay_name in names(sestats_all@stats_df)) {
      stats_df_all <- sestats_all@stats_df[[iassay_name]]
      stats_df_test <- sestats_test@stats_df[[iassay_name]]
      add_from_all <- setdiff(stats_df_all[[1]], stats_df_test[[1]])
      colnames_from_all <- colnames(stats_df_all)
      if (length(exclude_from_all) == 1) {
         colnames_from_all <- jamba::unvigrep(
            exclude_from_all,
            colnames(stats_df_all)
         )
      }
      stats_df_new <- jamba::mergeAllXY(
         stats_df_test,
         subset(
            stats_df_all[, colnames_from_all, drop = FALSE],
            stats_df_all[[1]] %in% add_from_all
         )
      )[, colnames(stats_df_test), drop = FALSE]
      rownames(stats_df_new) <- stats_df_new[[1]]
      # add back into sestats_test
      sestats_test@stats_df[[iassay_name]] <- stats_df_new
   }
   return(sestats_test)
}
