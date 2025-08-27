
#' Combine SEStats objects by contrast,signal
#'
#' Combine SEStats objects by contrast,signal
#'
#' This function concatenates two or more objects 'sestats' as
#' returned by `se_contrast_stats()`. It assumes that contrast:signal
#' are unique across all entries, however the same contrast can
#' be present if associated with different signal values.
#'
#' The cutoffs are not considered when determining which results can
#' be combined, in other words only the 'contrast' and 'signal' values
#' are used to define a "unique key". Only unique contrast:signal
#' combinations are valid as input.
#'
#' In other words, input will not accept multiple sestats objects
#' which contain the same contrast:signal and two different cutoffs.
#'
#' @family jamses utilities
#'
#' @return `list` equivalent to output from `se_contrast_stats()` with
#'    minor exceptions:
#'    * List elements 'idesign','icontrasts','normgroup' are each
#'    returned as a `list` with corresponding data from each sestats object.
#'    The reason is that these values cannot logically be combined into
#'    one new object, and therefore the next level of support is to
#'    retain the underlying data.
#'
#'    The 'hit_array' element will contain the unique set of values
#'    for each dimension: Cutoffs, Contrasts, Signal.
#'
#' @param sestats_list `list` with two or more `list` objects in sestats
#'    format, as returned by `se_contrast_stats()`.
#' @param create_stats_df `logical`, default TRUE, whether to create the
#'    'stats_df' `data.frame` object by merging corresponding
#'    'stats_dfs' `data.frame` objects together. This step assumes that
#'    columns of class "character" or "factor" will have identical values
#'    across all `data.frame` objects.
#'    Specifically, it assumes that at least one column has the same name,
#'    with row identifiers sufficient to merge each result.
#'    This assumption also implies that one should not combine sestats
#'    objects when the rows are not equivalent.
#' @param ... additional arguments are ignored.
#'
#' @export
combine_sestats <- function
(sestats_list,
 create_stats_df=TRUE,
 ...)
{
   # assign names if empty
   if (length(names(sestats_list)) == 0) {
      names(sestats_list) <- jamba::colNum2excelName(seq_along(sestats_list));
   }

   # get list of contrast names
   # dimnames(sestats_list[[2]]$hit_array)$Contrasts
   statlist_df <- jamba::rbindList(lapply(names(sestats_list), function(sename){
      sestats <- sestats_list[[sename]];
      ihit_array <- sestats$hit_array;
      idimnames <- dimnames(ihit_array);
      sedf <- jamba::rbindList(lapply(idimnames[[1]], function(icutoff){
         jamba::rbindList(lapply(idimnames[[2]], function(icontrast){
            jamba::rbindList(lapply(idimnames[[3]], function(isignal){
               ihits <- ihit_array[icutoff, icontrast, isignal][[1]];
               sedf <- data.frame(sename=sename,
                  signal=isignal,
                  contrast=icontrast,
                  cutoff=icutoff)
               if (!is.null(ihits) && is.numeric(ihits)) {
                  sedf
               } else {
                  sedf[0, , drop=FALSE]
               }
            }))
         }))
      }))
   }))

   # confirm no duplicated signal:contrast
   # Todo:
   # - consider ignoring duplicates
   # - consider renaming signal internally
   statlist_df$sc <- jamba::pasteByRow(statlist_df[, c("signal", "contrast")])
   statlist_check <- unique(statlist_df[, c("sename", "sc")]);
   if (any(duplicated(statlist_df$sc))) {
      cli::cli_abort(paste0("combine_sestats() cannot combine duplicated ",
         "{.var signal}:{.var contrast} from more than one {.var sestats} ",
         "object. Suggestion: Rename one {.var signal} and try again."))
      stop("signal:contrast cannot be present in multiple sestats.");
   }

   # create new dimnames for hit_array
   xcutoffs <- unique(statlist_df$cutoff)
   xcontrasts <- unique(statlist_df$contrast)
   xsignals <- unique(statlist_df$signal)

   # create new stats_dfs
   new_stats_dfs <- list();

   # create new hit_list
   new_hit_list <- list();

   # create new hit_array
   arrayDimnames <- list(
      Cutoffs=xcutoffs,
      Contrasts=xcontrasts,
      Signal=xsignals);
   arrayDim <- lengths(arrayDimnames);
   new_hit_array <- array(
      vector("list", prod(arrayDim)),
      dim=arrayDim,
      dimnames=arrayDimnames);
   for (irow in seq_len(nrow(statlist_df))) {
      isename <- statlist_df$sename[irow]
      isignal <- statlist_df$signal[irow]
      icontrast <- statlist_df$contrast[irow]
      icutoff <- statlist_df$cutoff[irow]
      ivalues <- (
         sestats_list[[isename]]$hit_array[icutoff, icontrast, isignal])
      new_hit_array[icutoff, icontrast, isignal][[1]] <- ivalues[[1]];
      # add to new_hit_list
      new_hit_list[[isignal]][[icontrast]][[icutoff]] <- ivalues[[1]];
      # add to new_stats_dfs
      new_stats_dfs[[isignal]][[icontrast]] <- (
         sestats_list[[isename]]$stats_dfs[[isignal]][[icontrast]]);
   }

   # create new stats_df
   # - note that some columns are repeated, and values differ by extremely
   #   small fraction, causing the merge to fail
   # - we take the first occurrence of each column
   new_stats_df <- lapply(new_stats_dfs, function(idfs){
      ucl <- lapply(idfs, colnames);
      acl <- unique(unlist(ucl));
      kcl <- character(0);
      for (iuc in seq_along(ucl)) {
         cc <- jamba::sclass(idfs[[iuc]])
         use_uc <- names(cc)[cc %in% c("character", "factor") |
            !names(cc) %in% kcl];
         ucl[[iuc]] <- use_uc;
         kcl <- unique(c(kcl, use_uc));
      }
      newdf <- jamba::mergeAllXY(lapply(seq_along(idfs), function(i){
         idfs[[i]][, ucl[[i]], drop=FALSE]
      }))
      # if (all(sort(colnames(newdf)) == sort(acl))) {
      #    newdf <- newdf[, acl, drop=FALSE]
      # }
      newdf
   })

   # Create new sestats list
   ret_list <- list();
   ret_list$stats_df <- new_stats_df;
   ret_list$stats_dfs <- new_stats_dfs;
   ret_list$hit_array <- new_hit_array;
   ret_list$hit_list <- new_hit_list

   # just concatenate idesign,icontrasts as list
   ret_list$idesign <- lapply(sestats_list, function(sestats){
      sestats$idesign
   })
   ret_list$icontrasts <- lapply(sestats_list, function(sestats){
      sestats$icontrasts
   })
   ret_list$normgroup <- lapply(sestats_list, function(sestats){
      sestats$normgroup
   })

   # Return the owl
   return(ret_list);

   # Todo:
   # - In future, refactor to use the new S4 SEStats object format.
   #
   # - Consider changing idesign,icontrasts to sedesign.
   # - Consider just using a list of sedesign objects to associate them.
   #
   # - Consider creating new merged form of stats_df

   # Assumptions / Rules:
   # - contrast,signal: cannot be duplicated in sestats_list
   # - cutoff: multiple can be encoded in one data.frame, which
   #      is why it is not included in the contrast,signal rule.
   #
   # Confirm no contrast is present more than once for the same signal.
   # Potential for contrast to be represented for two signals.
   #
   # Consider storing idesign,icontrast in metadata somewhere that
   # it can be reviewed later - but not done in a formal way.
   # (Do the easiest thing that works.)
}
