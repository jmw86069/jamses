
# jam_sestats.R
#
# placeholder for proper object type sestats

#' Convert sestats to table summary
#'
#' Convert sestats to table summary
#'
#' Note: This function is intended to provide a simple data.frame
#' summary with the number of hits for each contrast, signal, cutoff.
#' It is still being tested, and updated for usability.
#'
#' TODO: The order of `dimnames(hit_array)` should be user-customizable.
#' The series of dimnames in each lapply should use this order.
#'
#' @param sestats `list` output from `se_contrast_stats()`
#' @param style `character` string indicating what values to use:
#'    * `"text"`: number of hits (number up, number down)
#'    * `"integer"`: only the integer number of hits
#' @param rename_contrasts `logical` indicating whether to rename
#'    contrasts using `contrast2comp()`. The main benefit is reduction
#'    in length of the string that describes each contrast.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' if (FALSE) {
#' hitdf <- sestats_to_df(list(hit_array=hit_array));
#' hitdf_rowindex <- table(hitdf[[1]])[unique(hitdf[[1]])]
#' jamba::kable_coloring(
#'    hitdf[, -1, drop=FALSE],
#'    row.names=FALSE) %>%
#'    kableExtra::pack_rows(index=hitdf_rowindex)
#' }
#'
#' @export
sestats_to_df <- function
(sestats,
 style=c("text", "integer"),
 dimname_order=c(3, 1, 2),
 rename_contrasts=FALSE,
 ...)
{
   #
   style <- match.arg(style);

   #
   if ("list" %in% class(sestats) && "hit_array" %in% names(sestats)) {
      hit_array <- sestats$hit_array;
   } else {
      hit_array <- sestats;
   }

   #
   hit_array <- sestats$hit_array;
   idf <- jamba::rbindList(lapply(dimnames(hit_array)[[3]], function(icutoff){
      imatrix <- matrix(hit_array[, , icutoff],
         nrow=dim(hit_array)[1],
         ncol=dim(hit_array)[2],
         dimnames=dimnames(hit_array)[1:2]);
      ilist <- lapply(jamba::nameVector(rownames(imatrix)), function(isignal){
         ilist <- imatrix[isignal, ];
         # handle NA values, which occurs when interaction cutoffs are
         # used that differ from contrast cutoffs.
         # An NA value means the cutoff was not tested, it does not mean
         # the cutoff was tested and no hits were found.
         # Therefore we keep NA.
         # ilist <- sapply(ilist, jamba::rmNA);
         ilist_lengths <- jamba::rbindList(lapply(ilist, function(i){
            if (length(i) > 0 && all(is.na(i))) {
               c(hits=NA,
                  up=NA,
                  down=NA)
            } else {
               c(hits=length(i),
                  up=sum(i > 0),
                  down=sum(i < 0))
            }
         }))
         if ("integer" %in% style) {
            ilist_lengths[,1];
         } else {
            ifelse(is.na(ilist_lengths[,1]),
               "",
               paste0(
                  jamba::formatInt(ilist_lengths[,"hits"]),
                  ifelse(ilist_lengths[,1] == 0,
                     " hits",
                     paste0(
                        " hits (",
                        jamba::formatInt(ilist_lengths[,"up"]),
                        " up, ",
                        jamba::formatInt(ilist_lengths[,"down"]),
                        " down)")
                  )
               )
            )
         }
      });
      im <- do.call(cbind, ilist);
      rownames(im) <- dimnames(hit_array)[[2]];
      data.frame(check.names=FALSE,
         cutoff=icutoff,
         contrast=dimnames(hit_array)[[2]],
         im);
   }))

   if (rename_contrasts) {
      idf$contrast <- tryCatch({
         contrast2comp(idf$contrast);
      }, error=function(e){
         idf$contrast
      });
   }

   idf
}
