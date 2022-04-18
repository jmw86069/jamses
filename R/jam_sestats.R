
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
   # if (length(dimname_order) == 3 && all(c(1, 2, 3) %in% dimname_order)) {
   #    order1 <- dimname_order[1];
   #    order2 <- dimname_order[2];
   #    order3 <- dimname_order[3];
   #    names1 <- dimnames(hit_array)[[order1]];
   #    names2 <- dimnames(hit_array)[[order2]];
   #    names3 <- dimnames(hit_array)[[order3]];
   #    dim1 <- names(dimnames(hit_array))[order1];
   #    dim2 <- names(dimnames(hit_array))[order2];
   #    dim3 <- names(dimnames(hit_array))[order3];
   #    dimnames123 <- list(names1, names2, names3);
   #    names(dimnames123) <- c(dim1, dim2, dim3);
   # } else {
   #    stop("dimname_order must contain c(1, 2, 3) in any order, with length=3.");
   # }

   #
   hit_array <- sestats$hit_array;
   idf <- jamba::rbindList(lapply(dimnames(hit_array)[[3]], function(icutoff){
      imatrix <- matrix(hit_array[, , icutoff],
         nrow=dim(hit_array)[1],
         ncol=dim(hit_array)[2],
         dimnames=dimnames(hit_array)[1:2]);
      ilist <- lapply(jamba::nameVector(rownames(imatrix)), function(isignal){
         ilist <- imatrix[isignal, ];
         if ("integer" %in% style) {
            lengths(ilist)
         } else {
            paste0(
               jamba::formatInt(lengths(ilist)),
               " hits (",
               jamba::formatInt(sapply(ilist, function(i){sum(i > 0)})),
               " up, ",
               jamba::formatInt(sapply(ilist, function(i){sum(i < 0)})),
               " down)")
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
