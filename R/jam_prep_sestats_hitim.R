
#' Process sestats input into a hit incidence matrix
#'
#' Process sestats input into a hit incidence matrix
#'
#' @param sestats one of the following:
#'    * `list` object output from `se_contrast_stats()`, containing `"hit_array"`
#'    * `array` in format `"hit_array"` with dimnames
#'    `"Cutoffs","Contrasts","Signals"`.
#'    * `list` of `character` vectors representing `rownames(se)` for
#'    the parent `heatmap_se()` function.
#'    * `list` of `numeric` vectors named by `rownames(se)`.
#' @param cutoff_names,contrast_names,assay_names `character` or `numeric`
#'    passed to `hit_array_to_list()` when the input is `sestats` or
#'    `hit_array`.
#' @param contrast_suffix `character` optional suffix appended to the
#'    end of each contrast name.
#' @param rename_contrasts `logical` indicating whether to rename contrasts
#'    by calling `contrast2comp()`
#' @param rows `character` or `NULL` with optional fixed set of rownames
#'    expected in the output matrix.
#'    When `rows=NULL` all rows are returned using data from `sestats`.
#'    Otherwise, only rows defined by `rows` are returned.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to `contrast2comp()` if
#'    relevant.
#'
#' @export
process_sestats_to_hitim <- function
(sestats,
   cutoff_names=NULL,
   contrast_names=NULL,
   assay_names=NULL,
   contrast_suffix=NULL,
   rename_contrasts=FALSE,
   rows=NULL,
   verbose=FALSE,
   ...)
{
   #
   # if input is list, and does not contain name "hit_array"
   # it is not sestats
   gene_hits_im <- NULL;
   hit_array <- NULL;
   if ("list" %in% class(sestats) && !"hit_array" %in% names(sestats)) {
      # assume input is a hit list
      if (is.numeric(sestats[[1]]) && length(names(sestats[[1]])) > 0) {
         # assume named directional sestats
         if (verbose) {
            jamba::printDebug("process_sestats_to_hitim(): ",
               "converting sestats with venndir::list2im_value().");
         }
         gene_hits_im <- venndir::list2im_value(sestats,
            do_sparse=FALSE);
      } else {
         # assume list of character vectors with rownames(se)
         if (verbose) {
            jamba::printDebug("process_sestats_to_hitim(): ",
               "converting sestats with venndir::list2im_opt().");
         }
         gene_hits_im <- venndir::list2im_opt(sestats,
            do_sparse=FALSE);
      }
   } else if ("list" %in% class(sestats) && "hit_array" %in% names(sestats)) {
      # if input is list, and does contain name "hit_array"
      # it is sestats, so we grab hit_array
      hit_array <- sestats$hit_array;
   } else if ("matrix" %in% class(sestats)) {
      # if input is matrix, use directly as gene_hits_im
      gene_hits_im <- sestats;
      hit_array <- NULL;
   } else {
      # everything else assumed to be hit_array and
      # let it throw an error otherwise
      hit_array <- sestats;
   }
   if (length(hit_array) == 0) {
      if (verbose) {
         jamba::printDebug("process_sestats_to_hitim(): ",
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
         jamba::printDebug("process_sestats_to_hitim(): ",
            "sestats is generating an incidence matrix with hit_array_to_list().");
      }
      if (length(contrast_names) == 0) {
         contrast_names <- dimnames(hit_array)[[2]];
      }
      gene_hitlist <- hit_array_to_list(hit_array,
         cutoff_names=cutoff_names,
         contrast_names=contrast_names,
         assay_names=assay_names);
      gene_hits <- names(jamba::tcount(names(unlist(unname(
         gene_hitlist)))));
      # confirm all gene_hits are present in the data provided
      if (!all(gene_hits %in% rownames(se))) {
         gene_hits <- intersect(gene_hits, rownames(se));
      }
      gene_hits_im <- venndir::list2im_value(gene_hitlist,
         do_sparse=FALSE)[gene_hits, , drop=FALSE];
   }

   # optionally rename contrasts
   if (TRUE %in% rename_contrasts) {
      colnames(gene_hits_im) <- tryCatch({
         contrast2comp(colnames(gene_hits_im),
            ...)
      }, error=function(e){
         colnames(gene_hits_im)
      });
   }

   # optionally handle fixed rows
   if (length(rows) > 0) {
      rows_im <- matrix(nrow=length(rows),
         ncol=ncol(gene_hits_im),
         dimnames=list(rows, colnames(gene_hits_im)),
         data=0);
      rows_use <- intersect(rows, rownames(gene_hits_im));
      if (length(rows_use) > 0) {
         rows_use_match1 <- match(rows_use, rownames(rows_im));
         rows_use_match2 <- match(rows_use, rownames(gene_hits_im));
         rows_im[rows_use_match1, colnames(gene_hits_im)] <- gene_hits_im[rows_use_match2, colnames(gene_hits_im), drop=FALSE];
      }
   }

   # append optional contrast suffix to colnames
   if (length(contrast_suffix) > 0 && any(nchar(contrast_suffix)) > 0) {
      colnames(gene_hits_im) <- paste0(colnames(gene_hits_im),
         contrast_suffix);
   }
   return(gene_hits_im);
}
