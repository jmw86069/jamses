
#' Quick conversion of hit_array to hit_list
#'
#' Quick conversion of hit_array to hit_list
#'
#' This function is mainly useful when there are multiple dimensions
#' unresolved in a hit_array, in which case this function will combine
#' hits across the different cutoffs and signals contained
#' in the `hit_array` of `sestats` output from `se_contrast_stats()`.
#'
#' @return `list` named by `contrast_names`, that contains unique statistical
#'    hits by combining entries across the `cutoff_names` and `assay_names`
#'    for each contrast.
#'
#' @param hit_array `array` output from `se_contrast_stats()`, list element
#'    `"hit_array"`.
#' @param contrast_names `character` vector of contrasts.
#'
#' @export
hit_array_to_list <- function
(hit_array,
 contrast_names=NULL,
 cutoff_names=NULL,
 assay_names=NULL,
 ...)
{
   if (length(cutoff_names) == 0) {
      cutoff_names <- dimnames(hit_array)$Cutoffs;
   }
   if (is.numeric(cutoff_names)) {
      cutoff_names <- dimnames(hit_array)$Cutoffs[cutoff_names];
   }
   cutoff_names <- intersect(cutoff_names,
      dimnames(hit_array)$Cutoffs);
   if (length(contrast_names) == 0) {
      contrast_names <- dimnames(hit_array)$Contrasts;
   }
   if (is.numeric(contrast_names)) {
      contrast_names <- dimnames(hit_array)$Contrasts[contrast_names];
   }
   contrast_names <- intersect(contrast_names,
      dimnames(hit_array)$Contrasts);
   if (length(assay_names) == 0) {
      assay_names <- dimnames(hit_array)$Signal;
   }
   if (is.numeric(assay_names)) {
      assay_names <- dimnames(hit_array)$Signal[assay_names];
   }
   assay_names <- intersect(assay_names,
      dimnames(hit_array)$Signal);

   if (any(length(cutoff_names) == 0 ||
         length(contrast_names) == 0 ||
         length(assay_names) == 0)) {
      stop("cutoff_names, contrast_names, assay_names must match hit_array dimnames.");
   }
   hit_list <- apply(hit_array[cutoff_names, contrast_names, assay_names, drop=FALSE], 2, function(i){
      idf <- jamba::rbindList(lapply(i, function(j){
         j <- jamba::rmNA(j);
         if (length(j) == 0) {
            return(NULL)
         }
         data.frame(item=names(j), value=j)
      }))
      if (length(idf) == 0 || nrow(idf) == 0) {
         return(NULL)
      }
      jdf <- subset(idf, !duplicated(item));
      jamba::nameVector(jdf$value, jdf$item)
   })
}
