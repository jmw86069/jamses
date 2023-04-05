
#' Convert contrast names to sedesign object
#'
#' Convert contrast names to sedesign object
#'
#' This function is a convenience function intended only to convert
#' a vector of contrast names into `sedesign` format for use in other
#' functions. It assumes only one sample replicate per group for this
#' purpose.
#'
#' One utility of this function is to convert two-way contrast names
#' into a contrast matrix, to test whether the contrast defined
#' is equivalent for the two names.
#'
#' @return `sedesign` object, or if the input `contrast_names` contain
#'    mixed number of factors per contrast, the output is split into
#'    a `list` of `sedesign` objects based upon the number of factors.
#'
#' @param contrast_names `character` vector of contrast names where
#'    factors are separated by `factor_sep` and contrasts are separated
#'    by `contrast_sep`.
#' @param factor_sep,contrast_sep `character` strings used as delimiters
#'    between factors and contrasts, respectively.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' contrast_names_3fac <- c(
#'    "(CellA_Treated-CellA_Control)-(CellB_Treated-CellB_Control)",
#'    "(CellA_Treated-CellB_Treated)-(CellA_Control-CellB_Control)")
#' contrast_names_to_sedesign(contrast_names_3fac)
#'
#' contrast_names_3way <- c(
#'    paste0("((CellA_Treated_Mut-CellB_Treated_Mut)-",
#'       "(CellA_Control_Mut-CellB_Control_Mut))-",
#'       "((CellA_Treated_WT-CellB_Treated_WT)-",
#'       "(CellA_Control_WT-CellB_Control_WT))"),
#'    paste0("((CellA_Treated_Mut-CellA_Control_Mut)-",
#'       "(CellB_Treated_Mut-CellB_Control_Mut))-",
#'       "((CellA_Treated_WT-CellA_Control_WT)-",
#'       "(CellB_Treated_WT-CellB_Control_WT))"))
#' contrast_names_to_sedesign(contrast_names_3way)
#'
#' @export
contrast_names_to_sedesign <- function
(contrast_names,
 factor_sep="_",
 contrast_sep="-",
 ...)
{
   # extract group names
   contrast_groups <- strsplit(
      gsub("[()]", "", contrast_names),
      contrast_sep)
   contrast_factor_count <- sapply(contrast_groups, function(igroups){
      max(lengths(strsplit(igroups, factor_sep)))
   })
   contrast_names_list <- split(contrast_names, contrast_factor_count)

   group_names <- unique(unlist(
      strsplit(
         gsub("[()]", "", contrast_names),
         contrast_sep)))
   groups_df <- data.frame(check.names=FALSE,
      stringsAsFactors=FALSE,
      group_names,
      row.names=group_names,
      factor_count=lengths(strsplit(group_names, factor_sep)))

   # for now, treat different factor_count independently
   groups_df_list <- split(groups_df, groups_df$factor_count);
   sedesign_list <- lapply(jamba::nameVectorN(groups_df_list), function(iname){
      igroups_df <- groups_df_list[[iname]]
      groups_to_sedesign(igroups_df[,1, drop=FALSE],
         contrast_names=unique(contrast_names_list[[iname]]))
   })

   # return single sedesign when input has one factor depth
   if (length(sedesign_list) == 1) {
      return(sedesign_list[[1]])
   }
   return(sedesign_list)

   # optionally deconstruct the factors
   factor_df_list <- lapply(groups_df_list, function(igroups_df){
      idf <- data.frame(check.names=FALSE,
         stringsAsFactors=FALSE,
         jamba::rbindList(strsplit(igroups_df[,1], factor_sep)))
      colnames(idf) <- jamba::makeNames(rep("factor", ncol(idf)),
         suffix="",
         renameOnes=TRUE)
      idf
   })

   ## for display, model the layout using vcd::mosaic options
   i <- 1;
   vcdm <- vcd::mosaic(
      vcd::structable(factor_df_list[[i]],
         split_vertical=rep(
            c(FALSE, TRUE, TRUE),
            length.out=ncol(factor_df_list[[i]]))),
      shade=FALSE,
      labeling=vcd::labeling_values)

}
