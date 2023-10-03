
#' Convert contrast names to Venn setlists for visual comparison
#'
#' Convert contrast names to Venn setlists for visual comparison
#'
#' The motivation is to take a set of contrast names, and return
#' reasonable subsets of contrasts suitable for visual comparison
#' using Venn diagrams. Ultimately, the process is analogous to
#' defining contrasts themselves: keep experimental factors fixed
#' while varying one factor at a time. The difference is that
#' experimental "factors" may themselves involve a comparison.
#'
#' The process is currently being tested for two-factor design
#' scenarios, and will be extended to handle higher factor
#' designs in future.
#'
#' @returns `list` with contrast names suitable for use in Venn diagrams.
#'
#' @param contrast_names `character` vector of contrast names, used
#'    as the priority input when supplied.
#' @param sestats `list` output from `se_contrast_stats()`, used only when
#'    `contrast_names` is not supplied.
#' @param sedesign `SEDesign` object, used only when neither `contrast_names`
#'    nor `sestats` are supplied.
#' @param include_multifactor `logical` indicating whether to include
#'    twoway contrasts in the Venn diagram logic. Currently the logic
#'    includes comparisons only across compatible twoway comparisons,
#'    it does not (yet) include Venn diagrams that include oneway and
#'    corresponding twoway comparisons together.
#'    Note: One of `include_multifactor` and `include_singlefactor`
#'    must be true.
#' @param include_singlefactor `logical` indicating whether to include
#'    single-factor contrasts in the Venn diagram logic.
#'    Note: One of `include_multifactor` and `include_singlefactor`
#'    must be true.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to `contrasts_to_factors()`.
#'
#' @examples
#' group_names <- paste0(
#'    rep(c("UL3", "dH1A", "dH1B"), each=5), "_",
#'    c("Veh", "DEX", "PMA", "SF", "Ins"))
#' sedesign <- groups_to_sedesign(group_names)
#'
#' # by default it returns contrast names
#' venn_setlists <- contrasts_to_venn_setlists(sedesign=sedesign,
#'    factor_names=c("Genotype", "Treatment"))
#' venn_setlists
#'
#' venn_set_comps <- contrasts_to_venn_setlists(sedesign=sedesign,
#'    contrast_style="comp",
#'    factor_names=c("Genotype", "Treatment"))
#' venn_set_comps
#' data.frame(names(venn_set_comps))
#'
#' venn_set_factors <- contrasts_to_venn_setlists(sedesign=sedesign,
#'    contrast_style="factors",
#'    factor_names=c("Genotype", "Treatment"))
#' venn_set_factors
#' data.frame(names(venn_set_factors))
#'
#' @export
contrasts_to_venn_setlists <- function
(contrast_names=NULL,
 sestats=NULL,
 sedesign=NULL,
 include_multifactor=TRUE,
 include_singlefactor=TRUE,
 contrast_style=c("contrast",
    "comp",
    "factors"),
 verbose=FALSE,
 ...)
{
   contrast_style <- match.arg(contrast_style)
   # convert to table of factor comparisons
   if (length(contrast_names) == 0) {
      if (length(sestats) == 0) {
         if (length(sedesign) == 0) {
            stop("Must supply one of: contrast_names, sestats, sedesign")
         }
         contrast_names <- contrast_names(sedesign);
      } else {
         contrast_names <- dimnames(sestats$hit_array)$Contrasts;
      }
   }
   # ensure contrast_names are unique
   contrast_names <- unique(contrast_names);
   comps_df <- contrasts_to_factors(contrast_names,
      verbose=verbose,
      rowname="comp",
      ...)
   names(contrast_names) <- rownames(comps_df);

   # extract subsets of contrasts by shared factors
   venn_setlist_comps <- unlist(recursive=FALSE, jamba::rmNULL(
      lapply(jamba::nameVector(colnames(comps_df)), function(icol){
         if (verbose) {
            jamba::printDebug("icol: ", icol);# debug
         }
         xcols <- setdiff(colnames(comps_df), icol);
         subcomps_df <- subset(comps_df, grepl("-", comps_df[[icol]]))
         if (!TRUE %in% include_singlefactor) {
            for (xcol in xcols) {
               # remove contrasts that compare other factors
               subcomps_df <- subset(subcomps_df,
                  grepl("-", subcomps_df[[xcol]]))
            }
         }
         if (!TRUE %in% include_multifactor) {
            for (xcol in xcols) {
               # remove contrasts that compare other factors
               subcomps_df <- subset(subcomps_df,
                  !grepl("-", subcomps_df[[xcol]]))
            }
         }
         for (icol1 in colnames(subcomps_df)) {
            subcomps_df[[icol1]] <- factor(subcomps_df[[icol1]],
               levels=unique(subcomps_df[[icol1]]));
         }
         if (verbose) {
            jamba::printDebug("subcomps_df:");print(subcomps_df);# debug
         }
         if (length(xcols) == 0) {
            return(list(subcomps_df))
         }
         splitlist <- lapply(xcols, function(xcol){
            subcomps_df[[xcol]]
         })
         # split by oneway or multiway
         m1 <- jamba::pasteByRow(do.call(cbind, lapply(xcols, function(xcol){
            grepl("-", subcomps_df[[xcol]])*1 + 1
         })))
         splitlist <- c(splitlist,
            list(depth=m1))
         print(m1);
         splitcomps_dfs <- unlist(recursive=FALSE,
            lapply(splitlist, function(splitvector){
               # if (verbose) {
               #    jamba::printDebug("xcol: ", xcol);# debug
               # }
               splitcomps_df <- split(subcomps_df, splitvector);
               splitcomps_df <- splitcomps_df[
                  jamba::sdim(splitcomps_df)$rows > 1];
               if (verbose) {
                  print(splitcomps_df);# debug
               }
               # optionally sub-split again
               splitcomps_df2 <- jamba::rmNULL(unlist(recursive=FALSE,
                  lapply(splitcomps_df, function(splitcomp_df){
                     icolsplit <- rep(1, nrow(splitcomp_df))
                     if (nrow(splitcomp_df) >= 6) {
                        icolsplit <- jamba::gsubOrdered("^[^-]+-", "", splitcomp_df[[icol]]);
                     }
                     lapply(split(splitcomp_df, icolsplit), function(jdf){
                        if (nrow(jdf) < 2) {
                           NULL
                        } else {
                           jdf
                        }
                     })
                  })
               ))
            }))
         return(splitcomps_dfs)
      })))
   setlist_names <- lapply(venn_setlist_comps, function(idf){
      if ("contrast" %in% contrast_style) {
         contrast_names[rownames(idf)]
      } else if ("factors" %in% contrast_style) {
         idf;
      } else {
         rownames(idf)
      }
   })

   # remove duplicates
   setlist_names <- setlist_names[!duplicated(setlist_names)];

   return(setlist_names);
}
