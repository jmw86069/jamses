
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
#' ## The process
#'
#' * Contrasts are converted to `data.frame` with `contrasts_to_factors()`
#' * Each factor column is iterated to produce sets of contrasts as follows:
#'
#'    * The factor `data.frame` is subset for rows with a comparison
#'    in the factor column.
#'    * When `include_multifactor=FALSE` (not default) then data is
#'    filtered to remove rows with comparisons in any other factor columns.
#'    * When `include_singlefactor=FALSE` (not default) then data is
#'    filtered to remove rows with single values in any other factor columns.
#'    * The remaining rows are iteratively split using values in the other
#'    factor columns.
#'    * Remaining rows are also iteratively split by the depth of the contrast,
#'    oneway comparisons, and twoway comparisons.
#'    * If any subset contains more than `max_venn_size` rows, it is first
#'    split by the control factor level in the factor comparison, then
#'    it is split by the depth of the comparison.
#'    * In all cases, the resulting sets are split into subsets with
#'    size `max_venn_size`.
#'
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
#' @param contrast_style `character` string indicating how to return
#'    the resulting Venn set lists:
#'    * `"contrast"` (default) returns a `list` with contrast names
#'    * `"comp"` returns a `list` with comp names from `contrast2comp()`
#'    * `"factors"` returns a `data.frame` from `contrasts_to_factors()`
#'    with one column per design factor. In this case it may be useful
#'    to pass `factor_names` in order to assign column names to
#'    the factor columns.
#' @param max_venn_size `numeric` maximum number of groups to include
#'    per Venn diagram. When the number of contrasts to be included
#'    in a Venn set contains one extra contrast, the last two sets will
#'    be adjusted to accomodate the extra set:
#'    * when `max_venn_size=2` the last set will contain 3 members;
#'    * when `max_venn_size=3` (or higher) the last set will contain 2 members,
#'    and the previous set will contain `(max_venn_size - 1)` members.
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
#'    include_multifactor=FALSE,
#'    factor_names=c("Genotype", "Treatment"))
#' jamba::sdim(venn_setlists)
#'
#' # plot the contrasts included in one particular Venn setlist
#' par("mfrow"=c(2, 2));
#' for (n in 1:4) {
#' setest <- sedesign;
#' contrast_names(setest) <- venn_setlists[[n]];
#' plot_sedesign(setest, contrast_style="none")
#' }
#' par("mfrow"=c(1, 1))
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
 factor_names=NULL,
 contrast_style=c("contrast",
    "comp",
    "factors"),
 max_venn_size=4,
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
            jamba::printDebug("=========", "icol: ", icol);# debug
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

         # split by other factor column values
         splitlist <- lapply(jamba::nameVector(xcols), function(xcol){
            subcomps_df[[xcol]]
         })

         # split by oneway or multiway
         depth <- jamba::pasteByRow(do.call(cbind, lapply(xcols, function(xcol){
            grepl("-", subcomps_df[[xcol]])*1 + 1
         })))

         # split by current factor comparison, subset by depth
         splitlist1 <- list(depth=paste0(subcomps_df[[icol]], ".", depth))

         splitlist <- c(splitlist,
            splitlist1)
            # list(depth=depth)
         # decided not to split only by depth for now

         splitcomps_dfs <- unlist(recursive=FALSE,
            lapply(jamba::nameVectorN(splitlist), function(splitname){
               splitvector <- splitlist[[splitname]];
               if (verbose) {
                  jamba::printDebug("splitname: ", splitname);# debug
               }
               splitcomps_df <- split(subcomps_df, splitvector);
               splitcomps_df <- splitcomps_df[
                  jamba::sdim(splitcomps_df)$rows > 1];
               # optionally sub-split again
               splitcomps_df2 <- jamba::rmNULL(unlist(recursive=FALSE,
                  lapply(splitcomps_df, function(splitcomp_df){
                     icolsplit <- rep("x", nrow(splitcomp_df))
                     # if >5 rows, then split into smaller subsets
                     if (nrow(splitcomp_df) > 5) {
                        # try to split by common contol factor level
                        icolsplit <- jamba::gsubOrdered("^[^-]+-", "", splitcomp_df[[icol]]);
                     }
                     # adjust icolsplit to have 4 or fewer members
                     # - subsplit each group
                     icolsplit1 <- sub_split_vector(icolsplit,
                        max_size=max_venn_size);
                     if (verbose) {
                        jamba::printDebug("splitcomp_df:");
                        print(data.frame(splitcomp_df, icolsplit, icolsplit1));
                     }
                     lapply(split(splitcomp_df, icolsplit1), function(jdf){
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

   # assign user-friendly names
   names(setlist_names) <- tryCatch({
      venn_complists <- lapply(setlist_names, contrast2comp)
      venn_compnames <- sapply(venn_complists, function(i){
         idf <- data.frame(jamba::rbindList(strsplit(i, ":")));
         jamba::cPaste(
            jamba::cPasteU(
               as.list(idf)),
            sep=" : ")
      })
      jamba::makeNames(venn_compnames)
   }, error=function(e){
      names(setlist_names)
   })

   return(setlist_names);
}

#' Sub-split a split vector (Internal)
#'
#' Sub-split a split vector (Internal)
#'
#' This function is used to support `contrasts_to_venn_setlists()`
#' in order to limit the number of Venn sets included by subgroup.
#' This function adds another level of splitting, while keeping
#' the elements in the original order after splitting.
#'
#' @param icolsplit `numeric` vector used to split another vector
#' @param max_size `numeric` max number of identical values permitted
#'    in `icolsplit`.
#' @param ... additional arguments are ignored
#'
#' @family jamses utilities
#'
#' @returns `factor` vector suitable to use for splitting a vector
#'    into subsets, in order.
#'
#' @examples
#' icolsplit <- rep(c(1, 2, 3), c(6, 5, 4))
#' newsplit <- sub_split_vector(icolsplit)
#' split(icolsplit, newsplit)
#'
#' newsplit3 <- sub_split_vector(icolsplit, max_size=3)
#' split(icolsplit, newsplit3)
#'
#' @export
sub_split_vector <- function
(icolsplit,
 max_size=4,
 ...)
{
   # validate input
   if (max_size < 2) {
      max_size <- 2;
   }
   # convert to numeric order
   icolsplit <- match(icolsplit, unique(icolsplit))

   # split into subsets
   icolsplit <- unname(unlist(
      lapply(split(icolsplit, icolsplit), function(k){
         if (length(k) > max_size) {
            splitk <- head(
               rep(
                  seq_len(ceiling(length(k) / max_size)),
                  each=max_size),
               length(k))
            # if only one member left over,
            # flip the next to last member
            if ((length(k) %% max_size) == 1) {
               if (max_size == 2) {
                  splitk[length(splitk)] <- splitk[
                     length(splitk) - 1];
               } else {
                  splitk[length(splitk) - 1] <- splitk[
                     length(splitk)];
               }
            }
            paste0(k, ".", splitk)
         } else {
            paste0(k)
         }
      })
   ))
   icolsplit <- factor(icolsplit,
      levels=unique(icolsplit));
   icolsplit
}
