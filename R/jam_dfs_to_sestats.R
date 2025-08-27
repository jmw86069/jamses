
#' Prepare SEStats from a list of stat data.frame
#'
#' Prepare SEStats from a list of stat data.frame
#'
#' This function converts a list of stat `data.frame` objects
#' to a `list` object as returned by `se_contrast_stats()`.
#' In future this function will return a S4 `SEStats` object.
#'
#' ## Limitations
#'
#' This function only populates certain components:
#'
#' * `stats_dfs` as a `list` by assay_name (signal), with `list`
#' by contrast name with `data.frame` objects.
#' * `hit_array` as a 3-dimension array with 'Cutoffs', 'Contrasts',
#' and 'Signal' dimensions.
#' * `hit_list` as a `list` by assay name (signal), with `list` by
#' contrast names, containing `numeric` values named by measurement.
#'
#' This function does not return certain data which is not described
#' in a `list` of `data.frame` objects:
#'
#' * `idesign`, `icontrasts`, `normgroup` cannot be defined or inferred.
#'
#' Other limitations:
#'
#' * `stats_df` with the combined stats across contrasts as
#' one wide `data.frame`, is not created, partly because it does not
#' assume the contrast name is included in each corresponding
#' stat colname in each `data.frame`. Potential future enhancement.
#'
#' ## Rules Used
#'
#' In each `data.frame` certain colnames are recognized:
#'
#' * The hit columns is recognized by one or more colnames
#' that begin with `"hit "`. The hit column should have values
#' of `-1`, `0`, or `1` to indicate whether the row met the
#' statistical criteria. Typically the statistical criteria
#' are also included in the column header, followed by the
#' contrast name.
#'
#'    * If there is no hit column, the arguments '_cutoff' are applied
#'    to each row, and a new hit column is created accordingly.
#'    If there are multiple cutoff values, multiple hit columns are
#'    created.
#'
#' * Fold change column is recognized by a column beginning `"fold "`
#' (with space), typically followed by the contrast name. For example
#' `"fold Dex-Ctrl"`.
#'
#'    * If no fold change column exists, it searches for log fold column,
#'    specifically assumed to be log2 fold change even when `"logFC"`.
#'    By convention, edgeR and limma use `"logFC"` for log2 fold change.
#'    * If neither fold change, nor log fold change, columns are found,
#'    this function stops.
#'    * The log2 fold change is converted to normal fold change for
#'    filtering.
#'
#' * The raw P-value column should begin `"P-Value "` or any variation
#' of "pval", "p.val", "pvalue", "p.value", case-insensitive.
#' * The adjusted P-value column should begin `"adj-P-Val "` or any variation
#' of "adjp", "padj", "pvalue", "p.value", case-insensitive.
#'
#' @return `list` consistent with `se_contrast_stats()` output, except
#'    limited to elements: 'stats_df', 'hit_array', 'hit_list'.
#'    The content should be sufficient for use in `heatmap_se()` and
#'    other related functions.
#'
#' @family jamses utilities
#'
#' @param stats_dfs `list` with one of two formats:
#'    1. Containing a `list` named by `assay_names` (signal), which
#'    then contains a `list` of `data.frame` objects.
#'    2. Containing `data.frame` objects.
#'
#'    Note: The `list` of `data.frame` objects can be named by contrast,
#'    and these contrast names will be used if the contrast is not already
#'    encoded in the colnames the corresponding `data.frame`.
#' @param use_assay_name `character` string, default "norm", used when
#'    argument 'stats_dfs' is supplied as a `list` of `data.frames`
#'    (option 2 described for 'stats_dfs'). In this case, the assay name
#'    will be 'use_assay_name'.
#' @param p_cutoff,adjp_cutoff,fold_cutoff,mgm_cutoff `numeric` values used
#'    when there is no matching hit column in each `data.frame`, matching
#'    the `hit_pattern`.
#' @param hit_pattern,fold_pattern,lfc_pattern,p_pattern,adjp_pattern,mgm_pattern
#'    `character` string with regular expression pattern used to match
#'    each corresponding column type: Hit column, fold change, log fold change,
#'    P-value (unadjusted), adjusted P-value (or FDR or Q-value),
#'    max group mean (mgm).
#'    * The pattern is assumed to be at the start of each column name string.
#' @param contrast_pattern `character` string used to match the contrast
#'    encoded in any of the cutoff columns for the '_pattern' arguments.
#'    * The pattern is essentially any non-whitespace string that contains
#'    a hyphen '-' inside it. If it also contains a colon ':' the string
#'    is assumed to be a comp (see `comp2contrast()`) otherwise it is
#'    considered a contrast name.
#'    * When present, the contrast is extracted from the column name.
#'    * When not present, the name from the `list` of `data.frame` objects
#'    is assumed to be the contrast name.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
list_to_sestats <- function
(stats_dfs,
 use_assay_name="norm",
 p_cutoff=1,
 adjp_cutoff=0.05,
 fold_cutoff=1.5,
 mgm_cutoff=0,
 hit_pattern="hit",
 fold_pattern="fold",
 lfc_pattern="^logFC|^log2F|^logfold|^log.fold|^lfc",
 p_pattern="p.{0,1}val(ue|)",
 adjp_pattern="adj.{0,1}p.{0,1}val(ue|)|p.{0,1}adj|fdr|q.{0,1}val(ue|)",
 mgm_pattern="mgm|max.{0,1}group.{0,1}mean",
 contrast_pattern="^.*[ ]([^ ]+-[^ ]+)$",
 verbose=FALSE,
 ...)
{
   #
   ## start with stats_dfs and intuit the hit columns

   ## check whether input data is a list per assay_name
   if ("data.frame" %in% class(stats_dfs[[1]])) {
      stats_dfs <- list(assay=stats_dfs);
      names(stats_dfs) <- use_assay_name;
   }

   ## confirm stats_dfs
   stats_dfs <- lapply(stats_dfs, function(iDFs){
      if (length(names(iDFs)) == 0) {
         names(iDFs) <- paste0("contrast",
            jamba::padInteger(seq_along(iDFs)))
      }
      unlist(recursive=FALSE, lapply(names(iDFs), function(iname){
         iDF <- iDFs[[iname]];
         use_hit_pattern <- paste0("^(", hit_pattern, ")[^ ]*[ ]{0,1}");
         iHitCols <- jamba::nameVector(
            jamba::provigrep("^hit[ .]", colnames(iDF)));

         ## recognize stat colnames where possible
         # fold change
         use_fold_pattern <- paste0("^(", fold_pattern, ")[^ ]*[ ]{0,1}");
         fc_colname <- head(jamba::vigrep(use_fold_pattern,
            colnames(iDF)), 1)
         use_lfc_pattern <- paste0("^(", lfc_pattern, ")[^ ]*[ ]{0,1}");
         lfc_colname <- head(jamba::vigrep(use_lfc_pattern,
            colnames(iDF)), 1);
         # P-value
         use_p_pattern <- paste0("^(", p_pattern, ")[^ ]*[ ]{0,1}");
         p_colname <- head(vigrep(use_p_pattern, colnames(iDF)), 1);
         use_adjp_pattern <- paste0("^(", adjp_pattern, ")[^ ]*[ ]{0,1}");
         adjp_colname <- head(vigrep(use_adjp_pattern, colnames(iDF)), 1);
         # mgm - max group mean
         use_mgm_pattern <- paste0("^(", mgm_pattern, ")[^ ]*[ ]{0,1}");
         mgm_colname <- head(vigrep(use_mgm_pattern, colnames(iDF)), 1);

         # extract contrast name from stat columns if possible
         use_colnames <- c(p_colname,
            adjp_colname,
            fc_colname,
            lfc_colname,
            mgm_colname)
         use_contrast <- NULL;
         if (jamba::igrepHas(contrast_pattern, use_colnames)) {
            use_contrast <- gsub(contrast_pattern, "\\1",
               head(jamba::vigrep(contrast_pattern, use_colnames), 1));
         }
         if (length(use_contrast) == 0) {
            use_contrast <- iname;
         }
         # convert comp to contrast if necessary
         if (jamba::igrepHas(":", use_contrast)) {
            use_contrast <- jamses::contrast2comp(use_contrast);
         }

         # Prepare hit column if needed
         if (length(iHitCols) == 0) {
            iHitCol <- paste0("hit");
            if (length(mgm_cutoff) > 0) {
               iHitCol <- paste0(iHitCol, " mgm", mgm_cutoff);
            }
            if (length(p_cutoff) > 0 && p_cutoff > 0) {
               iHitCol <- paste0(iHitCol, " p", p_cutoff);
            }
            if (length(adjp_cutoff) > 0 && adjp_cutoff > 0) {
               iHitCol <- paste0(iHitCol, " adjp", adjp_cutoff);
            }

            ## Fold change
            if (length(fold_cutoff) > 0 && fold_cutoff >= 1) {
               iHitCol <- paste0(iHitCol, " fc", fold_cutoff);
            }
            if (length(fc_colname) == 0) {
               # prepare fold from log fold
               if (length(lfc_colname) == 0 &&
                     length(fold_cutoff) > 0 &&
                     any(fold_cutoff > 1)) {
                  cli::cli_abort(paste0(
                     "Cannot find column matching {.var fold_pattern} or ",
                     "{.var lfc_pattern}."))
                  stop("Cannot find fold,logFC column to apply fold_cutoff.")
               }
               # create fold column from lfc_colname
               fc_colname <- gsub("^[ ]+|[ ]+$", "",
                  gsub(uselfc_pattern, "fold ", lfc_colname));
               # convert log2 fold to fold
               iDF[[fc_colname]] <- jamses::fold_to_log2fold(
                  iDF[[lfc_colname]]);
            }

            iHitCol <- paste0(iHitCol, " ", use_contrast);
            is_hit <- rep(TRUE, nrow(iDF));
            if (length(fc_colname) == 1 &&
                  length(fold_cutoff) == 1 &&
                  fold_cutoff > 1) {
               is_hit <- is_hit & abs(iDF[[fc_colname]]) >= fold_cutoff;
            }
            if (length(p_colname) == 1 &&
                  length(p_cutoff) == 1 &&
                  p_cutoff < 1) {
               is_hit <- is_hit & iDF[[p_colname]] <= p_cutoff;
            }
            if (length(adjp_colname) == 1 &&
                  length(adjp_cutoff) == 1 &&
                  adjp_cutoff < 1) {
               is_hit <- is_hit & iDF[[adjp_colname]] <= adjp_cutoff;
            }
            if (length(mgm_colname) == 1 &&
                  length(mgm_cutoff) == 1) {
               is_hit <- is_hit & iDF[[mgm_colname]] >= mgm_cutoff;
            }
            iDF[[iHitCol]] <- is_hit * sign(iDF[[fc_colname]]);
            iHitCols <- c(iHitCols, iHitCol);
         }
         if (verbose > 1) {
            jamba::printDebug("prepare_sestats(): ",
               "added iHitCols:", iHitCols);
         }
         # assign rownames
         if (!any(duplicated(iDF[[1]]))) {
            rownames(iDF) <- iDF[[1]];
         }
         setNames(object=list(iDF),
            nm=use_contrast)
      }))
   })

   ## list of named lists
   stats_hits <- lapply(stats_dfs, function(iDFs){
      lapply(iDFs, function(iDF){
         iHitCols <- jamba::nameVector(
            jamba::provigrep("^hit[ .]", colnames(iDF)));
         if (verbose) {
            jamba::printDebug("prepare_sestats(): ",
               "iHitCols:", iHitCols);
         }
         lapply(iHitCols, function(iHitCol){
            iHitRows <- (!is.na(iDF[,iHitCol]) & iDF[,iHitCol] != 0);
            jamba::nameVector(
               iDF[iHitRows,iHitCol],
               rownames(iDF)[iHitRows]);
         });
      });
   });

   ## array of named lists
   all_cutoffs_list <- lapply(stats_hits[[1]], function(i){
      gsub("[ ]+[^ ]+$", "", names(i))
   })
   all_cutoffs <- unique(unlist(all_cutoffs_list));
   arrayDim <- c(length(all_cutoffs),
      length(stats_hits[[1]]),
      length(stats_hits));
   arrayDimnames <- list(all_cutoffs,
      names(stats_hits[[1]]),
      names(stats_hits));
   names(arrayDimnames) <- c("Cutoffs",
      "Contrasts",
      "Signal");

   ## assemble hit_array
   # - allow for missing entries
   # first define the array indices
   kji <- jamba::rbindList(lapply(names(stats_hits), function(i){
      jamba::rbindList(lapply(names(stats_hits[[i]]), function(j){
         jamba::rbindList(lapply(names(stats_hits[[i]][[j]]), function(k){
            k1 <- gsub("[ ][^ ]+$", "", k);
            kji <- cbind(match(k1, arrayDimnames[[1]]),
               match(j, arrayDimnames[[2]]),
               match(i, arrayDimnames[[3]]))
            kji
         }))
      }))
   }))
   hit_array <- array(dim=arrayDim,
      dimnames=arrayDimnames);
   # fill data into the array, which somehow converts it into a list
   hit_array[kji] <- unlist(recursive=FALSE,
      unlist(recursive=FALSE,
         stats_hits))
   # create the array again using the list data
   hit_array <- array(dim=arrayDim,
      data=hit_array,
      dimnames=arrayDimnames);


   ## return the owl
   list(
      stats_dfs=stats_dfs,
      hit_list=stats_hits,
      hit_array=hit_array)
}
