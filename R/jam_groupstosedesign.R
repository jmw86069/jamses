
#' Create SEDesign from experimental groups
#'
#' Create SEDesign from experimental groups
#'
#' This function creates `SEDesign` with appropriate design
#' and contrasts, based upon experimental groups.
#' This approach will use multiple experimental factors
#' to create appropriate one-way and n-way contrasts,
#' where each contrast represents a symmetric comparison
#' of each independent factor.
#'
#' Input can be provided in one of two ways:
#'
#' 1. `SummarizedExperiment` where experiment design is derived from
#' `SummarizedExperiment::colData()` of the `se` object, and
#' uses columns defined by `group_colnames`. This input should be
#' equivalent to providing a `data.frame` whose `rownames()` are
#' equal to `colnames(se)`.
#' 2. `data.frame` where each column represents a design factor.
#'
#'     * An example of `data.frame` input:
#'    ```R
#'    ifactors <- data.frame(
#'       treatment=c("Control", "Control", "Treated", "Treated"),
#'       genotype=c("Wildtype", "Knockout", "Wildtype", "Knockout"))
#'    ```
#'
#' 3. `character` vector, where design factor levels are separated
#' by a delimiter such as underscore `"_"`. This input will be
#' converted to `data.frame` before processing.
#'
#'    * An example of `character` input:
#'    ```R
#'    ifactors <- c(
#'       "Control_Wildtype",
#'       "Control_Knockout",
#'       "Treated_Wildtype",
#'       "Treated_Knockout")
#'    ```
#'
#' When rownames are provided in the `data.frame`, or names
#' are provided with a `character` vector, they are retained
#' and used as sample identifiers.
#'
#' Note:
#' This function will change any `"-"` in a factor name to
#' `"."` prior to detecting valid contrasts, in order to
#' prevent confusion and potential problems using the
#' contrast names in downstream analyses.
#' This step does not call `base::make.names()`, so that
#' step should be run beforehand if required.
#'
#' ## Troubleshooting
#'
#' * When this function returns no contrasts, or returns an unexpected
#' error during processing, it is most likely due to the limitation
#' of comparing one factor at a time. For example, the logic will
#' not define contrast `time1_treatment1-time2_treatment2`, because
#' this contrast changes two factors, it will only permit either
#' `time1_treatment1-time1_treatment2` or `time1_treatment1-time2_treatment1`.
#' * `max_depth` and `factor_order` are used to define the order in
#' which factors are compared, but do not affect the order of factors
#' used for things like group names.
#'
#' @return `SEDesign` object with the following slots:
#'    * `design`: `numeric` matrix with sample-to-group association
#'    * `contrasts`: `numeric` matrix with group-to-contrast association
#'    * `samples`: `character` vector that represents individual sample
#'    replicates, equivalent to `rownames()` of the `design` matrix.
#'
#'
#' @param ifactors `data.frame` or `character` vector.
#'    * When `data.frame` is supplied, each column is used as a
#'    design factor, and rownames are recognized as sample identifiers.
#'    * When `character` vector is supplied, it is converted to
#'    `data.frame` by splitting values with a delimiter
#'    `factor_sep`, and names are recognized as sample identifiers.
#' @param group_colnames `character` vector or `NULL`, used to
#'    define a subset of columns to use when `ifactors` is supplied
#'    as a `data.frame`. When `ifactors` is supplied as a `character`
#'    vector, this argument is used to define the `colnames`.
#' @param isamples `character` vector or `NULL`, optionally used to subset
#'    the sample identifiers used in subsequent steps. Note that only
#'    groups and contrasts that contain samples will be defined.
#' @param idesign `numeric` matrix or `NULL`, intended as an optional
#'    method to use an existing design matrix.
#' @param factor_order `integer` or `character` vector, used to define a
#'    specific order of factors when generating contrasts,  useful
#'    when there are multiple experimental factors.
#'    It can be helpful to force a secondary factor to be
#'    compared before a primary factor especially in two-way contrasts.
#'    Note that `factor_order` refers to the columns (factors) and not
#'    the factor levels (not column values).
#' @param omit_grep `character` regular expression pattern used to
#'    exclude secondary factors from contrasts.
#' @param max_depth `integer` value indicating the maximum depth of
#'    statistical contrasts to create. For example `max_depth=2` will
#'    allow two-way contrasts, and `max_depth=1` will only create
#'    one-way contrasts.
#' @param factor_sep `character` string used as a
#'    delimiter to separate experimental factors, when recognizing
#'    or creating experimental group names.
#' @param contrast_sep `character` string used as a
#'    delimiter to separate groups within each contrast name.
#' @param remove_pairs `list` of `character` vectors of factors
#'    that should not be compared. Each `character` vector should
#'    contain two factor levels for any given experimental factor,
#'    where those two factor levels should not be compared in
#'    the same pairwise contrast. For example, consider an experimental
#'    factor defined `treatment <- c("control", "dex", "compoundx")`.
#'    To prevent a direct comparison of `"dex"` to `"compoundx"`,
#'    use argument `remove_pairs=list(c("dex", "compoundx"))`.
#' @param make_unique `logical` indicating whether to make output
#'    contrasts unique.
#' @param pre_control_terms `character` vector used to
#'    place factor levels first in the order of levels, so these
#'    terms will be the denominator for contrasts. This approach
#'    is useful when the input `ifactors` does not already contain
#'    a `factor` with a specific order of factor levels.
#' @param add_contrastdf `data.frame` or `character` or `NULL`,
#'    intended to include a specific contrast in the output.
#'    This argument is typically used during iterative processing,
#'    and is not usually user-defined. It must contain
#' @param contrast_names `character` optional vector of specific
#'    contrasts to use when creating the contrast matrix. When
#'    `contrast_names=NULL` as default, the function defines contrasts
#'    using its internal logic. When `contrast_names` is supplied,
#'    only these `contrast_names` are used, with no other contrasts.
#' @param current_depth `integer` value used during iterative
#'    operations of this function.
#' @param rename_first_depth `logical` value used during iterative
#'    operations of this function.
#' @param return_sedesign `logical` used during iterative
#'    operations of this function. When `return_sedesign=FALSE`
#'    this function returns a `list`:
#'    * `"contrast_df"`: a `data.frame` as used in argument
#'    `add_contrastdf`, which describes each unique contrast.
#'    * `"contrast_names"`: a `character` vector of contrast names,
#'    which become `colnames()` of the contrast matrix.
#'    * `"idesign"`: a `numeric` design matrix as defined by the input data,
#'    suitable for debugging purposes for example.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @family jam experiment design
#'
#' @examples
#' # first define a vector of sample groups
#' igroups <- jamba::nameVector(paste(rep(c("WT", "KO"), each=6),
#'    rep(c("Control", "Treated"), each=3),
#'    sep="_"),
#'    suffix="_rep");
#' igroups <- factor(igroups, levels=unique(igroups));
#' igroups;
#'
#' condes <- groups_to_sedesign(igroups);
#' design(condes);
#' contrasts(condes);
#'
#' # now you can visualize the samples used in each contrast
#' iDesignL$idesign %*%  iDesignL$icontrasts;
#'
#' # you can adjust the order of factor levels per comparison
#' groups_to_sedesign(as.character(iGroups))$contrastName
#'
#' # make "WT" the first control term
#' groups_to_sedesign(as.character(iGroups), pre_control_terms=c("WT"), factor_order=2:1)$contrastName
#'
#' # prevent comparisons of WT to WT, or KO to KO
#' groups_to_sedesign(as.character(iGroups),
#'    remove_pairs=list(c("WT"), c("KO")))
#'
#' # input as a data.frame with ordered factor levels
#' ifactors <- data.frame(Genotype=factor(c("WT","WT","KO","KO"),
#'    levels=c("WT","KO")),
#'    Treatment=factor(c("Treated","Control"),
#'       levels=c("Control","Treated")))
#' ifactors;
#' groups_to_sedesign(ifactors)
#'
#'
#' # Again remove WT-WT and KO-KO contrasts
#' groups_to_sedesign(ifactors,
#'    remove_pairs=list(c("WT"), c("KO")))
#'
#' # demonstrating default_order "asis"
#' # contrasts show A-B, because B appears fist
#' # contrasts show Untreated-Treated because Treated appears first
#' df_test <- data.frame(
#'    set=c("B", "B", "A", "A"),
#'    treat=c("Treated", "Untreated"))
#' groups_to_sedesign(df_test)
#' groups_to_sedesign(jamba::pasteByRow(df_test))
#'
#' # demonstrating default_order "sort_samples"
#' # contrasts show B-A, because A is sorted first
#' # contrasts show Treated-Untreated because sort_samples()
#' #    determines "Untreated" is a preferred control term
#' groups_to_sedesign(df_test,
#'    default_order="sort_samples")
#' groups_to_sedesign(jamba::pasteByRow(df_test),
#'    default_order="sort_samples")
#'
#' # demonstrating default_order "mixedSort"
#' # contrasts show B-A, because A is sorted first
#' # contrasts show Untreated-Treated because Treated is sorted first
#' groups_to_sedesign(df_test,
#'    default_order="mixedSort")
#'
#' @export
groups_to_sedesign <- function
(ifactors,
 group_colnames=NULL,
 isamples=NULL,
 idesign=NULL,
 factor_order=NULL,
 omit_grep="[-,]",
 max_depth=2,
 factor_sep="_",
 contrast_sep="-",
 remove_pairs=NULL,
 pre_control_terms=NULL,
 add_contrastdf=NULL,
 contrast_names=NULL,
 current_depth=1,
 rename_first_depth=TRUE,
 return_sedesign=TRUE,
 default_order=c("asis",
    "sort_samples",
    "mixedSort"),
 verbose=FALSE,
 ...)
{
   ## Purpose is to take a data.frame, whose rows are groups,
   ## and whose columns are factors with factor levels as column values,
   ## and generate pairwise contrast names where only one factor changes
   ## at a time
   ##
   ## ifactors can be one of the following:
   ##
   ## - data.frame whose columns represent each statistical factor,
   ## whose values are either character, numeric, or factor, the latter
   ## can be ordered in order to provide preference to control groups.
   ##
   ## - vector of character strings representing each group,
   ## where the factors are separated by factor_sep, e.g. "WT_Dex", "NT_Veh"
   ##
   ## - idesign matrix whose colnames represent group names, and rownames
   ## represent samples present in those groups.
   ##
   ## - allNorm list object, with "targets" containing a data.frame of sample
   ## annotations, and group_colnames defines the columns to use for grouping.
   ##
   ## - remove_pairs is a list of vectors, where each vector is expected to
   ## contain two elements representing two factor levels not to be compared.
   ## For example, an experiment with Control, NTC, Vehicle, Dex, might not
   ## want to compare NTC-Control, Vehicle-Control, Dex-Control,
   ## remove_pairs <- list(c("NTC","Control"),c("Vehicle","Control"),c("Dex","Control"));
   ##
   ## TODO: enable remove_pairs to filter out contrasts after they are defined,
   ## for example c("NTC,Control", "d0") would remove the contrast NTC_d0-Control-d0
   ##
   ## makeUnique=TRUE will only return one entry for each set of factors compared,
   ## which will remove cases where factor 2 is tested, then factor 1 tested as an
   ## interaction; if factor 1 and factor 2 are already represented in another
   ## interaction contrast.
   ##
   ## Ultimately a table of experiment design is created, with number of columns
   ## equal to the number of factors. By default the contrasts are applied for
   ## each factor in order of colnames, but factor_order can be used to specify
   ## a custom order. This change can affect the way two-way contrasts are
   ## defined, by forcing the first/internal contrast to use a particular
   ## factor. In theory the result is simply aesthetic, as the underlying
   ## significance of the two-way comparison will be identical. But if not
   ## for aesthetics, what are we doing?
   ##
   ## TODO: fix issue when one column contains numeric values instead of
   ## character or factor, e.g. when "Time" contains c(15,45).
   ## One solution is convert to factor, then proceed.
   if (!suppressPackageStartupMessages(require(limma))) {
      stop("limma is required for groups_to_sedesign()).");
   }
   sample2group <- NULL;

   # validate default_order
   default_order <- match.arg(default_order);

   ## Handle remove_pairs by expanding to both orientations of contrast
   if (!is.null(remove_pairs)) {
      if (!is.list(remove_pairs)) {
         stop("remove_pairs must be a list of 1- or 2-member character vectors");
      }
      remove_pairsFull <- jamba::cPasteS(remove_pairs)
      if (verbose >= 2) {
         jamba::printDebug("groups_to_sedesign(): ",
            "remove_pairsFull:");
         print(remove_pairsFull);
      }
   }

   ## Special case where one data.frame column is sent, which is delimited.
   ## Mainly we treat as a vector, except that we keep the rownames
   ## so we can derive isamples.
   if ("SummarizedExperiment" %in% class(ifactors)) {
      se <- ifactors;
      if (length(isamples) == 0) {
         isamples <- colnames(se)
      }
      group_colnames <- intersect(group_colnames,
         colnames(SummarizedExperiment::colData(se)));
      if (length(group_colnames) > 0) {
         ifactors <- data.frame(check.names=FALSE,
            stringsAsFactors=FALSE,
            SummarizedExperiment::colData(se)[,group_colnames, drop=FALSE])
         rownames(ifactors) <- colnames(se);
         ifactors <- ifactors[match(isamples, rownames(ifactors)),,drop=FALSE];
      }
      if (verbose) {
         jamba::printDebug("groups_to_sedesign(): ",
            "ifactors from SummarizedExperiment input:");
         print(ifactors);
      }
   } else if (jamba::igrepHas("data.frame|matrix", class(ifactors)) &&
         ncol(ifactors) == 1) {
      ifactors <- jamba::nameVector(ifactors[,1], rownames(ifactors));
   }

   if (jamba::igrepHas("factor|character", class(ifactors))) {
      #####################################################
      ## Vector input
      ##
      if (verbose) {
         jamba::printDebug("groups_to_sedesign(): ",
            "splitting vector into groups");
      }
      if (length(names(ifactors)) == 0) {
         if (length(isamples) == 0) {
            ## Create isamples
            isamples <- jamba::makeNames(rep("sample", length(ifactors)));
            names(ifactors) <- isamples;
         } else if (length(isamples) != length(ifactors)) {
            stop(paste0("length(isamples) must be equal length(ifactors) ",
               "when there are no names(ifactors)."));
         }
         names(ifactors) <- isamples;
      } else if (length(isamples) == 0) {
         isamples <- names(ifactors);
      } else {
         if (!any(isamples %in% names(ifactors)) && length(isamples) == length(ifactors)) {
            ## Use isamples as-is
            names(ifactors) <- isamples;
         } else if (!all(isamples %in% names(ifactors))) {
            stop(paste0("isamples is present in some not not all names(ifactors). ",
               "isamples must either: all be present in names(ifactors); or ",
               "present in none of names(ifactors) and length(isamples) == length(ifactors)."))
         } else {
            ## Re-order ifactors to match isamples
            ifactors <- ifactors[match(isamples, names(ifactors))];
         }
      }
      if (jamba::igrepHas("factor", class(ifactors))) {
         ## Convert factor to a data.frame where each column
         ## is a factor with ordered levels that match the order
         ## the factor levels appear in the original factor.
         iFactorsL <- strsplitOrdered(ifactors, factor_sep);
         names(iFactorsL) <- names(ifactors);
         iFactorsLevels <- levels(iFactorsL[[1]]);
         ifactors <- data.frame(check.names=FALSE,
            stringsAsFactors=FALSE,
            jamba::rbindList(
               strsplit(as.character(ifactors),
                  factor_sep)));
         rownames(ifactors) <- names(iFactorsL);
         for (i in seq_len(ncol(ifactors))) {
            ifactors[,i] <- factor(ifactors[,i],
               levels=intersect(iFactorsLevels, ifactors[,i]));
         }
      } else {
         ## Convert to data.frame
         ifactors <- data.frame(check.names=FALSE,
            stringsAsFactors=FALSE,
            jamba::rbindList(strsplit(ifactors, factor_sep)));
         ## Convert each column to factor for proper sort order
         for (iCol in seq_len(ncol(ifactors))) {
            if ("asis" %in% default_order) {
               ifactors[,iCol] <- factor(ifactors[,iCol],
                  levels=unique(ifactors[,iCol]));
            } else if ("sort_samples" %in% default_order) {
               ifactors[,iCol] <- factor(ifactors[,iCol],
                  levels=sort_samples(unique(ifactors[,iCol]),
                     pre_control_terms=pre_control_terms,
                     ...));
            } else {
               ifactors[,iCol] <- factor(ifactors[,iCol],
                  levels=jamba::mixedSort(unique(ifactors[,iCol]),
                     ...));
            }
         }
      }
      if (length(group_colnames) > 0) {
         colnames(ifactors) <- jamba::makeNames(rep(group_colnames,
            length.out=ncol(ifactors)),
            renameFirst=FALSE);
      } else {
         colnames(ifactors) <- jamba::makeNames(
            rep("factor",
               length.out=ncol(ifactors)),
            renameOnes=TRUE);
      }
      if (length(rownames(ifactors)) == 0) {
         rownames(ifactors) <- jamba::makeNames(
            jamba::pasteByRow(ifactors, sep=factor_sep),
            suffix="_rep");
      }
      if (verbose) {
         jamba::printDebug("groups_to_sedesign(): ",
            "ifactors:");
         print(head(ifactors, 40));
      }

      # Assume sample rows and group columns
      rowGroups <- jamba::pasteByRowOrdered(ifactors, sep=factor_sep);
      sample2group <- split(rownames(ifactors), rowGroups);
      if (length(idesign) == 0) {
         idesign <- list2im_opt(sample2group, do_sparse=FALSE)[rownames(ifactors),levels(rowGroups),drop=FALSE];
      }
   } else if (jamba::igrepHas("data.frame|dataframe|matrix", class(ifactors))) {
      #####################################################
      ## data.frame input
      ##
      if (verbose) {
         jamba::printDebug("groups_to_sedesign(): ",
            "using existing data.frame");
      }
      if (length(rownames(ifactors)) == 0) {
         if (length(isamples) == 0) {
            ## Create isamples
            isamples <- jamba::makeNames(rep("sample", nrow(ifactors)));
         } else if (length(isamples) == nrow(ifactors)) {
            # use isamples as-is
         } else {
            stop(paste0("ifactors has no rownames, and ",
               "length(isamples) != nrow(ifactors). ",
               "Please make length(isamples) == nrow(iFactor)"));
         }
      } else {
         if (length(isamples) == 0) {
            isamples <- rownames(ifactors);
         } else {
            if (!any(isamples %in% ifactors) && length(isamples) == nrow(ifactors)) {
               ## use isamples as-is
            } else if (!all(isamples %in% rownames(ifactors))) {
               stop(paste0("isamples is not present in all rownames(ifactors). ",
                  "Either: all isamples must be present in rownames(ifactors); or ",
                  "no isamples are present in rownames(ifactors) and ",
                  "length(isamples) == nrow(ifactors)."));
            } else {
               ## Subset or re-order ifactors using matching isamples
               ifactors <- ifactors[match(isamples, rownames(ifactors)),,drop=FALSE];
               if (verbose) {
                  jamba::printDebug("groups_to_sedesign(): ",
                     "Specifying ifactors[isamples,]");
                  print(head(ifactors));
               }
            }
         }
         if (verbose >= 2) {
            jamba::printDebug("groups_to_sedesign(): ",
               "head(ifactors):");
            print(head(ifactors, 100));
         }
      }
      if (length(group_colnames) == 0) {
         if (length(colnames(ifactors)) == 0) {
            ## Create colnames
            group_colnames <- jamba::makeNames(
               renameOnes=TRUE,
               rep("factor",
                  length.out=ncol(ifactors)));
            colnames(ifactors) <- group_colnames;
         } else {
            group_colnames <- colnames(ifactors);
         }
      } else {
         if (!all(group_colnames %in% colnames(ifactors))) {
            stop(paste0("Not all group_colnames are in colnames(ifactors), please remedy."));
         }
         ## Use ifactors as-is
         #ifactors <- ifactors[,group_colnames,drop=FALSE];
      }
      if (verbose) {
         jamba::printDebug("groups_to_sedesign(): ",
            "Specifying ifactors[,group_colnames,drop=FALSE]");
         jamba::printDebug("groups_to_sedesign(): ",
            "group_colnames:",
            group_colnames);
      }

      # default_order == "asis" will convert character columns to factor
      #    using the observed order of terms as factor levels
      for (icol in group_colnames) {
         if (!"factor" %in% class(ifactors[,icol])) {
            if ("asis" %in% default_order) {
               if (verbose) {
                  jamba::printDebug("groups_to_sedesign(): ",
                     "Converting '", icol, "' to factor using default_order: ", "asis");
               }
               ifactors[,icol] <- factor(ifactors[,icol],
                  levels=unique(ifactors[,icol]),
                  exclude=NULL);
            } else if ("sort_samples" %in% default_order) {
               if (verbose) {
                  jamba::printDebug("groups_to_sedesign(): ",
                     "Converting '", icol, "' to factor using default_order: ", "sort_samples");
               }
               ifactors[,icol] <- factor(ifactors[,icol],
                  levels=sort_samples(unique(ifactors[,icol]),
                     pre_control_terms=pre_control_terms,
                     ...));
            } else {
               if (verbose) {
                  jamba::printDebug("groups_to_sedesign(): ",
                     "Converting '", icol, "' to factor using default_order: ", "mixedSort");
               }
               ifactors[,icol] <- factor(ifactors[,icol],
                  levels=jamba::mixedSort(unique(ifactors[,icol]),
                     ...));
            }
         }
      }
      # default_order == "mixedSort" will use alphanumeric sort
      # jamba::mixedSortDF() will honor factor level orders when present,
      # otherwise will use alphanumeric sort order.
      # To influence the sort order, use factors with ordered levels.
      ifactors <- jamba::mixedSortDF(ifactors,
         byCols=group_colnames);
      if (verbose >= 2) {
         jamba::printDebug("groups_to_sedesign(): ",
            "ifactors:");
         print(head(ifactors));
      }

      ## rowGroups is the unique set of group names, used to keep the original order
      #rowGroups <- jamba::pasteByRowOrdered(ifactors[,group_colnames,drop=FALSE],
      #   sep=factor_sep);
      ## Unclear whether to re-order columns to match group_colnames, for now we do not
      rowGroups <- jamba::pasteByRowOrdered(ifactors,
         sep=factor_sep);
      if (length(rownames(ifactors)) == 0) {
         iFactors_names <- jamba::makeNames(rowGroups,
            suffix="_rep");
         rownames(ifactors) <- iFactors_names;
      } else {
         iFactors_names <- rownames(ifactors);
      }
      ## Assume for now sample rows and group columns
      sample2group <- split(iFactors_names, rowGroups);
      if (length(idesign) == 0) {
         idesign <- list2im_opt(sample2group,
            do_sparse=FALSE)[iFactors_names, as.character(unique(rowGroups)), drop=FALSE];
         if (all(isamples %in% iFactors_names)) {
            idesign <- idesign[match(isamples, iFactors_names),,drop=FALSE];
         }
      } else {
         if (length(isamples) > 0) {
            idesign <- idesign[match(isamples, rownames(idesign)),,drop=FALSE];
         }
      }
   } else if (jamba::igrepHas("matrix", class(ifactors)) && all(c(0,1) %in% ifactors)) {
      ##################################
      ## idesign input
      ##
      if (verbose) {
         jamba::printDebug("groups_to_sedesign(): ",
            "converting idesign into ifactors data.frame");
      }
      ## Assume for now, idesign matrix with sample rows and group columns
      sample2group <- split(rownames(ifactors), sapply(seq_len(nrow(ifactors)), function(i){
         colnames(ifactors)[which(ifactors[i,] != 0)];
      }));
      idesign <- list2im_opt(sample2group, do_sparse=FALSE)[rownames(ifactors),names(sample2group)];
      iFactorsCols <- colnames(ifactors);
      ifactors <- jamba::rbindList(strsplit(iFactorsCols, factor_sep));
      if (!is.null(group_colnames)) {
         colnames(ifactors) <- jamba::makeNames(rep(group_colnames, length.out=ncol(ifactors)),
            renameFirst=FALSE);
      } else {
         colnames(ifactors) <- jamba::makeNames(
            rep("groupFactor",
               length.out=ncol(ifactors)),
            renameOnes=TRUE,
            suffix="_");
      }
      rownames(ifactors) <- unname(jamba::pasteByRow(ifactors, sep=factor_sep));
      jamba::printDebug("ifactors:");print(ifactors);
   }
   if (verbose >= 2) {
      jamba::printDebug("groups_to_sedesign(): ",
         "ifactors:");
      print(head(ifactors));
      if (!is.null(sample2group)) {
         jamba::printDebug("sample2group:");
         print(head(sample2group));
      }
   }

   ##########################################################
   ## Check to make sure no factor levels contain "-"
   for (i in colnames(ifactors)) {
      if (jamba::igrepHas("-", ifactors[,i])) {
         ifactors[,i] <- gsub("-", ".", ifactors[,i]);
      }
   }

   ##########################################################
   ## First check to make sure the ifactors values are unique
   ## and if not, use only unique entries
   iContrastGroupsUse <- colnames(ifactors);
   iFactorsV <- jamba::pasteByRow(ifactors, sep=factor_sep);
   iKeepRows <- match(unique(iFactorsV), iFactorsV);
   ifactors <- ifactors[iKeepRows,,drop=FALSE];
   if (rename_first_depth && current_depth==1) {
      rownames(ifactors) <- jamba::pasteByRow(ifactors, sep=factor_sep);
   }

   if (verbose >= 2) {
      jamba::printDebug("groups_to_sedesign(): ",
         "ifactors:");
      print(head(ifactors));
   }


   if (verbose) {
      jamba::printDebug("groups_to_sedesign(): ",
         "current_depth:",
         current_depth);
   }

   ##########################################################
   ## Iterate each factor in order, and create valid contrasts
   ## Note: we allow applying contrasts in a different order than the
   ## columns in iFactor, if !is.null(factor_order)
   ##
   # ensure factor_order only matches columns provided
   if (length(contrast_names) == 0) {
      factor_order <- factor_order[factor_order <= ncol(ifactors)]
      if (length(factor_order) == 0) {
         factor_order <- seq_along(colnames(ifactors));
      }
      # ensure max_depth is no larger than the number of factors
      max_depth <- min(c(max_depth, length(factor_order)))

      ##
      if (verbose) {
         jamba::printDebug("groups_to_sedesign(): ",
            "factor_order values:",
            colnames(ifactors)[factor_order]);
      }
      iContrastNames <- data.frame(check.names=FALSE,
         stringsAsFactors=FALSE,
         jamba::rbindList(lapply(factor_order, function(iChange){
            if (verbose) {
               jamba::printDebug("groups_to_sedesign(): ",
                  "factor_order iChange:",
                  colnames(ifactors)[iChange]);
            }
            iNoChange <- setdiff(seq_len(ncol(ifactors)), iChange);
            ## Optionally omit certain values from consideration,
            ## notably for "," or "-" which already contain changing factors
            iFactorUseRows <- jamba::unigrep(omit_grep, ifactors[,iChange]);

            if (length(iNoChange) == 0) {
               iSplit <- rep("", length(iFactorUseRows));
            } else {
               iSplit <- jamba::pasteByRowOrdered(ifactors[iFactorUseRows,iNoChange,drop=FALSE],
                  sep=factor_sep);
            }

            ## Split rows by constant values in non-changing factor columns
            iSplitL <- split(iFactorUseRows, iSplit);
            iSplitL <- iSplitL[lengths(iSplitL) > 1];
            ## Only consider contrasts when there are multiple rows
            if (length(iSplitL) > 0) {
               iDF <- jamba::rbindList(lapply(iSplitL, function(iSplitRows) {
                  use_factor_order <- unique(c(factor_order,
                     seq_len(ncol(ifactors))));
                  iFactorsSub <- ifactors[iSplitRows, use_factor_order, drop=FALSE];
                  if (verbose >= 2) {
                     jamba::printDebug("groups_to_sedesign(): ",
                        "   iSplitRows:",
                        iSplitRows,
                        ", use_factor_order:", use_factor_order);
                     jamba::printDebug("groups_to_sedesign(): ",
                        "   iFactorsSub:");
                     print(iFactorsSub);
                  }
                  iFactorVals <- iFactorsSub[,colnames(ifactors)[iChange]];
                  iMatch <- match(
                     sort_samples(iFactorVals,
                        pre_control_terms=pre_control_terms),
                     iFactorVals);
                  # 0.0.27.900: fix for one factor column input
                  if (length(iMatch) < 2) {
                     return(NULL)
                  }
                  iCombn <- combn(iMatch, 2);
                  iGrp1 <- ifelse(grepl("-", rownames(iFactorsSub)[iCombn[2,]]),
                     paste0("(", rownames(iFactorsSub)[iCombn[2,]], ")"),
                     rownames(iFactorsSub)[iCombn[2,]]);
                  iGrp2 <- ifelse(grepl("-", rownames(iFactorsSub)[iCombn[1,]]),
                     paste0("(", rownames(iFactorsSub)[iCombn[1,]], ")"),
                     rownames(iFactorsSub)[iCombn[1,]]);
                  iContrastName <- paste0(iGrp1, "-", iGrp2);
                  icondf <- iFactorsSub[intercalate(iCombn[2,], iCombn[1,]),,drop=FALSE];
                  iconfac <- factor(rep(iContrastName, each=2),
                     levels=unique(iContrastName));
                  iContrastDF <- data.frame(check.names=FALSE,
                     stringsAsFactors=FALSE,
                     lapply(jamba::nameVector(colnames(icondf)), function(i){
                        jamba::cPasteU(split(icondf[,i], iconfac))
                     }),
                     contrastName=iContrastName,
                     row.names=iContrastName);

                  # Create a string representing the combination of factors.
                  # which we will use to prevent re-creating the same contrasts.
                  #
                  # Modified the string to include colname, to ensure that two
                  # factors which may share some levels, will not be confused.
                  iContrastDF[,"contrastString"] <- jamba::pasteByRow(
                     iContrastDF[,colnames(iFactorsSub),drop=FALSE],
                     includeNames=TRUE,
                     sep=";",
                     sepName=":");
                  iContrastDF;
               }));
               rownames(iDF) <- iDF[,"contrastName"];
               if (verbose) {
                  jamba::printDebug("groups_to_sedesign(): ",
                     "   new contrasts:\n",
                     rownames(iDF),
                     sep=",\n");
               }
               iDF;
            } else {
               NULL;
            }
         })));

      ## Optionally spike in some pre-defined non-standard contrasts
      if (!is.null(add_contrastdf)) {
         if (verbose) {
            jamba::printDebug("groups_to_sedesign(): ",
               "Adding custom ",
               "add_contrastdf");
         }
         iContrastNames <- rbind(iContrastNames, add_contrastdf);
      }

      # Always make each row unique in terms of the factors compared.
      # Note: This step enforces order of comparison in two-way contrasts.
      # if (make_unique) {
      if (TRUE) {
         iDFcomponents <- jamba::pasteByRow(
            iContrastNames[,setdiff(colnames(iContrastNames), "contrastName"),drop=FALSE],
            sep="!");
         if (verbose >= 2) {
            jamba::printDebug("groups_to_sedesign(): ",
               "iDFcomponents:\n",
               iDFcomponents, sep="\n");
            jamba::printDebug("groups_to_sedesign(): ",
               "unique(iDFcomponents):\n",
               unique(iDFcomponents), sep="\n");
         }
         if (verbose && any(duplicated(iDFcomponents))) {
            dupe_comps <- iDFcomponents[duplicated(iDFcomponents)];
            dupe_kept_df <- data.frame(
               dupe_comp=iDFcomponents[iDFcomponents %in% dupe_comps],
               contrast=rownames(subset(iContrastNames, iDFcomponents %in% dupe_comps)),
               outcome=ifelse(!duplicated(iDFcomponents[iDFcomponents %in% dupe_comps]), "(kept)", "(removed)"))
            dupe_kept_df <- jamba::mixedSortDF(byCols=1, dupe_kept_df);
            jamba::printDebug("groups_to_sedesign(): ",
               "   removed duplicate (equivalent) contrasts:");
            print(dupe_kept_df[,-1, drop=FALSE]);
         }
         iContrastNames <- subset(iContrastNames, !duplicated(iDFcomponents));
      }

      if ("contrastName" %in% colnames(iContrastNames)) {
         if (verbose >= 2) {
            jamba::printDebug("groups_to_sedesign(): ",
               "tcount(iContrastNames$contrastName):")
            print(jamba::tcount(iContrastNames[,"contrastName"]));
         }
         rownames(iContrastNames) <- jamba::makeNames(iContrastNames[,"contrastName"]);
      }

      # Optionally remove contrasts with factor pairs in remove_pairs
      if (length(remove_pairs) > 0) {
         if (verbose) {
            jamba::printDebug("groups_to_sedesign(): ",
               "Processing any remove_pairs contrasts.");
         }
         for (iCol in setdiff(colnames(iContrastNames), "contrastName")) {
            if (verbose) {
               jamba::printDebug("groups_to_sedesign(): ",
                  "   Checking for remove_pairs in column:", iCol);
            }
            iColVals <- jamba::cPasteS(strsplit(as.character(iContrastNames[[iCol]]), ","));
            if (any(iColVals %in% remove_pairsFull)) {
               iWhich1 <- which(iColVals %in% remove_pairsFull);
               iWhich <- which(!iColVals %in% remove_pairsFull);
               if (verbose) {
                  jamba::printDebug("      removedPair with values:\n",
                     unique(iColVals[iWhich1]),
                     fgText=c("yellow", "purple"), sep="\n");
               }
               iContrastNames <- iContrastNames[iWhich,,drop=FALSE];
            }
         }
         if (nrow(iContrastNames) == 0) {
            warning("No contrasts remain after filtering remove_pairs.");
            return(NULL);
         }
      }

      if (verbose >= 2) {
         jamba::printDebug("groups_to_sedesign(): ",
            "iContrastNames:");
         print(head(iContrastNames, 100));
      }

      ##################################################
      # Interaction contrasts (iterative processing)
      if (length(setdiff(colnames(iContrastNames), "contrastName")) > 1 &&
            current_depth < max_depth) {
         iContrastNamesUse <- iContrastNames[,iContrastGroupsUse,drop=FALSE];
         for (i in iContrastGroupsUse) {
            j <- jamba::provigrep(c("^[^,]+$", "."), iContrastNamesUse[[i]]);
            iContrastNamesUse[[i]] <- factor(iContrastNamesUse[[i]],
               levels=unique(j));
         }
         if (verbose >= 2) {
            jamba::printDebug("groups_to_sedesign(): ",
               "   Defining interactions contrasts.");
            print(head(iContrastNamesUse[,iContrastGroupsUse,drop=FALSE], 100));
         }
         iContrastNamesInt <- groups_to_sedesign(iContrastNamesUse,
            omit_grep=omit_grep,
            current_depth=current_depth + 1,
            max_depth=max_depth,
            return_sedesign=FALSE,
            factor_sep=factor_sep,
            factor_order=rev(factor_order),
            contrast_sep=contrast_sep,
            rename_first_depth=rename_first_depth,
            remove_pairs=remove_pairs,
            pre_control_terms=pre_control_terms,
            verbose=verbose,
            ...);
         if (verbose >= 2) {
            jamba::printDebug("groups_to_sedesign(): ",
               "length(iContrastNamesInt):",
               length(iContrastNamesInt));
            print(iContrastNamesInt);
         }
         ## If length==0 then there are no valid interaction contrasts
         if (length(iContrastNamesInt) > 0 &&
               jamba::igrepHas("[(]", rownames(iContrastNamesInt[[1]]))) {
            return(iContrastNamesInt);
         }
         if (length(iContrastNamesInt) > 0 &&
               ncol(iContrastNamesInt) > 1 &&
               any(is.na(iContrastNamesInt[,1]))) {
            iContrastNamesInt <- iContrastNamesInt[!is.na(iContrastNamesInt[,1]),,drop=FALSE];
         }
         if (length(iContrastNamesInt) == 0 || ncol(iContrastNamesInt) > 1) {
            if (verbose >= 2) {
               jamba::printDebug("groups_to_sedesign(): ",
                  "begin iContrastNamesInt:");
               print(head(iContrastNamesInt));
               jamba::printDebug("  end iContrastNamesInt:");
            }
            iContrastNames <- jamba::rbindList(list(iContrastNames,
               iContrastNamesInt));
         }
      } else {
         if (verbose >= 2) {
            jamba::printDebug("groups_to_sedesign(): ",
               "   Skipping interactions");
            jamba::printDebug("      ncol(iContrastNames):",
               ncol(iContrastNames));
            jamba::printDebug("      head(iContrastNames):");
            print(head(iContrastNames));
         }
      }
      if ("contrastName" %in% colnames(iContrastNames)) {
         rownames(iContrastNames) <- jamba::makeNames(iContrastNames[["contrastName"]]);
         contrast_names <- unique(iContrastNames[["contrastName"]]);
      }
   }
   # end of automatic contrast definition
   ######################################################

   if (return_sedesign && current_depth == 1) {
      icontrasts <- NULL;
      if (!is.null(idesign) && length(contrast_names) > 0) {
         icontrasts <- limma::makeContrasts(contrasts=contrast_names,
            levels=idesign);
      }
      retvals <- validate_sedesign(
         new("SEDesign",
            design=idesign,
            contrasts=icontrasts));
   } else {
      retvals <- list();
      retvals$contrast_df <- iContrastNames;
      retvals$contrast_names <- contrast_names;
      retvals$idesign <- idesign;
   }
   return(retvals);
}



#' Sort biological sample labels for experimental design
#'
#' Sort biological sample labels for experimental design
#'
#' This function sorts a vector of sample labels using typical
#' heuristics that order typical control groups terms before
#' test groups. For example, `"Vehicle"` would be returned
#' before `"Treatment"` since `"Vehicle"` is a recognized control
#' term.
#'
#' It also employs `jamba::mixedSort()` for
#' proper alphanumeric sorting, for example so `"Time_5hr"` would
#' be sorted before `"Time_12hr"`.
#'
#' @return character vector ordered such that control terms are
#' preferentially first before non-control terms.
#'
#' @param x character vector or factor
#' @param control_terms vector of regular expression patterns used to
#'    determine control terms, where the patterns are matched and
#'    returned in order.
#' @param pre_control_terms vector or NULL, optional control
#'    terms or regular expressions to use before the `control_terms`
#'    above. This argument is used as a convenient prefix to the
#'    default terms.
#' @param post_control_terms vector or NULL, optional control
#'    terms or regular expressions to use after the `control_terms`
#'    above. This argument is used as a convenient suffix to the
#'    default terms.
#' @param ignore.case logical passed to `jamba::provigrep()` indicating
#'    whether to ignore case-sensitive matching.
#' @param boundary logical indicating whether to require a word
#'    boundary at either the start or end of the control terms.
#'    When TRUE, it uses `perl=TRUE` by default, and allows either
#'    perl boundary or an underscore `"_"`.
#' @param perl logical indicating whether to use Perl regular
#'    expression pattern matching.
#' @param keep_factor_order logical indicating whether to maintain
#'    factor level order, if `x` is supplied as a factor. If
#'    `keep_factor_order==TRUE` then only `sort(x)` is returned.
#' @param ... additional arguments are ignored.
#'
#' @family jam string functions
#' @family jam RNA-seq functions
#'
#' @examples
#' # the defaults perform well for clear descriptors
#' sort_samples(c("Trt_12h", "Trt_9h", "Trt_1h", "Trt_9h", "Vehicle"));
#'
#' # custom terms can be added before the usual control terms
#' sort_samples(c("Trt_12h", "Trt_9h", "Trt_1h", "Trt_9h", "Fixated", "Vehicle"),
#'    pre_control_terms="fixate");
#'
#' # custom terms can be added after the usual control terms
#' sort_samples(c("Trt_12h", "Trt_9h", "Trt_1h", "Trt_9h", "Fixated", "Vehicle"),
#'    post_control_terms="fixate");
#'
#' @export
sort_samples <- function
(x,
 control_terms=c(
    "WT|wildtype",
    "normal|healthy|healthycontrol|^hc$",
    "control|ctrl|ctl",
    "(^|[-_ ])(NT|NTC)($|[-_ ]|[0-9])",
    "none|empty|blank",
    "untreated|untrt|untreat",
    "Vehicle|veh",
    "ETOH|ethanol",
    #"(time|day|hour|min|minute)[s]{0,1}[0]",
    "scramble|mock|sham",
    "ttx|PBS",
    "knockout",
    "mutant"),
 sortFunc=jamba::mixedSort,
 pre_control_terms=NULL,
 post_control_terms=NULL,
 ignore.case=TRUE,
 boundary=TRUE,
 perl=boundary,
 keep_factor_order=TRUE,
 ...)
{
   ## Purpose is to order sample names by typical descriptions
   ## of control groups versus treatment groups
   ##
   ## Test set:
   ## sort_samples(c("Trt_12h", "Trt_9h", "Trt_1h", "Vehicle"))
   ## sort_samples(c("RA_Brg1", "EtOH_WT", "RA_WT", "EtOH_Brg1"))
   ## sort_samples(c("HCTWT_DXR6", "HCTWT_DXR12", "HCTWT_DXR24", "HCTWT_NT24"))
   #order1 <- jamba::proigrep(c(control_terms), x);
   #order2 <- jamba::proigrep(c(control_terms, "."), sortFunc=sortFunc, x);
   ##
   ## keep_factor_order=TRUE will keep factor levels unchanged, and use those levels in the sort
   ## instead of looking for control terms
   if (keep_factor_order && jamba::igrepHas("factor", class(x))) {
      sort(x);
   } else {
      control_terms <- unique(c(pre_control_terms,
         control_terms,
         post_control_terms));
      if (any(boundary)) {
         # Require regular expression boundary
         control_terms1 <- unlist(lapply(control_terms, function(i){
            paste0("(_|\\b)(", i, ")|(", i, ")(_|\\b)")
         }))
         if (any(!boundary)) {
            control_terms <- c(control_terms1,
               control_terms);
         } else {
            control_terms <- control_terms1;
         }
      }
      xU <- jamba::provigrep(c(control_terms, "."),
         sortFunc=sortFunc,
         perl=perl,
         ignore.case=ignore.case,
         x);
      xOrder <- order(match(x, xU));
      x <- x[xOrder];
      #attr(x, "control_terms") <- control_terms;
      x;
   }
}

#' Split the elements of an ordered factor vector
#'
#' Split the elements of an ordered factor vector
#'
#' This function performs `base::strsplit()` while trying to maintain
#' the order of factor levels in the output, based upon the order of
#' factor levels in the input data.
#'
#' @return list of factor vectors, where each factor shares the same
#'    global factor levels based upon the input data.
#'
#' @param x character or factor vector.
#' @param split character split value sent to `base::strsplit()`.
#' @param fixed,perl,useBytes additional arguments sent to `base::split()`.
#' @param sortFunc function used to sort character values when the input
#'    `x` is a character vector. The default `jamba::mixedSort()` applies
#'    alphanumeric sort.
#' @param keepOrder logical indicating whether to keep the order of values
#'    in the input data, for example with character input the values will
#'    be ordered by the first appearance of each term.
#' @param ... additional arguments are ignored.
#'
#' @family jamses utilities
#'
#' @examples
#' # first define a vector of sample groups
#' iGroups <- jamba::nameVector(paste(rep(c("WT", "KO"), each=6),
#'    rep(c("Control", "Treated"), each=3),
#'    sep="_"));
#' iGroups <- factor(iGroups, levels=unique(iGroups));
#' iGroups;
#' strsplitOrdered(iGroups, "_");
#'
#' @export
strsplitOrdered <- function
(x,
 split="_",
 fixed=FALSE,
 perl=FALSE,
 useBytes=FALSE,
 sortFunc=jamba::mixedSort,
 keepOrder=TRUE,
 ...)
{
   ## Purpose is to run strsplit() on factors, ordering the new factor
   ## levels consistent with the input
   if (!jamba::igrepHas("factor", class(x))) {
      if (keepOrder) {
         x <- factor(x,
            levels=unique(x));
      } else {
         x <- factor(x,
            levels=sortFunc(unique(x)));
      }
   }
   soL <- strsplit(x=levels(x),
      split=split,
      fixed=fixed,
      perl=perl,
      useBytes=useBytes);
   so1 <- jamba::rbindList(soL);

   ## Note: the setdiff() is there to remove "" values
   so1levels <- setdiff(unique(unlist(apply(so1, 2, unique))), "");
   soSplitL <- strsplit(as.character(x),
      split=split,
      fixed=fixed,
      perl=perl,
      useBytes=useBytes);
   soLordered <- lapply(soSplitL, function(i){
      factor(i,
         levels=so1levels);
   });
   return(soLordered);
}


#' Intercalate two or more vectors
#'
#' @export
intercalate <- function
(...)
{
   ## Purpose is to take a list of vectors, and intercalate their values, e.g.
   ## list1 <- paste("name1", letters[1:10], sep="");
   ## list2 <- paste("name2", letters[1:10], sep="");
   ## intercalate(list1, list2);
   ## name1a, name2a, name1b, name2b, name1c, name2c, etc.
   ##
   ## The special case where there are two lists, and the first has
   ## one element more than the second, then the second will only have
   ## its values in between the first, e.g.
   ## A B A B A B A
   ##
   ## Note: rmNULL() will remove empty lists
   aList <- jamba::rmNULL(list(...));
   if (length(aList) == 1 && class(aList[[1]]) %in% "list") {
      aList <- aList[[1]];
   }
   ## do.call will automatically repeat any vector to fill each row
   ## up to the maximum number of columns.
   if (length(unique(lengths(aList))) > 1) {
      ## Unequal lengths, to avoid warning should we expand them?
   }
   aMatrix <- do.call(rbind, aList);
   newVector <- as.vector(aMatrix);

   ## The special case where intercalating two vectors,
   ## where the second vector has one fewer entry, we
   ## will not repeat the last entry.
   ## E.g.
   ## c("A","A","A")
   ## c("B","B")
   ##
   ## desired output is
   ## c("A","B","A","B","A")
   if (length(aList) == 2 && length(aList[[1]]) == (length(aList[[2]]) + 1)) {
      newVector <- head(newVector, -1);
   }
   return(newVector);
}


#' Convert list to incidence matrix
#'
#' Convert list to incidence matrix
#'
#' @param setlist `list` of vectors
#' @param empty default single value used for empty/missing entries.
#' @param do_sparse `logical` indicating whether to convert output
#'    to `ngCMatrix` which is best for extremely large incidence
#'    matrix data.
#' @param ... additional arguments are ignored.
#'
#' @export
list2im_opt <- function
(setlist,
 empty=0,
 do_sparse=TRUE,
 ...)
{
   setnamesunion <- Reduce("union", setlist);
   if (length(empty) == 0) {
      empty <- NA;
   } else {
      empty <- head(empty, 1);
   }
   setlistim <- do.call(cbind, lapply(setlist, function(i){
      i_match <- match(i, setnamesunion);
      j <- rep(empty,
         length(setnamesunion));
      j[i_match] <- 1;
      j;
   }))
   rownames(setlistim) <- setnamesunion;
   if (do_sparse && suppressPackageStartupMessages(require(Matrix))) {
      setlistim <- as(setlistim, "ngCMatrix");
   }
   return(setlistim);
}
