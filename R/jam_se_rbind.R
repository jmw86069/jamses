
#' Combine SummarizedExperiment objects by row
#'
#' Combine SummarizedExperiment objects by row, using `rbind()` logic.
#'
#' This function is intended to help the process of calling
#' `SummarizedExperiment::rbind()`.
#'
#' The process:
#'
#' 1. Convert `colnames()` for each entry in `se_list` using `colnames_from`
#' and `colnames_to`. This step is useful when each object in `se_list`
#' may be using a different set of `colnames()`. For example `"sample_p_12"`
#' and `"sample_n_12"` might be equivalent, so renaming them
#' with `colnames_from=c("_[np]_")` and `colnames_to=c("_X_")` would convert
#' both values to `"sample_X_12"`.
#' 2. Subset each object in `se_list` using only shared `colnames()`.
#' 3. Determine how to handle `colData()` columns that are not identical:
#'
#'    * `colData_action="identical"`: will only keep columns whose values
#'    are identical across all objects in `se_list`.
#'    * `colData_action="all"`: will keep columns in `colData()`, however
#'    non-identical columns will be converted to `character` and values
#'    will be comma-delimited.
#'
#' 4. Perform `rbind()`.
#'
#' ## TODO:
#'
#' * Write equivalent `se_cbind()` - it will wait until there is a
#' driving use case.
#' * Consider retaining only shared `assayNames()` across `se_list`.
#' * Consider optionally retaining user-defined `assayNames()`.
#' (Alternatively, the user can subset the assayNames upfront, though
#' it might be tedious). The recommended pattern in that case:
#'
#' ```
#' se <- se_rbind(
#'    se_list=lapply(se_list, function(se){
#'       assays(se) <- assays(se)[assay_names];
#'       return(se)
#'    })
#' )
#' ```
#'
#' @family jamses utilities
#'
#' @returns `SummarizedExperiment` object whose `colData()` has been
#'    processed according to `colData_action` - either keeping only
#'    columns with identical values, or keeping all values delimited
#'    as a `character` string when values differ.
#'
#' @param se_list `list` of `SummarizedExperiment` objects.
#' @param colnames_from `character` vector of patterns used with `gsub()`
#'    to convert `colnames()` for each object in `se_list` to an identifier
#'    that will be shared across all objects in `se_list`.
#' @param colnames_to `character` vector of replacements used with `gsub()`
#'    alongside each entry in `colnames_from` to convert `colnames()`
#'    for each object in `se_list` to an identifier
#'    that will be shared across all objects in `se_list`.
#' @param colData_action `character` string indicating the action used to
#'    combine `colData()` across `se_list`:
#'    * `"identical"`: retain only those columns in `colData()` which are
#'    identical in all `se_list` objects.
#'    * `"all"`: retain all columns, but convert columns with mismatched
#'    values to store comma-delimited values.
#' @param colData_sep `character` string used as delimiter when
#'    `colData_action="all"` and when values in a column in `colData()`
#'    differs across objects in `se_list`. Only values that differ are
#'    delimited, to minimize redundancy.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' m1 <- matrix(rnorm(100), ncol=10);
#' colnames(m1) <- paste0("sample_p_", 1:10);
#' rownames(m1) <- paste0("row_", 1:10);
#' m2 <- matrix(rnorm(100), ncol=10);
#' colnames(m2) <- paste0("sample_n_", 1:10);
#' rownames(m2) <- paste0("row_", 11:20);
#' sample_id <- gsub("_[np]_", "_X_", colnames(m1));
#' m1
#' m2
#' se1 <- SummarizedExperiment::SummarizedExperiment(
#'    assays=list(counts=m1),
#'    rowData=data.frame(measurement=rownames(m1)),
#'    colData=data.frame(sample=colnames(m1),
#'       sample_id=sample_id))
#' se2 <- SummarizedExperiment::SummarizedExperiment(
#'    assays=list(counts=m2),
#'    rowData=data.frame(measurement=rownames(m2)),
#'    colData=data.frame(sample=colnames(m2),
#'       sample_id=sample_id))
#' # this step fails because colnames are not shared
#' # do.call(SummarizedExperiment::rbind, list(se1, se2))
#'
#' # keep only identical colData columns
#' se12 <- se_rbind(list(se1, se2))
#' SummarizedExperiment::colData(se12)
#'
#' # keep all colData columns
#' se12all <- se_rbind(list(se1, se2),
#'    colData_action="all")
#' SummarizedExperiment::colData(se12all)
#'
#' @export
se_rbind <- function
(se_list,
 colnames_from="_(n|p|neg|pos)_",
 colnames_to="_X_",
 colnames_keep=NULL,
 colData_action=c("identical",
    "all"),
 colData_sep=";",
 verbose=FALSE,
 ...)
{
   # edit colnames() if appropriate
   colData_sep <- head(colData_sep, 1);
   if (length(colData_sep) == 0) {
      colData_sep <- ";";
   }
   if (length(colnames_from) > 0 && any(nchar(colnames_from) > 0)) {
      colnames_from <- colnames_from[nchar(colnames_from) > 0];
      if (verbose) {
         jamba::printDebug("se_rbind(): ",
            "Processing colnames_from/colnames_to.")
      }
      se_list <- lapply(se_list, function(se){
         colnames(se) <- jamba::gsubs(colnames_from,
            colnames_to,
            colnames(se))
         se;
      })
   }

   # determine shared colnames
   se_colnames <- Reduce("intersect", lapply(se_list, colnames));
   if (length(colnames_keep) > 0) {
      se_keep_colnames <- intersect(colnames_keep,
         se_colnames);
      if (verbose) {
         jamba::printDebug("se_rbind(): ",
            "Subsetted ",
            jamba::formatInt(length(se_colnames)),
            " colnames to ",
            jamba::formatInt(length(se_keep_colnames)),
            " using colnames_keep.")
      }
      se_colnames <- se_keep_colnames;
   }
   if (length(se_colnames) == 0) {
      stop(paste0("No colnames were shared in se_list.\n",
         jamba::cPaste(sep="\n",
            sapply(se_list, function(se){
               jamba::cPaste(sep=", ", head(colnames(se)))
            })
         )));
   }

   # subset se_list by se_colnames
   se_list <- lapply(se_list, function(se){
      se[, se_colnames];
   })

   # validation check to ensure no rownames are shared
   se_rownames <- Reduce("intersect", lapply(se_list, rownames));
   if (length(se_rownames) > 0) {
      jamba::printDebug("se_rbind(): ",
         "There are ", jamba::formatInt(length(se_rownames)),
         " rownames shared across se_list, which is not permitted.");
      stop("rownames cannot be shared across se_list objects.");
   }

   # determine which colData colnames are identical
   se_colData_colnames <- Reduce("intersect", lapply(se_list, function(se){
      colnames(SummarizedExperiment::colData(se))
   }))
   if (length(se_colData_colnames) == 0) {
      jamba::printDebug("se_rbind(): ",
         "No colData colnames are shared across ", "se_list.")
   } else {
      # detect identity and process the data while it is convenient
      se_colData <- lapply(jamba::nameVector(se_colData_colnames),
         function(icol) {
            vals <- NULL;
            all_identical <- TRUE;
            # iterate each se_list object
            for (se_num in seq_along(se_list)) {
               if (length(vals) == 0) {
                  # note use unname() to prevent names from causing a mismatch
                  vals <- unname(
                     SummarizedExperiment::colData(
                        se_list[[se_num]])[[icol]]);
               } else {
                  all_identical <- all.equal(
                     unname(
                        SummarizedExperiment::colData(
                           se_list[[se_num]])[[icol]]),
                     vals)
                  if (!TRUE %in% all_identical) {
                     break;
                  }
               }
            }
            if (verbose) {
               jamba::printDebug("se_rbind(): ",
                  "colData column '", icol, "' identical values: ",
                  all_identical)
            }
            if (TRUE %in% all_identical) {
               return(vals);
            }
            if ("identical" %in% colData_action) {
               return(NULL);
            }
            if (verbose) {
               jamba::printDebug("se_rbind(): ",
                  "Merging values.", indent=6)
            }
            for (se_num in tail(seq_along(se_list), -1)) {
               newvals <- unname(
                  SummarizedExperiment::colData(
                     se_list[[se_num]])[[icol]]);
               diff_val <- (!vals == newvals)
               # delimit multiple values only when the new value differs
               if (any(diff_val)) {
                  if (!is.character(vals)) {
                     vals <- as.character(vals)
                  }
                  vals[diff_val] <- paste0(vals[diff_val],
                     colData_sep,
                     as.character(newvals[diff_val]))
               }
            }
            return(vals);
         })
      # now apply to each se_list entry
      se_list <- lapply(se_list, function(se){
         for (icol in names(se_colData)) {
            SummarizedExperiment::colData(se)[[icol]] <- se_colData[[icol]];
         }
         se;
      })
   }

   # now perform rbind
   se <- do.call(SummarizedExperiment::rbind, se_list);
   # now order colnames to match se_colnames
   se[, se_colnames]
}
