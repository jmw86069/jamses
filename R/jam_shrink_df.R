
#' Shrink data.frame by row groups
#'
#' Shrink data.frame by row groups
#'
#' Purpose is to shrink a `data.frame` to have one row per row grouping.
#' The row grouping can use a single column of identifiers, or multiple
#' columns. The challenge is to apply a relevant function to each column,
#' expecting there will be columns with `numeric`, `character`, or `factor`
#' types.
#'
#' The default behavior:
#'
#' * `numeric` columns are summarized with `mean(x, na.rm=TRUE)`, so that
#' NA values are ignored when there are non-NA values present.
#' * `character` columns are combined using unique, sorted `character` strings.
#'
#'    * This step uses `jamba::cPasteSU()` where the
#'    `S` activates sorting using `jamba::mixedSort()`, and
#'    `U` calls `unique()`.
#'    * To retain all values, remove the `U` and call `jamba::cPasteS()`
#'    * To skip the sort, remove the `S` and call `jamba::cPasteU()`
#'    * To keep all values, and skip sorting, call `jamba::cPaste()`
#'
#' @family jamses utilities
#'
#' @param x `data.frame` (or equivalent)
#' @param groupBy `character` vector with one of the following:
#'    * one or more columns in `colnames(x)`. The values in these columns
#'    will define the row groups used.
#'    * `character` or `factor` with length equal to `nrow(x)`. These
#'    values will define the row groups used.
#' @param string_func `function`, default uses `jamba::cPasteSU()`,
#'    used for `character` or `factor` columns. Note that string
#'    columns are handled differently than `numeric` columns by applying
#'    vectorized operations across the complete set of rows in one step,
#'    rather than calling `data.table` on each subgroup.
#' @param num_func `function`, default `function(x)mean(x, na.rm=TRUE)`,
#'    used for `numeric` columns. Note that this function is applied
#'    to each row group by `data.table`, and is typically very efficient
#'    for `numeric` values.
#' @param add_string_cols `character` with optional `numeric` columns that
#'    should be handled as if they were `character` columns. Default `NULL`.
#' @param num_to_string_func `function` used for `add_string_cols` when converting
#'    `numeric` columns to `character`. Default `as.character()` retains
#'    the full `numeric` value, however it may be useful to use something
#'    like `function(x)signif(x, digits=3)` to limit the output to
#'    only three significant digits, or `function(x)format(x, digits=3)`.
#' @param keep_na_groups `logical`, default TRUE, whether to convert `NA`
#'    values in row groups to `""` so they are retained in the output.
#'    * You may want to use `keep_na_groups=FALSE` when there are a large
#'    number of un-annotated rows that should not be aggregated together.
#'    This situation may occur if converting a probe to a gene symbol,
#'    where a subset of probes cannot be converted to a gene symbol
#'    and instead receive `NA`.
#' @param include_num_reps `logical` indicating whether to add a column
#'    `"num_reps"` to the output, with the `integer` number of rows
#'    in each row group.
#' @param collapse_method `integer` default 2, indicating the internal
#'    collapse method used. Experimental.
#'    * `1` collapses each `numeric` column independently.
#'    * `2` collapses each set of `numeric` columns that use the same
#'    numeric shrink function. When all `numeric` columns use the same
#'    shrink function, they are all calculated in a single step, which
#'    is typically much faster.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' testdf <- data.frame(check.names=FALSE,
#'    SYMBOL=rep(c("ACTB", "GAPDH", "PPIA"), c(2, 3, 1)),
#'    `logFC B-A`=c(1.4, 1.4, 2.3, NA, 2.5, 5.1),
#'    probe=paste0("probe", 1:6))
#' shrink_df(testdf, by="SYMBOL")
#'
#' shrink_df(testdf, by="SYMBOL", num_func=mean)
#'
#' shrink_df(testdf, by="SYMBOL", add_string_cols="logFC B-A")
#'
#' testdftall <- do.call(rbind, lapply(1:10000, function(i){
#'    idf <- testdf;
#'    idf$SYMBOL <- paste0(idf$SYMBOL, "_", i);
#'    idf;
#' }))
#' shrunk_tall <- shrink_df(testdftall,
#'    by="SYMBOL")
#' head(shrunk_tall, 6)
#'
#' shrunk_tall2 <- jamses::shrinkDataFrame(testdftall,
#'    groupBy="SYMBOL")
#' head(shrunk_tall2, 6)
#'
#' @export
shrinkDataFrame <- function
(x,
 groupBy,
 na.rm=TRUE,
 string_func=function(x)jamba::cPasteSU(x, na.rm=TRUE),
 num_func=function(x){mean(x, na.rm=TRUE)},
 add_string_cols=NULL,
 num_to_string_func=as.character,
 keep_na_groups=TRUE,
 include_num_reps=FALSE,
 collapse_method=2,
 verbose=FALSE,
 ...)
{
   #
   ## Recognize some legacy arguments for backward compatibility
   oldlist <- list(...);
   if (length(oldlist) > 0) {
      if ("stringShrinkFunc" %in% names(oldlist)) {
         string_func <- oldlist$stringShrinkFunc;
         if (verbose) {
            jamba::printDebug("shrinkDataFrame(): ",
               "recognized stringShrinkFunc to use for string_func");
         }
      }
      if ("numShrinkFunc" %in% names(oldlist)) {
         num_func <- oldlist$numShrinkFunc;
         if (verbose) {
            jamba::printDebug("shrinkDataFrame(): ",
               "recognized numShrinkFunc to use for num_func");
         }
      }
      if ("addStringCols" %in% names(oldlist)) {
         add_string_cols <- oldlist$addStringCols;
         if (verbose) {
            jamba::printDebug("shrinkDataFrame(): ",
               "recognized addStringCols to use for add_string_cols");
         }
      }
      if ("includeNumReps" %in% names(oldlist)) {
         include_num_reps <- oldlist$includeNumReps;
         if (verbose) {
            jamba::printDebug("shrinkDataFrame(): ",
               "recognized includeNumReps to use for include_num_reps");
         }
      }
   }

   x <- x;
   if (!"data.frame" %in% class(x)) {
      if (verbose) {
         jamba::printDebug("shrinkDataFrame(): ",
            "coerced with data.frame(check.names=FALSE, x).");
      }
      x <- data.frame(check.names=FALSE,
         x);
   }

   #####################################
   ## Define groupBy
   groupByCols <- NULL;
   if (length(groupBy) == nrow(x)) {
      ## groupBy use as-is
      ## If factor, make sure NA is a level if there are NAs
      if (jamba::igrepHas("factor", class(groupBy))) {
         groupBy <- addNA(groupBy, ifany=TRUE);
      }
      if (!keep_na_groups && any(is.na(as.character(groupBy)))) {
         if (verbose) {
            jamba::printDebug("shrinkDataFrame(): ",
               "dropping NA groups.");
         }
         nonNaGroupBy <- !is.na(as.character(groupBy));
         groupBy <- droplevels(groupBy[nonNaGroupBy]);
         x <- x[nonNaGroupBy, , drop=FALSE];
      }
      groupByDF <- data.frame(groupBy=unique(groupBy),
         row.names=jamba::rmNA(naValue="NA", as.character(unique(groupBy))));
      if (verbose) {
         jamba::printDebug("shrinkDataFrame(): ",
            "head(groupBy, 50):");
         print(head(groupBy, 50));
      }
   } else {
      if (any(groupBy %in% colnames(x))) {
         groupByCols <- intersect(groupBy, colnames(x));
         if (FALSE %in% keep_na_groups) {
            # check for rows that are entirely NA
            all_na <- Reduce("&", lapply(groupByCols, function(icol){
               is.na(x[[icol]])
            }));
            if (any(all_na)) {
               x <- subset(x, !all_na)
               if (verbose) {
                  jamba::printDebug("shrinkDataFrame(): ",
                     "removed ", jamba::formatInt(sum(all_na)),
                     " rows with NA group rows due to",
                     "keep_na_groups=FALSE");
               }
            }
         }
         groupBy <- jamba::pasteByRow(
            x[, groupByCols, drop=FALSE], sep="::");
         if (verbose) {
            jamba::printDebug("shrinkDataFrame(): ",
               "head(groupByCols):",
               head(groupByCols));
            jamba::printDebug("shrinkDataFrame(): ",
               "head(groupBy, 50):");
            print(head(groupBy, 50));
         }
         matchg <- match(unique(groupBy), groupBy);
         groupByDF <- data.frame(check.names=FALSE,
            stringsAsFactors=FALSE,
            x[matchg, groupByCols, drop=FALSE],
            row.names=unique(groupBy));
      } else {
         stop(paste0("groupBy must either be a vector equal to nrow(x),",
            " or a vector of colnames contained in x."));
      }
   }
   if (any(is.na(as.character(groupBy)))) {
      groupBy <- jamba::rmNA(groupBy, naValue="NA");
   }

   #########################################################
   ## iCols contains all columns which are not groupByCols,
   iCols <- setdiff(colnames(x), groupByCols);
   if (verbose) {
      jamba::printDebug("shrinkDataFrame(): ",
         "head(iCols):",
         head(iCols));
   }

   #########################################################
   ## Determine the column classes
   colClasses <- jamba::sclass(x);
   #colClasses <- sapply(jamba::nameVector(iCols), function(iCol){
   #   class(x[[iCol]]);
   #});

   ## Detect numeric columns
   numCols <- names(colClasses)[sapply(colClasses, function(iClass){
      jamba::igrepHas("integer|numeric|float", iClass);
   })];
   ## Detect string/character
   stringCols <- names(colClasses)[sapply(colClasses, function(iClass){
      jamba::igrepHas("character|factor|ordered|logical", iClass);
   })];
   ## All other classes are ignored
   if (verbose) {
      jamba::printDebug("shrinkDataFrame(): ",
         "head(numCols, 20):",
         head(numCols, 20));
      jamba::printDebug("shrinkDataFrame(): ",
         "head(stringCols, 20):",
         head(stringCols, 20));
   }

   #########################################################
   ## Optionally move numeric columns to be handled as string columns
   if (any(add_string_cols %in% numCols)) {
      moveNumCols <- intersect(add_string_cols, numCols);
      numCols <- setdiff(numCols, moveNumCols);
      stringCols <- c(stringCols, moveNumCols);
      if (verbose) {
         jamba::printDebug("shrinkDataFrame(): ",
            sep="",
            c("Moved columns: {",
            jamba::cPaste(moveNumCols),
            "} from numCols to stringCols."))
      }
      if (length(num_to_string_func) == 1 && is.function(num_to_string_func)) {
         for (icol in moveNumCols) {
            x[[icol]] <- num_to_string_func(x[[icol]])
         }
         if (verbose) {
            jamba::printDebug("shrinkDataFrame(): ",
               "Applied num_to_string_func().")
         }
      }
   }

   #########################################################
   ## Define the num_func list
   if (length(numCols) > 0) {
      if (!is.list(num_func)) {
         num_func <- jamba::nameVector(rep(list(num_func), length(numCols)),
            numCols);
      } else {
         if (!all(numCols %in% names(num_func))) {
            ## Ensure we keep direct assignments to named columns,
            ## then reuse num_func for subsequent columns
            ## whose names are not in names(num_func).
            ## After all, we subset for the observed numCols to
            ## allow for num_func to contain more entries than
            ## we have numCols.
            iNewCols <- setdiff(numCols, names(num_func));
            num_func <- c(num_func,
               jamba::nameVector(
                  rep(num_func, length.out=length(iNewCols)),
                  iNewCols));
         }
         num_func <- num_func[numCols];
      }
      if (verbose) {
         jamba::printDebug("shrinkDataFrame(): ",
            "numCols:",
            numCols,
            ", names(num_func):",
            names(num_func));
      }
   }

   #########################################################
   ## Define the string_func list
   if (length(stringCols) > 0) {
      if (!is.list(string_func)) {
         string_func <- jamba::nameVector(
            rep(list(string_func),
               length(stringCols)),
            stringCols);
      } else if (!all(stringCols %in% names(string_func))) {
         iNewCols <- setdiff(stringCols,
            names(string_func));
         string_func <- c(string_func,
            jamba::nameVector(
               rep(string_func,
                  length.out=length(iNewCols)),
               iNewCols));
      }
      if (verbose) {
         jamba::printDebug("shrinkDataFrame(): ",
            "stringCols:",
            stringCols,
            ", names(string_func):",
            names(string_func));
      }
   }

   #########################################################
   ## Count the number of entries grouped (optional)
   if (include_num_reps) {
      if (verbose) {
         jamba::printDebug("shrinkDataFrame(): ",
            "including num_reps");
      }
      groupByT <- table(groupBy);
      iMatch <- match(unique(groupBy), names(groupByT));
      num_repsDF <- data.frame(num_reps=as.vector(groupByT[iMatch]));
      groupByDF <- cbind(groupByDF,
         num_repsDF);
   }

   #####################################
   ## Shrink numeric columns
   if (length(numCols) > 0) {
      if (verbose) {
         jamba::printDebug("shrinkDataFrame(): ",
            "shrinking num columns:",
            numCols);
      }
      ## Collapse them all at once if all num_func are identical
      ## TODO: split the num_func by digest::digest() has value
      ## so sets of identical functions will be performed together.
      if (collapse_method == 2) {
         #num_func
         numShrinkFuncD <- sapply(num_func,
            digest::digest,
            algo="md5");
         numShrinkFuncDL <- split(names(numShrinkFuncD), numShrinkFuncD);
         if (verbose > 1) {
            jamba::printDebug("shrinkDataFrame(): ",
               "collapse_method:", collapse_method,
               " running on sets of identical num_func functions.");
            jamba::printDebug("shrinkDataFrame(): ",
               "num_func: ");
            print(num_func);
            jamba::printDebug("shrinkDataFrame(): ",
               "numShrinkFuncD: ");
            print(numShrinkFuncD);
            jamba::printDebug("shrinkDataFrame(): ",
               "numShrinkFuncDL: ");
            print(numShrinkFuncDL);
         }
         numShrunk <- do.call(cbind, lapply(numShrinkFuncDL, function(iCol){
            if (verbose > 1) {
               jamba::printDebug("shrinkDataFrame(): ",
                  indent=5,
                  "iCol:",
                  iCol,
                  fgText=c("orange", "lightgreen"));
               # jamba::printDebug("   shrinkDataFrame(): ",
               #    "colnames(x): ",
               #    colnames(x));
               jamba::printDebug("shrinkDataFrame(): ",
                  indent=5,
                  "head(data.frame):");
               print(head(
                  data.frame(x[, iCol, drop=FALSE],
                     check.names=FALSE,
                     stringsAsFactors=FALSE)
               ));
            }
            grOLi <- shrink_matrix(
               data.frame(x[, iCol, drop=FALSE],
                  check.names=FALSE,
                  stringsAsFactors=FALSE),
               groupBy=groupBy,
               shrink_func=num_func[[head(iCol, 1)]],
               return_class="matrix");
         }));
      } else {
         allNumIdentical <- all(sapply(num_func[-1], function(x){
            identical(num_func[[1]], x);
         }));
         if (allNumIdentical) {
            if (verbose > 1) {
               jamba::printDebug("shrinkDataFrame(): ",
                  "collapse_method:", collapse_method,
                  " running only once for all numeric values.");
            }
            numShrunk <- shrink_matrix(
               as.data.frame(x[, numCols, drop=FALSE]),
               groupBy=jamba::rmNA(naValue="NA", groupBy),
               shrink_func=num_func[[1]],
               return_class="matrix");
         } else {
            if (verbose > 1) {
               jamba::printDebug("shrinkDataFrame(): ",
                  "collapse_method:", collapse_method,
                  " running on each numeric column.");
            }
            numShrunk <- do.call(cbind,
               lapply(jamba::nameVector(numCols), function(iCol){
                  if (verbose) {
                     jamba::printDebug("  - ", iCol,
                        fgText=c("orange", "lightgreen"));
                     print(num_func[[iCol]]);
                  }
                  grOLi <- shrink_matrix(
                     as.data.frame(x[, iCol, drop=FALSE]),
                     groupBy=groupBy,
                     shrink_func=num_func[[iCol]], return_class="matrix");
               }));
         }
      }
      numShrunkDF <- data.frame(check.names=FALSE,
         stringsAsFactors=FALSE,
         numShrunk);
      if (verbose > 1) {
         jamba::printDebug("shrinkDataFrame(): ",
            "head(numShrunkDF):");
         print(head(numShrunkDF));
      }
      iMatch <- match(unique(groupBy), rownames(numShrunkDF));
      icolnames1 <- colnames(groupByDF);
      groupByDF <- cbind(groupByDF,
         numShrunkDF[iMatch, , drop=FALSE]);
      colnames(groupByDF) <- c(icolnames1,
         colnames(numShrunkDF));
   }

   #####################################
   ## Shrink string columns
   if (length(stringCols) > 0) {
      if (verbose) {
         jamba::printDebug("shrinkDataFrame(): ",
            "Shrinking string columns:",
            stringCols);
      }
      options("stringsAsFactor"=FALSE);
      stringShrunkDF <- as.data.frame(
         lapply(jamba::nameVector(stringCols), function(iCol){
            if (verbose > 1) {
               jamba::printDebug("   shrinkDataFrame(): ",
                  "string_func(x[['", iCol, "']])",
                  fgText=c("orange", "lightgreen"));
            }
            ## split vector into a list, then run the function on the list
            ## which is intended to generate a vector of values
            f1 <- x[[iCol]];
            iValsL <- string_func[[iCol]](split(as.character(f1),
               groupBy));
            if (jamba::igrepHas("factor", class(x[[iCol]]))) {
               if (verbose > 1) {
                  jamba::printDebug("      shrinkDataFrame(): ",
                     "shrinking factor column values.");
               }
               ## Generate new factor levels in the original order, followed by
               ## sorted values, then keeping only those levels contained in
               ## the resulting data
               iLevels <- intersect(unique(c(levels(f1), mixedSort(iValsL))),
                  iValsL);
               if (!NA %in% levels(x[[iCol]])) {
                  ## Bonus points: if NA was not in the levels previously, keep it out
                  ## otherwise NA will be a formal factor level
                  iLevels <- jamba::rmNA(iLevels);
               }
               iValsL <- factor(iValsL, levels=iLevels);
            }
            iValsL;
         }));
      colnames(stringShrunkDF) <- stringCols;
      iMatch <- match(unique(groupBy),
         rownames(stringShrunkDF));
      groupByDF <- cbind(groupByDF,
         stringShrunkDF[iMatch, , drop=FALSE]);
   }
   ## keep_cols is used to protect against columns whose types are
   ## not considered num or string, and may be lost above
   keep_cols <- intersect(colnames(x),
      colnames(groupByDF));
   # if (apply_keep_cols) {
   # } else {
   #    keep_cols <- colnames(groupByDF);
   # }
   if (include_num_reps) {
      keep_cols <- unique(c(keep_cols, "num_reps"));
   }
   return(groupByDF[, keep_cols, drop=FALSE]);
}
