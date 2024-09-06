
#' Shrink data.frame by row groups
#'
#' Shrink data.frame by row groups
#'
#' Purpose is to take a data.frame with multiple column classes,
#' and a groupBy key, and to collapse multiple rows per key
#' into one row per key.
#'
#' @family jamses utilities
#'
#' @param x `data.frame` (or equivalent)
#' @param groupBy `character` vector with one of the following:
#'    * one or more columns in `colnames(x)`. The values in these columns
#'    will define the row groups used.
#'    * `character` or `factor` with length equal to `nrow(x)`. These
#'    values will define the row groups used.
#' @param stringShrinkFunc `function`, default uses `jamba::cPasteSU()`,
#'    used for `character` or `factor` columns. Note that string
#'    columns are handled differently than `numeric` columns by applying
#'    vectorized operations across the complete set of rows in one step,
#'    rather than calling `data.table` on each subgroup.
#' @param numShrinkFunc `function`, default `function(x)mean(x, na.rm=TRUE)`,
#'    used for `numeric` columns. Note that this function is applied
#'    to each row group by `data.table`, and is typically very efficient
#'    for `numeric` values.
#' @param addStringCols `character` with optional `numeric` columns that
#'    should be handled as if they were `character` columns. Default `NULL`.
#' @param stringToNumFunc `function` used for `addStringCols` when converting
#'    `numeric` columns to `character`. Default `as.character()` retains
#'    the full `numeric` value, however it may be useful to use something
#'    like `function(x)signif(x, digits=3)` to limit the output to
#'    only three significant digits, or `function(x)format(x, digits=3)`.
#' @param keepNAgroupBy `logical`, default TRUE, whether to convert `NA`
#'    values in row groups to `""` so they are retained in the output.
#' @param includeNumReps `logical` indicating whether to add a column
#'    `"numReps"` to the output, with the `integer` number of rows
#'    in each row group.
#' @param collapseMethod `integer` default 2, indicating the internal
#'    collapse method used. Experimental.
#'    * `1` collapses each `numeric` column independently.
#'    * `2` collapses each set of `numeric` columns that use the same
#'    numeric shrink function. When all `numeric` columns use the same
#'    shrink function, they are all calculated in a single step, which
#'    is typically much faster.
#' @param applyKeepCols `logical` whether to reorder columns to match the
#'    input `colnames(x)`.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' testdf <- data.frame(check.names=FALSE,
#'    SYMBOL=rep(c("ACTB", "GAPDH", "PPIA"), c(2, 3, 1)),
#'    `logFC B-A`=c(1.4, 1.4, 2.3, NA, 2.3, 5.1),
#'    probe=paste0("probe", 1:6))
#' shrink_df(testdf, by="SYMBOL", num_func=function(x){mean(x, na.rm=TRUE)})
#'
#' testdftall <- do.call(rbind, lapply(1:10000, function(i){
#'    idf <- testdf;
#'    idf$SYMBOL <- paste0(idf$SYMBOL, "_", i);
#'    idf;
#' }))
#' head(testdftall, 12)
#' shrunk_tall <- shrink_df(testdftall,
#'    by="SYMBOL",
#'    num_func=function(x){mean(x, na.rm=TRUE)})
#'
#' shrunk_tall2 <- jamses::shrinkDataFrame(testdftall,
#'    groupBy="SYMBOL")
#' head(shrunk_tall2, 12)
#'
#' @export
shrinkDataFrame <- function
(x,
 groupBy,
 na.rm=TRUE,
 stringShrinkFunc=function(x)jamba::cPasteSU(x, na.rm=TRUE),
 numShrinkFunc=function(x){mean(x, na.rm=TRUE)},
 addStringCols=NULL,
 stringToNumFunc=as.character,
 keepNAgroupBy=TRUE,
 includeNumReps=FALSE,
 collapseMethod=2,
 applyKeepCols=TRUE,
 verbose=FALSE,
 ...)
{
   ## Purpose is to take a data.frame with multiple column classes,
   ## and a groupBy key, and to collapse multiple rows per key
   ## into one row per key.
   ##
   ## groupBy must be a factor or vector representing entries to be grouped
   ## together.  It can also be a character vector of colnames.
   ##
   ## numShrinkFunc is a function, or optionally a named list of functions,
   ## whose names match the numeric columns as detected.
   ## If a numeric colname is not present in the numShrinkFunc names,
   ## numShrinkFunc is recycled among the unmatched numeric colnames.
   ##
   ## stringShrinkFunc is a function, or optionally a named list of functions,
   ## whose names match the numeric columns as detected.
   ## If a string colname is not present in the stringShrinkFunc names,
   ## stringShrinkFunc is recycled among the unmatched string colnames.
   ##
   ## Note that numShrinkFunc is used in shrinkMatrix(..., shrinkFunc=numShrinkFunc[[1]])
   ##
   ## Note that stringShrinkFunc is used directly, stringShrinkFunc[[1]](...)
   ##
   ## addStringCols optionally specifies numeric columns which will be
   ## treated as text columns, e.g. the values will not be averaged, but
   ## instead will be concatenated together.
   ##
   ## only for addStringCols, where numeric is converted to character,
   ## this function is used for that coersion.  By default as.character(x)
   ## simply converts directly, but format(x, digits=2, trim=TRUE) can be used
   ## to minimize the number of significant digits.
   ##
   ## keepNAgroupBy=TRUE will convert NA to "" in the groupBy values,
   ## otherwise these entries will simply be dropped.
   ##
   ## includeNumReps=TRUE will include a column "numReps" which will contain
   ## an integer count of the rows grouped together
   ##
   ## collapseMethod==1 will run shrinkMatrix on each column
   ## collapseMethod==2 will run shrinkMatrix on sets of columns which
   ##    share the same shrink function

   ## iTopTableEx is a legacy variable to be migrated to use x in future...
   iTopTableEx <- x;
   if (!"data.frame" %in% class(iTopTableEx)) {
      if (verbose) {
         jamba::printDebug("shrinkDataFrame(): ",
            "coerced with data.frame(check.names=FALSE, x).");
      }
      iTopTableEx <- data.frame(check.names=FALSE,
         iTopTableEx);
   }

   #####################################
   ## Define groupBy
   groupByCols <- NULL;
   if (length(groupBy) == nrow(iTopTableEx)) {
      ## groupBy use as-is
      ## If factor, make sure NA is a level if there are NAs
      if (jamba::igrepHas("factor", class(groupBy))) {
         groupBy <- addNA(groupBy, ifany=TRUE);
      }
      if (!keepNAgroupBy && any(is.na(as.character(groupBy)))) {
         nonNaGroupBy <- !is.na(as.character(groupBy));
         groupBy <- droplevels(groupBy[nonNaGroupBy]);
         iTopTableEx <- iTopTableEx[nonNaGroupBy,,drop=FALSE];
      }
      groupByDF <- data.frame(groupBy=unique(groupBy),
         row.names=rmNA(naValue="NA", as.character(unique(groupBy))));
      if (verbose) {
         jamba::printDebug("shrinkDataFrame(): ",
            "head(groupBy):",
            head(groupBy));
      }
   } else {
      if (any(groupBy %in% colnames(iTopTableEx))) {
         groupByCols <- intersect(groupBy, colnames(iTopTableEx));
         groupBy <- jamba::pasteByRow(
            iTopTableEx[, groupByCols, drop=FALSE], sep="::");
         if (verbose) {
            jamba::printDebug("shrinkDataFrame(): ",
               "head(groupByCols):",
               head(groupByCols));
            jamba::printDebug("shrinkDataFrame(): ",
               "head(groupBy):",
               head(groupBy));
         }
         matchg <- match(unique(groupBy), groupBy);
         groupByDF <- data.frame(check.names=FALSE,
            stringsAsFactors=FALSE,
            iTopTableEx[matchg, groupByCols, drop=FALSE],
            row.names=unique(groupBy));
      } else {
         stop(paste0("groupBy must either be a vector equal to nrow(x),",
            " or a vector of colnames contained in x."));
      }
   }
   if (any(is.na(as.character(groupBy)))) {
      groupBy <- rmNA(groupBy, naValue="NA");
   }

   #########################################################
   ## iCols contains all columns which are not groupByCols,
   iCols <- setdiff(colnames(iTopTableEx), groupByCols);
   if (verbose) {
      jamba::printDebug("shrinkDataFrame(): ",
         "head(iCols):",
         head(iCols));
   }

   #########################################################
   ## Determine the column classes
   colClasses <- jamba::sclass(iTopTableEx);
   #colClasses <- sapply(jamba::nameVector(iCols), function(iCol){
   #   class(iTopTableEx[[iCol]]);
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
   if (!is.null(addStringCols) && any(addStringCols %in% numCols)) {
      moveNumCols <- intersect(addStringCols, numCols);
      numCols <- setdiff(numCols, moveNumCols);
      stringCols <- c(stringCols, moveNumCols);
      if (verbose) {
         jamba::printDebug("Moving columns: ",
            jamba::cPaste(moveNumCols),
            " from numCols to stringCols.")
      }
   }

   #########################################################
   ## Define the numShrinkFunc list
   if (length(numCols) > 0) {
      if (!is.list(numShrinkFunc)) {
         numShrinkFunc <- jamba::nameVector(rep(list(numShrinkFunc), length(numCols)),
            numCols);
      } else {
         if (!all(numCols %in% names(numShrinkFunc))) {
            ## Ensure we keep direct assignments to named columns,
            ## then reuse numShrinkFunc for subsequent columns
            ## whose names are not in names(numShrinkFunc).
            ## After all, we subset for the observed numCols to
            ## allow for numShrinkFunc to contain more entries than
            ## we have numCols.
            iNewCols <- setdiff(numCols, names(numShrinkFunc));
            numShrinkFunc <- c(numShrinkFunc,
               jamba::nameVector(
                  rep(numShrinkFunc, length.out=length(iNewCols)),
                  iNewCols));
         }
         numShrinkFunc <- numShrinkFunc[numCols];
      }
      if (verbose) {
         jamba::printDebug("shrinkDataFrame(): ",
            "numCols:",
            numCols,
            ", names(numShrinkFunc):",
            names(numShrinkFunc));
      }
   }

   #########################################################
   ## Define the stringShrinkFunc list
   if (length(stringCols) > 0) {
      if (!is.list(stringShrinkFunc)) {
         stringShrinkFunc <- jamba::nameVector(
            rep(list(stringShrinkFunc),
               length(stringCols)),
            stringCols);
      } else if (!all(stringCols %in% names(stringShrinkFunc))) {
         iNewCols <- setdiff(stringCols,
            names(stringShrinkFunc));
         stringShrinkFunc <- c(stringShrinkFunc,
            jamba::nameVector(
               rep(stringShrinkFunc,
                  length.out=length(iNewCols)),
               iNewCols));
      }
      if (verbose) {
         jamba::printDebug("shrinkDataFrame(): ",
            "stringCols:",
            stringCols,
            ", names(stringShrinkFunc):",
            names(stringShrinkFunc));
      }
   }

   #########################################################
   ## Count the number of entries grouped (optional)
   if (includeNumReps) {
      if (verbose) {
         jamba::printDebug("shrinkDataFrame(): ",
            "including numReps");
      }
      groupByT <- table(groupBy);
      iMatch <- match(unique(groupBy), names(groupByT));
      numRepsDF <- data.frame(numReps=as.vector(groupByT[iMatch]));
      groupByDF <- cbind(groupByDF,
         numRepsDF);
   }

   #####################################
   ## Shrink numeric columns
   if (length(numCols) > 0) {
      if (verbose) {
         jamba::printDebug("shrinkDataFrame(): ",
            "shrinking num columns:",
            numCols);
      }
      ## Collapse them all at once if all numShrinkFunc are identical
      ## TODO: split the numShrinkFunc by digest::digest() has value
      ## so sets of identical functions will be performed together.
      if (collapseMethod == 2) {
         #numShrinkFunc
         numShrinkFuncD <- sapply(numShrinkFunc,
            digest::digest,
            algo="md5");
         numShrinkFuncDL <- split(names(numShrinkFuncD), numShrinkFuncD);
         if (verbose) {
            jamba::printDebug("shrinkDataFrame(): ",
               "collapseMethod:", collapseMethod,
               " running on sets of identical numShrinkFunc functions.");
            jamba::printDebug("shrinkDataFrame(): ",
               "numShrinkFunc: ");
            print(numShrinkFunc);
            jamba::printDebug("shrinkDataFrame(): ",
               "numShrinkFuncD: ");
            print(numShrinkFuncD);
            jamba::printDebug("shrinkDataFrame(): ",
               "numShrinkFuncDL: ");
            print(numShrinkFuncDL);
         }
         numShrunk <- do.call(cbind, lapply(numShrinkFuncDL, function(iCol){
            if (verbose) {
               jamba::printDebug("   shrinkDataFrame(): ",
                  "iCol:",
                  iCol,
                  fgText=c("orange", "lightgreen"));
               # jamba::printDebug("   shrinkDataFrame(): ",
               #    "colnames(iTopTableEx): ",
               #    colnames(iTopTableEx));
               jamba::printDebug("   shrinkDataFrame(): ",
                  "head(data.frame):");
               print(head(
                  data.frame(iTopTableEx[, iCol, drop=FALSE],
                     check.names=FALSE,
                     stringsAsFactors=FALSE)
               ));
            }
            grOLi <- shrink_matrix(
               data.frame(iTopTableEx[, iCol, drop=FALSE],
                  check.names=FALSE,
                  stringsAsFactors=FALSE),
               groupBy=groupBy,
               shrinkFunc=numShrinkFunc[[head(iCol, 1)]],
               returnClass="matrix");
         }));
      } else {
         allNumIdentical <- all(sapply(numShrinkFunc[-1], function(x){
            identical(numShrinkFunc[[1]], x);
         }));
         if (allNumIdentical) {
            if (verbose) {
               jamba::printDebug("shrinkDataFrame(): ",
                  "collapseMethod:", collapseMethod,
                  " running only once for all numeric values.");
            }
            numShrunk <- shrink_matrix(
               as.data.frame(iTopTableEx[, numCols, drop=FALSE]),
               groupBy=rmNA(naValue="NA", groupBy),
               shrinkFunc=numShrinkFunc[[1]],
               returnClass="matrix");
         } else {
            if (verbose) {
               jamba::printDebug("shrinkDataFrame(): ",
                  "collapseMethod:", collapseMethod,
                  " running on each numeric column.");
            }
            numShrunk <- do.call(cbind,
               lapply(jamba::nameVector(numCols), function(iCol){
                  if (verbose) {
                     jamba::printDebug("  - ", iCol,
                        fgText=c("orange", "lightgreen"));
                     print(numShrinkFunc[[iCol]]);
                  }
                  grOLi <- shrink_matrix(
                     as.data.frame(iTopTableEx[, iCol, drop=FALSE]),
                     groupBy=groupBy,
                     shrinkFunc=numShrinkFunc[[iCol]], returnClass="matrix");
               }));
         }
      }
      numShrunkDF <- data.frame(check.names=FALSE,
         stringsAsFactors=FALSE,
         numShrunk);
      if (verbose) {
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
            if (verbose) {
               jamba::printDebug("   shrinkDataFrame(): ",
                  "stringShrinkFunc[['", iCol, "']]",
                  fgText=c("orange", "lightgreen"));
            }
            ## split vector into a list, then run the function on the list
            ## which is intended to generate a vector of values
            f1 <- iTopTableEx[[iCol]];
            iValsL <- stringShrinkFunc[[iCol]](split(as.character(f1),
               groupBy));
            if (jamba::igrepHas("factor", class(iTopTableEx[[iCol]]))) {
               if (verbose) {
                  jamba::printDebug("      shrinkDataFrame(): ",
                     "shrinking factor column values.");
               }
               ## Generate new factor levels in the original order, followed by
               ## sorted values, then keeping only those levels contained in
               ## the resulting data
               iLevels <- intersect(unique(c(levels(f1), mixedSort(iValsL))),
                  iValsL);
               if (!NA %in% levels(iTopTableEx[[iCol]])) {
                  ## Bonus points: if NA was not in the levels previously, keep it out
                  ## otherwise NA will be a formal factor level
                  iLevels <- rmNA(iLevels);
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
   ## keepCols is used to protect against columns whose types are
   ## not considered num or string, and may be lost above
   if (applyKeepCols) {
      keepCols <- intersect(colnames(iTopTableEx),
         colnames(groupByDF));
   } else {
      keepCols <- colnames(groupByDF);
   }
   if (includeNumReps) {
      keepCols <- unique(c(keepCols, "numReps"));
   }
   return(groupByDF[, keepCols, drop=FALSE]);
}
