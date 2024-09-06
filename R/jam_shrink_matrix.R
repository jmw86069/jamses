
#' Shrink a numeric matrix by groups of rows
#'
#' Shrink a numeric matrix across groups of rows, by applying a summary
#' function.
#'
#' This function is mainly a wrapper to use the amazingly fast
#' `data.table` package, with the ability to
#' provide a custom function to shrink row values.
#'
#' The default function uses `mean(x, na.rm=TRUE)` so that `NA`
#' values are ignored where possible.
#'
#' This function applies the same `shrink_func` to all columns, and
#' it optimal for `numeric` values. For more control over which
#' function to apply to specific columns, see `shrinkDataFrame()`.
#'
#' Trivia:
#' This function is identical to `splicejam::shrinkDataFrame()` except
#' that the default `shrink_func` includes `na.rm=TRUE` and no
#' longer calls the `.Internal()` function, since that is not
#' permitted by CRAN package guidelines.
#'
#' @param x `numeric` matrix
#' @param groupBy `character` or `factor` vector of group labels,
#'    whose length equals `nrow(x)`.
#'    These values will become rownames in the output data.
#' @param shrink_func `function` that takes vector input and returns
#'    **single value** output. The vector class can be checked, in order to
#'    call a function on numeric or character data separately, as
#'    needed.
#' @param return_class `character` string indicating the return data type.
#'    * `"data.frame"` returns a `data.frame` whose first column contains
#'    entries from `groups`.
#'    * `"matrix"` returns a numeric matrix whose rownames are entries
#'    from `groups`.
#' @param verbose logical indicating whether to print verbose output.
#'
#' @returns `data.frame` or `matrix` based upon argument `return_class`.
#'
#' @import data.table
#'
#' @family jamses utilities
#'
#' @export
shrink_matrix <- function
(x,
 groupBy,
 shrink_func=function(x){mean(x, na.rm=TRUE)},
 return_class=c("data.frame",
    "matrix"),
 verbose=FALSE,
 ...)
{
   ## Note: using .Internal(mean(x)) is (used to be) 5x faster than mean(x)
   ##
   ## Timings which show data.table and sqldf are fastest at apply() functions
   ## on groups:
   ## http://zvfak.blogspot.com/2011/03/applying-functions-on-groups-sqldf-plyr.html
   ##
   ## Another alternative is package sqldf, example syntax:
   ## n <- 100000;
   ## grp1 <- sample(1:750, n, replace=TRUE);
   ## grp2 <- sample(1:750, n, replace=TRUE);
   ## d <- data.frame(x=rnorm(n), y=rnorm(n), grp1=grp1, grp2=grp2, n, replace=TRUE);
   ## rsqldf <- system.time(sqldf("select grp1, grp2, avg(x), avg(y) from dgroup by grp1, grp2"));
   ##
   ## DT <- data.table(d);
   ## rdataT <- system.time(DT[,list(.Internal(mean(x)), .Internal(mean(y))), by=list(grp1,grp2)]);
   # if (!suppressPackageStartupMessages(require(data.table))) {
   #    stop("This method requires the data.table package.");
   # }
   return_class <- match.arg(return_class);

   ## Create DT object
   if (verbose) {
      t1 <- Sys.time();
   }
   DT <- data.table(
      data.frame(check.names=FALSE,
         stringsAsFactors=FALSE,
         x,
         groupBy=groupBy),
      key="groupBy");
   if (verbose) {
      t2 <- Sys.time();
   }

   ## Operate on the DT object
   byDT <- DT[, lapply(.SD, shrink_func),
      by="groupBy"];
   if (verbose) {
      t3 <- Sys.time();
   }

   if (verbose) {
      jamba::printDebug("shrink_matrix(): ",
         "Duration for data.table DT creation: ",
         format(t2-t1));
      jamba::printDebug("shrink_matrix(): ",
         "Duration for data.table shrink_matrix: ",
         format(t3-t2));
   }
   retData <- as(byDT, "data.frame");
   if (return_class %in% "matrix") {
      retData <- matrix(ncol=ncol(retData)-1,
         data=(as.matrix(retData[,-1,drop=FALSE])),
         dimnames=list(retData[,"groupBy"], colnames(retData)[-1]));
   }
   return(retData);
}
