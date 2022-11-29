
#' Mock-up shrink data.frame function
#'
#' Mock-up shrink data.frame function
#'
#' This function is a simplified version of `shrinkDataFrame()`
#' which is intended to use more modern methods from the R package
#' `data.table`.
#'
#' The general idea is to collapse `numeric` columns using `num_func`,
#' and collapse `character` and all other columns using `string_func`.
#'
#' Any exceptions, where a different function should be applied, are
#' passed via argument `extra_funcs` which is a `list` of functions
#' named by values in `colnames(df)`.
#'
#' @param df `data.frame` or compatible input class.
#' @param by `character` vector of one or more `colnames(df)`, used to define
#'    the row grouping.
#' @param string_func `function` used for `character` and other non-numeric
#'    column types.
#' @param num_func `function` used for `numeric` column types.
#' @param extra_funcs `list` of functions with list names that match values
#'    in `colnames(df)`, to be applied to specific columns in `df`. These
#'    functions will therefore override the default functions defined
#'    by `string_func` and `num_func`.
#' @param do_test `logical` indicating whether to perform an internal test
#'    with internally-generated argument values.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @family jamses utilities
#'
#' @export
shrink_df <- function
(df,
 by,
 string_func=jamba::cPasteU,
 num_func=mean,
 extra_funcs=NULL,
 do_test=FALSE,
 verbose=FALSE,
 ...)
{
   if (!suppressPackageStartupMessages(require(data.table))) {
      stop("The data.table package is required.");
   }
   if (do_test) {
      df <- data.frame(A=rep(LETTERS[1:3], c(1,2,3)),
         B=1:6,
         C=rep(LETTERS[4:6], c(3,2,1)));
      by <- "C";
   }
   by <- intersect(by, colnames(df));
   if (length(by) == 0) {
      stop("'by' not found colnames(df).");
   }
   use_names <- jamba::nameVector(setdiff(colnames(df), by));
   func_set <- lapply(use_names, function(i){
      if (is.numeric(df[[i]])) {
         num_func
      } else {
         string_func
      }
   });
   extra_names <- intersect(names(extra_funcs), use_names);
   if (length(extra_names) > 0) {
      func_set[extra_names] <- extra_funcs[extra_names];
   }

   compare_func_list <- function
   (l)
   {
      func_i <- rep(NA, length(func_set));
      names(func_i) <- names(func_set);
      i_seq <- seq_along(func_i);
      for (i in i_seq) {
         if (is.na(func_i[i])) {
            func_i[i] <- i;
            j_seq <- tail(i_seq, -i);
            for (j in j_seq) {
               k <- identical(func_set[[i]], func_set[[j]]);
               if (k) {
                  func_i[j] <- i;
               }
            }
         }
      }
      func_name_l <- split(names(func_set), func_i);
      func_l <- func_set[match(unique(func_i), func_i)];
      names(func_l) <- unique(func_i);
      return(list(names=func_name_l, fn=func_l));
   }

   func_sets <- compare_func_list(func_set);

   dt <- tryCatch({
      data.table::data.table(df, key=by);
   }, error=function(e){
      data.table::data.table(as.data.frame(df), key=by);
   })
   if (verbose) {
      jamba::printDebug("shrink_df(): ",
         "Running each data.table function set.");
   }
   dts <- lapply(seq_along(func_sets[[1]]), function(i){
      i_names <- func_sets[[1]][[i]];
      i_func <- func_sets[[2]][[i]];
      id1 <- which(names(dt) %in% i_names);
      dt1 <- dt[, lapply(.SD, i_func), by=by, .SDcols=id1];
   });
   if (verbose) {
      jamba::printDebug("shrink_df(): ",
         "Merging data.table function sets.");
   }
   dt2 <- do.call(`[`, dts);
   dt3 <- dt2[,colnames(dt), with=FALSE];
   if (verbose) {
      jamba::printDebug("shrink_df(): ",
         "Applying original object class.");
   }
   df3 <- as(dt3, head(class(df), 1));
   return(df3);
}
