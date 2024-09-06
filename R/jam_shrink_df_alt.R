
#' Shrink data.frame by row groups
#'
#' Shrink data.frame by row groups
#'
#' This function is currently a wrapper for `shrinkDataFrame()`,
#' it was formerly a simplified version of `shrinkDataFrame()`
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
#'    For efficiency, `string_func` by default is applied to the entire
#'    column, with `list` input, expecting `vector` output. It is not
#'    applied using `data.table`.
#' @param num_func `function` used for `numeric` column types. This function
#'    is applied using `data.table` and should expect a `vector` input,
#'    and provide a single atomic value output.
#' @param extra_funcs `list`, default `NULL`, containing `function` objects.
#'    The list names should match  `colnames(x)`, in order to apply a
#'    function to a specific column in `x`.
#'    These functions will therefore override the default functions defined
#'    by `string_func` and `num_func`.
#'    Only one function is applied per column.
#' @param do_test `logical`, default FALSE, indicating whether to perform an internal test
#'    with internally-generated argument values.
#' @param use_new_method `logical` default FALSE, whether to call newer
#'    tidy/data.table methods (TRUE), or call `shrinkDataFrame()` (FALSE).
#'    Currently `shrinkDataFrame()` is remarkably faster.
#'    More research necessary.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @family jamses utilities
#'
#' @export
shrink_df <- function
(x,
 by,
 string_func=function(x)jamba::cPasteSU(x, na.rm=TRUE),
 num_func=function(x)mean(x, na.rm=TRUE),
 add_string_cols=NULL,
 num_to_string_func=as.character,
 keep_na_groups=TRUE,
 include_num_reps=FALSE,
 extra_funcs=NULL,
 do_test=FALSE,
 use_new_method=FALSE,
 verbose=FALSE,
 ...)
{
   #
   if (do_test) {
      jamba::printDebug("shrink_df(): ",
         c("Running do_test=TRUE"))
      x <- data.frame(A=rep(LETTERS[1:3], c(1,2,3)),
         B=1:6,
         C=rep(LETTERS[4:6], c(3,2,1)));
      by <- "C";
   }
   if (nrow(x) == 0) {
      return(x);
   }
   if (FALSE %in% use_new_method) {
      return(shrinkDataFrame(x=x,
         groupBy=by,
         string_func=string_func,
         num_func=num_func,
         add_string_cols=add_string_cols,
         num_to_string_func=num_to_string_func,
         keep_na_groups=keep_na_groups,
         include_num_reps=include_num_reps,
         verbose=verbose,
         ...))
   }

   if (all(by %in% colnames(x))) {
      by <- intersect(by, colnames(x));
      use_names <- jamba::nameVector(setdiff(colnames(x), by));
   } else if (length(by) ==nrow(x)) {
      use_names <- colnames(x);
   }
   if (length(by) == 0) {
      stop("'by' not found colnames(df).");
   }
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
