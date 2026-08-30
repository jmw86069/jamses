
#' Internal incidence matrix to list conversion
#' 
#' Internal incidence matrix to list conversion, for internal use
#' 
#' @param x `matrix`
#' @param empty `numeric` value for empty values
#' @param ... additional arguments are ignored
#' 
#' @returns `list`
#' 
#' @keywords internal
#' @noRd
im2list_internal <- function
(
   x,
   empty = 0,
   ...
) {
   # the reciprocal of list2im()
   x_rows <- rownames(x)
   x_cols <- colnames(x)
   l <- lapply(jamba::nameVector(x_cols), function(i) {
      # i_empty <- as(empty, class(x[, i]))
      # has_value <- (!x[, i] %in% i_empty)
      has_value <- (!x[, i] %in% empty)
      x_rows[has_value]
   })
   return(l)
}


#' Convert list to incidence matrix
#'
#' Convert list to incidence matrix
#'
#' @param setlist `list` of vectors
#' @param empty default `0`, single value used for missing entries.
#' @param do_sparse `logical` indicating whether to convert output
#'    to sparse matrix, not currently implemented.
#' @param ... additional arguments are ignored.
#'
#' @family jamses utilities
#' 
#' @returns `numeric` matrix
#'
#' @export
list2im_opt <- function(
   setlist,
   empty=0,
   do_sparse=FALSE,
   ...
) {
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
   if (TRUE %in% do_sparse) {
      setlistim <- as(setlistim, "ngCMatrix");
   }
   return(setlistim);
}

#' List to value incidence matrix
#' 
#' List to value incidence matrix, intended as an internal function
#' 
#' @family jamses utilities
#' 
#' @inheritParams list2im_opt
#' @param force_sign `logical` default FALSE, using TRUE
#'    will convert `numeric` to `integer` sign using `sign()`.
#' 
#' @export
list2im_value_internal <- function(
   setlist,
   empty=NULL,
   do_sparse=FALSE,
   force_sign=FALSE,
   ...
) {
   setnames <- lapply(setlist, names);
   setnamesunion <- Reduce("union", setnames);
   
   # check for any character or factor input
   # so the resulting im will use consistent empty values
   setlist_hascharacter <- any(sapply(setlist, function(i){
      # is.character(i) | is.factor(i)
      inherits(i, c("character", "factor"))
   }))
   
   # define empty when not defined
   if (length(empty) == 0) {
      if (TRUE %in% setlist_hascharacter) {
         empty <- ""
      } else {
         empty <- 0
      }
   } else {
      empty <- head(empty, 1);
   }
   
   setlistim <- do.call(cbind, lapply(setlist, function(i){
      i_match <- match(names(i), setnamesunion);
      j <- rep(NA, length(setnamesunion));
      if (isTRUE(force_sign)) {
         if (!is.numeric(i)) {
            cli::cli_abort(
               "{.code force_sign=TRUE} but data is not {.cls numeric}.")
         }
         j[i_match] <- sign(i);
      } else {
         if (is.factor(i)) {
            cli::cli_warn(c("{.fun list2im_value}:",
               "coerced some input values from {.cls factor} to {.cls character}."))
            j[i_match] <- as.character(i);
         } else {
            j[i_match] <- i;
         }
      }
      if (anyNA(j)) {
         j[is.na(j)] <- empty;
      }
      j;
   }))
   rownames(setlistim) <- setnamesunion;
   return(setlistim);
}
