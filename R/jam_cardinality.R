
#' Determine cardinality between two vectors
#' 
#' Determine cardinality between two vectors, defined in terms
#' of "one to one", "one to many", "many to one", or
#' "many to many".
#' 
#' This function is called by `heatmap_se()` to make a reasonable
#' guess for which columns from `colData()` to include in
#' `top_colnames` for a heatmap.
#' 
#' @returns `integer` vector with max `x`, max `y` repeated values.
#' @keywords internal
#' @noRd
cardinality <- function(
   x,
   y=NULL,
   verbose=FALSE,
   ...
) {
   #
   df_classes <- c("matrix",
      "data.frame",
      "DataFrame",
      "DFrame",
      "tbl");
   if (length(y) > 0 && inherits(y, df_classes)) {
      if (verbose) {
         jamba::printDebug("cardinality(): ",
            c("Converting y to vector with ", "pasteByRow()"),
            sep="");
      }
      y <- jamba::pasteByRow(y);
   }
   if (length(x) > 0) {
      if (inherits(x, df_classes)) {
         if (length(y) > 0) {
            if (verbose) {
               jamba::printDebug("cardinality(): ",
                  c("Converting x to vector with ", "pasteByRow()"),
                  sep="");
            }
            x <- jamba::pasteByRow(x);
            x <- data.frame(x=x,
               y=y);
         } else {
            if (ncol(x) != 2) {
               cli::cli_abort(
                  "When only x is supplied it must have two columns.");
            }
            x <- data.frame(
               x=x[, 1],
               y=x[, 2]
            );
         }
      } else {
         if (length(y) == 0) {
            stop("When x is a vector, y must be supplied.");
         }
         x <- data.frame(
            x=x,
            y=y
         );
      }
   }

   # x is a data.frame with two columns
   x_uniq <- unique(x);

   if (verbose > 1) {
      jamba::printDebug("cardinality(): ",
         "head(x, 20):")
      print(head(x, 20));
      jamba::printDebug("cardinality(): ",
         "head(unique(x), 20):")
      print(head(x_uniq, 20));
   }

   x_tc <- jamba::tcount(x_uniq[,2]);
   y_tc <- jamba::tcount(x_uniq[,1]);
   c(`from`=max(x_tc),
      `to`=max(y_tc))
}
