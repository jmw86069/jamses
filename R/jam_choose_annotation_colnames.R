
#' Choose interesting annotation colnames from a data.frame
#'
#' Choose interesting annotation colnames from a data.frame
#'
#' @param df `data.frame` with annotations that could be interesting
#'    to display at the top or side of a heatmap.
#' @param min_reps `numeric` minimum number of replicates required for
#'    a column to be considered interesting. For example, `min_reps=3`
#'    would require any value in a column to be repeated at least
#'    `3` times for that column to be interesting. This filter is
#'    intended to remove columns whose values are all unique, such
#'    as row identifiers.
#' @param min_values `numeric` minimum number of unique values
#'    required for a column to be considered interesting.
#' @param max_values `numeric` maximum number of unique values
#'    required for a column to be considered interesting. Too
#'    many values and the interest is lost. Also, too many values,
#'    and the color key becomes unbearable with too many labels.
#' @param keep_numeric `logical` indicating whether to keep columns
#'    with `numeric` values. When `keep_numeric == TRUE` it will
#'    override the rules above.
#' @param simplify `logical` indicating whether to filter out columns
#'    whose data already matches another column with 1:1 cardinality.
#'    This step requires `platjam::cardinality()` until that function
#'    is moved into the `jamba` package.
#' @param max_colnames `numeric` maximum number of colnames to return.
#'    Note that columns are not sorted for priority, so they will be
#'    returned in the order they appear in `df` after applying the
#'    relevant criteria.
#' @param ... additional arguments are ignored.
#'
#' @family jamses utilities
#'
#' @returns `character` vector of colnames in `df` that meet the criteria.
#'    If no colnames meet the criteria, this function returns `NULL`.
#'
#' @examples
#' df <- data.frame(
#'    threereps=paste0("threereps_", letters[c(1,1,1,3,5,7,7)]),
#'    time=paste0("time_", letters[c(1:7)]),
#'    tworeps=paste0("tworeps_", letters[c(12,12,14,14,15,15,16)]),
#'    num=sample(1:7),
#'    class=paste0("class_", LETTERS[c(1,1,1,3,5,7,7)]),
#'    blah=rep("blah", 7),
#'    maxvalues=c("one", "two", "three", "four", "five", "six", "six"))
#' df
#'
#' choose_annotation_colnames(df)
#' df[,choose_annotation_colnames(df)]
#'
#' choose_annotation_colnames(df, max_values=5)
#' df[,choose_annotation_colnames(df, max_values=5)]
#'
#' choose_annotation_colnames(df, simplify=FALSE)
#' df[,choose_annotation_colnames(df, simplify=FALSE)]
#'
#' choose_annotation_colnames(df, min_reps=3)
#'
#' choose_annotation_colnames(df, min_reps=1)
#'
#' choose_annotation_colnames(df, keep_numeric=TRUE)
#'
#' choose_annotation_colnames(df, min_reps=1)
#'
#' choose_annotation_colnames(df, min_reps=1, keep_numeric=TRUE)
#'
#' @export
choose_annotation_colnames <- function
(df,
 min_reps=2,
 min_values=2,
 max_values=Inf,
 keep_numeric=FALSE,
 simplify=TRUE,
 max_colnames=20,
   ...)
{
   #
   df_colnames <- jamba::nameVector(colnames(df));
   col_max_reps <- lapply(df_colnames, function(df_colname){
      unname(head(jamba::tcount(df[[df_colname]]), 1))
   })
   col_num_values <- lapply(df_colnames, function(df_colname){
      length(unique(df[[df_colname]]))
   })
   keep_flag <- (col_max_reps >= min_reps &
         col_num_values <= max_values &
         col_num_values >= min_values);
   if (keep_numeric) {
      col_classes <- jamba::sclass(x=df);
      keep_flag <- (keep_flag |
            col_classes %in% c("numeric", "integer"))
   }
   keep_colnames <- df_colnames[keep_flag];

   if (simplify %in% TRUE) {
      test_colnames <- keep_colnames;
      keep_colnames <- NULL;
      for (test_colname in test_colnames) {
         if (keep_numeric && col_classes[test_colname] %in% c("numeric", "integer")) {
            keep_colnames <- c(keep_colnames,
               test_colname);
            next;
         }
         if (length(keep_colnames) == 0) {
            keep_colnames <- test_colname;
         } else {
            min_cardinality <- min(sapply(keep_colnames, function(keep_colname){
               max(platjam::cardinality(
                  x=df[[test_colname]],
                  y=df[[keep_colname]]))
            }))
            if (min_cardinality > 1) {
               keep_colnames <- c(keep_colnames,
                  test_colname);
            }
         }
      }
   }
   if (length(max_colnames) && is.numeric(max_colnames)) {
      keep_colnames <- head(keep_colnames,
         max_colnames);
   }
   return(keep_colnames);
}
