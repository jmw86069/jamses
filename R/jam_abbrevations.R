
#' Find the shortest abbrevation to retain unique values
#'
#' Find the shortest abbrevation to retain unique values
#'
#' This function is intended to abbreviate factor levels used in
#' statistical contrasts to the smallest substring that uniquely
#' represents the unique entries provided in `x`.
#'
#' For example, `c("one", "two", "three", "four")` would be converted
#' to `c("on", "tw", "th", "fo")`.
#'
#' The default `retain_contig_numbers=TRUE` will attempt to retain
#' numeric values at the end of a string, to avoid splitting the number
#' at an intermediate position. This option only applies when the
#' character substring is not already unique before encountering
#' the numeric substring.
#'
#'     * For this input:
#'    `c("a", "p6", "p12", "p21")`
#'    the output keeps the contiguous numbers together:
#'    `c("a", "p6", "p12", "p21")`
#'
#'    * For this input:
#'    `c("a", "b6", "c12", "d21")`
#'    only the first character is retained,
#'    because it is already unique:
#'    `c("a", "b", "c", "d")`
#'
#' # Todo
#'
#' * Consider some method to retain contiguous numbers at the end
#' of a long string, while abbreviating the long string.
#'
#'    * For this input:
#'    `c("adult", "prenatal6", "prenatal12", "prenatal21")`
#'    the ideal output would be:
#'    `c("a", "p6", "p12", "p21")`
#'
#'    * To be fair, I do not know how to describe this logic.
#'    It may required breaking into words by character/non-character
#'    breakpoints, then applying substring to each?
#'
#' @family jamses utilities
#'
#' @returns `character` vector named using unique values in `x`, and
#'    whose values are the shortest abbreviated substrings which
#'    maintain consistent uniqueness.
#'
#' @param x `character` vector
#' @param retain_contig_numbers `logical`, default `TRUE`, whether numbers
#'    at the end of an abbreviated string should remain contiguous.
#'    * When `TRUE`, the goal is not to split a numeric value in the middle
#'    of the number.
#'    * When `FALSE` the string will be abbreviated at the first position
#'    of uniqueness.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' x <- c("a", "p6", "p12", "p21");
#' shortest_unique_abbreviation(x)
#'
#' shortest_unique_abbreviation(x, retain_contig_numbers=TRUE)
#'
#' x1 <- c("male", "female");
#' shortest_unique_abbreviation(x1)
#'
#' x2 <- c("Control", "Nicotine");
#' shortest_unique_abbreviation(x2)
#'
#' x3 <- c("Control", "Nicotine10", "Nicotine12", "Nicotine20");
#' shortest_unique_abbreviation(x3)
#'
#' x4 <- c("one", "two", "three", "four");
#' shortest_unique_abbreviation(x4)
#'
#' @export
shortest_unique_abbreviation <- function
(x,
 retain_contig_numbers=TRUE,
 verbose=FALSE,
 ...)
{
   #
   uniquex <- unique(x);
   max_nchar <- max(nchar(unique(x)), na.rm=TRUE)
   seq_nchar <- seq_len(max_nchar);
   for (use_nchar in seq_nchar) {
      if (use_nchar == max_nchar) {
         # we cannot abbreviate any further
         return(setNames(uniquex, uniquex))
      }
      use_uniquex <- substr(uniquex, 1, use_nchar);
      if (TRUE %in% retain_contig_numbers) {
         # experimental method to retain contiguous numbers
         seq_addnchar <- seq(from=use_nchar + 1, to=max_nchar);
         isnum_eligible <- (nchar(use_uniquex) %in%  use_nchar &
            grepl("[0-9]", gsub("^.*(.)$", "\\1", use_uniquex)))
         if (verbose) jamba::printDebug("use_nchar: ", use_nchar);# debug
         if (verbose) jamba::printDebug("use_uniquex: ", use_uniquex);# debug
         if (verbose) jamba::printDebug("isnum_eligible: ", isnum_eligible);# debug
         if (verbose) jamba::printDebug("seq_addnchar: ", seq_addnchar);# debug
         if (any(isnum_eligible)) {
            for (add_nchar in seq_addnchar) {
               new_uniquex <- substr(uniquex, add_nchar, add_nchar);
               new_isnum <- grepl("[0-9]", new_uniquex) & isnum_eligible;
               if (any(new_isnum)) {
                  # add only numeric digits
                  if (verbose) jamba::printDebug("Adding: ", new_uniquex);# debug
                  use_uniquex[new_isnum] <- paste0(use_uniquex[new_isnum],
                     new_uniquex[new_isnum]);
                  isnum_eligible <- new_isnum;
               }
            }
         }
         if (verbose) jamba::printDebug("\n");# debug
      }
      if (length(unique(use_uniquex)) == length(uniquex)) {
         break;
      }
   }
   return(setNames(use_uniquex, uniquex))
}
