
# jam_sestats.R
#
# placeholder for proper object type sestats

#' Convert sestats to table summary
#'
#' Convert sestats to table summary
#'
#' Note: This function is intended to provide a simple data.frame
#' summary with the number of hits for each contrast, signal, cutoff.
#' It is still being tested, and updated for usability.
#'
#' TODO: The order of `dimnames(hit_array)` should be user-customizable.
#' The series of dimnames in each lapply should use this order.
#'
#' @param sestats `list` output from `se_contrast_stats()`
#' @param style `character` string indicating what values to use:
#'    * `"text"`: number of hits (number up, number down)
#'    * `"integer"`: only the integer number of hits
#' @param rename_contrasts `logical` indicating whether to rename
#'    contrasts using `contrast2comp()`. The main benefit is reduction
#'    in length of the string that describes each contrast.
#' @param ... additional arguments are passed to `contrast2comp()`, which
#'    is relevant when `rename_contrasts=TRUE`, relevant arguments include:
#'    * `contrast_delim="-"` to
#'    * `contrast_factor_delim="_"` to customize the delimiter between factors
#'    in the input group names, for example `"Treatment1_Time1"` uses "_".
#'    * `comp_factor_delim=":"` to customize the delimiter between factors
#'    * `factor_order=NULL` to customize the order of factor comparisons
#'
#' @examples
#' if (FALSE) {
#' hitdf <- sestats_to_df(list(hit_array=hit_array));
#' hitdf_rowindex <- table(hitdf[[1]])[unique(hitdf[[1]])]
#' jamba::kable_coloring(
#'    hitdf[, -1, drop=FALSE],
#'    row.names=FALSE) %>%
#'    kableExtra::pack_rows(index=hitdf_rowindex)
#' }
#'
#' @export
sestats_to_df <- function
(sestats,
 style=c("text", "integer"),
 dimname_order=c(3, 2, 1),
 rename_contrasts=FALSE,
 ...)
{
   #
   style <- match.arg(style);

   #
   if ("list" %in% class(sestats) && "hit_array" %in% names(sestats)) {
      hit_array <- sestats$hit_array;
   } else {
      hit_array <- sestats;
   }

   # validate dimname_order
   if (!all(seq_along(dimnames(hit_array)) %in% dimname_order) &&
         !all(names(dimnames(hit_array)) %in% dimname_order)) {
      stop(paste0(
         "dimname_order must contain all integer index values for dimnames(hit_array),",
         " or must contain all names(dimnames(hit_array))"));
   }

   # optionally rename contrasts upfront
   if (rename_contrasts) {
      ha_dimnames <- dimnames(hit_array)
      contrast_dimname <- head(jamba::vigrep("contrast",
         names(ha_dimnames)), 1);
      if (length(contrast_dimname) == 1) {
         ha_dimnames[[contrast_dimname]] <- tryCatch({
            contrast2comp(ha_dimnames[[contrast_dimname]],
               ...);
         }, error=function(e){
            ha_dimnames[[contrast_dimname]]
         });
      }
      dimnames(hit_array) <- ha_dimnames;
   }

   # internal helper function
   # converts numeric vector of c(-1, 0, 1, NA) named by measurement
   # into character label


   # redefine hit_array so dimensions are in order, simplifying the process
   hit_array_new <- aperm(hit_array, dimname_order)
   idf <- jamba::rbindList(lapply(dimnames(hit_array_new)[[1]], function(dim1value){
      im <- format_hits(hit_array_new[dim1value, ,])
      if ("matrix" %in% class(im)) {
         if (all(colnames(im) %in% dimnames(hit_array_new)[[3]])) {
            im <- do.call(cbind, lapply(jamba::nameVector(colnames(im)), function(j){
               unlist(im[,j])
            }))
         }
      } else if ("list" %in% class(im)) {
         if (all(names(im) %in% dimnames(hit_array_new)[[2]])) {
            im <- jamba::rbindList(im);
         } else if (all(names(im) %in% dimnames(hit_array_new)[[3]])) {
            im <- do.call(cbind, im);
         }
         rownames(im) <- dimnames(hit_array_new)[[2]];
         colnames(im) <- dimnames(hit_array_new)[[3]];
      }
      idf1 <- data.frame(
         check.names=FALSE,
         stringsAsFactors=FALSE,
         dim1=dim1value,
         dim2=rownames(im),
         im)
      colnames(idf1)[1:2] <- names(dimnames(hit_array_new))[c(1, 2)];
      idf1
   }))

   if (rename_contrasts) {
      if ("Contrasts" %in% colnames(idf)) {
         idf$Contrasts <- tryCatch({
            contrast2comp(idf$Contrasts);
         }, error=function(e){
            idf$Contrasts
         });
      } else {
         contrast_colnums <- which(!colnames(idf) %in% names(dimnames(hit_array_new)))
      }
   }

   idf
}

#' Format list of hit vectors into summary counts
#'
#' Format list of hit vectors into summary counts
#'
#' This function is used by `sestats_to_df()`,
#' to summarize hits for each contrast into one of these formats:
#' * string summary of statistical hits with `"hits,up,down"`
#' when `style="text"`, for example: `"623 hits (267 up, 379 down)"`
#' * integer count of statistical hits when `style="integer"`,
#' for example: `623`.
#' * `vector` of integer counts for `c("hits", "up", "down")`,
#' for example: `c(hits=623, up=267, down=379)`.
#'
#' The function may be useful outside of `sestats_to_df()` so it
#' is exported as a convenience function.
#'
#' @param hits `list`, `array`, or `vector` with values `c(-1, 1)`
#'    indicating hits down, or hits up, respectively. When any vector
#'    contains only `NA` values, then `NA` is also returned. This
#'    distinction is as follows:
#'    * `NA` values indicate the statistical contrast was not performed,
#'    which happens when `se_contrast_stats()` arguments for interaction
#'    contrasts differ from pairwise contrasts, for example `int_adjp_cutoff`,
#'    `int_p_cutoff`, or `int_fold_cutoff`.
#'    * `length==0` indicates the statistical contrast was performed, and
#'    there were no statistical hits for the given cutoffs.
#' @param style `character` string indicating the output format:
#'    * `"text"`: `character` string with `"hits(hits up, hits down)"`.
#'    * `"integer"`: `integer` number of hits, or `NA` when the test was
#'    not performed.
#'    * `"vector"`: `integer` vector with names `("hit", "up", "down")`.
#' @param ... additional arguments are ignored.
#'
#' @return the same data type as input, where the hit vector is replaced
#'    with a single value summarizing the hits. The data types have
#'    three expected options:
#'    * `array`: when used with `hit_array` data from `se_contrast_stats()`
#'    * `list`: when used on a particular subset of `hit_array` data
#'    * `numeric`: when used with a single contrast
#'
#' @examples
#' set.seed(123)
#' hitlist <- list(
#'    `groupA-groupB`=sample(c(-1, 1), size=25, replace=TRUE),
#'    `groupA-groupC`=sample(c(-1, 1), size=50, replace=TRUE))
#' format_hits(hitlist, style="text")
#' format_hits(hitlist, style="integer")
#'
#' @export
format_hits <- function
(hits,
 style=c("text", "integer", "vector"),
 ...) {
   # convert hit vector into integer values by hits,up,down
   # using NA values when all value are NA, which indicates the test
   # was not performed, distinct from length==0 which indicates the
   # test was performed but produced no statistical hits.
   style <- match.arg(style)

   ilist <- hits;
   if (!is.list(hits)) {
      ilist <- list(hits);
   }
   ilist_lengths <- jamba::rbindList(lapply(ilist, function(i){
      if (length(i) > 0 && all(is.na(i))) {
         c(hits=NA,
            up=NA,
            down=NA)
      } else {
         c(hits=length(i),
            up=sum(i > 0),
            down=sum(i < 0))
      }
   }))
   # create the appropriate style
   if ("integer" %in% style) {
      iout <- ilist_lengths[,1];
   } else if ("vector" %in% style) {
      iout <- lapply(seq_len(nrow(ilist_lengths)), function(j){
         i <- unname(unlist(ilist_lengths[j,]))
         if (is.na(i[1])) {
            c(hits=NA, up=NA, down=NA)
         } else {
            c(hits=i[1], up=i[2], down=i[3])
         }
      })
      #iout <- ifelse(is.na(ilist_lengths[,1]),
   } else {
      iout <- ifelse(is.na(ilist_lengths[,1]),
         "",
         paste0(
            jamba::formatInt(ilist_lengths[,"hits"]),
            ifelse(ilist_lengths[,1] == 0,
               " hits",
               paste0(
                  " hits (",
                  jamba::formatInt(ilist_lengths[,"up"]),
                  " up, ",
                  jamba::formatInt(ilist_lengths[,"down"]),
                  " down)")
            )
         )
      )
   }
   if (is.list(hits)) {
      hits[] <- iout;
   } else {
      hits <- unlist(iout);
   }
   hits;
}
