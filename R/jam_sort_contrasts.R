
#' Sort contrasts by factor and level
#' 
#' Sort contrasts by factor and level, given contrast names, or
#' `SEDesign`, or `data.frame` of factors.
#' 
#' # Sort order
#' 
#' * Contrasts are sorted by depth: one-factor comparisons,
#' two-factor comparisons, etc.
#' * Contrasts are sorted for factor level contrasts in each
#' factor in `factor_order`. If `factor_order=1` then
#' contrasts appear first in the first factor column.
#' * Contrasts are then sorted by that factor column
#' using the observed contrast order.
#' * Contrasts are then sorted by other factor columns
#' as they appear in `factor_order`, using the
#' observed level order in each factor column.
#' 
#' # Sort order, described another way:
#' 
#' * All one-way contrasts will appear at the top.
#' 
#'    * The first one-way contrasts will be comparisons using
#'    the first value in `factor_order`, sorted by the
#'    comparison, then sorted by each remaining column
#'    in `factor_order`.
#'    * The next set of one-way contrasts will be the
#'    second value in `factor_order`, sorted by that
#'    comparison, then sorted by each remaining column
#'    in `factor_order`.
#' 
#' * All two-way contrasts will appear at the end,
#' since each two-way contrast involves two factors.
#' 
#'    * They will be sorted to involve contrasts in the
#'    same order as as `factor_order`.
#'    * The first two-way contrasts will involve comparisons
#'    involving the first two values in `factor_order`.
#' 
#' @family jam experiment design
#' 
#' @returns `data.frame` with factors in each column, and
#'    factor levels or factor level contrasts as column values.
#'    The `rownames()` contain contrast names.
#' 
#' @param x one of:
#'    * `character` vector of contrast names
#'    * `SEDesign`
#'    * `data.frame` with factors
#' @param factor_order `integer`, default NULL uses all factors,
#'    or specify the factor ordering to use for sorting.
#' @param ... additional arguments are passed to internal
#'    functions. For example 'factor_names' is passed to
#'    `contrasts_to_factors()`.
#' 
#' @examples
#' isamples_1 <- paste0(
#'    rep(c("DMSO", "Etop", "DMSO", "Etop"), each=6),
#'    "_",
#'    rep(c("NF", "Flag"), each=12),
#'    "_",
#'    rep(c("WT", "KO", "WT", "KO", "WT", "D955N", "WT", "D955N"), each=3),
#'    "_",
#'    LETTERS[1:3])
#' # simple data.frame with group information
#' idf <- data.frame(jamba::rbindList(strsplit(isamples_1, "_")))[,1:3]
#' rownames(idf) <- isamples_1;
#' colnames(idf) <- c("Treatment", "Flag", "Genotype")
#' # convert to sedesign
#' sedesign <- groups_to_sedesign(idf)
#' sort_contrasts(sedesign)
#' 
#' @export
sort_contrasts <- function(
   x,
   factor_order = NULL,
   ...
) {
   sedesign <- NULL;
   contrasts_df <- NULL;
   design_df <- NULL;
   if (inherits(x, "SEDesign")) {
      contrasts_df <- x@contrasts_df
      design_df <- x@design_df;
   } else {
      if (inherits(x, "data.frame")) {
         contrasts_df <- x
      } else if (is.atomic(x)) {
         contrasts_df <- contrasts_to_factors(x, ...)
      }
      design_df <- apply(contrasts_df, 2, simplify=FALSE,
         function(i){
            unique(jamba::unvigrep("-", i))
         });
   }
   n <- ncol(contrasts_df)
   if (length(factor_order) == 0) {
      factor_order <- seq_len(n)
   }
   
   cdf2 <- apply(contrasts_df, 2, function(i) {
      1 - grepl("-", i)
   })
   cdf3 <- n - rowSums(cdf2)
   cdf <- cbind(contrasts_df, cdf2, cdf3);
   for (i in seq_len(n)) {
      cdf[, i] <- factor(cdf[, i], 
         unique(c(
            design_df[[i]],
         cdf[, i])))
   }
   new_df <- jamba::rbindList(lapply(unname(split(cdf, cdf3)), function(cdf1) {
      unique(jamba::rbindList(lapply(factor_order, function(i1) {
         jamba::mixedSortDF(
            subset(cdf1, cdf1[, i1 + n] %in% 0),
            byCols = unique(c(i1, factor_order + n, factor_order))
         )
      })))
   }))
   new_df[, seq_len(n), drop=FALSE]
}
