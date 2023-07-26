
#' Convert SEDesign to data.frame of design factors
#'
#' Convert SEDesign to data.frame of design factors
#'
#' This function is a utility function intended to convert `SEDesign`
#' objects, which contain a design matrix, into a `data.frame`
#' representing the experimental design factors used to produce
#' the design matrix.
#'
#' @family jam experiment design
#'
#' @returns `data.frame` with one column for each experiment design factor,
#'    where each column is converted to `factor` for non-numeric columns,
#'    and where colnames are defined as follows:
#'    * use `factor_colnames` when supplied
#'    * use the corresponding column in `colData(se)` when `se` is supplied
#'    * use `factor#` format for all other columns
#'
#' @param sedesign `SEDesign` object
#' @param se `SummarizedExperiment` object which contains
#'    `SummarizedExperiment::colData()` with columns of design annotations.
#'    When present, values in the factor `data.frame` will be matched
#'    to values in these columns in order to re-use factor levels
#'    where possible. When `factor_names` is not defined, this step
#'    will also re-use the associated colnames from `colData(se)`
#'    when the values are matched.
#' @param factor_names `character` vector of colnames to use for the resulting
#'    design factor colnames.
#' @param default_order `character` string indicating how to order factor
#'    levels when a colname in the factor `data.frame` is not a `factor`
#'    in the matching `colData()` column.
#'    * `"appearance"`: factor levels are defined in the order they
#'    appear in the sample data
#'    * `"mixedSort"`: factor levels are defined by sorting data using
#'    `jamba::mixedSort()` which performs proper alphanumeric sorting,
#'    such that "Treat2" appears before "Treat10" for example.
#' @param verbose `logical` indicating whether to print verbose output.
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
#' # convert to sedesign
#' sedesign <- groups_to_sedesign(idf)
#'
#' # prepare colData data.frame
#' cdf <- data.frame(check.names=FALSE, stringsAsFactors=FALSE,
#'    jamba::rbindList(strsplit(isamples_1, "_")))
#' colnames(cdf) <- c("Treatment", "Flag", "Genotype", "Rep")
#' rownames(cdf) <- samples(sedesign);
#' cdf
#' # prepare assay matrix
#' imatrix <- matrix(data=seq_len(nrow(cdf) * 10), ncol=nrow(cdf));
#' colnames(imatrix) <- rownames(cdf);
#' rownames(imatrix) <- paste0("row", 1:10);
#' # prepare se
#' se <- SummarizedExperiment::SummarizedExperiment(
#'    assays=list(raw=imatrix),
#'    colData=cdf)
#'
#' sedesign_to_factors(sedesign, se=se)
#' sedesign_to_factors(sedesign, se=se)[,1]
#'
#' SummarizedExperiment::colData(se)$Treatment <- factor(
#'    SummarizedExperiment::colData(se)$Treatment,
#'    levels=c("Etop", "DMSO"))
#' sedesign_to_factors(sedesign, se=se)[,1]
#'
#' SummarizedExperiment::colData(se)$Treatment <- factor(
#'    SummarizedExperiment::colData(se)$Treatment,
#'    levels=c("DMSO", "Etop"))
#' sedesign_to_factors(sedesign, se=se)[,1]
#'
#' sedesign_to_factors(sedesign, factor_names=c("Treatment", "Flag", "Genotyoe"), se=se)
#' sedesign_to_factors(sedesign, factor_names=c("Treatment", "Flag", "Genotyoe"), se=se)[,1]
#'
#' sedesign_to_factors(sedesign, factor_names=c("Treat", "Flag", "Genotyoe"), se=se)[,1]
#'
#' @export
sedesign_to_factors <- function
(sedesign,
 se=NULL,
 factor_names=NULL,
 default_order=c("appearance",
    "mixedSort"),
 verbose=FALSE,
 ...)
{
   # validate arguments
   default_order <- match.arg(default_order);

   # convert sedesign to vector of group names
   group_names <- groups(sedesign);

   # associated samples to groups
   im2list_internal <- function
   (x,
      empty=0,
      ...)
   {
      # the reciprocal of list2im()
      x_rows <- rownames(x);
      x_cols <- colnames(x);
      l <- lapply(jamba::nameVector(x_cols), function(i){
         i_empty <- as(empty, class(x[,i]));
         has_value <- (!x[,i] %in% i_empty);
         x_rows[has_value];
      });
      return(l);
   }
   groups_list <- im2list_internal(jamses::design(sedesign));
   isamples <- unlist(unname(groups_list));
   sample_groups <- jamba::nameVector(
      rep(names(groups_list), lengths(groups_list)),
      isamples)

   # split group names by "_" into data.frame of factor levels
   factors_df <- data.frame(check.names=FALSE,
      stringsAsFactors=FALSE,
      jamba::rbindList(strsplit(group_names, factor_sep)))
   rownames(factors_df) <- group_names;

   # assemble samples_df with one row per sample
   sample_match <- match(sample_groups,
      group_names)
   samples_df <- factors_df[sample_match, , drop=FALSE];
   rownames(samples_df) <- isamples;

   # assign colnames using factor_names
   if (length(factor_names) == 0) {
      # check colData(se)
      colnames(samples_df) <- paste0("factor", seq_len(ncol(factors_df)));
      icols <- jamba::nameVector(colnames(samples_df));
      if (length(se) > 0 && inherits(se, "SummarizedExperiment")) {
         colData_colnames <- colnames(SummarizedExperiment::colData(se));
         imatch <- match(isamples, colnames(se));
         # iterate to compare column values
         for (icol in icols) {
            if (verbose) {
               jamba::printDebug("Using values in '",
                  icol, "'",
                  indent=0)
            }
            ivals <- samples_df[[icol]];
            icol_name <- NULL;
            for (colData_colname in colData_colnames) {
               if (verbose) {
                  jamba::printDebug("Testing values in '",
                     colData_colname, "'",
                     indent=4)
               }
               colData_vals <- SummarizedExperiment::colData(
                  se[, imatch])[[colData_colname]];
               if (all(as.character(ivals) == as.character(colData_vals))) {
                  if (is.factor(colData_vals)) {
                     samples_df[[icol]] <- factor(ivals,
                        levels=levels(colData_vals))
                  } else if (is.numeric(colData_vals)) {
                     samples_df[[icol]] <- colData_vals;
                  } else {
                     if ("mixedSort" %in% default_order) {
                        samples_df[[icol]] <- factor(ivals,
                           levels=jamba::mixedSort(unique(colData_vals)))
                     } else {
                        samples_df[[icol]] <- factor(ivals,
                           levels=unique(colData_vals))
                     }
                  }
                  # rename colname
                  if (verbose) {
                     jamba::printDebug("Matched '",
                        indent=8,
                        icol, "' to '", colData_colname, "'.");
                  }
                  colData_colnames <- setdiff(colData_colnames,
                     colData_colname);
                  samples_df <- jamba::renameColumn(samples_df,
                     from=icol,
                     to=colData_colname)
                  break;
               }
            }
         }
      }
      factor_names <- colnames(samples_df);
   } else {
      n <- min(c(length(factor_names), ncol(samples_df)));
      colnames(samples_df)[seq_len(n)] <- head(factor_names, n);
      # assign given factor_names
      factor_names <- colnames(samples_df);
      print(head(samples_df))

      # use se if supplied
      if (length(se) > 0 && inherits(se, "SummarizedExperiment")) {
         colData_colnames <- colnames(SummarizedExperiment::colData(se));
         if (any(factor_names %in% colData_colnames)) {
            imatch <- match(isamples, colnames(se));
            for (icol in intersect(factor_names, colData_colnames)) {
               if (verbose) {
                  jamba::printDebug("Comparing values in '",
                     icol, "'",
                     indent=4)
               }
               ivals <- samples_df[[icol]];
               colData_vals <- SummarizedExperiment::colData(
                  se[, imatch])[[icol]];
               if (all(as.character(ivals) == as.character(colData_vals))) {
                  if (is.factor(colData_vals)) {
                     samples_df[[icol]] <- factor(ivals,
                        levels=levels(colData_vals))
                  } else if (is.numeric(colData_vals)) {
                     samples_df[[icol]] <- colData_vals;
                  } else {
                     if ("mixedSort" %in% default_order) {
                        samples_df[[icol]] <- factor(ivals,
                           levels=jamba::mixedSort(unique(colData_vals)))
                     } else {
                        samples_df[[icol]] <- factor(ivals,
                           levels=unique(colData_vals))
                     }
                  }
                  if (verbose) {
                     jamba::printDebug("Matched '",
                        indent=8,
                        " values in '", icol, "'.");
                  }
               }
            }
         }
      }
   }

   # one more pass to assign factor order as needed
   for (icol in factor_names) {
      if (!is.factor(samples_df[[icol]])) {
         if ("mixedSort" %in% default_order) {
            samples_df[[icol]] <- factor(samples_df[[icol]],
               levels=jamba::mixedSort(unique(samples_df[[icol]])));
         } else {
            samples_df[[icol]] <- factor(samples_df[[icol]],
               levels=unique(samples_df[[icol]]));
         }
      }
   }

   # colnames(samples_df) <- factor_names;
   return(samples_df)
}
