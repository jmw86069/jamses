
#' Get SE colData and rowData
#'
#' Get SE colData and rowData consistently for SummarizedExperiment,
#' SingleCellExperiment, Seurat data, and generic Biobase `eSet`
#' compatible objects that provide `featureData()` and `phenoData()`.
#'
#' This function provides a straightforward way to return
#' the equivalent of `data.frame(check.names=FALSE, rowData(se))`
#' and `data.frame(check.names=FALSE, colData(se))` for several
#' types of object types.
#' * It also defines `rownames()` and `colnames()` if either are missing.
#' * When `rowData` has no annotation columns, it defines one column `"rows"`
#' using `rownames(se)`.
#' * If slot name `"rowRanges"` exists, and `"rowData"` either does not
#' exist or has zero columns, it will use `rowRanges()`.
#' * When `colData` has no annotation columns, it defines one column
#' `"columns"` using `colnames(se)`.
#' * When input `se` class is `"Seurat"`, it converts the object with
#' `Seurat::as.SingleCellExperiment()`
#' * Any other object uses `Biobase::featureData()` for `"rowData_se"`,
#' and `Biobase::phenoData()` for `"colData_se"`.
#'
#' To verify the logic used at each step, set `verbose=TRUE`.
#'
#' For Class `"NanoStringGeoMxSet"` defined in `GeomxTools`, the
#' `colData_se` is defined using `GeomxTools::sData()` in order to
#' return the combined `data.frame` with `protocolData()` and `pData()`
#' together. Accordingly, it is possible to have duplicated colnames,
#' which becomes a problem for many downstream tools, so the second
#' instance of any duplicated colnames have `"_v1"` added by using
#' `jamba::makeNames(x, renameFirst=FALSE)`. The first instance of
#' each duplicated colname is not renamed. Use `verbose=TRUE` to confirm
#' when duplicated columns are detected and renamed.
#'
#' @family jamses utilities
#'
#' @returns `list` with two components:
#' * `"colData_se"` as a `data.frame` with column metadata from `se`
#' * `"rowData_se"` as a `data.frame` with row metadata from `se`
#'
#' @param se recognized object:
#' * `SummarizedExperiment` or any object that inherits from this class
#' * `SingleCellExperiment` because it inherits from `SummarizedExperiment`
#' * `Seurat` which is converted to `SingleCellExperiment`
#' * `ExpressionSet` and any object that inherits from this class, using
#' `Biobase::featureData()` and `Biobase::phenoData()`.
#' * `NanoStringGeoMxSet` which uses `GeomxTools::sData()` or
#' `NanoStringNCTools::sData()` to define `colData_se`, otherwise is
#' handled equivalent to `ExpressionSet`.
#' @param ... additional arguments are ignored
#'
#' @export
se_to_rowcoldata <- function
(se,
 verbose=FALSE,
 ...)
{
   # Experimental convenience function to generate colData_df and rowData_df
   #
   if (inherits(se, "Seurat")) {
      if (verbose) {
         jamba::printDebug("se_to_colrowdata(): ",
            "Converting Seurat to SingleCellExperiment");
      }
      se <- Seurat::as.SingleCellExperiment(se,
         # assay=assay_name,
         ...);
   }

   # Validate rownames/colnames
   if (length(rownames(se)) == 0) {
      if (verbose) {
         jamba::printDebug("se_to_colrowdata(): ",
            "Defining rownames which were empty, using 'row00#' format.");
      }
      rownames(se) <- paste0("row",
         jamba::padInteger(seq_len(nrow(se))));
   }
   if (length(colnames(se)) == 0) {
      if (verbose) {
         jamba::printDebug("se_to_colrowdata(): ",
            "Defining colnames which were empty, using 'column00#' format.");
      }
      colnames(se) <- paste0("column",
         jamba::padInteger(seq_len(ncol(se))));
   }
   # Todo: Add GeoMx object type

   # Generate colData_df, rowData_df
   if (inherits(se, "SummarizedExperiment") ||
         inherits(se, "SingleCellExperiment")) {
      if (verbose) {
         jamba::printDebug("se_to_colrowdata(): ",
            "Handling SingleCellExperiment");
      }
      # colData
      if (length(SummarizedExperiment::colData(se)) > 0) {
         if (verbose) {
            jamba::printDebug("se_to_colrowdata(): ",
               "Used colData() for colData_se");
         }
         colData_se <- data.frame(check.names=FALSE,
            SummarizedExperiment::colData(se))
      } else {
         if (verbose) {
            jamba::printDebug("se_to_colrowdata(): ",
               "Defined new colData_se.");
         }
         colData_se <- data.frame(check.names=FALSE,
            row.names=colnames(se),
            columns=colnames(se))
      }
      if (ncol(colData_se) == 0) {
         if (verbose) {
            jamba::printDebug("se_to_colrowdata(): ",
               "Added 'columns' to zero-column colData_se using colnames()");
         }
         colData_se$columns <- rownames(colData_se);
      }
      # rowData
      # 0.0.69.900 - check using rowData() accessor
      if ((length(SummarizedExperiment::rowData(se)) > 0 &&
            ncol(SummarizedExperiment::rowData(se)) > 0) ||
            !"rowRanges" %in% slotNames(se)) {
         if (verbose) {
            jamba::printDebug("se_to_colrowdata(): ",
               "Used rowData() for rowData_se");
         }
         rowData_se <- data.frame(check.names=FALSE,
            SummarizedExperiment::rowData(se));
      } else if (length(SummarizedExperiment::rowRanges(se)) > 0) {
         # 0.0.69.900 - check using rowRanges() accessor
         if (verbose) {
            jamba::printDebug("se_to_colrowdata(): ",
               "Used rowRanges() for rowData_se");
         }
         rowData_se <- data.frame(check.names=FALSE,
            row.names=rownames(se),
            SummarizedExperiment::values(
               SummarizedExperiment::rowRanges(se)))
      } else {
         if (verbose) {
            jamba::printDebug("se_to_colrowdata(): ",
               "Defined new rowData_se.");
         }
         rowData_se <- data.frame(check.names=FALSE,
            row.names=rownames(se),
            rows=rownames(se))
      }
      if (ncol(rowData_se) == 0) {
         if (verbose) {
            jamba::printDebug("se_to_colrowdata(): ",
               "Added 'rows' to zero-column rowData_se using rownames()");
         }
         rowData_se$rows <- rownames(rowData_se);
      }
   } else {
      ## Todo
      # - handle non-standard column classes in AnnotatedDataFrame

      # rowData using tryCatch()
      rowData_se <- tryCatch({
         rowData_se1 <- as(Biobase::featureData(se), "data.frame");
         rownames(rowData_se1) <- rownames(se);
         if (verbose) {
            jamba::printDebug("se_to_colrowdata(): ",
               "Defined rowData_se using featureData().");
         }
         rowData_se1;
      }, error=function(e){
         if (verbose) {
            jamba::printDebug("se_to_colrowdata(): ",
               "Defined new rowData_se using rownames.");
         }
         data.frame(check.names=FALSE,
            row.names=rownames(se),
            rows=rownames(se))
      });
      if (ncol(rowData_se) == 0) {
         if (verbose) {
            jamba::printDebug("se_to_colrowdata(): ",
               "Added 'rows' to zero-column rowData_se using rownames()");
         }
         rowData_se$rows <- rownames(rowData_se);
      }

      # colData using tryCatch()
      colData_se <- NULL;
      if (inherits(se, "NanoStringGeoMxSet")) {
         # for NanoStringGeoMxSet data, try to use GeomxTools::sData()
         # or NanoStringNCTools::sData()
         if (jamba::check_pkg_installed("GeomxTools")) {
            colData_se <- GeomxTools::sData(se)
            rownames(colData_se) <- colnames(se);
            if (verbose) {
               jamba::printDebug("se_to_colrowdata(): ",
                  "Defined colData_se using GeomxTools::sData().");
            }
         } else if (jamba::check_pkg_installed("NanoStringNCTools")) {
            colData_se <- NanoStringNCTools::sData(se)
            rownames(colData_se) <- colnames(se);
            if (verbose) {
               jamba::printDebug("se_to_colrowdata(): ",
                  "Defined colData_se using NanoStringNCTools::sData().");
            }
         }
         ## data.frame is returned, consider checking then coercing in future

         # check for duplicated colnames(colData_se)
         # - duplicated colnames arise when sData()
         #   combines pData() and protocolData() and they may contain
         #   the same colnames, unintentionally
         # - the first instance retains the colname, the second adds "_v1"
         if (length(colData_se) > 0) {
            if (any(duplicated(colnames(colData_se)))) {
               if (verbose) {
                  jamba::printDebug("se_to_colrowdata(): ",
                     "Renaming duplicated colnames(colData_se) using ",
                     "jamba::makeNames(x, renameFirst=FALSE).");
               }
               colnames(colData_se) <- jamba::makeNames(colnames(colData_se),
                  renameFirst=FALSE,
                  ...)
            }
         }
      }
      if (length(colData_se) == 0) {
         colData_se <- tryCatch({
            colData_se1 <- as(Biobase::phenoData(se), "data.frame");
            rownames(colData_se1) <- colnames(se);
            if (verbose) {
               jamba::printDebug("se_to_colrowdata(): ",
                  "Defined colData_se using phenoData().");
            }
            colData_se1;
         }, error=function(e){
            if (verbose) {
               jamba::printDebug("se_to_colrowdata(): ",
                  "Defined new rowData_se using colnames.");
            }
            data.frame(check.names=FALSE,
               row.names=colnames(se),
               columns=colnames(se))
         })
      }
      # check for empty columns
      if (ncol(colData_se) == 0) {
         if (verbose) {
            jamba::printDebug("se_to_colrowdata(): ",
               "Added 'columns' to zero-column colData_se using colnames()");
         }
         colData_se$columns <- rownames(colData_se);
      }
   }
   return(list(
      rowData_se=rowData_se,
      colData_se=colData_se));
}
