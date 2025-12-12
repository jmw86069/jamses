
#' Get SE assay names
#'
#' Get SE assay names consistently for `SummarizedExperiment`,
#' `SingleCellExperiment`, `Seurat` data, `NanoStringGeoMxSet`,
#' and generic Biobase `eSet` compatible objects that provide
#' `assayData()`.
#'
#' @family jamses SE utilities
#'
#' @returns `character` vector with corresponding assay names
#'
#' @param se `SummarizedExperiment` or other recognized data type
#'    inheriting from either `SummarizedExperiment`, `ExpressionSet`,
#'    `eSet`, or `Seurat`.
#' @param ... additional arguments are ignored.
#'
#' @export
se_to_assay_names <- function
(se,
 ...)
{
   # Biobase
   biobase_types <- c("ExpressionSet",
      "eSet",
      "NanoStringGeoMxSet")
   # SummarizedExperiment
   se_types <- c("SummarizedExperiment")
   # Seurat
   seurat_types <- c("Seurat")

   if (inherits(se, biobase_types)) {
      #
      names(Biobase::assayData(se))
   } else if (inherits(se, se_types)) {
      SummarizedExperiment::assayNames(se)
   } else if (inherits(se, seurat_types)) {
      if (!requireNamespace("Seurat", quietly=TRUE)) {
         cli::cli_abort("Seurat must be installed to accept Seurat objects.")
         stop("Seurat must be installed to accept Seurat objects.")
      }
      # currently heatmap_se() converts to SingleCellExperiment
      names(SummarizedExperiment::assays(
         Seurat::as.SingleCellExperiment(se, ...)))
   } else {
      # try Biobase
      names(Biobase::assayData(se))
   }
}

#' Get SE assay data
#'
#' Get SE assay data consistently for `SummarizedExperiment`,
#' `SingleCellExperiment`, `Seurat` data, `NanoStringGeoMxSet`,
#' and generic Biobase `eSet` compatible objects that provide
#' `assayData()`.
#'
#' @family jamses SE utilities
#'
#' @returns `character` vector with corresponding assay names
#'
#' @param se `SummarizedExperiment` or other recognized data type
#'    inheriting from either `SummarizedExperiment`, `ExpressionSet`,
#'    `eSet`, or `Seurat`.
#' @param ... additional arguments are ignored.
#'
#' @export
se_to_assay_data <- function
(se,
 assay_name=NULL,
 ...)
{
   # Biobase
   biobase_types <- c("ExpressionSet",
      "eSet",
      "NanoStringGeoMxSet")
   # SummarizedExperiment
   se_types <- c("SummarizedExperiment")
   # Seurat
   seurat_types <- c("Seurat")

   # confirm assay_name
   if (length(assay_name) == 0) {
      se_assay_names <- se_to_assay_names(se);
      if (length(se_assay_names) == 1) {
         assay_name <- head(se_assay_names, 1);
      } else {
         cli::cli_abort(paste0(
            "{.var assay_name} must be supplied when there are ",
            "multiple available assays in {.var se}"));
         stop("assay_name must be supplied when there are multiple assays(se)")
      }
   }

   # get data
   if (inherits(se, biobase_types)) {
      #
      Biobase::assayData(se)[[assay_name]]
   } else if (inherits(se, se_types)) {
      SummarizedExperiment::assays(se)[[assay_name]]
   } else if (inherits(se, seurat_types)) {
      if (!requireNamespace("Seurat", quietly=TRUE)) {
         cli::cli_abort("Seurat must be installed to accept Seurat objects.")
         stop("Seurat must be installed to accept Seurat objects.")
      }
      # currently heatmap_se() converts to SingleCellExperiment
      SummarizedExperiment::assays(
         Seurat::as.SingleCellExperiment(se, ...))[[assay_name]]
   } else {
      # try Biobase
      Biobase::assayData(se)
   }
}
