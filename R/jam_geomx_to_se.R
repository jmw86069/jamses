
#' Convert NanoStringGeoMxSet to SummarizedExperiment
#'
#' Convert NanoStringGeoMxSet to SummarizedExperiment, combining sData
#' and phenoData into colData
#'
#' The main purpose of this function is to convert
#' `GeomxTools::NanoStringGeoMxSet-class` to `SummarizedExperiment`
#' while also handling the oddity that the GeoMx data stores
#' some sample (column) annotation in slot `sData` and others in
#' slot `phenoData`. These two tables are combined into one
#' `colData` entry in the output object.
#'
#' Each data matrix in `assayData` from GeoMx is stored into the
#' corresponding `assays` entry in the output object.
#'
#' @returns `SummarizedExperiment` with `rowData`, `colData`, and `assays`.
#'
#' @family jamses SE utilities
#'
#' @param x `NanoStringGeoMxSet` object
#' @param assay_names `character` vector with one or more assay names.
#'    The default NULL will use all available assay names recognized
#'    by `se_to_assay_names()`.
#' @param ... additional arguments are ignored.
#'
#' @export
geomx_to_se <- function
(x,
 assay_names=NULL,
 verbose=FALSE,
 ...)
{
   #
   x_assay_names <- jamses::se_to_assay_names(x);
   if (verbose) {
      jamba::printDebug("geomx_to_se(): ",
         "Available assay names:",
         x_assay_names)
   }
   if (length(assay_names) == 0) {
      use_assay_names <- x_assay_names;
   } else {
      use_assay_names <- intersect(assay_names, x_assay_names);
   }
   if (length(use_assay_names) == 0) {
      cli::cli_abort(paste0(
         "{.var assay_names} ({assay_names}) did not match the input ",
         "{.var NanoStringGeoMxSet} object values ({x_assay_names})."));
      stop("assay_names did not match the NanoStringGeoMxSet object.")
   }
   rowcoldata <- jamses::se_to_rowcoldata(x)

   use_se <- SummarizedExperiment::SummarizedExperiment(
      assays=lapply(jamba::nameVector(use_assay_names), function(i){
         jamses::se_to_assay_data(x, assay_name=i)
      }),
      rowData=rowcoldata$rowData_se,
      colData=rowcoldata$colData_se)
   use_se
}
