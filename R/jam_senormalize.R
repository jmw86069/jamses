
#' Normalize SummarizedExperiment data
#'
#' Normalize SummarizedExperiment data
#'
#' This function applies one or more data normalization methods
#' to an input `SummarizedExperiment` object. The normalization is
#' applied to one or more matrix data stored in `assays(se)`,
#' each one is run independently.
#'
#' Note that supplying `genes` and `samples` will apply normalization
#' to only those `genes` and `samples`, and this data will be
#' stored in the full `SummarizedExperiment` object `se` with
#' `NA` values used to fill any values not present in `genes`
#' or `samples`.
#'
#' For example if `assay_names` contains two assay names,
#' and `method` contains two methods, the output will include
#' four normalizations, where each assay name is normalized two ways.
#' The output assay names will be something like `"assay1_method1"`,
#' `"assay1_method2"`, `"assay2_method1"`, `"assay2_method2"`.
#' It is not always necessary to normalize data by multiple different
#' methods, however when two methods are similar and need to be
#' compared, the `SummarizedExperiment` object is a convenient
#' place to store different normalization results for downstream
#' comparison. Further, the method `se_contrast_stats()` is able
#' to apply equivalent statistical contrasts to each normalization,
#' and returns an array of statistical hits which is convenient
#' for direct comparison of results.
#'
#' This method calls `matrix_normalize()` to perform each normalization
#' step, see that function description for details on each method.
#'
#' @family jamses SE utilities
#'
#' @return `SummarizedExperiment` object where the normalized output
#'    is added to `assays(se)` using the naming format `method_assayname`.
#'
#' @param se `SummarizedExperiment` object
#' @param method `character` vector indicating which normalization method(s)
#'    to apply.
#'    * `"quantile"`: quantile normalization via `limma::normalizeQuantiles()`
#'    * `"jammanorm"`: log-ratio normalization via `jamma::jammanorm()`
#'    * `"limma_batch_adjust"`: batch adjustment via
#'    `limma::removeBatchEffect()`, recommended for data visualization,
#'    but not recommended for downstream statistical comparisons.
#'    * `"TMM"`: trimmed mean of M-values via `edgeR::calcNormFactors()`
#'    * `"TMMwsp"`: TMM with singleton pairing via `edgeR::calcNormFactors()`
#'    * `"RLE"`: relative log expression via `edgeR::calcNormFactors()`
#' @param assay_names `character` vector or one or more `names(assays(se))`
#'    that indicates which numeric matrix to use during normalization. When
#'    multiple values are provided, each matrix is normalized independently
#'    by each `method`.
#' @param output_method_prefix `character` vector (optional) with custom
#'    method prefix values to use when creating the new `assay_name` for
#'    each normalization. It must have length equal to `length(method)`,
#'    to be applied to each method in order.
#'    Note that `output_assay_names` takes priority, and when it is defined
#'    the `output_method_prefix` entries are ignored.
#'
#'    Consider these arguments:
#'    ```R
#'    assay_name="counts",
#'    method="limma_batch_adjust",
#'    output_method_prefix="lba"
#'    ```
#'    The assay_name created during normalization will be `"lba_counts"`.
#' @param output_assay_names `character` vector (optional) which overrides
#'    the default method for defining assay names for normalized data.
#'    This vector length must equal `length(method) * length(assay_names)`,
#'    and will be applied in the order data is normalized:
#'    1. `assay_names` are iterated.
#'    2. For each value in `assay_names`, each normalization in `method`
#'    is applied.
#'
#'    Therefore the order of `output_assay_names` could follow this order:
#'    `method1_assay1`, `method1_assay2`, `method2_assay1`, `method2_assay2`.
#' @param genes `character` vector (optional) used to define a subset of
#'    gene rows in `se` to use for normalization.
#'    Values must match `rownames(se)`.
#' @param samples `character` vector (optional) used to define a subset of
#'    sample columns in `se` to use for normalization.
#'    Values must match `colnames(se)`.
#' @param params `list` (optional) parameters specific to each
#'    normalization method, passed to `matrix_normalize()`. Any
#'    value which is not defined in the `params` provided will use
#'    the default value in `matrix_normalize()`, for example
#'    `params=list(jammanorm=list(minimum_mean=2))` will use
#'    `minimum_mean=2` then use other default values relevant
#'    to the `jammanorm` normalization method.
#' @param normgroup `character` or equivalent vector that defines subgroups
#'    of `samples` to be normalized indendently of each normgroup, or
#'    `character` vector with one or more colnames in `colData(se)`.
#'    * When `normgroup` is `NULL` all data is normalized together by default.
#'    * When supplying colnames in `colData(se)`, values are combined using
#'    `jamba::pasteByRow()` which uses `_` underscore delimiter between
#'    non-empty fields in each column.
#'    * When supplying a vector of normgroup values, it should be equal
#'    to `length(samples)` or should be named using values in `samples`
#'    such that all `samples` are present in `names(normgroup)`.
#'    * Note that when data are normalized in independent normgroups,
#'    these normgroups are not normalized relative to each other.
#'    Data are only normalized within each independent normgroup.
#' @param output_sep `character` string used as a delimited between the
#'    `method` and the `assay_names` to define the output assay name,
#'    for example when `assay_name="counts"`, `method="quantile"`,
#'    and `output_sep="_"` the new assay name will be `"quantile_counts"`.
#' @param override `logical` indicating whether to override any pre-existing
#'    matrix values with the same output assay name. When `override=FALSE`
#'    and the output assay name already exists, the normalization will
#'    not be performed.
#' @param populate_mcols `logical`, default TRUE, whether to populate
#'    normalization details into `mcols(assays(se))`, including
#'    the normalization `method`, the source `assay_name` used during
#'    normalization, and values from `params`.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to `matrix_normalize()`.
#'
#' @examples
#' if (jamba::check_pkg_installed("farrisdata")) {
#'
#'    # se_normalize
#'    # suppressPackageStartupMessages(library(SummarizedExperiment))
#'    GeneSE <- farrisdata::farrisGeneSE;
#'    samples <- colnames(GeneSE);
#'    genes <- rownames(GeneSE);
#'
#'    GeneSE <- se_normalize(GeneSE,
#'       genes=genes,
#'       samples=samples,
#'       assay_names=c("raw_counts", "counts"),
#'       method="jammanorm",
#'       params=list(jammanorm=list(minimum_mean=5)))
#'    SummarizedExperiment::mcols(SummarizedExperiment::assays(GeneSE))
#'    names(SummarizedExperiment::assays(GeneSE))
#'
#'    # review normalization factor values
#'    round(digits=3, attr(
#'       SummarizedExperiment::assays(GeneSE)$jammanorm_raw_counts, "nf"))
#'
#'    # the data in "counts" was already normalized
#'    # so the normalization factors are very near 0 as expected
#'    round(digits=3,
#'       attr(SummarizedExperiment::assays(GeneSE)$jammanorm_counts, "nf"))
#'
#'
#'    # note that housekeeper genes are supplied in params
#'    # also this demonstrates output_method_prefix
#'    set.seed(123);
#'    hkgenes <- sample(rownames(GeneSE), 1000)
#'    GeneSE <- se_normalize(GeneSE,
#'       genes=genes,
#'       samples=samples,
#'       assay_names=c("raw_counts"),
#'       method="jammanorm",
#'       output_method_prefix="hkjammanorm",
#'       params=list(jammanorm=list(minimum_mean=5,
#'          controlGenes=hkgenes)))
#'    SummarizedExperiment::mcols(SummarizedExperiment::assays(GeneSE))
#'
#'    # example showing quantile normalization
#'    GeneSE <- se_normalize(GeneSE,
#'       assay_names=c("raw_counts"),
#'       method="quantile")
#'    SummarizedExperiment::mcols(SummarizedExperiment::assays(GeneSE))
#'
#'    # example showing quantile normalization with custom output_assay_names
#'    GeneSE <- se_normalize(GeneSE,
#'       assay_names=c("raw_counts"),
#'       method="quantile",
#'       output_assay_names="newquantile_raw_counts")
#'    SummarizedExperiment::mcols(SummarizedExperiment::assays(GeneSE))
#' }
#'
#'
#' @export
se_normalize <- function
(se,
 method=c("quantile",
    "jammanorm",
    "limma_batch_adjust",
    "TMM",
    "TMMwsp",
    "RLE"),
 assay_names=NULL,
 output_method_prefix=NULL,
 output_assay_names=NULL,
 genes=NULL,
 samples=NULL,
 params=list(
    `quantile`=list(
       ties=TRUE),
    `jammanorm`=list(controlGenes=NULL,
       minimum_mean=0,
       controlSamples=NULL,
       centerGroups=NULL,
       useMedian=FALSE,
       noise_floor=NULL,
       noise_floor_value=NULL),
    `limma_batch_adjust`=list(
       batch=NULL,
       group=NULL),
    `TMM`=list(
       refColumn=NULL,
       logratioTrim=0.3,
       sumTrim=0.05,
       doWeighting=TRUE,
       Acutoff=NULL),
    `TMMwsp`=list(
       refColumn=NULL,
       logratioTrim=0.3,
       sumTrim=0.05,
       doWeighting=TRUE,
       Acutoff=NULL),
    `RLE`=list(
       refColumn=NULL,
       logratioTrim=0.3,
       sumTrim=0.05,
       doWeighting=TRUE,
       Acutoff=NULL)),
 normgroup=NULL,
 floor=0,
 enforce_norm_floor=TRUE,
 output_sep="_",
 override=TRUE,
 populate_mcols=TRUE,
 verbose=FALSE,
 ...)
{
   assay_names <- intersect(assay_names,
      names(SummarizedExperiment::assays(se)));
   if (length(assay_names) == 0) {
      cli::cli_abort(paste0(
         "{.var assay_names} must be supplied."));
      stop("assay_names must be supplied.");
   }
   method <- match.arg(method,
      several.ok=TRUE);

   if (length(output_method_prefix) > 0) {
      if (!length(output_method_prefix) == length(method)) {
         cli::cli_abort(paste0(
            "{.var length(method)}={length(method)}",
            " must equal {.var length(output_method_prefix)}=",
            "{length(output_method_prefix)}. ",
            "Cannot continue."))
      }
   } else {
      output_method_prefix <- method;
   }

   # validate output_sep
   if (length(output_sep) == 0) {
      output_sep <- "";
   }
   output_sep <- head(output_sep, 1);

   allgenes <- rownames(se);
   allsamples <- colnames(se);

   # validate genes and samples using rownames(se), colnames(se)
   if (length(genes) == 0) {
      genes <- allgenes;
   } else if (inherits(genes, "list")) {
      genes <- lapply(genes, function(igenes){
         intersect(igenes, allgenes)
      })
   } else {
      genes <- intersect(genes, allgenes);
   }
   if (length(samples) == 0) {
      samples <- allsamples;
   } else {
      samples <- intersect(samples, allsamples);
   }

   # check if any entries have length zero
   if (length(genes) == 0 || length(samples) == 0 ||
         (inherits(genes, "list") && any(lengths(genes) == 0))) {
      cli::cli_abort(paste0(
         "Only recognized {.var length(genes)}={length(genes)}",
         " and {.var length(samples)}={length(samples)} ",
         "in the input {.var se}. ",
         "Cannot continue."))
      stop(paste0("Only recognized ",
         length(genes),
         " genes, and ",
         length(samples),
         " samples in the input se."));
   }

   # optional normgroup
   if (length(normgroup) > 0) {
      if (all(normgroup %in% colnames(SummarizedExperiment::colData(se)))) {
         # 0.0.72.900 - recognize colnames in colData
         normgroup <- jamba::nameVector(
            jamba::pasteByRow(
               SummarizedExperiment::colData(
                  se[, samples])[, normgroup, drop=FALSE]),
            samples)
      } else if (length(names(normgroup)) > 0) {
         if (!all(samples %in% names(normgroup))) {
            cli::cli_abort(paste0(
               "When {.var names(normgroup)} is defined ",
               "all {.var samples} must be present in ",
               "{.var names(normgroup)}, ",
               "and this requirement was not met."));
            stop(paste(
               "When names(normgroup) is defined,",
               "all samples must be present in names(normgroup).",
               "This requirement was not met."));
         }
         normgroup <- normgroup[match(samples, names(normgroup))];
      } else {
         if (length(normgroup) != length(samples)) {
            cli::cli_abort(paste0(
               "When {.var names(normgroup)} is not defined, ",
               "values should match {.var colnames(colData(se))}, or ",
               "length(normgroup)=",
               "{length(normgroup)}",
               ", must equal ",
               "length(samples)=",
               "{length(samples)}."))
            stop(paste(
               "When names(normgroup) are not defined,",
               "length(normgroup) must equal length(samples).",
               "This requirement was not met."));
         }
      }
   }

   # define output_assay_names
   expected_output_length <- length(assay_names) * length(method);
   # define output_assay_names when none are supplied
   if (length(output_assay_names) == 0) {
      output_assay_names <- paste0(
         rep(output_method_prefix, length(assay_names)),
         output_sep,
         rep(assay_names, each=length(method)));
   }
   # define names(output_assay_names)
   if (length(output_assay_names) == expected_output_length) {
      names(output_assay_names) <- paste0(
         rep(assay_names, each=length(method)),
         "_",
         rep(method, length(assay_names)));
   } else {
      if (length(output_assay_names) > 0) {
         cli::cli_abort(paste0(
            "length(output_assay_names)=",
            "{length(output_assay_names)}",
            ", must equal ",
            "( length(assay_names) * length(method) )=",
            "{length(assay_names) * length(method)}."))
         stop(paste0("length(output_assay_names) must equal ",
            "length(assay_names) * length(method)"))
      }
   }

   # handle genes input as list
   normgroup_rows <- NULL;
   if (length(genes) > 0 && inherits(genes, "list")) {
      if (length(normgroup) == 0) {
         if (length(genes) == 1) {
            # use as vector
            genes <- unique(genes[[1]]);
         } else {
            # error
            cli::cli_abort(paste0(
               "normgroup is empty, genes should not be a list ",
               "with multiple entries."))
            stop(paste0("normgroup is empty, genes should not be a list ",
               "with multiple entries."))
         }
      } else {
         use_normgroups <- unique(as.character(normgroup));
         if (!all(use_normgroups %in% names(genes))) {
            # error, not all normgroup are defined in names(genes)
            cli::cli_abort(paste0(
               "normgroup values must all be present in {.var names(genes)}."))
            stop(paste0("normgroup values must all be present in names(genes)}."))
         } else {
            # define list for normgroup_rows
            normgroup_rows <- genes[use_normgroups];
            # use all rows in se
            genes <- rownames(se);
         }
      }
   }


   # iterate assay_names, within which iterate method
   for (assay_name in assay_names) {
      for (imethod in method) {
         output_assay_name <- unname(output_assay_names[
            paste0(assay_name, "_", imethod)]);
         if (output_assay_name %in% names(SummarizedExperiment::assays(se)) &&
               !TRUE %in% override) {
            if (verbose) {
               jamba::printDebug("se_normalize(): ",
                  sep="",
                  c("Skipped, normalized data already exists for '",
                     output_assay_name,
                     "'"));
            }
            next;
         }
         if (verbose) {
            jamba::printDebug("se_normalize(): ",
               c("Applying '", imethod, "' to '",
                  assay_name, "' to store in '",
                  output_assay_name,
                  "'"),
               sep="");
         }
         imatrix <- SummarizedExperiment::assays(
            se[genes, samples])[[assay_name]];
         inorm <- matrix_normalize(imatrix,
            method=imethod,
            params=params,
            normgroup=normgroup,
            normgroup_rows=normgroup_rows,
            verbose=(verbose - 1) > 0,
            floor=floor,
            enforce_norm_floor=enforce_norm_floor,
            ...);

         # generate matrix of NA values to fill for normalized genes, samples
         # this way we can add attributes to the full matrix
         # when normalizing a smaller matrix
         namatrix <- matrix(data=NA,
            ncol=length(allsamples),
            nrow=length(allgenes),
            dimnames=list(allgenes,
               allsamples));
         # assign normalized data for genes, samples
         namatrix[genes, samples] <- inorm;

         # re-assign any missing attributes
         # potentially useful information from the normalization method
         na_attrnames <- names(attributes(namatrix));
         new_attrs <- setdiff(names(attributes(inorm)), na_attrnames);
         if (length(new_attrs) > 0) {
            if (verbose > 1) {
               jamba::printDebug("se_normalize(): ",
                  "re-assigning new_attrs: ",
                  sep=", ",
                  new_attrs);
            }
            for (new_attr in new_attrs) {
               attr(namatrix, new_attr) <- attr(inorm, new_attr);
            }
         }

         # assign namatrix to se
         SummarizedExperiment::assays(se)[[output_assay_name]] <- namatrix;

         # Bonus points: populate mcols() with normalization details
         if (TRUE %in% populate_mcols) {
            # create default DataFrame if it does not exist
            if (length(SummarizedExperiment::mcols(
               SummarizedExperiment::assays(se))) == 0) {
               SummarizedExperiment::mcols(
                  SummarizedExperiment::assays(se)) <- S4Vectors::DataFrame(
                  assay_name=SummarizedExperiment::assayNames(se),
                  row.names=SummarizedExperiment::assayNames(se))
            }
            # define list to update mcols()
            mcol_list <- jamba::rmNULL(c(
               list(
                  assay_name=output_assay_name,
                  normalization_method=imethod,
                  source_assay_name=assay_name),
               params[[imethod]]));
            if (verbose > 1) {
               jamba::printDebug("se_normalize(): ",
                  "Populating mcol_list into mcols(assays(se)):");
               for (mname in names(mcol_list)) {
                  jamba::printDebug(paste0(mname, ": "),
                     indent=6,
                     mcol_list[[mname]])
               }
            }
            mcolnames <- colnames(SummarizedExperiment::mcols(
               SummarizedExperiment::assays(se)));
            for (newcolname in names(mcol_list)) {
               mvalues <- mcol_list[[newcolname]];
               mna <- head(c(NA, mvalues), 1);
               if (!newcolname %in% mcolnames) {
                  SummarizedExperiment::mcols(
                     SummarizedExperiment::assays(se))[,newcolname] <- mna;
               }
               if (length(mvalues) == 1) {
                  SummarizedExperiment::mcols(
                     SummarizedExperiment::assays(se))[
                        output_assay_name, newcolname] <- mvalues;
               } else if (length(mvalues) > 1) {
                  SummarizedExperiment::mcols(
                     SummarizedExperiment::assays(se))[
                        output_assay_name, newcolname] <- I(list(mvalues));
               }
            }
            # end names(mcol_list)
         }
         # end mcols(assays(se))
      }
   }
   return(se);
}

#' Normalize a numeric data matrix
#'
#' Normalize a numeric data matrix
#'
#' This function is a wrapper for several relevant normalization
#' methods that operate on a numeric matrix.
#'
#' # Normalization Methods Implemented:
#'
#' ## method='quantile'
#'
#' Quantile-normalization performed by
#' `limma::normalizeQuantiles()`. This method has one
#' parameter `"ties"` passed to `limma::normalizeQuantiles()`,
#' the default here `ties=TRUE` which handles tied numeric
#' expression values in a robust way to avoid unpredictability
#' otherwise. This option is especially relevant with expression
#' count data, where integer counts cause a large number
#' of values to be represented multiple times.
#'
#' ## method='jammanorm'
#'
#' Median-normalization performed by
#' `jamma::jammanorm()`. This method shifts expression
#' data as shown on MA-plots, so the median expression
#' is zero across all samples, using only the rows that
#' meet the relevant criteria.
#'
#' Some relevant criteria to
#' define rows used for normalization:
#'
#' * `controlGenes` defines specific genes to use for
#' normalization, such as housekeeper genes. It may also
#' be useful to use detected genes here, so the normalization
#' only considers those genes defined as detected by
#' the protocol.
#' * `minimum_mean` sets a numeric threshold and requires
#' the mean expression (shown on the x-axis of the MA-plot)
#' to be at least this value.
#'
#' Note that when both `controlGenes` and `minimum_mean`
#' are defined, both criteria are enforced. So the `controlGenes`
#' are also required to have expression of at least `minimum_mean`.
#'
#' Also note that all rows of data are normalized by this method,
#' only the subset of rows defined by `controlGenes` and `minimum_mean`
#' are used to compute the normalization factor.
#'
#'
#' ## method='limma_batch_adjust'
#'
#' Batch adjustment performed by
#' `limma::removeBatchEffect()` which is intended to apply
#' batch-adjustment as a form of normalization, but which
#' does not represent full normalization itself. There are
#' two relevant parameters: `"batch"` which is a vector of
#' batch values in order of `colnames(x)`, and `"group"`
#' which is a vector of sample groups in order of `colnames(x)`.
#'
#' # Normalization groups via `normgroup`
#'
#' The `normgroup` argument is intended as a convenient method to
#' apply a normalization method to each independent `normgroup`.
#' This situation is especially useful when a study contains
#' multiple tissue types, or multiple data types, that may not be
#' appropriate to normalize directly relative to one another.
#'
#' For example, one could normalize total RNA-seq and nascent 4sU-seq
#' data independently, without expectation that the two would ever
#' have a common frame of reference to normalize one relative to another.
#' However, both may be amenable to `"quantile"` or
#' `"jammanorm"` median normalization.
#'
#' Similarly, one could normalize each tissue type independently,
#' which may be appropriate when analyzing data that contains very
#' different mammalian tissue organ samples, such as muscle and brain.
#'
#' It would generally not be appropriate to use quantile normalization
#' across muscle and brain samples, since the overall pattern and
#' distribution of expression values is not expected to be similar.
#' Quantile normalize assumes (and imposes) a common distribution,
#' by adjusting mean expression signal at each quantile to a common
#' mean expression across all samples.
#'
#' For a rough approximation of cross-tissue normalization, one
#' could apply `"quantile"` normalization within each `normgroup` defined
#' by tissue type, then apply `"jammanorm"` median normalization to
#' apply a linear adjustment of signal across tissue types. The median
#' normalization does not affect distribution, thus will not affect
#' intra-tissue contrasts, except by adjusting its overall signal
#' which may change downstream assumptions regarding signal thresholds.
#'
#' It is recommended not to compare directly across tissue types.
#' In some cases a two-way contrast may be appropriate, where
#' fold change within one tissue type is compared to the
#' fold change within another tissue type. However, even in that
#' case the two tissue types do not need to be normalized relative to
#' each other upfront - the within-tissue fold change serves as
#' one method of normalizing the observations across tissue types.
#'
#' # Other useful parameters
#'
#' Note the `floor` and `enforce_norm_floor` have recommended
#' default values `floor=0` and `enforce_norm_floor=TRUE`.
#'
#' * `floor` is applied prior to normalization, typically
#' to minimize effects of low, noisy signal on the normalization
#' process itself. Specifically, this floor is used to remove negative
#' values, which may be by-products of upstream signal processing.
#' "A measured signal at or below the noise floor of a platform
#' technology is effectively the same as a signal at the noise floor."
#' * `enforce_norm_floor` is applied after normalization, typically
#' as a convenience, also to prevent low, noisy signal from
#' contributing to downstream analysis steps.
#'
#' These defaults will set any assay value at or below `0` to `0`,
#' and after normalization any values whose input values were
#' at or below `0` will also be set to `0` to prevent normalizing
#' a value of `0` to non-zero. Any normalized value at or
#' below `0` will also be set to `0` to prevent results from
#' containing negative normalized values.
#'
#' The assumption for this default is that a value of zero
#' is not a measurement but represents the lack of a measurement.
#' Similarly, the intent of `floor` is a numeric threshold at or
#' below there is no confidence in the reported measurement, therefore
#' values at or below this threshold are treated as equivalent
#' to the threshold for the purpose of downstream analyses.
#'
#' Some platforms like QPCR for example, have substantially lower
#' confidence at high CT values, where expression values
#' using the equation `2^(40-CT)` might impose a noise threshold
#' at expression 32 or lower. This noise threshold for QPCR
#' means any expression measurement of 32 or lower is as likely
#' to be `32` as it is to be `2`, and therefore any differences
#' between reported expression of `32` and `2` should not be
#' considered relevant. Applying `floor=32` in this case
#' accomplishes this goal by setting all values at or below
#' `32` to `32`. Of course when using this method `matrix_normalize()`
#' the data should be log2 transformed, which means the `floor`
#' should also be log2 transformed, e.g. `floor=log2(32)`
#' which is `floor=5`.
#'
#' One alternative might be to set values at or below zero to `NA`
#' prior to normalization, and before calling `matrix_normalize()`.
#' In this case, only non-NA values will be used during
#' normalization according to the `method` being used.
#'
#' @returns `numeric` matrix with the same dimensions as the
#'    input matrix `x`.
#'    Additional information may be returned as `attributes(x)`:
#'    * `"norm_method"`: a `character` string with the method used.
#'    * `"nf"`: a `numeric` vector with normalization factors, returned
#'    only by `"jammanorm"`, `"TMM"`, and `"TMMwsp"`.
#'    * `"hk"`: a `character` vector of `rownames(x)` used as housekeeper
#'    `controlGenes` by `"jammanorm"`.
#'
#' @family jamses utilities
#'
#' @param x `numeric` matrix with sample columns, and typically
#'    gene rows, but any measured assay row will meet the assumptions
#'    of the method.
#' @param method `character` string indicating which normalization
#'    method to apply.
#'    * `"quantile"`: quantile normalization via `limma::normalizeQuantiles()`
#'    * `"jammanorm"`: log-ratio normalization via `jamma::jammanorm()`
#'    * `"limma_batch_adjust"`: batch adjustment via
#'    `limma::removeBatchEffect()`, recommended for data visualization,
#'    but not recommended for downstream statistical comparisons.
#'    * `"TMM"`: trimmed mean of M-values via `edgeR::calcNormFactors()`
#'    * `"TMMwsp"`: TMM with singleton pairing via `edgeR::calcNormFactors()`
#'    * `"RLE"`: relative log expression via `edgeR::calcNormFactors()`
#' @param apply_log2 `character` string indicating whether to apply
#'    log2 transformation: `"ifneeded"` will apply log2 transform
#'    when any absolute value is greater than 40; `"no"` will not
#'    apply log2 transformation; `"always"` will apply log2 transform.
#'    Note the log2 transform is applied with `jamba::log2signed(x, offset=1)`
#'    which is equivalent to `log(1 + x)` except that negative values
#'    are also transformed using the absolute value, then multiplied
#'    by their original sign.
#' @param floor `numeric` value indicating the lowest accepted numeric
#'    value, below which values are assigned to this floor. The default
#'    `floor=0` requires all values are `0`, and any values below `0` are
#'    assigned `0`. Note that the `floor` is applied after log2 transform,
#'    when the log2 transform is performed.
#' @param floor_value `numeric` or `NA` used to replace values in `x`
#'    when they are at or below `floor`, default `floor_value=floor`,
#'    however it can be useful to assign `NA` to replace zero in
#'    circumstances when that is preferable.
#' @param enforce_norm_floor `logical` indicating whether to enforce the
#'    `floor` for the normalized results, default is `TRUE`. For example,
#'    when `floor=0` any values at or below `0` are set to `0` before
#'    normalization. After normalization some of these values will be
#'    above or below `0`. When `enforce_norm_floor=TRUE` these values
#'    will again be set to `0` because they are considered to be
#'    below the noise threshold of the protocol, and adjustments
#'    are not relevant; also any normalized values below the `floor`
#'    will also be set to `floor`.
#' @param params `list` of parameters relevant to the `method` of
#'    normalization. The `params` should be a `list` named by the `method`,
#'    whose values are a list named by the relevant method parameter.
#'    See examples.
#' @param normgroup `character` or equivalent vector that defines subgroups
#'    of `samples` to be normalized indendently of each normgroup.
#'    * When `normgroup` is `NULL` all data is normalized together as default.
#'    * The `normgroup` vector is expected to be in the order of
#'    `colnames(x)`.
#'    * Note that when data are normalized in independent normgroups,
#'    these normgroups are not normalized relative to each other.
#'    Data are only normalized within each independent normgroup.
#' @param normgroup_rows `list` with `character` values which must match
#'    `rownames(x)`, and where `names(normgroup_rows)` must also match
#'    the values in `normgroup`.
#'    * The purpose is to permit using independent rows for each normgroup,
#'    for example if the detected genes (rows) in each matrix differ
#'    for each normgroup, the normalization should only use each
#'    appropriate set of detected genes (rows).
#'    * When `normgroup_rows` are applied, data for rows not present in each
#'    particular normgroup will be set to `NA`, since those values were
#'    not included in the normalization for that normgroup.
#' @param subset_columns `integer` intended for internal use when
#'    `normgroup` is provided. This argument is used to instruct
#'    each normalization method to use an appropriate subset of
#'    `params` based upon the subset of columns being analyzed.
#' @param verbose `logical` indicating whether to print verbose output.
#'
#' @examples
#' # use farrisdata real world data if available
#' if (jamba::check_pkg_installed("farrisdata")) {
#'
#'    suppressPackageStartupMessages(library(SummarizedExperiment))
#'
#'    # test matrix_normalize()
#'    GeneSE <- farrisdata::farrisGeneSE;
#'    imatrix <- assays(GeneSE)$raw_counts;
#'    genes <- rownames(imatrix);
#'    samples <- colnames(imatrix);
#'    head(imatrix);
#'
#'    # matrix_normalize()
#'    # normalize the numeric matrix directly
#'    imatrix_norm <- matrix_normalize(imatrix,
#'       genes=genes,
#'       samples=samples,
#'       method="jammanorm",
#'       params=list(minimum_mean=5))
#'    names(attributes(imatrix_norm))
#'
#'    # review normalization factors
#'    round(digits=3, attr(imatrix_norm, "nf"));
#'
#'    # example for quantile normalization
#'    imatrix_quant <- matrix_normalize(imatrix,
#'       genes=genes,
#'       samples=samples,
#'       method="quantile")
#'    names(attributes(imatrix_quant))
#' }
#'
#'
#' # simulate reasonably common expression matrix
#' set.seed(123);
#' x <- matrix(rnorm(9000)/4, ncol=9);
#' colnames(x) <- paste0("sample", LETTERS[1:9]);
#' rownames(x) <- paste0("gene", jamba::padInteger(seq_len(nrow(x))))
#' rowmeans <- rbeta(nrow(x), shape1=2, shape2=5)*14+2;
#' x <- x + rowmeans;
#' for (i in 1:9) {
#'    x[,i] <- x[,i] + rnorm(1);
#' }
#'
#' # display MA-plot with jamma::jammaplot()
#' jamma::jammaplot(x)
#'
#' # normalize by jammanorm
#' xnorm <- matrix_normalize(x, method="jammanorm")
#' jamma::jammaplot(xnorm, maintitle="method='jammanorm'")
#'
#' # normalize by jammanorm with housekeeper genes
#' hk_genes <- sample(rownames(x), 10);
#' xnormhk <- matrix_normalize(x,
#'    method="jammanorm",
#'    params=list(jammanorm=list(controlGenes=hk_genes)))
#'
#' jamma::jammaplot(xnormhk,
#'    maintitle="method='jammanorm' with housekeeper genes",
#'    highlightPoints=list(housekeepers=hk_genes),
#'    highlightColor="green");
#'
#' xnormhk6 <- matrix_normalize(x,
#'    method="jammanorm",
#'    params=list(jammanorm=list(
#'       controlGenes=hk_genes,
#'       minimum_mean=6)))
#' hk_used <- attr(xnormhk6, "hk")[[1]];
#' jamma::jammaplot(xnormhk6,
#'    maintitle="method='jammanorm' with housekeeper genes, minimum_mean=6",
#'    highlightPoints=list(housekeepers=hk_genes,
#'       hk_used=hk_used),
#'    highlightColor=c("red", "green"));
#'
#' # normalize by quantile
#' xquant <- matrix_normalize(x, method="quantile")
#' jamma::jammaplot(xquant,
#'    maintitle="method='quantile'")
#'
#' # simulate higher noise for lower signal
#' rownoise <- rnorm(prod(dim(x))) * (3 / ((rowmeans*1.5) - 1.5));
#' xnoise <- x;
#' xnoise <- xnoise + rownoise;
#' jamma::jammaplot(xnoise,
#'    maintitle="simulated higher noise at lower signal");
#'
#' # simulate non-linearity across signal
#' # sin(seq(from=pi*4/10, to=pi*7/10, length.out=100))-0.8
#' rowadjust <- (sin(pi * jamba::normScale(rowmeans, from=3.5/10, to=5.5/10)) -0.9) * 20;
#' xwarp <- xnoise;
#' xwarp[,3] <- xnoise[,3] + rowadjust;
#' jamma::jammaplot(xwarp,
#'    maintitle="signal-dependent noise and non-linear effects");
#'
#' # quantile-normalization is indicated for this scenario
#' xwarpnorm <- matrix_normalize(xwarp,
#'    method="quantile");
#' jp <- jamma::jammaplot(xwarpnorm,
#'    maintitle="quantile-normalized: signal-dependent noise and non-linear effects");
#'
#' @export
matrix_normalize <- function
(x,
 method=c("quantile",
    "jammanorm",
    "limma_batch_adjust",
    "TMM",
    "TMMwsp",
    "RLE"),
 apply_log2=c("ifneeded",
    "no",
    "always"),
 floor=0,
 floor_value=floor,
 enforce_norm_floor=TRUE,
 params=list(
    #`vsn`=list(lts.quantile=0.5),
    #`rsn`=list(excludeFold=2,
    #   span=0.03),
    #cyclicLoess=list(method="default",
    #   span=0.8),
    `quantile`=list(
       ties=TRUE),
    `jammanorm`=list(controlGenes=NULL,
       minimum_mean=0,
       controlSamples=NULL,
       centerGroups=NULL,
       useMedian=FALSE,
       noise_floor=NULL,
       noise_floor_value=NULL),
    `limma_batch_adjust`=list(
       batch=NULL,
       group=NULL),
    `TMM`=list(
       refColumn=NULL,
       logratioTrim=0.3,
       sumTrim=0.05,
       doWeighting=TRUE,
       Acutoff=NULL),
    `TMMwsp`=list(
       refColumn=NULL,
       logratioTrim=0.3,
       sumTrim=0.05,
       doWeighting=TRUE,
       Acutoff=NULL)),
 normgroup=NULL,
 normgroup_rows=NULL,
 subset_columns=NULL,
 debug=FALSE,
 verbose=TRUE,
 ...)
{
   method <- match.arg(method);
   apply_log2 <- match.arg(apply_log2);

   ## Update the default methodParams without losing the original defaults if not overriden
   params <- update_function_params("matrix_normalize",
      "params",
      params,
      verbose=FALSE);

   # optional normgroup
   if (length(normgroup) > 0) {
      if (length(normgroup) != ncol(x)) {
         stop(paste("length(normgroup) was not equal to ncol(x)."));
      }
      x_colnames <- seq_len(ncol(x));
      x_split <- split(x_colnames, normgroup);
      x_split_v <- unlist(unname(x_split));
      x_split_order <- order(x_split_v);

      if (length(normgroup_rows) > 0) {
         # validate the names(normgroup_rows) match unique(normgroup)
         if (!all(unique(as.character(normgroup)) %in% names(normgroup_rows))) {
            cli::cli_abort(paste0(
               "Not all {.var normgroup} values are present in ",
               "{.var names(normgroup_rows)}."));
            stop(paste0("All normgroup values must be ",
               "provided in names(normgroup_rows)."));
         }
      }

      normgroup_out <- lapply(x_split, function(i_colnames){
         # 0.0.72.900 - optional normgroup_rows
         use_rows <- rownames(x);
         if (length(normgroup_rows) > 0) {
            use_normgroup <- as.character(normgroup[head(i_colnames, 1)]);
            use_rows <- normgroup_rows[[use_normgroup]];
         }
         # create NA matrix
         xnorm <- matrix(data=NA_real_,
            ncol=length(i_colnames),
            nrow=nrow(x));
         rownames(xnorm) <- rownames(x);
         colnames(xnorm) <- colnames(x)[i_colnames];
         # fill rows used for normalization
         xnorm[use_rows, i_colnames] <- matrix_normalize(
            x=x[use_rows, i_colnames, drop=FALSE],
            method=method,
            apply_log2=apply_log2,
            floor=floor,
            floor_value=floor_value,
            enforce_norm_floor=enforce_norm_floor,
            params=params,
            normgroup=NULL,
            normgroup_rows=NULL,
            subset_columns=i_colnames,
            verbose=verbose,
            ...);
         xnorm;
      });
      #return(normgroup_out);
      # combine normalized matrices
      normgroup_x <- do.call(cbind, normgroup_out)[, x_split_order, drop=FALSE];
      # combine attributes
      normgroup_attrnames <- unique(unlist(lapply(normgroup_out, function(inorm){
         attr_names <- setdiff(names(attributes(inorm)),
            c("dim", "names", "dimnames"))
      })))
      if (length(normgroup_attrnames) > 0) {
         normgroup_attrs <- lapply(jamba::nameVector(normgroup_attrnames), function(iattr){
            unlist(recursive=FALSE,
               lapply(unname(normgroup_out), function(inorm){
                  attributes(inorm)[[iattr]]
               }))[x_split_order]
         })
         for (iattr in normgroup_attrnames) {
            attr(normgroup_x, iattr) <- normgroup_attrs[[iattr]];
         }
      }
      if (debug) {
         return(list(
            normgroup_out=normgroup_out,
            normgroup_attrnames=normgroup_attrnames,
            normgroup_x=normgroup_x));
      }
      return(normgroup_x);
   }

   # apply log2 transform if needed
   if ("ifneeded" %in% apply_log2) {
      if (any(abs(x) > 40 & !is.na(x))) {
         x <- jamba::log2signed(x,
            offset=1);
         if (verbose) {
            jamba::printDebug("matrix_normalize(): ",
               c("Applied ",
                  "jamba::log2signed(x, offset=1)",
                  " because ",
                  "any(abs(x) > 40)"),
               sep="");
         }
      }
   } else if ("always" %in% apply_log2) {
      x <- jamba::log2signed(x,
         offset=1);
      if (verbose) {
         jamba::printDebug("matrix_normalize(): ",
            c("Applied ", "jamba::log2signed(x, offset=1)"),
            sep="");
      }
   }

   # apply optional floor
   if (length(floor) > 0 && any(x <= floor & !is.na(x))) {
      if (length(floor_value) == 0) {
         floor_value <- floor;
      }
      x_floored <- (x <= floor & !is.na(x));
      x[x_floored] <- floor_value;
      if (verbose) {
         jamba::printDebug("matrix_normalize(): ",
            c("Applied floor:", floor,
               ", replacing with:", floor_value),
            sep="");
      }
   } else {
      x_floored <- FALSE
   }

   if ("quantile" %in% method) {
      ties <- params$quantile$ties;
      if (length(ties) == 0) {
         ties <- TRUE;
      }
      inorm <- limma::normalizeQuantiles(A=x,
         ties=ties);
      attr(inorm, "norm_method") <- "quantile";
   } else if ("jammanorm" %in% method) {
      controlGenes <- params$jammanorm$controlGenes;
      minimum_mean <- params$jammanorm$minimum_mean;
      controlSamples <- params$jammanorm$controlSamples;
      centerGroups <- params$jammanorm$centerGroups;
      useMedian <- params$jammanorm$useMedian;
      noise_floor <- params$jammanorm$noise_floor;
      noise_floor_value <- params$jammanorm$noise_floor_value;
      if (length(noise_floor) == 0) {
         noise_floor <- -Inf;
      }
      if (length(noise_floor_value) == 0) {
         noise_floor_value <- noise_floor;
      }
      if (verbose) {
         jamba::printDebug("matrix_normalize(): ",
            "Calling ", "jammanorm():",
            c("\n      minimum_mean:", minimum_mean,
               "\n      useMedian:", useMedian,
               if (length(noise_floor) > 0 && noise_floor > -Inf) {
                  c("\n      noise_floor:", noise_floor)},
               if ((length(jamba::rmNA(noise_floor_value)) > 0 && jamba::rmNA(noise_floor_value) > -Inf) ||
                     any(is.na(noise_floor_value))) {
                  c("\n      noise_floor_value:", noise_floor_value)}),
            sep="");
      }
      # if no rownames(x) exist, use index integers
      x_rownames <- rownames(x);
      if (length(x_rownames) == 0) {
         rownames(x) <- seq_len(nrow(x));
      }
      inorm <- jamma::jammanorm(x,
         controlGenes=controlGenes,
         minimum_mean=minimum_mean,
         controlSamples=controlSamples,
         centerGroups=centerGroups,
         useMedian=useMedian,
         useMean=NULL,
         noise_floor=noise_floor,
         noise_floor_value=noise_floor_value,
         verbose=verbose);
      # if (nrow(inorm) == nrow(x) && length(x_rownames) == 0) {
      #    rownames(inorm) <- x_rownames;
      # }
      attr(inorm, "norm_method") <- "jammanorm";
   } else if ("limma_batch_adjust" %in% method) {
      batch <- params$limma_batch_adjust$batch;
      group <- params$limma_batch_adjust$group;
      if (length(subset_columns) > 0) {
         batch <- batch[subset_columns];
         group <- group[subset_columns];
      }

      designBatchGroup <- model.matrix(~0+group);
      rownames(designBatchGroup) <- colnames(x);

      ## limma::removeBatchEffect()
      inorm <- limma::removeBatchEffect(x,
         batch=batch,
         design=designBatchGroup);
      attr(inorm, "norm_method") <- "limma_batch_adjust";
   } else if (any(c("TMMwsp", "TMM", "RLE") %in% method)) {
      if (!jamba::check_pkg_installed("edgeR")) {
         stop("The edgeR package is required to normalize by 'TMM'");
      }
      # method parameters
      refColumn <- params$TMM$refColumn;
      logratioTrim <- params$TMM$logratioTrim;
      if (length(logratioTrim) == 0) {
         logratioTrim <- 0.3;
      }
      sumTrim <- params$TMM$sumTrim;
      if (length(sumTrim) == 0) {
         sumTrim <- 0.05;
      }
      doWeighting <- params$TMM$doWeighting;
      if (length(doWeighting) == 0) {
         doWeighting <- TRUE;
      }
      Acutoff <- params$TMM$Acutoff;
      if (length(Acutoff) == 0) {
         Acutoff <- -1e10;
      }
      # TMM disallows NA values, so we must replace with zero
      if (any(is.na(x))) {
         x[is.na(x)] <- 0;
      }
      # obtain normalization factors
      normFactors <- edgeR::calcNormFactors(
         object=2^x - 1,
         method=method,
         refColumn=refColumn,
         logratioTrim=logratioTrim,
         sumTrim=sumTrim,
         doWeighting=doWeighting,
         Acutoff=Acutoff);
      # apply normalization factors
      nfs <- log2(normFactors)
      inorm <- x + rep(nfs, each=nrow(x));
      # inorm <- x / rep(normFactors, each=nrow(x));

      # store normalization factors as attributes
      attr(inorm, "nf") <- nfs;
      attr(inorm, "normFactors") <- normFactors;
      attr(inorm, "logratioTrim") <- logratioTrim;
      attr(inorm, "sumTrim") <- sumTrim;
      attr(inorm, "doWeighting") <- doWeighting;
      attr(inorm, "Acutoff") <- Acutoff;
      attr(inorm, "norm_method") <- method;
   }

   # enforce the floor for output matrix
   if (enforce_norm_floor && any(x_floored)) {
      inorm[(x_floored | inorm <= floor) & !is.na(inorm)] <- floor_value;
   }

   return(inorm);
}

#' Update function default parameters
#'
#' Update function default parameters
#'
#' This function is a minor extension to `update_list_elements()`
#' intended to help update function parameters which are defined
#' as a nested list.
#'
#' @family jamses utilities
#'
#' @export
update_function_params <- function
(function_name=NULL,
   param_name=NULL,
   new_values=NULL,
   verbose=FALSE,
   ...)
{
   ## Purpose is to facilitate updating default parameters which are present in a list,
   ## but where defaults are defined in a function name.  So if someone wants to change
   ## one of the default values, but keep the rest of the list of defaults, this function
   ## does it.
   ##
   ## The default values are taken from the function formals, using eval(formals(functionName)).
   default_params <- eval(formals(function_name)[[param_name]]);
   if (verbose) {
      jamba::printDebug("update_function_params(): ",
         "str(default_params)");
      print(str(default_params));
   }
   default_params <- update_list_elements(default_params,
      new_values,
      verbose=verbose);
   return(default_params);
}

#' Update a subset of list elements
#'
#' Update a subset of list elements
#'
#' This function is intended to help update a nested `source_list`,
#' a subset of whose values should be replaced with entries
#' in `update_list`, leaving any original entries in `source_list`
#' which were not defined in `update_list`.
#'
#' This function may be useful when manipulating lattice or ggplot2
#' graphical parameters, which are often stored in a nested
#' list structure.
#'
#' @family jamses utilities
#'
#' @export
update_list_elements <- function
(source_list,
   update_list,
   list_layer_num=1,
   verbose=FALSE,
   ...)
{
   ## Purpose is to update elements in a list, allowing for multi-layered lists
   ## of lists. In case of a list-of-list, it will call this function again with
   ## each successive layer of listedness.
   ##
   ## Handy for updating lattice graphics settings, which are impossibly nested
   ## tangled ball of textual yarn.  An example:
   ## tp2 <- updateListElements(trellis.par.get(), list(fontsize=list(points=6)));
   ## trellis.par.set(tp2);
   ##
   ## Or in one line:
   ## trellis.par.set(updateListElements(trellis.par.get(), list(fontsize=list(points=6))));
   if (length(update_list) == 0) {
      return(source_list);
   }
   if ("list" %in% class(update_list) && "list" %in% class(update_list[[1]])) {
      for (update_list_name in names(update_list)) {
         if (update_list_name %in% names(source_list)) {
            ## If the name already exists, we must update items within the list
            source_list[[update_list_name]] <- update_list_elements(
               source_list=source_list[[update_list_name]],
               update_list=update_list[[update_list_name]],
               list_layer_num=list_layer_num+1);
         } else {
            ## If the name does not already exist, we can simply add it.
            source_list[[update_list_name]] <- update_list[[update_list_name]];
         }
      }
   } else {
      if (!is.null(names(update_list))) {
         source_list[names(update_list)] <- update_list;
      } else {
         source_list <- update_list;
      }
   }
   return(source_list);
}
