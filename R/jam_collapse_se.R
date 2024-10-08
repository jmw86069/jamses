
#se_collapse_by_row <- function

#' Collapse SummarizedExperiment data by row
#'
#' Collapse SummarizedExperiment data by row
#'
#' Purpose is to collapse rows of a `SummarizedExperiment` object,
#' where measurements for a given entity, usually a gene, are split
#' across multiple rows in the source data. The output of this function
#' should be measurements appropriately summarized to the gene level.
#'
#' The key arguments are `group_func_name`, and `data_transform`.
#' Note that data is inverse-transformed based upon `data_transform`,
#' prior to calculating group summary values defined by `group_func_name`.
#' The reason is to enable using `group_func_name="sum"` on normal
#' space abundance values, when input data has already been
#' transformed with `log2(1 + x)` for example. In this case it is most
#' appropriate to take the `sum` of normal space abundance values,
#' then to re-apply the transformation afterwards.
#'
#' However, when using `group_func_name="mean"` it is usually
#' recommended to use `data_transform="none"` so that data is maintained
#' in appropriately transformed state.
#'
#' The driving use case is proteomics mass spectrometry data,
#' where measurements are described in terms of peptide sequences,
#' with or without optional post-translational modification (PTM),
#' and the peptide sequences are annotated to a source protein or gene.
#' This function can be used to:
#'
#' * collapse peptide-PTM data to the peptide level
#' * collapse peptide data to the protein level
#'
#' In future it may be used to collapse multiple microarray probe
#' measurements to the gene level, although that process is more likely
#' to be useful and recommended after performing probe-level statistical
#' analysis.
#'
#' ### Proteomics mass spectrometry analysis
#'
#' For proteomics mass spectrometry data, proteins are inconsistently
#' fragmented into smaller peptides of varying sizes. The peptides are
#' usually separated on a chromatography column, from which aliquot
#' fractions are taken and measured by mass spectrometry. The total
#' signal derived from the original protein is therefore some combination
#' of the measured peptide parts.
#'
#' In some upstream data processing tools, such as Proteomics Discoverer,
#' and PEAKS, the peptide data may be annotated with observed
#' modification events (PTM). In this scenario, peptide measurements
#' are split across multiple rows of data, where each
#' row represents an observed combination of peptide and PTMs.
#'
#' ### Collapse methods
#'
#' It is fairly straightforward to observe peptide-PTM measurement
#' data is correlated with overall protein quantification, and that
#' the specific combination of peptide fragments may be inconsistent
#' across samples. That is, one may observe five peptides of protein A
#' in one sample, and may observe seven peptides of protein A in
#' another sample. The quantities of each peptide may be inconsistent,
#' due to variability in protein fragmentation across samples. However,
#' the general sum of peptide measurements is typically fairly stable
#' across samples, especially for proteins of moderate to high abundance
#' which are known to have stable abundance per cell.
#'
#' Choice of method to collapse measurements is not trivial, and is
#' therefore configurable. In general, proteomics abundances are
#' analyzed after `log2( 1 + x )` transformation. However, measurements
#' cannot be summed in log2 form, which would be equivalent to
#' multiplying measurements in normal form. Measurements can be summed
#' but only after exponentiating the data, for example the reciprocal
#' `( 2 ^ x ) - 1` is sufficient.
#'
#' @family jamses utilities
#'
#' @return `SummarizedExperiment` object with these changes:
#'    * rows will be collapsed by `row_groups`, for each `assays(se)`
#'    `numeric` matrix defined by `assay_names`. The collapse may
#'    optionally apply a data transformation defined in
#'    `data_transform` in order to apply an appropriate `numeric` summary
#'    calculation.
#'    * `rowData(se)` will also be collapsed by `shrinkDataFrame()` to
#'    combine unique values from each row annotation.
#'
#' @param se `SummarizedExperiment`
#' @param rows `character` vector of `rows(se)` to use for analysis. When
#'    `rows=NULL` the default is to use all `rows(se)`.
#' @param row_groups `character` vector representing groups of rows to
#'    be combined.
#' @param assay_names `character` vector of `names(assays(se))` to use for
#'    the collapse operation. When `assay_names=NULL` the default is to
#'    use all `assays(se)`.
#' @param group_func_name `character` name of function used to aggregate
#'    measurement data within `row_groups`.
#'    * `sum` - takes the `sum()` of each value in the group. This option
#'    should be used together with `data_transform` when there has been
#'    any data transformation, so that the data is inverse-transformed
#'    prior to calculating the `sum()`, after which data is re-transformed
#'    to its original state. This method is appropriate for log2p
#'    `log2(1 + x)` transformed abundance measurements for example.
#'    * `mean` - calculates the mean value per group. Note that in this
#'    case is it usually recommended not to define `data_transform`
#'    so that values are averaged in the appropriately transformed
#'    numeric space.
#'    * `weighted.mean` - calculates `weighted.mean()` where weights `w`
#'    are defined by the values used. This method may be appropriate and
#'    effective with normal space abundance values derived from
#'    proteomics mass spec quantitation.
#'    * `geomean` - calculates geometric mean of values in each group.
#'    * `none` -
#' @param rowStatsFunc `function` optional function used instead of
#'    `group_func_name`.
#' @param rowDataColnames `character` subset of colnames in `rowData(se)`
#'    to be retained in the output data. Multiple values are combined
#'    usually by comma-delimited concatenation within `row_groups`,
#'    therefore it may be beneficial to include only relevant columns
#'    in that output.
#' @param keepNULLlevels `logical` indicating whether to drop unused
#'    factor levels in `row_groups`, this argument is passed to
#'    `jamba::rowGroupMeans()`.
#' @param delim `character` string indicating a delimiter.
#' @param data_transform `character` string indicating which
#'    transformation was used when preparing the assay data. The
#'    assumption is that all assays were transformed by this method.
#'    During processing, data is inverse-transformed prior to applying
#'    the `group_func_name` or `rowStatsFunc` if supplied. After that
#'    function is applied, data is transformed using this function.
#'    The purpose is to enable taking the `sum()` in proper measured
#'    absolute units (in normal space for example) where relevant,
#'    after which is original numeric transformation is re-applied.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to `jamba::rowGroupMeans()`.
#'
#' @export
se_collapse_by_row <- function
(se,
 rows=rownames(se),
 row_groups,
 assay_names=NULL,
 group_func_name=c("sum",
    "mean",
    "weighted.mean",
    "geomean",
    "none"),
 rowStatsFunc=NULL,
 rowDataColnames=NULL,
 keepNULLlevels=FALSE,
 delim="[ ]*[;,]+[ ]*",
 data_transform=c("none",
    "log2p+sqrt",
    "log2+sqrt",
    "log2p",
    "log2"),
 verbose=TRUE,
 ...)
{
   ## Purpose is to collapse rows of a SummarizedExperiment object, intended
   ## for proteomics data to combine per-peptide rows into per-protein rows
   ##
   ## SE is a SummarizedExperiment object
   ##
   ## rows is a vector of rownames(SE) or a subset of rows to include
   ##
   ## row_groups is a vector of length(rows) whose names are rownames(SE)
   ## that represents the row grouping. For example:
   ##    rows <- c("Apoe_pep1", "Apoe_pep2", "Gapdh_pep1")
   ##    row_groups <- c(Apoe_pep1="Apoe", Apoe_pep2="Apoe", Gapdh_pep1="Gapdh")
   ##
   ## signal a vector of names(assays(SE)) or subset to include in the
   ## operation, any names not included will be dropped from the resulting
   ## SummarizedExperiment object.
   ##
   ## group_func_name is a vector of character names of row stats functions,
   ## provided for convenience. If group_func_name is length=1, it is used
   ## for all signal values. Otherwise, group_func_name is expected to be
   ## a vector named by signal, allowing each signal to use a different
   ## row stats function appropriate to that signal type.
   ##
   ## rowStatsFunc is an optional custom function, or list of custom
   ## functions named by signal, intended to produce a single numeric
   ## value per row, for example rowMeans() or rowMedians(). It
   ## is sent to rowGroupMeans(...,rowStatsFunc=rowStatsFunc) and must
   ## accept matrix input, and perform row-wise operation to produce a
   ## single numeric value per row.
   ## Note that if a list is supplied, any names matching signal will
   ## cause that function to be used instead of the corresponding
   ## group_func_name function.
   ##
   ## keepNULLlevels=TRUE, if groups is a factor, then all factor levels
   ## will be maintained even if no values exist with that factor level.
   ## keepNULLlevels=FALSE, if groups is a factor, then only factor levels
   ## containing values will be returned.
   ##
   ## rowDataColnames is a vector of colnames(rowData(SE)) or subset, which
   ## should be included in the resulting collapsed rowData data.frame. It
   ## is useful to omit any column no longer necessary because it saves
   ## substantial space and time in this operation.
   ##
   ## data_transform describes how input expression values
   ## were transformed, so that they can be un-transformed
   ## before applying the math operation. Use "none" to perform
   ## no additional transformations before applying the math.
   ##
   ## data_transform "log2+sqrt" means data was transformed log2(1+sqrt(x))
   ## which means data was sqrt() transformed, then log2(1+x)
   ## transformed. Therefore, data will be exponentiated with
   ## (2^(x)-1)^2 which reverses the log2 transform with "(2^x)-1"
   ## then reverses the sqrt by raising to the power of 2 with "x^2".
   ##
   if (!jamba::check_pkg_installed("SummarizedExperiment")) {
      stop("se_collapse_by_row() requires SummarizedExperiment");
   }

   ## Parameter validation
   ##
   if (length(assay_names) == 0) {
      assay_names <- names(SummarizedExperiment::assays(se));
   } else {
      assay_names <- intersect(assay_names, names(SummarizedExperiment::assays(se)));
   }
   if (length(assay_names) == 0) {
      stop("se_collapse_by_row() needs assay_names in names(assays(se)).");
   }

   if (length(rowDataColnames) == 0) {
      rowDataColnames <- colnames(SummarizedExperiment::rowData(se));
   }

   ## group_func_name will become a vector of text or NA values, named by signal
   group_func_args <- c("sum",
      "mean",
      "weighted.mean",
      "geomean",
      "none");
   if (identical(group_func_args, group_func_name)) {
      group_func_name <- head(group_func_name, 1);
   }
   if (length(group_func_name) > 0) {
      if (length(names(group_func_name)) == 0) {
         ## If given group_func_name without names, assign names by signal
         group_func_name <- rep(group_func_name, length.out=length(assay_names));
         names(group_func_name) <- assay_names;
      }
   } else {
      group_func_name <- NA;
   }
   group_func_name <- jamba::nameVector(group_func_name[assay_names],
      assay_names);
   group_func_name[group_func_name %in% c("NA","","none","custom")] <- NA;

   ## rowStatsFunc will become a list named by signal after this section
   if (length(rowStatsFunc) > 0) {
      if (jamba::igrepHas("function", class(rowStatsFunc))) {
         if (length(group_func_name) == 0) {
            customSignal <- assay_names;
         } else {
            customSignal <- assay_names[group_func_name %in%
                  c(NA, "", "none", "NA")];
         }
         rowStatsFunc <- lapply(jamba::nameVector(customSignal), function(i){
            rowStatsFunc
         });
      } else if (!jamba::igrepHas("list", class(rowStatsFunc))) {
         stop(paste0("rowStatsFunc must be NULL, a function,",
            " or a list of functions"));
      }
      rowStatsFunc <- rowStatsFunc[assay_names];
      names(rowStatsFunc) <- assay_names;
      rowStatsFunc;
   } else {
      rowStatsFunc <- lapply(jamba::nameVector(assay_names), function(i)NULL);
   }

   ## data.frame summary of which stats function to use per signal
   statsFuncDF <- as.data.frame(jamba::rmNULL(
      list(group_func_name=group_func_name,
         rowStatsFunc=sapply(rowStatsFunc[assay_names], class))));
   if (verbose) {
      jamba::printDebug("se_collapse_by_row(): ",
         "Performing these operations by assay_names:");
      print(statsFuncDF);
   }

   ###############################################################
   ## if supplied group_func_name, define the relevant function
   ## sum
   rowSum <- function(x, ..., na.rm=TRUE){
      rowSums(x,
         ...,
         na.rm=na.rm)
   }
   ## mean
   rowMean <- function(x, ..., na.rm=TRUE){
      rowMeans(x,
         ...,
         na.rm=na.rm)
   }
   ## weighted.mean
   rowWeightedMean <- function(x, w=x, ..., na.rm=TRUE){
      #isNA <- is.na(x);
      #x <- x[!isNA];
      #w[is.na(w)] <- 0;
      apply(x, 1, function(x1){
         weighted.mean(x=x1,
            w=x1,
            na.rm=TRUE,
            ...);
      })
      #rowWeightedMeans(x=x,
      #   w=w,
      #   ...,
      #   na.rm=TRUE);
   }
   # custom function geomean()
   geomean <- function
   (x,
    na.rm=TRUE,
    naValue=NA,
    offset=0,
    ...)
   {
      ## Purpose is to calculate the classical geometric mean
      if (na.rm && any(is.na(x))) {
         x <- jamba::rmNA(x,
            naValue=naValue);
      }
      2^mean(log2(x + offset)) - offset;
   }
   # custom function rowGeomeans()
   rowGeomeans <- function
   (x,
    na.rm=TRUE,
    verbose=FALSE,
    ...)
   {
      ## Purpose is to provide a wrapper to geomean(x, matrixBy="row", ...)
      geomean(x=x,
         matrixBy="row",
         na.rm=na.rm,
         verbose=verbose,
         ...)
   }

   ###############################################################
   ## populate rowStatsFunc
   for (iSignal in assay_names) {
      if (statsFuncDF[iSignal,2] %in% "NULL") {
         if (is.na(statsFuncDF[iSignal,1])) {
            if (verbose) {
               jamba::printDebug("se_collapse_by_row(): ",
                  "Dropping assay_names:",
                  iSignal,
                  " which has no stats function defined.",
                  fgText=c("yellow","red"));
            }
            assay_names <- setdiff(assay_names, iSignal);
         } else {
            if ("sum" %in% statsFuncDF[iSignal,1]) {
               rowStatsFunc[[iSignal]] <- rowSum;
            } else if ("mean" %in% statsFuncDF[iSignal,1]) {
               rowStatsFunc[[iSignal]] <- rowMean;
            } else if ("weighted.mean" %in% statsFuncDF[iSignal,1]) {
               rowStatsFunc[[iSignal]] <- rowWeightedMean;
            } else if ("geomean" %in% statsFuncDF[iSignal,1]) {
               rowStatsFunc[[iSignal]] <- rowGeomeans;
            } else {
               if (verbose) {
                  jamba::printDebug("se_collapse_by_row(): ",
                     "Dropping assay_names:",
                     iSignal,
                     " unrecognized group_func_name:",
                     statsFuncDF[iSignal,1],
                     fgText=c("yellow","red"));
               }
               assay_names <- setdiff(assay_names, iSignal);
            }
         }
      }
   }
   if (length(assay_names) == 0) {
      stop("se_collapse_by_row() had no assay_names matching group_func_name or rowStatsFunc.");
   }

   ###############################################################
   # Ensure we use only rows which are present in rownames(se)
   # and in names(row_groups)
   if (length(names(row_groups)) == 0) {
      if (length(row_groups) == length(rows)) {
         names(row_groups) <- rows;
      } else if (length(row_groups) == nrow(se)) {
         names(row_groups) <- rownames(se);
      }
   }
   rows <- intersect(rows, rownames(se));
   rows <- intersect(rows, names(row_groups));
   if (length(rows) == 0) {
      stop("se_collapse_by_row() needs 'rows' present in rownames(se) and names(row_groups)");
   }
   row_groups <- row_groups[rows];
   # any row_group with no value will not be collapsed,
   # it will retain its original value
   if (any(nchar(row_groups) == 0)) {
      iRG <- (is.na(row_groups) | nchar(row_groups) == 0);
      if (any(iRG)) {
         row_groups[iRG] <- rows[iRG];
         if (verbose) {
            jamba::printDebug("se_collapse_by_row(): ",
               "Note: Replaced ",
               "row_groups[nchar(row_groups) == 0]",
               " with the corresponding value from rows",
               fgText=c("orange","red"));
         }
         #stop("se_collapse_by_row() row_groups contains some empty values, any(nchar(row_groups)=0)");
      }
   }

   ###############################################################
   ## rowDataColnames the colnames to keep in the new rowData
   se_rowDataColnames <- colnames(SummarizedExperiment::rowData(se));
   rowDataColnames <- intersect(rowDataColnames, se_rowDataColnames);
   if (length(rowDataColnames) == 0 && length(se_rowDataColnames) > 0) {
      rowDataColnames <- colnames(SummarizedExperiment::rowData(se));
   }

   #################################################################
   ## Iterate each signal and perform the collapse
   ## Any signals not included will be dropped from the resulting object
   assaysGroupedL <- lapply(jamba::nameVector(assay_names), function(iSignal){
      if (verbose) {
         jamba::printDebug("se_collapse_by_row(): ",
            "Collapsing assay_names:",
            iSignal);
      }
      useStatsFunc <- rowStatsFunc[[iSignal]];

      iMatrix <- SummarizedExperiment::assays(se[rows,])[[iSignal]];
      if (data_transform %in% "log2p+sqrt") {
         if (verbose) {
            jamba::printDebug("se_collapse_by_row():",
               "Applying de-transform: ",
               "( 2 ^ x - 1 ) ^2");
         }
         ## log2(1+sqrt(x))
         iMatrix <- (2 ^ iMatrix - 1) ^ 2;
      } else if (data_transform %in% "log2+sqrt") {
         if (verbose) {
            jamba::printDebug("se_collapse_by_row():",
               "Applying de-transform: ",
               "( 2 ^ x ) ^ 2");
         }
         ## log2(sqrt(x))
         iMatrix <- (2 ^ iMatrix) ^ 2;
      } else if (data_transform %in% "log2p") {
         ## log2( 1 + x )
         if (verbose) {
            jamba::printDebug("se_collapse_by_row():",
               "Applying de-transform: ",
               "( 2 ^ x ) - 1");
         }
         iMatrix <- 2 ^ iMatrix - 1;
      } else if (data_transform %in% "log2") {
         ## log2( x )
         if (verbose) {
            jamba::printDebug("se_collapse_by_row():",
               "Applying de-transform: ",
               "( 2 ^ x )");
         }
         iMatrix <- 2 ^ iMatrix;
      } else {
         if (verbose) {
            jamba::printDebug("se_collapse_by_row():",
               "No de-transform was applied.");
         }
      }
      iMatrixGrp <- t(jamba::rowGroupMeans(
         x=t(iMatrix),
         groups=row_groups,
         rowStatsFunc=useStatsFunc,
         verbose=FALSE,
         keepNULLlevels=keepNULLlevels,
         ...
      ));
      # re-apply the appropriate transformation
      if (data_transform %in% "log2p+sqrt") {
         ## log2(1+sqrt(x))
         if (verbose) {
            jamba::printDebug("se_collapse_by_row():",
               "Applying re-transform: ",
               "log2( 1 + sqrt(x) )");
         }
         iMatrixGrp <- log2( 1 + sqrt(iMatrixGrp));
      } else if (data_transform %in% "log2+sqrt") {
         ## log2(sqrt(x))
         if (verbose) {
            jamba::printDebug("se_collapse_by_row():",
               "Applying re-transform: ",
               "log2( sqrt(x) )");
         }
         iMatrixGrp <- log2( 1 + sqrt(iMatrixGrp) );
      } else if (data_transform %in% "log2p") {
         ## log2( 1 + x )
         if (verbose) {
            jamba::printDebug("se_collapse_by_row():",
               "Applying re-transform: ",
               "log2( 1 + x )");
         }
         iMatrixGrp <- log2( 1 + iMatrixGrp );
      } else if (data_transform %in% "log2") {
         ## log2( x )
         if (verbose) {
            jamba::printDebug("se_collapse_by_row():",
               "Applying re-transform: ",
               "log2( x )");
         }
         iMatrixGrp <- log2( iMatrixGrp );
      } else {
         if (verbose) {
            jamba::printDebug("se_collapse_by_row():",
               "No transform was re-applied.");
         }
      }
      # Revert zero to NA
      # - review whether this option should be configurable
      iMatrixGrp[iMatrixGrp == 0] <- NA;
      iMatrixGrp;
   });

   # Now try to be clever and create a new rowData() which indicates
   # the row groupings
   if (verbose) {
      jamba::printDebug("se_collapse_by_row(): ",
         "Collapsing rowData(se)");
   }
   rowDataShrunk <- shrinkDataFrame(
      x=SummarizedExperiment::rowData(se[rows,])[, rowDataColnames, drop=FALSE],
      groupBy=row_groups,
      include_num_reps=TRUE,
      verbose=(verbose - 1) > 0,
      ...);

   ## Optional polish step to re-delimit column values to make them unique
   if (length(delim) > 0) {
      for (iCol in colnames(rowDataShrunk)) {
         if (jamba::igrepHas("character", class(rowDataShrunk[[iCol]])) &&
               jamba::igrepHas(delim, rowDataShrunk[[iCol]])) {
            if (verbose) {
               jamba::printDebug("se_collapse_by_row(): ",
                  "Re-delimiting unique values in column:",
                  iCol);
            }
            ## split fields using delim, then call cPasteUnique()
            ## to paste only unique values back together
            rowDataShrunk[[iCol]] <- jamba::cPasteU(
               strsplit(rowDataShrunk[[iCol]], delim),
               na.rm=TRUE);
         }
      }
   }
   if (verbose) {
      jamba::printDebug("se_collapse_by_row(): ",
         "Re-creating SummarizedExperiment");
   }

   se_shrunk <- SummarizedExperiment::SummarizedExperiment(
      assays=assaysGroupedL,
      rowData=rowDataShrunk,
      colData=SummarizedExperiment::colData(se),
      metadata=list(
         genes=rownames(rowDataShrunk),
         samples=colnames(se)));
   #rownames(SEshrunk) <- rownames(rowDataShrunk);
   se_shrunk;
}
