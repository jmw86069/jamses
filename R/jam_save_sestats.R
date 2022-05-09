
# jam_save_sestats.R

#' Save SE contrast stats output
#'
#' Save SE contrast stats output
#'
#' @param sestats `list` object output from `se_contrast_stats()`
#' @param file `character` string indicating the filename to save.
#' @param assay_names `character` string indicating which assay names
#'    to save, stored in `dimnames(sestats$hit_array)$Signal`.
#'    When `NULL` then all assay names are saved.
#' @param contrast_names `character` string indicating which contrasts
#'    to save, stored in `dimnames(sestats$hit_array)$Contrasts`.
#'    When `NULL` then all contrasts are saved.
#' @param type `character` string indicating the type of file to save.
#'    * `"xlsx"`: saves an Excel xlsx file using `jamba::writeOpenxlsx()`.
#' @param max_nchar_sheetname `integer` number of characters allowed in
#'    MS Excel worksheet names, currently 31 characters.
#' @param review_output `logical` indicating whether a summary of output
#'    should be returned as a `data.frame` without exporting data. This
#'    summary will indicate all worksheets to be saved, in addition
#'    to the sheetName for each worksheet.
#' @param sheet_prefix `character` string with optional character prefix
#'    to use when creating worksheet names.
#' @param use_assay_suffix `logical` indicating whether to include
#'    `assay_names` as suffix when forming sheet names, when there is more
#'    than one unique assay name to be saved. This step will
#'    attempt to abbreviate `assay_names` by taking up to 4 characters
#'    from each word in the assay name, where each word is
#'    delimited by `"[-.:_ ]+"`.
#'    Otherwise, sheet names are forced to be unique by taking a substring
#'    of the contrast name of up to `max_nchar_sheetname`, passing any
#'    duplicate strings to `jamba::makeNames()` with suffix `"_v"`
#'    followed by an integer number.
#' @param width_factor `numeric` used to adjust relative column widths
#'    in the output Excel worksheets.
#' @param colorSub `character` vector of colors, optional, used to define
#'    categorical background colors for text string fields in Excel.
#'    The `names(colorSub)` are matched to character strings to assign
#'    colors.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to `jamba::writeOpenxlsx()`
#'
#' @export
save_sestats <- function
(sestats,
 file,
 assay_names=NULL,
 contrast_names=NULL,
 cutoff_names=NULL,
 type=c("xlsx"),
 max_nchar_sheetname=31,
 review_output=TRUE,
 sheet_prefix=NULL,
 use_assay_suffix=TRUE,
 width_factor=1,
 max_rows=NULL,
 colorSub=NULL,
 rename_contrasts=TRUE,
 se=NULL,
 rowData_colnames=NULL,
 row_type="gene_name",
 verbose=FALSE,
 ...)
{
   #
   type <- match.arg(type);

   # determine which results to save
   # validate assay_names
   assay_names <- intersect(assay_names,
      dimnames(sestats$hit_array)$Signal);
   if (length(assay_names) == 0) {
      assay_names <- dimnames(sestats$hit_array)$Signal;
   }

   # validate contrast_names
   contrast_names <- intersect(contrast_names,
      dimnames(sestats$hit_array)$Contrasts);
   if (length(contrast_names) == 0) {
      contrast_names <- dimnames(sestats$hit_array)$Contrasts;
   }

   # validate cutoff_names
   cutoff_names <- intersect(cutoff_names,
      dimnames(sestats$hit_array)$Cutoffs);
   if (length(cutoff_names) == 0) {
      cutoff_names <- head(dimnames(sestats$hit_array)$Cutoffs, 1);
   }

   # assembly export_df to describe what data will be exported
   # assay_names <- paste0("assay", 1);
   # cutoff_names <- paste0("cutoff", 1);
   # contrast_names <- paste0("contrast", 1:5);
   export_df <- data.frame(check.names=FALSE,
      assay_names=rep(assay_names,
         each=length(cutoff_names) * length(contrast_names)),
      cutoff_names=rep(
         rep(cutoff_names,
            each=length(contrast_names)), length(assay_names)),
      contrast_names=rep(
         rep(contrast_names,
            length(assay_names) * length(cutoff_names)))
      )
   # create sheet names
   export_df$sheetName <- export_df$contrast_names;

   # optionally rename contrasts
   if (rename_contrasts) {
      export_df$sheetName <- tryCatch({
         gsub("[:]+", ";",
            contrast2comp(export_df$contrast_names))
      }, error=function(e){
         export_df$contrast_names
      });
   }

   if (length(sheet_prefix) > 0) {
      export_df$sheetName <- paste0(sheet_prefix,
         export_df$sheetName);
   }

   # abbreviate assay_names
   short_assay_names <- "";
   if (use_assay_suffix &&
         length(unique(export_df$assay_names)) > 1) {
      for (n in 4:1) {
         short_assay_names <- export_df$assay_names
         ipattern <- paste0("([a-zA-Z]{1,",
            n,
            "})[a-zA-Z. ]*(_|$)")
         short_assay_names <- gsub(ipattern,
            "\\1",
            short_assay_names)
         test_nchar <- max(nchar(export_df$sheetName)) +
            max(nchar(short_assay_names)) + 4;
         if (test_nchar <= max_nchar_sheetname) {
            break;
         }
      }
      export_df$sheetName <- paste0(export_df$sheetName,
         "_",
         short_assay_names);
   }

   # quick function to ensure unique sheetName,
   # up to max_nchar_sheetname nchar length
   validate_sheetNames <- function(sheetNames, max_nchar_sheetname=31, trim_n=3, ...) {
      sheetNames <- gsub("^[-.:_ ]+|[-.:_ ]+$", "", sheetNames);
      is_toolong <- (nchar(sheetNames) > max_nchar_sheetname);
      if (any(is_toolong)) {
         sheetNames[is_toolong] <- substr(sheetNames, 1, max_nchar_sheetname)
      }
      sheetNames <- gsub("^[-.:_ ]+|[-.:_ ]+$", "", sheetNames);
      is_dupe <- duplicated(sheetNames);
      if (any(is_dupe)) {
         sheetNames[is_dupe] <- jamba::makeNames(
            substr(sheetNames[is_dupe], 1, max_nchar_sheetname - trim_n),
            renameOnes=TRUE,
            suffix="_v");
      }
      add_n <- 0;
      while (any(nchar(sheetNames) > max_nchar_sheetname)) {
         add_n <- add_n + 1;
         sheetNames[is_dupe] <- jamba::makeNames(
            substr(sheetNames[is_dupe], 1, max_nchar_sheetname - (trim_n + add_n)),
            renameOnes=TRUE,
            suffix="_v");
         if (add_n > 3) {
            break;
         }
      }
      if (any(duplicated(sheetNames))) {
         sheetNames <- validate_sheetNames(sheetNames,
            max_nchar_sheetname,
            trim_n + 1);
      }
      sheetNames
   }
   # ensure each sheet name is unique and no greater than max_nchar_sheetname
   export_df$sheetName <- validate_sheetNames(export_df$sheetName,
      max_nchar_sheetname);

   if (length(max_rows) > 0) {
      max_rows <- rep(max_rows, nrow(export_df));
      export_df$max_rows <- max_rows;
   }
   if (review_output) {
      return(export_df);
   }


   # Iterate each stat table and save to Excel worksheet
   append <- FALSE;
   wb <- NULL;
   for (irow in seq_len(nrow(export_df))) {
      assay_name <- export_df$assay_names[irow];
      cutoff_name <- export_df$cutoff_names[irow];
      contrast_name <- export_df$contrast_names[irow];
      sheetName <- export_df$sheetName[irow];

      if (verbose) {
         jamba::printDebug("save_sestats(): ",
            c("Saving sheet (", irow, "): "),
            sheetName,
            sep="");
      }
      iDF <- sestats$stats_dfs[[assay_name]][[contrast_name]];
      if (length(max_rows) > 0) {
         iDF <- head(iDF,
            max_rows[irow]);
      }
      if (length(iDF) == 0 || nrow(iDF) == 0) {
         if (verbose) {
            jamba::printDebug("save_sestats(): ",
               "      No data availale to save.");
         }
         next;
      }
      iDF <- jamba::renameColumn(iDF,
         from="probes",
         to=row_type);
      icols <- jamba::provigrep(c("probes|symbol|gene|protein",
         "^hit ",
         "."),
         colnames(iDF));
      iDF <- iDF[,icols, drop=FALSE];

      # optionally add rowData_colnames
      if (length(se) > 0 && length(rowData_colnames) > 0) {
         rowData_colnames <- intersect(rowData_colnames,
            colnames(SummarizedExperiment::rowData(se)));
         if (length(rowData_colnames) > 0) {
            imatch <- match(iDF[,1], rownames(se));
            if (any(is.na(imatch))) {
               jmatch <- match(rownames(iDF), rownames(se));
               if (!any(is.na(jmatch))) {
                  imatch <- jmatch;
               }
            }
            if (!any(is.na(imatch))) {
               jDF <- data.frame(check.names=FALSE,
                  SummarizedExperiment::rowData(se[imatch,])[, rowData_colnames, drop=FALSE]);
               iDF <- data.frame(check.names=FALSE,
                  iDF[, 1, drop=FALSE],
                  jDF,
                  iDF[, -1, drop=FALSE]);
            }
         }
      }

      iDF$assay_name <- assay_name;

      # detect column types
      highlightColumns <- jamba::igrep("symbol|gene|descr", colnames(iDF));
      hitColumns <- jamba::igrep("^hit ", colnames(iDF));
      pvalueColumns <- setdiff(
         jamba::igrep("p.val|pval|adjp", colnames(iDF)),
         hitColumns);
      fcColumns <- jamba::igrep("^fold ", colnames(iDF));
      lfcColumns <- jamba::igrep("^logFC ", colnames(iDF));
      intColumns <- jamba::igrep(" mean$|^mgm ", colnames(iDF));

      # adjust column widths
      smallColumns <- unique(c(pvalueColumns,
         fcColumns,
         lfcColumns,
         intColumns));
      save_widths <- rep(10, ncol(iDF));
      save_widths[length(save_widths)] <- 40;
      save_widths[hitColumns] <- 20;
      save_widths[smallColumns] <- 15;
      save_widths <- save_widths * width_factor;

      # save to Excel
      wb <- jamba::writeOpenxlsx(
         file=NULL,
         wb=wb,
         append=append,
         sheetName=sheetName,
         x=iDF,
         verbose=FALSE,
         autoWidth=FALSE,
         highlightColumns=highlightColumns,
         hitColumns=hitColumns,
         hitRule=c(-1, 0, 1),
         pvalueColumns=pvalueColumns,
         pvalueRule=c(1e-05, 0.05, 1),
         fcColumns=fcColumns,
         lfcColumns=lfcColumns,
         intColumns=intColumns,
         freezePaneColumn=2,
         colWidths=save_widths,
         colorSub=colorSub,
         ...);
      append <- TRUE;
   }
   openxlsx::saveWorkbook(wb=wb,
      file=file,
      overwrite=TRUE);
   return(invisible(export_df));
}



