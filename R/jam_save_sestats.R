
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
 type=c("xlsx"),
 max_nchar_sheetname=31,
 review_output=TRUE,
 sheet_prefix=NULL,
 width_factor=1,
 max_rows=NULL,
 colorSub=NULL,
 verbose=FALSE,
 ...)
{
   #
   type <- match.arg(type);

   # determine which results to save
   if (length(assay_names) == 0) {
      assay_names <- dimnames(sestats$hit_array)$Signal;
   } else {
      assay_names <- intersect(assay_names,
         dimnames(sestats$hit_array)$Signal);
   }
   # if (length(cutoff_names) == 0) {
   #    cutoff_names <- dimnames(sestats$hit_array)$Cutoffs;
   # } else {
   #    cutoff_names <- intersect(cutoff_names,
   #       dimnames(sestats$hit_array)$Cutoffs);
   # }
   cutoff_names <- "spacer";
   if (length(contrast_names) == 0) {
      contrast_names <- dimnames(sestats$hit_array)$Contrasts;
   } else {
      contrast_names <- intersect(contrast_names,
         dimnames(sestats$hit_array)$Contrasts);
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
   if (length(sheet_prefix) > 0) {
      export_df$sheetName <- paste0(sheet_prefix,
         export_df$sheetName);
   }
   if (any(nchar(export_df$sheetName) > max_nchar_sheetname)) {
      export_df$sheetName <- jamba::makeNames(
         substr(export_df$sheetName,
            1, max_nchar_sheetname - 3));
      if (any(nchar(export_df$sheetName) > max_nchar_sheetname)) {
         export_df$sheetName <- jamba::makeNames(
            substr(export_df$sheetName,
               1, max_nchar_sheetname - 4));
      }
   }
   if (length(max_rows) > 0) {
      max_rows <- rep(max_rows, nrow(export_df));
      export_df$max_rows <- max_rows;
   }
   if (review_output) {
      return(export_df);
   }


   # Iterate each stat table and save to Excel worksheet
   append <- FALSE;
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
         to="gene_name");
      icols <- jamba::provigrep(c("probes|symbol|gene|protein",
         "^hit ",
         "."),
         colnames(iDF));
      iDF <- iDF[,icols, drop=FALSE];
      iDF$assay_name <- assay_name;

      # detect column types
      highlighColumns <- jamba::igrep("symbol|gene|descr", colnames(iDF));
      hitColumns <- jamba::igrep("^hit ", colnames(iDF));
      pvalueColumns <- jamba::igrep("p.val", colnames(iDF));
      fcColumns <- jamba::igrep("^fold ", colnames(iDF));
      lfcColumns <- jamba::igrep("^logFC ", colnames(iDF));
      intColumns <- jamba::igrep(" mean$|^mgm ", colnames(iDF));

      # save to Excel
      jamba::writeOpenxlsx(
         file=file,
         append=append,
         sheetName=sheetName,
         x=iDF,
         verbose=FALSE,
         autoWidth=FALSE,
         highlighColumns=highlighColumns,
         hitColumns=hitColumns,
         hitRule=c(-1, 0, 1),
         pvalueColumns=pvalueColumns,
         pvalueRule=c(1e-05, 0.05, 1),
         fcColumns=fcColumns,
         lfcColumns=lfcColumns,
         intColumns=intColumns,
         freezePaneColumn=2,
         colorSub=colorSub,
         ...);
      append <- TRUE;

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
      if (verbose) {
         jamba::printDebug("save_sestats(): ",
            "      Adjusting colwidths: ", save_widths);
      }
      jamba::set_xlsx_colwidths(xlsxFile=file,
         sheet=sheetName,
         widths=save_widths);
   }
   return(invisible(export_df));
}
