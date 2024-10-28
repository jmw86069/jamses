
# jam_save_sestats.R

#' Save SE contrast stats output
#'
#' Save SE contrast stats output
#'
#' This function is intended as a convenient method to export
#' a series of statistical tables into organized, formatted Excel
#' worksheets.
#'
#' The output will generally contain two types of worksheets:
#'
#' 1. Each contrast in its own worksheet.
#' This is option is enabled by including `"contrasts"` in
#' argument `data_content`, which is default.
#'
#'    * If there are multiple "Signals" (e.g. multiple `assay_name`)
#'    then each contrast/signal combination will be saved to its
#'    own worksheet.
#'
#' 2. One table will be created with one column for each contrast,
#' using values `c(1, 0, -1)` to indicate whether the row met the
#' statistical criteria.
#' This is option is enabled by including `"hits"` in
#' argument `data_content`, which is default.
#'
#'    * If there are multiple "Signals" (e.g. multiple `assay_name`)
#'    then one table for each signal will be saved to its
#'    own worksheet.
#'
#' ## Output to Excel xlsx or as `list`
#'
#' * To output a `list` of `data.frame` objects, use:
#' `type="list", review_output=FALSE`
#'
#' ## Worksheet names
#'
#' Because the `xlsx` worksheet name is limited to 31 characters, the
#' contrast names are abbreviated by default. The output of this
#' function is a `data.frame` which includes the relationship of
#' the contrast and signal, to the abbreviated worksheet name.
#'
#' @returns `data.frame` when `type="xlsx"` or `list` when `type="list"`.
#'
#' @family jamses stats
#'
#' @param sestats `list` object output from `se_contrast_stats()`
#' @param file `character` string indicating the filename to save.
#'    When `file` is `NULL`, output is returned as a `list`, equivalent
#'    to `type="list"`.
#' @param assay_names `character` string indicating which assay names
#'    to save, stored in `dimnames(sestats$hit_array)$Signal`.
#'    When `NULL` then all assay names are saved.
#' @param contrast_names `character` string indicating which contrasts
#'    to save, stored in `dimnames(sestats$hit_array)$Contrasts`.
#'    The default `NULL` will save all contrasts.
#' @param type `character` string indicating the type of file to save.
#'    * `"xlsx"` - saves an Excel xlsx file using `jamba::writeOpenxlsx()`.
#'    Each worksheet is renamed so the string length does not exceed
#'    `max_nchar_sheetname`, whose default is 31 characters.
#'    * `"list"` - returns a `list` of `data.frame` objects, equivalent
#'    to the data to be stored in an output file.
#'    This option will **not** save data to `file`.
#' @param data_content `character` string describing the data content
#'    to include:
#'    * `"contrasts","hits"` - include worksheets per `contrast_names`,
#'    then assemble one `"hit sheet"` across all contrasts.
#'    One hit sheet is created for each value in `assay_names`.
#'    * `"contrasts"` - (default) include worksheets per `contrast_names`
#'    * `"hits"` - include only one `"hit sheet"` per value in
#'    `assay_names`.
#' @param hits_use_lfc `logical` default FALSE, indicating whether values
#'    in `"hits"` columns should use the log2 fold change.
#'    * `FALSE` (default) assigns `c(-1, 0, 1)` to indicate directionality
#'    after applying stat thresholds.
#'    * `TRUE` assigns the actual log2 fold change *only for hits* as defined
#'    by the stat thresholds.
#' @param max_nchar_sheetname `integer` number of characters allowed in
#'    Microsoft Excel worksheet names, default 31 characters.
#' @param abbreviate `logical` indicating whether to abbreviate factor
#'    levels using `shortest_unique_abbreviation()`.
#'    This option is `FALSE` by default, but may become preferred
#'    after more testing.
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
#' @param rename_contrasts `logical` indicating whetheer to apply
#'    `contrasts2comp()` to shorten long contrast names.
#' @param se `SummarizedExperiment`, default NULL, used when
#'    `rowData_colnames` is defined.
#' @param rowData_colnames `character`, default NULL, with optional colnames
#'    used only when `se` is also provided. When defined, it provides
#'    additional annotations for each row as defined by `rowData(se)`.
#' @param row_type `character` with custom column name to use for the
#'    primary row identifier. The default `"probes"` is often not accurate,
#'    though this may not be problematic in practice.
#'    When defined, the first column is renamed to `row_type`.
#' @param hitRule,hitFormat,freezePaneColumn arguments passed to
#'    `jamba::writeOpenxlsx()`, used only to define the color thresholds
#'    used with conditional formatting. It changes none of the data.
#'    The `freezePaneColumn` defines the first non-fixed column when
#'    viewed in Excel, and by default keeps only the first column fixed
#'    when scrolling to the right. Use a higher value if columns added
#'    by `rowData_colnames` should also be fixed columns.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to `jamba::writeOpenxlsx()`
#'
#' @examples
#' se <- make_se_test();
#'
#' # create SEDesign
#' sedesign <- groups_to_sedesign(se, group_colnames="group")
#'
#' # limma contrasts
#' sestats <- se_contrast_stats(se=se,
#'    sedesign=sedesign,
#'    assay_names="counts")
#'
#' # review_output=TRUE
#' info_df <- save_sestats(sestats, review_output=TRUE)
#' info_df
#'
#' # review_output=FALSE
#' stat_dfs1 <- save_sestats(sestats, review_output=FALSE, type="list")
#' head(stat_dfs1[[1]])
#'
#' # review_output=FALSE, hits_use_lfc=TRUE
#' stat_dfs <- save_sestats(sestats, review_output=FALSE, type="list", hits_use_lfc=TRUE)
#' head(stat_dfs[[1]])
#'
#' set.seed(12)
#' heatmap_se(se, sestats=sestats)
#'
#' set.seed(12)
#' heatmap_se(stat_dfs[[2]], column_names_rot=80,
#'    column_cex=0.2, row_cex=0.5) +
#' heatmap_se(se, sestats=sestats, rows=rownames(se))
#'
#' set.seed(12)
#' heatmap_se(stat_dfs[[2]], column_names_rot=80,
#'    column_cex=0.2, row_cex=0.5) +
#' heatmap_se(stat_dfs1[[2]], column_names_rot=80,
#'    column_cex=0.2, row_cex=0.5) +
#' heatmap_se(se, sestats=sestats, rows=rownames(se))
#'
#' @export
save_sestats <- function
(sestats,
 file=NULL,
 assay_names=NULL,
 contrast_names=NULL,
 cutoff_names=NULL,
 type=c("xlsx",
    "list"),
 data_content=c("data",
    "hits"),
 hits_use_lfc=FALSE,
 max_nchar_sheetname=31,
 abbreviate=FALSE,
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
 hitRule=c(-1, 0, 1),
 hitFormat="#,##0",
 freezePaneColumn=2,
 verbose=FALSE,
 ...)
{
   # validate arguments
   type <- match.arg(type);
   data_content <- match.arg(data_content,
      choices=c("data", "hits"),
      several.ok=TRUE);
   if (length(file) == 0) {
      type <- "list";
   }

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
   if (verbose) {
      jamba::printDebug("save_sestats(): ",
         "contrast_names:", contrast_names, sep="; ");
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
   export_df$saved <- ifelse("data" %in% data_content, "Yes", "No")
   if (verbose) {
      jamba::printDebug("save_sestats(): ",
         "export_df (input):");
      print(export_df);
   }

   # optionally rename contrasts
   if (rename_contrasts) {
      export_df$sheetName <- tryCatch({
         gsub("[:]+", ";",
            contrast2comp(export_df$contrast_names,
               abbreviate=abbreviate,
               ...))
      }, error=function(e){
         export_df$contrast_names
      });
   }
   # if (verbose) {
   #    jamba::printDebug("save_sestats(): ",
   #       "export_df (edited):");
   #    print(export_df);
   # }

   if (length(sheet_prefix) > 0) {
      if (verbose) {
         jamba::printDebug("save_sestats(): ",
            "applying sheet_prefix: ", sheet_prefix);
      }
      export_df$sheetName <- paste0(sheet_prefix,
         export_df$sheetName);
   }

   # abbreviate assay_names
   short_assay_names <- "";
   if (use_assay_suffix &&
         length(unique(export_df$assay_names)) > 1) {
      if (verbose) {
         jamba::printDebug("save_sestats(): ",
            "abbreviating assay_names.");
      }
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

   # optionally add hit sheets per assay_name
   if ("hits" %in% data_content) {
      hit_df <- data.frame(assay_names=assay_names,
         cutoff_names=NA,
         contrast_names="hits",
         sheetName=paste("hit", assay_names),
         saved="Yes")
      export_df <- rbind(export_df, hit_df)
   }

   # quick function to ensure unique sheetName,
   # up to max_nchar_sheetname nchar length
   validate_sheetNames <- function
   (sheetNames,
    max_nchar_sheetname=31,
    trim_n=3,
    ...)
   {
      sheetNames <- gsub("^[-.:_ ]+|[-.:_ ]+$", "", sheetNames);
      is_toolong <- (nchar(sheetNames) > max_nchar_sheetname);
      if (any(is_toolong)) {
         sheetNames[is_toolong] <- substr(sheetNames[is_toolong],
            1, max_nchar_sheetname);
      }
      # remove leading/trailing whitespace/punctuation
      sheetNames <- gsub("^[-.:_ ]+|[-.:_ ]+$", "", sheetNames);

      # look for duplicates
      is_dupe <- duplicated(sheetNames);
      if (any(is_dupe)) {
         sheetNames[is_dupe] <- jamba::makeNames(
            substr(sheetNames[is_dupe],
               1, max_nchar_sheetname - trim_n),
            renameOnes=TRUE,
            suffix="_v");
      }
      add_n <- 0;
      while (any(nchar(sheetNames) > max_nchar_sheetname)) {
         add_n <- add_n + 1;
         sheetNames[is_dupe] <- jamba::makeNames(
            substr(sheetNames[is_dupe],
               1, max_nchar_sheetname - (trim_n + add_n)),
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
   if (verbose) {
      jamba::printDebug("save_sestats(): ",
         "export_df (edited):");
      print(export_df);
   }
   # ensure each sheet name is unique and no greater than max_nchar_sheetname
   export_df$sheetName <- validate_sheetNames(export_df$sheetName,
      max_nchar_sheetname);
   if (verbose) {
      jamba::printDebug("save_sestats(): ",
         "export_df (validated):");
      print(export_df);
   }

   if (length(max_rows) > 0) {
      max_rows <- rep(max_rows, nrow(export_df));
      export_df$max_rows <- max_rows;
   }
   if (TRUE %in% review_output) {
      return(export_df);
   }


   # Iterate each stat table and save to Excel worksheet
   append <- FALSE;
   wb <- NULL;
   export_dfs <- list();
   hit_list_dfs <- list();
   for (irow in seq_len(nrow(export_df))) {
      assay_name <- export_df$assay_names[irow];
      cutoff_name <- export_df$cutoff_names[irow];
      contrast_name <- export_df$contrast_names[irow];
      sheetName <- export_df$sheetName[irow];

      if (verbose) {
         jamba::printDebug("save_sestats(): ",
            c("Preparing sheet (", irow, "): "),
            sheetName,
            sep="");
      }

      # assemble iDF to be saved to Excel
      name_colnames <- NULL;
      if ("hits" %in% export_df$contrast_names[irow]) {
         if (TRUE %in% verbose) {
            jamba::printDebug("save_sestats(): ",
               "Assembling hit sheet for assay_name '", assay_name, "'");
         }
         # assemble one hit data.frame
         name_colnames <- unique(unlist(
            lapply(hit_list_dfs[[assay_name]], function(xdf){
               attr(xdf, "nameColumns")
            })));
         # make first column a factor to help order rows consistently
         iDF <- jamba::mergeAllXY(
            lapply(hit_list_dfs[[assay_name]], function(xdf){
               xdf[, 1] <- factor(xdf[, 1], levels=unique(xdf[, 1]))
               xdf
            }));
         iDF <- jamba::mixedSortDF(iDF, byCols=1);
         # optionally keep only max_rows rows
         if (length(max_rows) > 0) {
            iDF <- head(iDF,
               max_rows[irow]);
         }
         # revert to character for compatibility with openxlsx
         iDF[, 1] <- as.character(iDF[, 1]);
         attr(iDF, "nameColumns") <- name_colnames;
      } else {
         iDF <- sestats$stats_dfs[[assay_name]][[contrast_name]];
         # optionally keep only max_rows rows
         if (length(max_rows) > 0) {
            iDF <- head(iDF,
               max_rows[irow]);
         }
         if (length(iDF) == 0 || nrow(iDF) == 0) {
            if (verbose) {
               jamba::printDebug("save_sestats(): ",
                  "      No data available to save.");
            }
            next;
         }
         iDF <- jamba::renameColumn(iDF,
            from="probes",
            to=row_type);
         name_grep <- unique(c(
            gsub("[()]", ".", row_type),
            "probes",
            "symbol",
            "gene",
            "protein"));
         name_grep <- jamba::cPaste(sep="|",
            name_grep[nchar(name_grep) > 0]);
         name_colnames <- jamba::provigrep(name_grep, colnames(iDF));
         icols <- jamba::provigrep(c(
            name_grep,
            "^hit ",
            "."),
            colnames(iDF));
         iDF <- iDF[,icols, drop=FALSE];

         # optionally convert hits to use log2 fold change
         if (TRUE %in% hits_use_lfc) {
            lfccol <- head(jamba::vigrep("^logFC|^log2f", icols), 1);
            if (length(lfccol) == 0) {
               # warn here and skip
            } else {
               if (verbose) {
                  jamba::printDebug("save_sestats(): ",
                     "Converted hitcols to use log2 fold change.");
               }
               hitcols <- jamba::vigrep("^hit ", icols);
               for (hitcol in hitcols) {
                  iDF[[hitcol]] <- (abs(sign(iDF[[hitcol]])) * iDF[[lfccol]])
               }
            }
         }

         # optionally add rowData_colnames
         if (length(se) > 0 && length(rowData_colnames) > 0) {
            rowData_colnames <- setdiff(intersect(rowData_colnames,
               colnames(SummarizedExperiment::rowData(se))),
               name_colnames);
            nonname_colnames <- setdiff(icols, name_colnames);
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
                     SummarizedExperiment::rowData(se[imatch,])[,
                        rowData_colnames, drop=FALSE]);
                  iDF <- data.frame(check.names=FALSE,
                     iDF[, name_colnames, drop=FALSE],
                     jDF,
                     iDF[, nonname_colnames, drop=FALSE]);
                  # add to name_colnames
                  name_colnames <- c(name_colnames,
                     rowData_colnames);
                  # name_colnames <- c(name_colnames[1],
                  #    colnames(jDF),
                  #    name_colnames[-1]);
               }
            }
         }
      }

      iDF$assay_name <- assay_name;

      # detect column types
      highlightColumns <- jamba::igrep("symbol|gene|descr", colnames(iDF));
      hitColumns <- jamba::igrep("^hit ", colnames(iDF));
      pvalueColumns <- setdiff(
         jamba::igrep("p.val|pval|adjp|adj.p", colnames(iDF)),
         hitColumns);
      fcColumns <- jamba::igrep("^fold ", colnames(iDF));
      lfcColumns <- jamba::igrep("^logFC ", colnames(iDF));
      numColumns <- jamba::igrep(" mean$|^mgm ", colnames(iDF));

      # adjust column widths
      smallColumns <- unique(c(pvalueColumns,
         fcColumns,
         lfcColumns,
         numColumns));
      save_widths <- rep(10, ncol(iDF));
      save_widths[length(save_widths)] <- 40;
      save_widths[hitColumns] <- 20;
      save_widths[smallColumns] <- 15;
      save_widths <- save_widths * width_factor;

      # Add attributes with Excel column types
      attr(iDF, "colWidths") <- save_widths;
      attr(iDF, "highlightColumns") <- highlightColumns;
      attr(iDF, "hitColumns") <- hitColumns;
      attr(iDF, "pvalueColumns") <- pvalueColumns;
      attr(iDF, "fcColumns") <- fcColumns;
      attr(iDF, "lfcColumns") <- lfcColumns;
      attr(iDF, "numColumns") <- numColumns;
      attr(iDF, "nameColumns") <- name_colnames;

      # optionally assemble hit sheets
      if ("hits" %in% data_content) {
         # add the hit data.frame to a list to assemble later
         # jamba::printDebug("head(iDF, 5):");print(head(iDF, 5));# debug
         # jamba::printDebug("name_colnames:");print(name_colnames);# debug
         # jamba::printDebug("hitColumns:");print(hitColumns);# debug
         hit_list_dfs[[assay_name]][[sheetName]] <- iDF[, c(name_colnames,
            colnames(iDF)[hitColumns]), drop=FALSE];
         attr(hit_list_dfs[[assay_name]][[sheetName]],
            "nameColumns") <- name_colnames;
      }
      # save to Excel
      if ("xlsx" %in% type && "Yes" %in% export_df[irow, "saved"]) {
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
            hitRule=hitRule,
            hitFormat=hitFormat,
            pvalueColumns=pvalueColumns,
            pvalueRule=c(1e-05, 0.05, 1),
            fcColumns=fcColumns,
            lfcColumns=lfcColumns,
            numColumns=numColumns,
            freezePaneColumn=freezePaneColumn,
            colWidths=save_widths,
            colorSub=colorSub,
            ...);
         append <- TRUE;
      }
      # add to list
      if ("list" %in% type && "Yes" %in% export_df[irow, "saved"]) {
         export_dfs[[sheetName]] <- iDF;
      }
   }
   if ("xlsx" %in% type) {
      openxlsx::saveWorkbook(wb=wb,
         file=file,
         overwrite=TRUE);
   }
   return(invisible(c(export_dfs,
      list(key=export_df))));
}



