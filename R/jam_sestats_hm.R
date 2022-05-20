
#' Expression heatmap of sestats hits
#'
#' Expression heatmap of sestats hits
#'
#' Note: Still a work in progress.
#'
#' This function is a bold attempt to automate the intricate task
#' of creating an expression heatmap, by default using the `sestats`
#' output from `se_contrast_stats()`, specifically the `hit_array`.
#'
#' The intent is to display expression values from `assays(se)`,
#' centered across all columns, or with customization defined by
#' `centerby_colnames` and `normgroup_colnames`.
#'
#' The top heatmap annotations use `colData(se)` with user-supplied
#' `top_colnames` or by auto-detecting those colnames that apply
#' to multiple `colnames(se)`.
#'
#' The `hit_array` data is used to define an incidence matrix of up/down
#' hits, which is displayed to the left of the heatmap. The contrasts
#' can optionally be subset with `contrast_names`.
#'
#' The heatmap title is returned as an `attr(hm, "hm_title")` that
#' describes the assay data used, data centering, and total rows
#' displayed.
#'
#' For comparison across other `sestats` results, argument `alt_sestats`
#' allows supplying an alternative hit array. These hit arrays are placed
#' as `left_annotation`, alongside optional data defined by `rowData_colnames`.
#'
#' When `rowData_colnames` is supplied, data in the corresponding colnames
#' of `rowData(se)` are also displayed in `left_annotation`. Colors can
#' be defined in `sample_color_list`.
#'
#' Argument `sample_color_list` is intended to be the output
#' from `platjam::design2colors()` which will very likely be moved
#' into this package. It is a named list of named color vectors.
#' The list names should contain all `top_colnames`, and each vector
#' name should correspond to a value in the relevant column of
#' `colData(se)`. If no colors are defined, `ComplexHeatmap::Heatmap()`
#' will determine its own colors.
#'
#' A custom `left_annotation` can be supplied, but this method currently
#' prevents the other annotations described above from being displayed.
#' Currently the best way to supply custom row annotations in addition
#' to those described above, supply `right_annotation` to be displayed
#' on the right side of the heatmap.
#'
#' Data provided in `se` is expected to be `SummarizedExperiment`, however
#' it also accepts other Bioconductor data types that provide
#' accessor functions `featureData()`, `phenoData()`, and `assayData()`,
#' including for example `"MethyLumiSet"` class.
#'
#' @param se `SummarizedExperiment` object with accessor functions:
#'    `rowData()`, `colData()`, and `assays()`;
#'    or another suitable Bioconductor object with accessor functions:
#'    `featureData()`, `phenoData()`, and `assayData()`.
#' @param sestats `list` output from `se_contrast_stats()`, which
#'    specifically contains `hit_array` as a 3-dimensional array of hits,
#'    with dimensions "Cutoffs", "Contrasts", "Signal". When `sestats`
#'    is supplied, then all genes with a non-zero entry in `hit_array`,
#'    for the corresponding `contrast_names`, will be included in the
#'    heatmap. When `rows` is also supplied, then the intersection
#'    with `rows` is displayed.
#' @param rows `character` vector of `rownames(se)` to define a specific
#'    set of rows to display. When `sestats` is supplied, then the
#'    intersection of `rows` with genes defined by `sestats` is displayed.
#' @param row_type `character` string used in the title of the heatmap
#'    which indicates how many rows are displayed. For example
#'    `"1,234 genes detected above background"` or
#'    `"1,234 DEGs by limma-voom"`.
#' @param correlation `logical` indicating whether to calculate sample
#'    correlation, and plot a sample-by-sample correlation heatmap.
#'    This option is included here since many of the same arguments
#'    are required for data centering, and sample annotations.
#'    Note that `color_max` is forced to a maximum value of `1.0`,
#'    representing the maximum correlation value.
#' @param assay_name `character` string indicating the name in
#'    `assays(se)` to use for data to be displayed in the heatmap.
#' @param contrast_names `character` vector of contrasts in
#'    `sestats$hit_array` to use for the heatmap. When `contrast_names=NULL`
#'    then all contrasts are displayed, which is the default.
#' @param contrast_suffix `character` string with optional suffix to append
#'    to the end of each contrast name label for `sestats` hit incidence
#'    matrix beside the heatmap. This suffix may be useful when comparing
#'    two methods for the same set of contrast names, with `sestats` and
#'    `alt_sestats`.
#' @param cutoff_name `character` or `integer` index used to define the
#'    specific statistical cutoffs to use from `sestats$hit_array`.
#' @param alt_sestats,alt_assay_name,alt_contrast_names,alt_contrast_suffix
#'    arguments analogous to those described above for `sestats` which
#'    are used when `alt_sestats` is supplied.
#' @param isamples `character` vector of `colnames(se)` used to provide a
#'    specific subset, or specific order of columns displayed in the heatmap.
#' @param normgroup_colname `character` vector of colnames in `colData(se)`
#'    used during data centering. When supplied, samples are centered
#'    independently within each normgroup grouping.
#' @param centerby_colnames `character` vector of colnames in `colData(se)`
#'    used during data centering. When supplied, samples are centered
#'    independently within each centerby grouping. It is typically used
#'    for things like cell lines, to center each cell line by a time
#'    point control, or untreated control.
#' @param controlSamples `character` vector of samples to use as the
#'    reference during data centering. Note that samples are still
#'    centered within each normgroup and centerby grouping, and within
#'    that grouping samples are centered to the `controlSamples`
#'    which are present in that grouping. In absence of `controlSamples`
#'    defined, or within the grouping, samples are centered relative
#'    to the median value of the grouping.
#' @param control_name `character` string used to describe the control
#'    used during data centering, displayed in the heatmap title.
#' @param top_colnames `character` vector of colnames to use from
#'    `colData(se)` as annotations to display in `top_annotation` above
#'    the heatmap. When not supplied, reasonable colnames are detected
#'    internally:
#'    * columns with more than one unique value
#'    * columns with at least one duplicated value
#' @param top_annotation specific heatmap annotation as defined by
#'    `ComplexHeatmap::HeatmapAnnotation()`. When supplied, the `top_colnames`
#'    described above is not used.
#' @param rowData_colnames `character` vector of colnames in `rowData(se)`
#'    to use for heatmap annotations displayed on the left side of
#'    the heatmap. Specific colors can be included in `sample_color_list`
#'    as a named `list` of color vectors or color functions. The names
#'    of this list must match colnames to be displayed, otherwise
#'    `ComplexHeatmap::Heatmap()` will define its own color function.
#' @param left_annotation specific heatmap annotation as defined by
#'    `ComplexHeatmap::rowAnnotation()`. When supplied, the `rowData_colnames`
#'    and `sestats` row annotations are not displayed. In order to supply
#'    custom row annotations and not lose `left_annotation` defined above,
#'    supply the row annotations as `right_annotation`.
#' @param alt_sestats `list` output from `se_contrast_stats()`, which
#'    specifically contains `hit_array` as a 3-dimensional array of hits,
#'    with dimensions "Cutoffs", "Contrasts", "Signal". This data is
#'    only used when `sestats` is also supplied.
#'    Note that the rows displayed is defined by `rows` and `sestats` above,
#'    and is not defined here.
#' @param `row_split` is used to define heatmap split by row, ultimately
#'    passed to `ComplexHeatmap::Heatmap()` argument `row_split`. However,
#'    the input type can vary:
#'    * `integer` number of row splits based upon row clustering
#'    * `character` value or values in colnames of `rowData(se)` to split
#'    using row annotation in `se`.
#'    * `character` or `factor` vector named by `rownames(se)` with another
#'    custom row split, passed directly to `ComplexHeatmap::Heatmap()`
#'    argument `row_split`, with proper order for rows being displayed.
#' @param row_title_rot `numeric` value indicating text rotation in degrees
#'    to use for row titles.
#' @param sample_color_list named `list` of color vectors or color functions,
#'    where names correspond to colnames in either `colData(se)` or
#'    `rowData(se)`, and which are passed to corresponding left or top
#'    annotation functions. When colors are not defined,
#'    `ComplexHeatmap::Heatmap()` will define colors using its own internal
#'    function.
#' @param row_cex,column_cex `numeric` values used to adjust the row and
#'    column name font size, relative to the automatic adjustment that
#'    is already done based upon the number of rows and columns being
#'    displayed.
#' @param useMedian `logical` passed to `jamma::centerGeneData()` during
#'    data centering.
#' @param show_row_names,show_row_dend `logical` indicating whether to
#'    display row names, and row dendrogram, respectively. With more than
#'    2,000 rows this step can become somewhat slow.
#' @param row_label_colname `character` string used as a row label, where
#'    this value is a colname in `rowData(se)`. It is useful when rownames
#'    are some identifier that is not user-friendly, and where another column
#'    in the data may provide a more helpful label, for example `"SYMBOL"`
#'    to display gene symbol instead of accession number.
#' @param cluster_colnames `logical` indicating whether to cluster columns
#'    by hierarchical clustering.
#' @param column_split `character` or `integer` vector used to define
#'    heatmap column split.
#' @param color_max `numeric` value passed to `colorjam::col_div_xf()`
#'    which defines the upper limit of color gradient used in the heatmap.
#' @param lens `numeric` value passed to `colorjam::col_div_xf()` to control
#'    the intensity of color gradient applied to the numeric range.
#' @param rename_contrasts,rename_alt_contrasts `logical` indicating
#'    whether to rename long contrast names in `sestats` and `alt_sestats`
#'    using `contrast2comp()`.
#' @param debug `logical` indicating debug mode, data is returned in a `list`:
#'    * `hm` object `ComplexHeatmap::Heatmap`
#'    * `top_annotation` object `ComplexHeatmap::HeatmapAnnotation` for columns
#'    * `left_annotation` object `ComplexHeatmap::HeatmapAnnotation` for rows
#'    * `hm_title` object `character` string with the heatmap title.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to supporting functions.
#'
#' @export
heatmap_se <- function
(se,
   sestats=NULL,
   rows=NULL,
   row_type="rows",
   correlation=FALSE,
   assay_name=NULL,
   contrast_names=NULL,
   contrast_suffix="",
   cutoff_name=1,
   alt_sestats=NULL,
   alt_assay_name=assay_name,
   alt_contrast_names=NULL,
   alt_contrast_suffix="",
   alt_cutoff_name=1,
   isamples=colnames(se),
   normgroup_colname=NULL,
   centerby_colnames=NULL,
   controlSamples=NULL,
   control_label="",
   top_colnames=NULL,
   top_annotation=NULL,
   rowData_colnames=NULL,
   left_annotation=NULL,
   row_split=NULL,
   row_title_rot=0,
   sample_color_list=NULL,
   row_cex=0.8,
   column_cex=1,
   useMedian=FALSE,
   show_row_names=length(rows) < 2000,
   show_row_dend=length(rows) < 2000,
   row_label_colname=NULL,
   cluster_columns=FALSE,
   cluster_rows=function(x, ...){
      amap::hcluster(rmNA(naValue=0, x),
         ...,
         method="euclidean",
         link="ward")},
   column_split=NULL,
   color_max=3,
   lens=2,
   rename_contrasts=TRUE,
   rename_alt_contrasts=TRUE,
   verbose=FALSE,
   debug=FALSE,
   ...)
{
   #
   if (!jamba::check_pkg_installed("ComplexHeatmap")) {
      stop("This function requires Bioconductor package ComplexHeatmap.");
   }
   if (!jamba::check_pkg_installed("venndir")) {
      stop("This function requires Github package venndir from 'jmw86069/venndir'");
   }
   if (!suppressPackageStartupMessages(require(SummarizedExperiment))) {
      stop("This function requires Bioconductor package SummarizedExperiment.");
   }
   if (length(correlation) == 0) {
      correlation <- FALSE;
   }

   # define rows to use
   gene_hitlist <- NULL;
   alt_gene_hitlist <- NULL;
   gene_hits_im <- NULL;
   alt_gene_hits_im <- NULL;
   if (length(sestats) > 0) {
      if ("list" %in% class(sestats) && "hit_array" %in% names(sestats)) {
         hit_array <- sestats$hit_array;
      } else {
         hit_array <- sestats;
      }
      if (length(contrast_names) == 0) {
         contrast_names <- dimnames(sestats$hit_array)[[2]];
      }
      # contrast_names1 <- contrast_names;
      # if (length(contrast_names) == 1) {
      #    contrast_names1 <- rep(contrast_names, 2);
      # }
      gene_hitlist <- hit_array_to_list(hit_array,
         cutoff_names=cutoff_name,
         contrast_names=contrast_names,
         assay_names=assay_name);
      gene_hits <- names(jamba::tcount(names(unlist(unname(
         gene_hitlist)))));
      gene_hits_im <- venndir::list2im_value(gene_hitlist,
         do_sparse=FALSE)[gene_hits,,drop=FALSE];

      # optionally rename contrasts
      if (rename_contrasts) {
         colnames(gene_hits_im) <- tryCatch({
            contrast2comp(colnames(gene_hits_im));
         }, error=function(e){
            colnames(gene_hits_im)
         });
      }

      if (length(contrast_suffix) > 0 && any(nchar(contrast_suffix)) > 0) {
         colnames(gene_hits_im) <- paste0(colnames(gene_hits_im),
            contrast_suffix);
      }
   }

   # rows is user-defined
   rows <- intersect(rows, rownames(se));
   if (length(rows) > 0) {
      if (length(sestats) > 0) {
         rows_im <- (gene_hits_im * 0)[rep(1, length(rows)), , drop=FALSE];
         rownames(rows_im) <- rows;
         gene_hits_rows <- intersect(rows, gene_hits);
         if (length(gene_hits_rows) > 0) {
            rows_im[match(gene_hits_rows, rows),] <-
               gene_hits_im[match(gene_hits_rows, rownames(gene_hits_im)),,drop=FALSE];
         }
         gene_hits_im <- rows_im;
      } else {
         gene_hits_im <- NULL;
      }
      gene_hits <- rows;
      rows <- gene_hits;
   } else {
      rows <- rownames(se);
   }

   # alt gene_hitlist
   if (length(sestats) > 0 && length(alt_sestats) > 0) {
      if ("list" %in% class(alt_sestats) && "hit_array" %in% names(alt_sestats)) {
         alt_hit_array <- alt_sestats$hit_array;
      } else {
         alt_hit_array <- alt_sestats;
      }
      if (length(alt_contrast_names) == 0) {
         alt_contrast_names <- dimnames(alt_hit_array)[[2]];
      }
      # alt_contrast_names1 <- alt_contrast_names;
      # if (length(alt_contrast_names) == 1) {
      #    alt_contrast_names1 <- rep(alt_contrast_names, 2);
      # }
      gene_hitlist_alt <- hit_array_to_list(hit_array,
         cutoff_names=alt_cutoff_name,
         contrast_names=alt_contrast_names,
         assay_names=alt_assay_name);
      gene_hits_alt <- names(tcount(names(unlist(unname(
         gene_hitlist_alt)))));
      gene_hits_im_alt1 <- venndir::list2im_value(gene_hitlist_alt,
         do_sparse=FALSE)[gene_hits_alt,,drop=FALSE];
      gene_hits_im_alt <- (gene_hits_im * 0)[,rep(1, ncol(gene_hits_im_alt1)), drop=FALSE];
      colnames(gene_hits_im_alt) <- colnames(gene_hits_im_alt1);
      genes_shared <- intersect(gene_hits_alt,
         gene_hits);
      if (length(genes_shared) > 0) {
         gene_hits_im_alt[genes_shared,] <- gene_hits_im_alt1[genes_shared,];
      }

      # optionally rename contrasts
      if (rename_contrasts) {
         colnames(gene_hits_im_alt) <- tryCatch({
            contrast2comp(colnames(gene_hits_im_alt));
         }, error=function(e){
            colnames(gene_hits_im_alt)
         })
      }

      if (length(alt_contrast_suffix) > 0 && any(nchar(alt_contrast_suffix)) > 0) {
         colnames(gene_hits_im_alt) <- paste0(colnames(gene_hits_im_alt),
            alt_contrast_suffix);
      }
   }

   # validate sample_color_list
   if (length(sample_color_list) > 0) {
      sample_color_list <- lapply(sample_color_list, function(i){
         if (is.function(i)) {
            i
         } else {
            jamba::rmNA(i)
         }
      })
   }

   # pull colData and rowData as data.frame
   # to be tolerant of other data types
   if (grepl("SummarizedExperiment", ignore.case=TRUE, class(se))) {
      rowData_se <- data.frame(check.names=FALSE,
         rowData(se));
      colData_se <- data.frame(check.names=FALSE,
         colData(se[,isamples]))
   } else {
      if (verbose) {
         jamba::printDebug("heatmap_se(): ",
            "using accessor functions: ",
            c("featureData()", "phenoData()"))
      }
      rowData_se <- as(featureData(se[gene_hits,]), "data.frame");
      rownames(rowData_se) <- gene_hits;
      colData_se <- as(phenoData(se[,isamples]), "data.frame")
      rownames(colData_se) <- isamples;
   }

   # normgroup for column split
   normgroup_colname <- intersect(normgroup_colname,
      colnames(colData_se));
   if (length(column_split) == 0) {
      if (length(normgroup_colname) > 0 &&
            length(unique(colData_se[[normgroup_colname]])) > 0) {
         column_split <- colData_se[[normgroup_colname]];
      } else {
         column_split <- NULL;
      }
   } else {
      if (any(c("factor", "character") %in% class(column_split))) {
         if (all(column_split %in% colnames(colData_se))) {
            column_split <- jamba::pasteByRowOrdered(
               data.frame(check.names=FALSE,
                  colData_se[, column_split, drop=FALSE]),
               keepOrder=TRUE);
         } else if (all(names(column_split) %in% isamples)) {
            column_split <- column_split[isamples];
         } else if (length(column_split) == length(isamples)) {
            # leave as-is
         } else {
            column_split <- NULL;
         }
      } else if (length(column_split) == 1 && is.numeric(column_split)) {
         # leave as-is
      } else {
         column_split <- NULL;
      }
   }

   # column font size
   column_fontsize <- jamba::noiseFloor(
      column_cex * 60/(length(isamples))^(1/2),
      ceiling=20,
      minimum=2);

   # row font size
   if (correlation) {
      row_fontsize <- jamba::noiseFloor(
         row_cex * (60*(14 / 10))/(length(isamples))^(1/2),
         minimum=1,
         ceiling=20);
   } else {
      row_fontsize <- jamba::noiseFloor(
         row_cex * (60*(14 / 10))/(length(gene_hits))^(1/2),
         minimum=1,
         ceiling=20);
   }

   # determine annotations atop samples
   if (length(top_colnames) == 0) {
      top_colnames <- names(which(
         sapply(colnames(colData_se), function(i){
            any(duplicated(colData_se[[i]])) &&
               length(unique(colData_se[[i]])) > 1
         })));
      if (length(top_colnames) > 0 && verbose) {
         jamba::printDebug("heatmap_se(): ",
            "derived top_colnames: ",
            top_colnames);
      }
   }
   if (length(top_annotation) == 0 && length(top_colnames) > 0) {
      top_annotation <- ComplexHeatmap::HeatmapAnnotation(
         border=TRUE,
         df=data.frame(check.names=FALSE,
            colData_se[,top_colnames, drop=FALSE]),
         annotation_legend_param=list(
            border=TRUE,
            color_bar="discrete"
         ),
         col=sample_color_list);
   }

   # left_annotation
   if (length(left_annotation) == 0 && !correlation) {
      column_anno_fontsize <- jamba::noiseFloor(
         column_cex * 12,
         ceiling=24,
         minimum=2);
      left_anno_list <- list();
      left_color_list <- list();
      left_param_list <- list();
      show_left_legend <- logical(0);
      # sestats annotations
      if (length(sestats) > 0) {
         show_left_legend <- c(TRUE,
            show_left_legend);
         left_anno_list <- c(list(
            hits=gene_hits_im[gene_hits, , drop=FALSE]),
            left_anno_list);
         left_color_list <- c(list(
            hits=colorjam::col_div_xf(1.5)),
            left_color_list);
         left_param_list <- c(list(
            hits=list(
               at=c(-1, 0, 1),
               color_bar="discrete",
               border=TRUE,
               labels=c("down", "no change", "up"))),
            left_param_list);
      }
      # alt_sestats annotations
      if (length(alt_sestats) > 0) {
         show_left_legend <- c(FALSE,
            show_left_legend);
         left_anno_list <- c(list(
            hits_alt=gene_hits_im_alt[gene_hits, , drop=FALSE]),
            left_anno_list);
         left_color_list <- c(list(
            hits_alt=colorjam::col_div_xf(1.5)),
            left_color_list);
         left_param_list <- c(list(
            hits_alt=list(
               at=c(-1, 0, 1),
               color_bar="discrete",
               border=TRUE,
               labels=c("down", "no change", "up"))),
            left_param_list);
      }
      # rowData annotations
      if (length(rowData_colnames) > 0) {
         if (verbose) {
            jamba::printDebug("heatmap_se(): ",
               "rowData_colnames: ",
               rowData_colnames);
         }
         show_left_legend <- c(TRUE,
            show_left_legend);
         left_anno_list <- c(list(
            df=data.frame(check.names=FALSE,
               rowData_se[gene_hits, ][,rowData_colnames, drop=FALSE])),
            left_anno_list);
         use_color_list_names <- intersect(rowData_colnames,
            names(sample_color_list));
         left_color_list <- c(
            jamba::rmNULL(
               sample_color_list[use_color_list_names]),
            left_color_list);
         left_param_list <- c(
            lapply(jamba::nameVector(rowData_colnames), function(iname){
               if (iname %in% names(sample_color_list)) {
                  if (is.function(sample_color_list[[iname]])) {
                     list(border=TRUE,
                        color_bar="discrete",
                        at=attr(sample_color_list[[iname]], "breaks"))
                  } else {
                     list(border=TRUE,
                        color_bar="discrete",
                        at=jamba::rmNA(names(sample_color_list[[iname]])))
                  }
               } else {
                  list(border=TRUE)
               }
            }),
            left_param_list);
      }

      # put it all together
      if (length(left_anno_list) > 0) {
         left_alist <- alist(
            col=left_color_list,
            annotation_legend_param=left_param_list,
            annotation_name_gp=grid::gpar(fontsize=column_anno_fontsize),
            #show_legend=show_left_legend,
            border=TRUE);
         if (debug > 1) {
            jamba::printDebug("heatmap_se(): ",
               "left_alist:");
            print(sdim(left_alist));
            print(left_alist);
            jamba::printDebug("heatmap_se(): ",
               "left_anno_list:");
            print(sdim(left_anno_list));
            print(left_anno_list);
         }
         left_arglist <- c(
            left_alist,
            left_anno_list);
         left_annotation <- do.call(ComplexHeatmap::rowAnnotation,
            left_arglist);
      }
   }

   # optional row_split
   if (length(row_split) > 0) {
      if (correlation) {
         # correlation uses colData for split
         if (any(c("factor", "character") %in% class(row_split))) {
            if (all(row_split %in% colnames(colData_se))) {
               row_split <- data.frame(check.names=FALSE,
                  colData_se[gene_hits, row_split, drop=FALSE]);
            } else if (all(names(row_split) %in% isamples)) {
               row_split <- row_split[isamples];
            } else {
               row_split <- NULL;
            }
         } else if (length(row_split) == 1 && is.numeric(row_split)) {
            # leave as-is
         } else {
            row_split <- NULL;
         }
      } else {
         # non-correlation uses rowData for split
         if (any(c("factor", "character") %in% class(row_split))) {
            if (all(row_split %in% colnames(rowData_se))) {
               row_split <- data.frame(check.names=FALSE,
                  rowData_se[gene_hits, row_split, drop=FALSE]);
            } else if (all(names(row_split) %in% gene_hits)) {
               row_split <- row_split[gene_hits];
            } else {
               row_split <- NULL;
            }
         } else if (length(row_split) == 1 && is.numeric(row_split)) {
            # leave as-is
         } else {
            row_split <- NULL;
         }
      }
   } else if (correlation) {
      row_split <- column_split;
   }

   norm_label <- paste0(assay_name, " data");

   centerby_colnames <- intersect(centerby_colnames,
      colnames(colData_se));
   if (length(centerby_colnames) > 0) {
      centerby_label <- paste0("centered by ",
         jamba::cPaste(centerby_colnames,
            sep="/"));
      centerGroups <- jamba::pasteByRow(
         colData_se[,centerby_colnames]);
   } else {
      centerGroups <- NULL;
      centerby_label <- "global-centered";
   }
   if (length(control_label) > 0 && any(nchar(control_label)) > 0) {
      centerby_label <- paste(centerby_label,
         control_label);
   }

   # row_labels
   if (correlation) {
      if (length(row_label_colname) == 0) {
         row_labels <- isamples;
      } else {
         row_labels <- colData_se[isamples, , drop=FALSE][[row_label_colname]];
      }
   } else {
      if (length(row_label_colname) == 0) {
         row_labels <- gene_hits;
      } else {
         row_labels <- rowData_se[gene_hits, , drop=FALSE][[row_label_colname]];
      }
   }

   # heatmap legend labels
   if (correlation) {
      if (abs(color_max) > 1) {
         color_max <- 1;
      }
      legend_at <- seq(-ceiling(color_max),
         to=ceiling(color_max),
         by=0.25);
      legend_labels <- legend_at;
   } else {
      legend_at <- seq(-ceiling(color_max),
         to=ceiling(color_max));
      legend_labels <- round(jamba::exp2signed(legend_at+0.001, offset=0))
   }

   # pull assay data separately so we can tolerate other object types
   if (grepl("SummarizedExperiment", ignore.case=TRUE, class(se))) {
      se_matrix <- assays(se[gene_hits, isamples])[[assay_name]];
   } else {
      se_matrix <- assayData(se[gene_hits, isamples])[[assay_name]];
   }

   # cluster_columns
   if (cluster_columns %in% TRUE) {
      cluster_columns <- function(x, ...) {
         amap::hcluster(rmNA(naValue=0, x),
            ...,
            method="euclidean",
            link="ward")}
   }

   # define heatmap matrix
   se_matrix <- jamma::centerGeneData(
      useMedian=useMedian,
      centerGroups=centerGroups,
      x=se_matrix,
      controlSamples=controlSamples,
      ...);
   hm_name <- "centered\nexpression";
   if (correlation) {
      # call correlation function cor()
      se_matrix <- multienrichjam::call_fn_ellipsis(cor,
         x=se_matrix,
         use="pairwise.complete.obs",
         ...);
      cluster_rows <- cluster_columns;
      hm_name <- "centered\ncorrelation";
   }

   # pre-calculate row clusters
   # This step is required to enable row_split as integer number of clusters,
   # which is not accepted when supplying a function.
   # This step does not work with character or data.frame row_split
   if (length(row_split) == 1 &&
         is.numeric(row_split) &&
         is.function(cluster_rows)) {
      cluster_rows <- cluster_rows(se_matrix);
   }
   if (length(column_split) == 1 &&
         is.numeric(column_split) &&
         is.function(cluster_columns)) {
      cluster_columns <- cluster_columns(se_matrix);
   }

   # define heatmap
   hm_hits <- multienrichjam::call_fn_ellipsis(ComplexHeatmap::Heatmap,
      matrix=se_matrix,
      use_raster=TRUE,
      top_annotation=top_annotation,
      left_annotation=left_annotation,
      heatmap_legend_param=list(
         border=TRUE,
         color_bar="discrete",
         at=legend_at,
         labels=legend_labels
      ),
      clustering_method_rows="ward.D",
      column_split=column_split,
      row_split=row_split,
      row_title_rot=row_title_rot,
      cluster_column_slices=FALSE,
      border=TRUE,
      name="centered\nexpression",
      show_row_names=show_row_names,
      show_row_dend=show_row_dend,
      row_labels=row_labels,
      row_names_gp=grid::gpar(fontsize=row_fontsize),
      column_names_gp=grid::gpar(fontsize=column_fontsize),
      col=colorjam::col_div_xf(color_max,
         lens=lens,
         ...),
      cluster_columns=cluster_columns,
      cluster_rows=cluster_rows,
      ...)
   hm_title <- paste0(
      formatInt(length(gene_hits)),
      " ", row_type,
      "\n", norm_label, ", ",
      "\n", centerby_label)
   attr(hm_hits, "hm_title") <- hm_title;
   # hm_hits <- draw(hm_hits,
   #    column_title=paste0("Gene counts of ",
   #       formatInt(length(gene_hits)),
   #       " total RNA-seq DGEs",
   #       norm_label,
   #       "\ncentered by ", centerby_label))
   if (debug) {
      ret_list <- list(
         hm=hm_hits,
         top_annotation=top_annotation,
         left_annotation=left_annotation,
         hm_title=hm_title
      );
      ret_list$gene_hits_im <- gene_hits_im;
      ret_list$alt_gene_hits_im <- alt_gene_hits_im;
      ret_list$gene_hitlist <- gene_hitlist;
      ret_list$alt_gene_hitlist <- alt_gene_hitlist;
      ret_list$hit_array <- hit_array;
      return(ret_list)
   }
   hm_hits
}


#' Call function using safe ellipsis arguments
#'
#' Call function using safe ellipsis arguments
#'
#' This function is a wrapper function intended to help
#' pass ellipsis arguments `...` from a parent function
#' to an external function in a safe way.
#' It will only include arguments from `...` that are
#' recognized by the external function.
#'
#' When the external function FUN arguments `formals()`
#' includes ellipsis `...`, then the `...` will be passed
#' as-is without change.
#'
#' When the external function FUN arguments `formals()`
#' does not include ellipsis `...`, then only named
#' arguments in `...` that are recognized by FUN
#' will be passed, as defined by `names(formals(FUN))`.
#'
#' Note that arguments must be named.
#'
#' @param FUN `function` that should be called with arguments
#'    in `...`
#' @param ... arguments are passed to `FUN()` in safe manner.
#'
#' @examples
#' new_mean <- function(x, trim=0, na.rm=FALSE) {
#'    mean(x, trim=trim, na.rm=na.rm)
#' }
#' x <- c(1, 3, 5, NA);
#' new_mean(x, na.rm=TRUE);
#' tryCatch({
#'    new_mean(x, na.rm=TRUE, color="red");
#' }, error=function(e){
#'    print(e);
#' })
#'
#' call_fn_ellipsis(new_mean, x=x, na.rm=TRUE, color="red")
#' call_fn_ellipsis(new_mean, x=x, color="red")
#'
#' @export
call_fn_ellipsis <- function
(FUN,
 ...)
{
   FUN_argnames <- names(formals(FUN));
   if ("..." %in% FUN_argnames) {
      FUN(...)
   } else {
      arglist <- list(...)
      argkeep <- which(names(arglist) %in% FUN_argnames);
      arguse <- arglist[argkeep]
      do.call(FUN, arguse)
   }
}
