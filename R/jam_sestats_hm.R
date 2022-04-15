
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
#' allows supplying an alternative hit array.
#'
#' Argument `sample_color_list` is intended to be the output
#' from `platjam::design2colors()` which will very likely be moved
#' into this package. It is a named list of named color vectors.
#' The list names should contain all `top_colnames`, and each vector
#' name should correspond to a value in the relevant column of
#' `colData(se)`. If no colors are defined, `ComplexHeatmap::Heatmap()`
#' will determine its own colors.
#'
#' @export
heatmap_se <- function
(se,
   sestats,
   rows=NULL,
   assay_name=NULL,
   contrast_names=NULL,
   cutoff_name=1,
   alt_sestats=NULL,
   alt_assay_name=assay_name,
   alt_contrast_names=NULL,
   alt_cutoff_name=1,
   isamples=colnames(se),
   normgroup_colname="Type",
   centerby_colnames=c(normgroup_colname,
      "Run"),
   controlSamples=NULL,
   top_colnames=NULL,
   top_annotation=NULL,
   sample_color_list=NULL,
   row_cex=1,
   useMedian=FALSE,
   show_row_names=TRUE,
   row_label_colname=NULL,
   cluster_columns=FALSE,
   color_max=3,
   lens=2,
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

   #
   if (length(contrast_names) == 0) {
      contrast_names <- dimnames(sestats$hit_array)[[2]];
   }

   # define rows to use
   if (length(sestats) > 0) {
      if ("list" %in% class(sestats) && "hit_array" %in% names(sestats)) {
         hit_array <- sestats$hit_array;
      } else {
         hit_array <- sestats;
      }
      gene_hitlist <- hit_array[cutoff, contrast_names, assay_name];
      gene_hits <- names(jamba::tcount(names(unlist(unname(
         gene_hitlist)))));
      gene_hits_im <- venndir::list2im_value(gene_hitlist,
         do_sparse=FALSE)[gene_hits,,drop=FALSE];
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
      gene_hitlist_alt <- alt_hit_array[alt_cutoff_name, alt_contrast_names, alt_assay_name];
      gene_hits_alt <- names(tcount(names(unlist(unname(
         gene_hitlist_alt)))));
      gene_hits_im_alt1 <- venndir::list2im_value(gene_hitlist_alt,
         do_sparse=FALSE)[gene_hits_alt,,drop=FALSE];
      gene_hits_im_alt <- (gene_hits_im * 0)[,rep(1, ncol(gene_hits_im_alt1)), drop=FALSE];
      colnames(gene_hits_im_alt) <- colnames(gene_hits_im_alt1);
      genes_shared <- intersect(gene_hits_alt,
         gene_hits);
      gene_hits_im_alt[genes_shared,] <- gene_hits_im_alt1[genes_shared,];
   }

   # normgroup for column split
   normgroup_colname <- intersect(normgroup_colname,
      colnames(colData(se)));
   if (length(normgroup_colname) > 0 &&
         length(unique(colData(se[,isamples])[[normgroup_colname]])) > 0) {
      column_split <- colData(se[,isamples])[[normgroup_colname]]
   } else {
      column_split <- NULL;
   }

   # determine annotations atop samples
   if (length(top_colnames) == 0) {
      top_colnames <- names(which(
         sapply(colnames(colData(se)), function(i){
            any(duplicated(colData(se)[[i]]))
         })))
   }
   if (length(top_annotation) == 0) {
      top_annotation <- ComplexHeatmap::HeatmapAnnotation(
         border=TRUE,
         df=data.frame(check.names=FALSE,
            colData(se[,isamples])[,top_colnames, drop=FALSE]),
         annotation_legend_param=list(
            border=TRUE
         ),
         col=sample_color_list);
   }

   # left_annotation
   if (length(sestats) > 0) {
      if (length(alt_sestats) > 0) {
         left_annotation <- ComplexHeatmap::rowAnnotation(
            border=TRUE,
            hits_alt=gene_hits_im_alt[gene_hits,,drop=FALSE],
            hits=gene_hits_im[gene_hits,,drop=FALSE],
            col=list(
               hits_alt=colorjam::col_div_xf(1.5),
               hits=colorjam::col_div_xf(1.5)),
            annotation_legend_param=list(
               hits=list(
                  at=c(-1, 0, 1),
                  color_bar="discrete",
                  border=TRUE,
                  labels=c("down", "no change", "up")),
               hits_alt=FALSE)
         )
      } else {
         left_annotation <- ComplexHeatmap::rowAnnotation(
            border=TRUE,
            hits=gene_hits_im[gene_hits, , drop=FALSE],
            col=list(
               hits=colorjam::col_div_xf(1.5)),
            annotation_legend_param=list(
               hits=list(
                  at=c(-1, 0, 1),
                  color_bar="discrete",
                  border=TRUE,
                  labels=c("down", "no change", "up")))
         )
      }
   } else {
      left_annotation <- NULL;
   }

   norm_label <- paste0(assay_name, " data");

   centerby_colnames <- intersect(centerby_colnames,
      colnames(colData(se)));
   if (length(centerby_colnames) > 0) {
      centerby_label <- paste0("centered by ",
         jamba::cPaste(centerby_colnames,
            sep="/"));
      centerGroups <- jamba::pasteByRow(
         colData(se[,isamples])[,centerby_colnames]);
   } else {
      centerGroups <- NULL;
      centerby_label <- "global-centered";
   }

   # row font size
   row_fontsize <- jamba::noiseFloor(
      row_cex * (60*(14 / 10))/(length(gene_hits))^(1/2),
      minimum=1,
      ceiling=18);

   # row_labels
   if (length(row_label_colname) == 0) {
      row_labels <- gene_hits;
   } else {
      row_labels <- rowData(se[gene_hits,])[[row_label_colname]]
   }

   # heatmap legend labels
   legend_at <- seq(-ceiling(color_max), to=ceiling(color_max));
   legend_labels <- round(jamba::exp2signed(legend_at+0.001, offset=0))

   # define heatmap
   #hm_hits <- ComplexHeatmap::Heatmap(
   hm_hits <- multienrichjam::call_fn_ellipsis(ComplexHeatmap::Heatmap,
      matrix=jamma::centerGeneData(
         useMedian=useMedian,
         centerGroups=centerGroups,
         x=assays(se[gene_hits, isamples])[[assay_name]],
         controlSamples=controlSamples,
         ...),
      use_raster=TRUE,
      top_annotation=top_annotation,
      left_annotation=left_annotation,
      heatmap_legend_param=list(
         border=TRUE,
         color_bar="discrete",
         at=legend_at,
         labels=legend_labels
      ),
      clustering_method_rows="ward.D2",
      column_split=column_split,
      cluster_column_slices=FALSE,
      border=TRUE,
      name="centered\nexpression",
      show_row_names=show_row_names,
      row_labels=row_labels,
      row_names_gp=grid::gpar(fontsize=row_fontsize),
      col=col_div_xf(color_max,
         lens=lens),
      cluster_columns=cluster_columns,
      ...)
   hm_title <- paste0(
      formatInt(length(gene_hits)),
      " rows",
      "\n", norm_label, ", ",
      centerby_label)
   attr(hm_hits, "hm_title") <- hm_title;
   # hm_hits <- draw(hm_hits,
   #    column_title=paste0("Gene counts of ",
   #       formatInt(length(gene_hits)),
   #       " total RNA-seq DGEs",
   #       norm_label,
   #       "\ncentered by ", centerby_label))
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
