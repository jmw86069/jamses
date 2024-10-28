
#' Extract stats as data.frame from SEStats results
#'
#' Extract stat as data.frame from SEStats results
#'
#' The purpose is to output `data.frame` or `list` of `data.frame`
#' from SEStats results. This function is essentially a wrapper
#' to `save_sestats()` with `type="list"`.
#'
#' @returns `list` of `data.frame` based upon `data_content`:
#'    * `data_content="contrasts"` will return `list` with one `data.frame`
#'    for each contrast, and each assay_name.
#'    Each list element is named by the contrast name.
#'    * `data_content="hits"` will return `list` with one `data.frame`
#'    for each `assay_name`, including one column for each contrast.
#'    Each list element is named by "hits" plus the assay_name, for
#'    example `"hits quantile_counts"`.
#'    * `data_content=c("contrasts", "hits")` will return `list` which
#'    includes both the options above.
#'
#' @family jamses stats
#'
#' @param sestats `list` object output from `se_contrast_stats()`
#' @param assay_names `character` string indicating which assay names
#'    to save, stored in `dimnames(sestats$hit_array)$Signal`.
#'    When `NULL` then all assay names are saved.
#' @param contrast_names `character` string indicating which contrasts
#'    to save, stored in `dimnames(sestats$hit_array)$Contrasts`.
#'    The default `NULL` will save all contrasts.
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
#' @param rename_contrasts `logical` indicating whetheer to apply
#'    `contrasts2comp()` to shorten long contrast names.
#'    Currently this option only applies to the list element names,
#'    not the column headers.
#' @param se `SummarizedExperiment`, default NULL, used when
#'    `rowData_colnames` is defined.
#' @param rowData_colnames `character`, default NULL, with optional colnames
#'    used only when `se` is also provided. When defined, it provides
#'    additional annotations for each row as defined by `rowData(se)`.
#' @param row_type `character` with custom column name to use for the
#'    primary row identifier. The default `"probes"` is often not accurate,
#'    though this may not be problematic in practice.
#'    When defined, the first column is renamed to `row_type`.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to `save_sestats()`.
#'
sestats_to_dfs <- function
(sestats,
 assay_names,
 contrast_names=NULL,
 data_content=c("contrasts",
    "hits"),
 hits_use_lfc=FALSE,
 rename_contrasts=TRUE,
 se=NULL,
 rowData_colnames=NULL,
 row_type="gene_name",
 verbose=FALSE,
 ...)
{
   #
   df_list <- save_sestats(sestats=sestats,
      assay_names=assay_names,
      contrast_namescontrast_names,
      data_content=data_content,
      type="list",
      hits_use_lfc=hits_use_lfc,
      rename_contrasts=rename_contrasts,
      se=se,
      rowData_colnames=rowData_colnames,
      row_type=row_type,
      verbose=verbose,
      ...)

   return(df_list);
}


#' Make SummarizedExperiment test data
#'
#' Make SummarizedExperiment test data
#'
#' @returns `SummarizedExperiment` object
#'
#' @family jamses SE utilities
#'
#' @param ngroups `integer` number of experimental groups
#' @param mreps `integer` number of replicates per group, can be used
#'    to provide the number of replicates for each group in order.
#' @param nrow `integer` number of rows (measurements)
#' @param multiplier `numeric` value multiplied by `rnorm()` to adjust
#'    the magnitude of values produced, default 1.
#' @param offset `numeric` value added to the output of `rnorm()`, default 7.
#' @param hit_fraction `numeric` value between 0 and 1 indicating the
#'    fraction of rows to simulate as having a fold change, default 1/2.
#' @param hit_max `numeric` maximum value for a simulated fold change,
#'    default 2.8, intended to be interpreted as log2 of a 7-fold change.
#' @param noise_factor `numeric` multiplied by `rnorm()` to add additional
#'    noise.
#' @param seed `numeric` passed to `set.seed()` when provided, for
#'    reproducible random output.
#' @param assay_name `character` name to use for the assay name in the output.
#' @param sparsity `numeric` value from 0 to 1, default 0, indicating the
#'    fraction of values that are converted to `NA` to simulate sparse
#'    data measurements. It can be provided as a vector and applied to
#'    each group in order.
#'    In some proteomics datasets, control samples may be substantially
#'    more sparse than other groups, for example if the control is
#'    non-targeted IP, or negative control. In this case, data can
#'    be simulated by using `sparsity=c(0.6, 0, 0)` for three groups.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' se <- make_se_test();
#' heatmap_se(se)
#' heatmap_se(se, controlSamples=1:3)
#'
#' se <- make_se_test(ngroups=3, assay_name="expression");
#' heatmap_se(se, apply_hm_column_title=TRUE)
#' heatmap_se(se, controlSamples=1:3, control_label="versus A", apply_hm_column_title=TRUE)
#'
#' se <- make_se_test(ngroups=3, assay_name="expression", nreps=c(3, 5, 4));
#' heatmap_se(se, controlSamples=1:3, control_label="versus A", apply_hm_column_title=TRUE)
#'
#' se <- make_se_test(ngroups=3, assay_name="expression", nreps=c(3, 5, 4), hit_max=4);
#' heatmap_se(se, controlSamples=1:3, control_label="versus A", apply_hm_column_title=TRUE)
#'
#' sedesign <- groups_to_sedesign(SummarizedExperiment::colData(se)[, "group", drop=FALSE])
#' plot_sedesign(sedesign)
#' sestats <- se_contrast_stats(se=se, sedesign=sedesign, assay_names="expression")
#' plot_sedesign(sedesign, sestats=sestats,
#'    contrast_style="none", sestats_style="label")
#' heatmap_se(se, controlSamples=1:3, control_label="versus A",
#'    column_split="group",
#'    apply_hm_column_title=TRUE, sestats=sestats, color_max=5)
#'
#' if (jamba::check_pkg_installed("venndir")) {
#'    venndir::venndir(hit_array_to_list(sestats), overlap_type="each")
#'    venndir::venndir(hit_array_to_list(sestats), overlap_type="each", proportional=TRUE)
#'    venndir::venndir(hit_array_to_list(sestats), overlap_type="each", show_labels="ncs")
#' }
#'
#' # demonstrate sparsity
#' se2 <- make_se_test(sparsity=c(0.5, 0));
#' hm2a <- heatmap_se(se2,
#'    column_title="global centered\nall values shown")
#' # when centering versus a sparse control group, some values can be lost:
#' hm2b <- heatmap_se(se2, controlSamples=1:3,
#'    column_title="centered vs A\nsome values become NA")
#' hm2a + hm2b
#' # use naControlFloor
#' hm2c <- heatmap_se(se2,
#'    column_title="centered vs A\nnaControlFloor=7",
#'    naControlFloor=7, naControlAction="floor", controlSamples=1:3)
#' hm2a + hm2b + hm2c
#'
#' @export
make_se_test <- function
(ngroups=2,
 nreps=3,
 nrow=50,
 multiplier=1,
 offset=7,
 hit_fraction=1/2,
 hit_max=2.8,
 noise_factor=1,
 seed=123,
 assay_name="counts",
 sparsity=0,
 verbose=FALSE,
 ...)
{
   #
   # define ncol
   if (length(ngroups) != 1 || ngroups < 1) {
      ngroups <- 1
   }
   nreps <- rep(nreps, length.out=ngroups)
   hit_fraction <- rep(hit_fraction, length.out=ngroups)
   if (any(hit_fraction > 1 | is.na(hit_fraction))) {
      hit_fraction[hit_fraction > 1] <- 1;
   }
   group_names <- rep(paste0("group", head(LETTERS, ngroups)), nreps)
   names(nreps) <- unique(group_names);
   names(hit_fraction) <- unique(group_names);
   sample_names <- jamba::makeNames(group_names);
   ncol <- length(sample_names);

   # set random seed if provided
   if (length(seed) == 1) {
      set.seed(seed)
   }

   # generate data
   n <- ncol * nrow;
   expr <- rnorm(nrow) * multiplier + offset;
   noise <- rnorm(n) * multiplier * 0.2 * noise_factor;
   m <- matrix(data=expr + noise,
      ncol=ncol);
   colnames(m) <- sample_names;
   rownames(m) <- paste0("row_", jamba::padInteger(seq_len(nrow)));

   # sparsity
   if (length(sparsity) >= 1 && any(sparsity > 0)) {
      sparsity <- rep(sparsity, length.out=ngroups);
      group_list <- split(sample_names, group_names);
      for (i in seq_along(sparsity)) {
         ivals <- m[, group_list[[i]], drop=FALSE];
         nvals <- length(ivals);
         nblank <- ceiling(nvals * sparsity[[i]]);
         toblank <- sample(seq_len(nvals), size=nblank)
         ivals[toblank] <- NA;
         m[, group_list[[i]]] <- ivals;
      }
   }

   # define hit values
   hit_vals <- c(-2.8, -1.3, -0.585, -0.485, -0.38,
      0.38, 0.485, 0.585, 1.3, 2.8) / 2.8 * hit_max;

   # add known fold changes to group 2 and higher
   if (ngroups > 1) {
      group_list <- tail(split(sample_names, group_names), -1);
      for (iname in names(group_list)) {
         i <- group_list[[iname]];
         irows <- ceiling(nrow * hit_fraction[iname])
         irowsk <- sample(rownames(m), size=irows);
         ilfc <- sample(hit_vals,
            size=irows,
            replace=TRUE)
         if (verbose) {
            jamba::printDebug("make_se_test(): ",
               "iname:", iname,
               ", irows:", irows,
               ", ilfc:", ilfc);
         }
         # fold <- rnorm(length(i) * irows) * multiplier / 2.5;
         m[irowsk, i] <- m[irowsk, i] + ilfc;
      }
   }

   # simulate some "hits"
   if (FALSE) {
      m[1, 4:6] <- m[1, 4:6] + 2
      m[2, 4:6] <- m[2, 4:6] + 1.5
      m[3, 4:6] <- m[3, 4:6] - 1.3
   }

   # create SummarizedExperiment
   se <- SummarizedExperiment::SummarizedExperiment(
      assays=list(counts=m),
      colData=data.frame(sample=colnames(m),
         row.names=colnames(m),
         group=factor(group_names,
            levels=unique(group_names))),
      rowData=data.frame(measurement=rownames(m),
         row.names=rownames(m)))
   SummarizedExperiment::assayNames(se) <- head(assay_name, 1);

   se;
}

