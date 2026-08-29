
#' Simulate SummarizedExperiment tests
#'
#' Simulate SummarizedExperiment tests, with adjustable batch effect.
#' 
#' This function extends `make_se_test()` in two ways:
#' 
#' 1. For four-group scenario, it will generate two-factor changes.
#' 2. Batch effects are generated, adjustable with `batch_ex`.
#' Use `batch_ex=0` for no batch effect. See examples for
#' demonstration of centering within batch as visual cue that
#' batch effect is well-defined.
#'
#' @returns `SummarizedExperiment` object
#'
#' @family jamses SE utilities
#'
#' @param ngroups `integer` number of experimental groups
#' @param nreps `integer` number of replicates per group, can be used
#'    to provide the number of replicates for each group in order.
#' @param nrow `integer` number of rows (measurements)
#' @param multiplier `numeric` value multiplied by `rnorm()` to adjust
#'    the magnitude of values produced, default 1.
#' @param offset `numeric` value added to the output of `rnorm()`, default 7.
#' @param batch_ex `numeric` adjustment for batch effect, default `0.15`
#'    is a moderate batch effect, it would require adjustment or
#'    blocking or covariate for effective analysis.
#' @param hit_exp `numeric` default 2, adjust the "sharpness" of
#'    differential changes, higher values are sharper, with fewer
#'    strong hits.
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
#' # batch effects
#' seb <- simulate_se_test();
#' hmb <- heatmap_se(seb,
#'    controlSamples=colnames(seb)[c(1:5, 11:15)],
#'    column_title="With batch effect\nglobal-centered")
#' hmb
#' 
#' # when centering versus a sparse control group, some values can be lost:
#' hmbg <- heatmap_se(seb,
#'    controlSamples=colnames(seb)[c(1:5, 11:15)],
#'    centerby_colnames="batch",
#'    column_title="With batch effect\ncentered by batch")
#' hmb + hmbg
#' 
#' hmbc <- heatmap_se(seb,
#'    correlation=TRUE,
#'    controlSamples=colnames(seb)[c(1:5, 11:15)],
#'    column_title="With batch effect\nglobal-centered")
#' hmbc
#' 
#' hmbgc <- heatmap_se(seb,
#'    correlation=TRUE,
#'    controlSamples=colnames(seb)[c(1:5, 11:15)],
#'    centerby_colnames="batch",
#'    column_title="With batch effect\ncentered by batch")
#' hmbgc
#' hmbc + hmbgc
#'
#' @export
simulate_se_test <- function(
   ngroups = 4,
   nreps = 5,
   nrow = 250,
   multiplier = 1.5,
   offset = 7,
   batch_ex = 0.15,
   hit_exp = 2,
   hit_fraction = 1 / 2,
   hit_max = 2.8,
   noise_factor = 1,
   seed = 123,
   assay_name = "counts",
   sparsity = 0,
   verbose = FALSE,
   ...
) {
   #
   # define ncol
   if (length(ngroups) != 1 || ngroups < 1) {
      ngroups <- 1
   }
   nreps <- rep_len(nreps, ngroups)
   hit_fraction <- rep_len(hit_fraction, ngroups)
   if (any(hit_fraction > 1 | is.na(hit_fraction))) {
      hit_fraction[hit_fraction > 1] <- 1
   }
   if (ngroups == 4) {
      #
      fac1 <- rep(c("Ctl", "Dex"), each=nreps * 2)
      fac2 <- rep(rep(c("WT", "KO"), each=nreps), 2)
      group_names <- paste0(fac1, "_", fac2)
   } else {
      group_names <- rep(paste0("group", head(LETTERS, ngroups)), nreps)
   }
   names(nreps) <- unique(group_names)
   names(hit_fraction) <- unique(group_names)
   sample_names <- jamba::makeNames(group_names)
   ncol <- length(sample_names)

   # set random seed if provided
   if (length(seed) == 1) {
      set.seed(seed)
   }

   # generate data
   n <- ncol * nrow
   expr <- rnorm(nrow) * multiplier + offset
   noise <- rnorm(n) * multiplier * 0.1 * noise_factor
   m <- matrix(data = expr + noise, ncol = ncol)
   colnames(m) <- sample_names
   rownames(m) <- paste0("row_", jamba::padInteger(seq_len(nrow)))

   # sparsity
   if (length(sparsity) >= 1 && any(sparsity > 0)) {
      sparsity <- rep_len(sparsity, ngroups)
      group_list <- split(sample_names, group_names)
      for (i in seq_along(sparsity)) {
         ivals <- m[, group_list[[i]], drop = FALSE]
         nvals <- length(ivals)
         nblank <- ceiling(nvals * sparsity[[i]])
         toblank <- sample(seq_len(nvals), size = nblank)
         ivals[toblank] <- NA
         m[, group_list[[i]]] <- ivals
      }
   }

   # define batch effect values
   repnames <- paste0("v", seq_len(max(nreps)));
   if (batch_ex > 0) {
      # column-shared signal
      batch_offset <- rnorm(max(nreps)) * multiplier * batch_ex / 1.2;
      # row-shared signal
      batch_noise <- matrix(
         nrow=nrow(m),
         ncol=max(nreps),
         data=rnorm(max(nreps) * nrow(m)) * multiplier * batch_ex * 1.2);
      names(batch_offset) <- repnames;
      batch_m <- matrix(
         nrow=nrow(m),
         ncol=ncol(m),
         data=rep(batch_offset, each=nrow(m)) + batch_noise
      )
      colnames(batch_m) <- colnames(m);
      rownames(batch_m) <- rownames(m);
      m <- m + batch_m;
   }

   # define hit values
   # simulate log2fold
   lfc1 <- rnorm(nrow(m) * ngroups - 1) / 2;
   lfc <- sign(lfc1) * (abs(lfc1) ^ hit_exp);
   lfc_group <- matrix(
      ncol=ngroups,
      nrow=nrow(m),
      data=c(
         rep(0, nrow(m)),
         lfc
      )
   )
   colnames(lfc_group) <- unique(group_names);

   # Add lfc to each group
   if (ngroups == 4) {
      # two-factor style
      usegroup <- unique(group_names)[2];
      for (groupnum in unique(group_names)[c(2, 4)]) {
         for (repnum in seq_len(nreps)) {
            use_sample <- paste0(groupnum, "_v", repnum);
            m[, use_sample] <- (m[, use_sample] + lfc_group[, usegroup]);
         }
      }
      usegroup <- unique(group_names)[3];
      for (groupnum in unique(group_names)[c(3, 4)]) {
         for (repnum in seq_len(nreps)) {
            use_sample <- paste0(groupnum, "_v", repnum);
            m[, use_sample] <- (m[, use_sample] + lfc_group[, usegroup]);
         }
      }
   } else {
      # independent groups
      for (groupnum in tail(unique(group_names) -1)) {
         for (repnum in seq_len(nreps)) {
            use_sample <- paste0(groupnum, "_v", repnum);
            m[, use_sample] <- (m[, use_sample] + lfc_group[, groupnum]);
         }
      }
   }

   # create SummarizedExperiment
   se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = m),
      colData = data.frame(
         sample = colnames(m),
         row.names = colnames(m),
         group = factor(group_names, levels = unique(group_names))
      ),
      rowData = data.frame(measurement = rownames(m), row.names = rownames(m))
   )
   if (ngroups == 4) {
      SummarizedExperiment::colData(se)$Trt <- gsub(
         "^.+_", "", SummarizedExperiment::colData(se)$group);
      SummarizedExperiment::colData(se)$Geno <- gsub(
         "_.+$", "", SummarizedExperiment::colData(se)$group);
      SummarizedExperiment::colData(se)$group <- NULL;
   }
   SummarizedExperiment::colData(se)$batch <- gsub("^.+_", "", colnames(m))

   SummarizedExperiment::assayNames(se) <- head(assay_name, 1)

   se
}
