

#' Compute contrast statistics on SummarizedExperiment data
#'
#' Compute contrast statistics on SummarizedExperiment data
#'
#' This function is essentially a wrapper around statistical methods
#' in the `limma` package, with additional steps to apply statistical
#' thresholds to define "statistical hits" by three main criteria:
#'
#' * P-value or adjusted P-value
#' * fold change
#' * max group mean
#'
#' This function is unique in that it applies the statistical methods
#' to one or more "signals" in the input `SummarizedExperiment` assays,
#' specifically intended to compare things like normalization methods.
#'
#' If multiple statistical thresholds are defined, each one is applied
#' in order, which is specifically designed to compare the effect of
#' applying different statistical thresholds. For example one may want
#' to pre-compute "statistical hits" using adjusted P-value 0.05, and 0.01;
#' or using fold change >= 1.5, or fold change >= 2.0. The underlying
#' statistics are the same, but a column indicating hits is created
#' for each set of thresholds.
#'
#' Hits are annotated:
#'
#'  * `-1` for down-regulation
#'  * `0` for un-changed (by the given criteria)
#'  * `1` for up-regulation
#'
#'  The results are therefore intended to feed directional Venn
#'  diagrams, which display the overlapping gene hits, and whether the
#'  directions are shared or opposed.
#'
#' This function can optionally apply the limma-voom workflow,
#' which involves calculating matrix weights using `limma::voom()`,
#' then applying those weights during the model fit.
#'
#' The output is intended to include several convenient formats:
#'
#' * `stats_dfs` - list of `data.frame` stats one per contrast
#' * `stats_df` - one `data.frame` with all stats together
#' * `hit_array` - `array` with three dimensions: signal; contrast; threshold
#' whose cells contain hit flags (-1, 0, 1) named by `rownames(se)`.
#'
#' Design and contrast matrices can be defined using the function
#' `jamses::groups_to_sedesign()`. That function assigns each sample
#' to a sample group, then assembles all relevant group contrasts
#' which involve only one-factor contrast at a time. It optionally
#' defines two-factor contrasts (contrast of contrasts) where
#' applicable.
#'
#' A subset of genes (`rownames(se)`) or samples (`colnames(se)`) can
#' be defined, to restrict calculations to use only the subset data.
#'
#' @param se `SummarizedExperiment` object
#' @param assay_names `character` vector with one or more assay names
#'    from `names(assays(se))`.
#' @param normgroup `character` or `factor` vector with length
#'    `length(isamples)`, whose values define independent normalization
#'    groups, intended so that samples in `isamples` are analyzed only
#'    with samples that have the same `normgroup`. During limma model
#'    fit, all samples in all groups are used by default, however this
#'    technique allows subsets of samples to be analyzed independently.
#'    When `normgroup=NULL` the default is to assume all samples are in
#'    the same `normgroup="bulk"`.
#'    Each subset of samples begins with the same `sestats`, `idesign`,
#'    `icontrast`, however they are fed into `validate_sestats()` to
#'    subset the appropriate contrasts to use only samples within the
#'    current normgroup. As a result, any contrasts that span two
#'    normgroups will be removed, and will not appear in the output.
#' @param do_normgroups `logical` whether to enable normgroup processing,
#'    or to use the previous technique that kept all samples together.
#'    Note that when `normgroup=NULL` the output should be identical
#'    with `do_normgroups=TRUE` and `do_normgroups=FALSE`.
#' @param ... additional arguments are ignored.
#'
#' @family jamses stats
#'
#' @export
se_contrast_stats <- function
(se,
 assay_names,
 adjp_cutoff=0.05,
 p_cutoff=NULL,
 fold_cutoff=1.5,
 int_adjp_cutoff=adjp_cutoff,
 int_p_cutoff=p_cutoff,
 int_fold_cutoff=fold_cutoff,
 mgm_cutoff=NULL,
 ave_cutoff=NULL,
 confint=FALSE,
 floor_min=NULL,
 floor_value=NULL,
 sedesign=NULL,
 icontrasts=NULL,
 idesign=NULL,
 igenes=NULL,
 isamples=NULL,
 enforce_design=TRUE,
 use_voom=FALSE,
 posthoc_test=c("none",
    "DEqMS"),
 posthoc_args=list(DEqMS=list(
    PSM_counts=NULL,
    fit.method="loess")),
 weights=NULL,
 robust=FALSE,
 handle_na=c("partial", "full", "none"),
 collapse_by_gene=FALSE,
 block=NULL,
 correlation=NULL,
 normgroup=NULL,
 do_normgroups=TRUE,
 verbose=FALSE,
 ...)
{
   handle_na <- match.arg(handle_na);
   posthoc_test <- match.arg(posthoc_test);

   ## Validate input parameters
   if (length(isamples) == 0) {
      isamples <- colnames(se);
   }
   isamples <- intersect(isamples, colnames(se));
   if (length(igenes) == 0) {
      igenes <- rownames(se);
   }
   igenes <- intersect(igenes, rownames(se));

   if (length(sedesign) > 0 && "SEDesign" %in% class(sedesign)) {
      if (length(idesign) > 0 || length(icontrasts) > 0) {
         warning(paste0("Note when supplying sedesign,",
            "idesign and icontrasts are ignored."))
      }
      if (length(isamples) > 0) {
         sedesign <- validate_sedesign(sedesign,
            samples=isamples,
            verbose=verbose);
      }
      isamples <- sedesign@samples;
      idesign <- sedesign@design;
      icontrasts <- sedesign@contrasts;
   } else {
      if (length(idesign) == 0) {
         stop("idesign must be defined.");
      }
      isamples <- intersect(isamples, rownames(idesign));
      idesign <- idesign[match(isamples, rownames(idesign)),,drop=FALSE];
      if (length(rownames(idesign)) == 0) {
         stop("rownames(idesign) must contain values matching colnames(se).");
      }
      icontrasts <- icontrasts[match(rownames(icontrasts), colnames(idesign)),,drop=FALSE];
      if (length(rownames(icontrasts)) == 0) {
         stop("rownames(icontrasts) must match values in colnames(idesign).");
      }
   }

   ## normgroup implementation
   if (length(normgroup) == 0) {
      normgroup <- jamba::nameVector(rep("bulk", length(isamples)),
         isamples);
   }
   if (length(names(normgroup)) == 0 && length(normgroup) == length(isamples)) {
      names(normgroup) <- isamples;
   }
   # isamples_normgroup_list is a list separated by normgroup,
   # each vector will be analyzed independently
   isamples_normgroup_list <- split(isamples, normgroup[isamples]);

   ## Iterate each assay_name
   ## Run statistical tests for gene data
   stats_hits_dfs1 <- lapply(jamba::nameVector(assay_names), function(signalSet) {
      retVals <- list();
      imatrix <- SummarizedExperiment::assays(se[igenes, isamples])[[signalSet]];

      # iterate each normgroup independently
      if (TRUE %in% do_normgroups) {
         normgroup_stats <- lapply(jamba::nameVectorN(isamples_normgroup_list), function(normgroup_name) {
            normgroup_samples <- isamples_normgroup_list[[normgroup_name]];
            if (verbose) {
               jamba::printDebug("se_contrast_stats(): ",
                  "   Analyzing normgroup_name: ", normgroup_name);
            }
            # new imatrix_ng only uses samples in this normgroup
            imatrix_ng <- imatrix[,normgroup_samples, drop=FALSE];
            sestats_ng <- validate_sedesign(
               new("SEDesign",
                  design=idesign,
                  contrasts=icontrasts),
               samples=normgroup_samples);
            icontrasts_ng <- sestats_ng@contrasts;
            idesign_ng <- sestats_ng@design;
            if (length(icontrasts_ng) == 0 || ncol(icontrasts_ng) == 0 || nrow(icontrasts_ng) == 0) {
               if (verbose) {
                  jamba::printDebug("se_contrast_stats(): ",
                     "   No contrasts were defined for this normgroup_name.", fgText=c("darkorange2", "red"));
               }
               return(NULL)
            }

            if (!"none" %in% handle_na && any(is.na(imatrix_ng))) {
               imatrix_ng <- handle_na_values(imatrix_ng,
                  idesign=idesign_ng,
                  handle_na=handle_na);
            }
            ## Optionally determine voom weights prior to running limma
            if (use_voom) {
               if (verbose) {
                  jamba::printDebug("se_contrast_stats(): ",
                     "   Determining Voom weight matrix (within probe reps).");
               }
               imatrix_v <- voom_jam((2^imatrix_ng)-1,
                  design=idesign_ng,
                  normalize.method="none",
                  plot=FALSE,
                  verbose=verbose,
                  ...);
               weights <- imatrix_v$weights;
               rownames(weights) <- rownames(imatrix_ng);
               colnames(weights) <- colnames(imatrix_ng);
               if (verbose) {
                  jamba::printDebug("se_contrast_stats(): ",
                     "   Determined Voom weight matrix.");
               }
            }

            #######################################################
            ## Optionally convert zero (or less than zero) to NA
            if (length(floor_min) == 1 && !is.na(floor_min) && any(imatrix_ng <= floor_min)) {
               if (verbose) {
                  jamba::printDebug("se_contrast_stats(): ",
                     c("Applying floor_min:",
                        floor_min,
                        ", replacing with floor_value:",
                        floor_value),
                     sep="");
               }
               imatrix_ng[imatrix_ng <= floor_min] <- floor_value;
            }

            ## Run limma
            rlr_result_ng <- run_limma_replicate(imatrix=imatrix_ng,
               idesign=idesign_ng,
               icontrasts=icontrasts_ng,
               weights=weights,
               robust=robust,
               verbose=verbose,
               adjp_cutoff=adjp_cutoff,
               p_cutoff=p_cutoff,
               fold_cutoff=fold_cutoff,
               mgm_cutoff=mgm_cutoff,
               int_adjp_cutoff=int_adjp_cutoff,
               int_p_cutoff=int_p_cutoff,
               int_fold_cutoff=int_fold_cutoff,
               ave_cutoff=ave_cutoff,
               collapse_by_gene=collapse_by_gene,
               block=block,
               correlation=correlation,
               posthoc_test=posthoc_test,
               posthoc_args=posthoc_args,
               ...);
            return(rlr_result_ng);
         });
         # end normgroup processing
         # now assemble normgroup_stats into rlr_result as before
         rlr_result_stats_dfs <- unlist(recursive=FALSE,
            lapply(names(normgroup_stats), function(rlr_name){
               rlr <- normgroup_stats[[rlr_name]];
               rlr$stats_dfs;
            }));
         rlr_result_stats_df_colnames <- unique(unlist(lapply(rlr_result_stats_dfs, colnames)));
         rlr_result_stats_df <- jamba::mergeAllXY(rlr_result_stats_dfs)
         # rlr_result_stats_df <- rlr_result_stats_df[,rlr_result_stats_df_colnames, drop=FALSE];
         rlr_result <- list(
            stats_df=rlr_result_stats_df,
            stats_dfs=rlr_result_stats_dfs,
            rep_fits=lapply(normgroup_stats, function(rlr){rlr$rep_fits})
         )
      } else {
         if (!"none" %in% handle_na && any(is.na(imatrix))) {
            imatrix <- handle_na_values(imatrix,
               idesign=idesign,
               handle_na=handle_na);
         }
         ## Optionally determine voom weights prior to running limma
         if (use_voom) {
            if (verbose) {
               jamba::printDebug("se_contrast_stats(): ",
                  "   Determining Voom weight matrix (within probe reps).");
            }
            imatrix_v <- voom_jam((2^imatrix)-1,
               design=idesign,
               normalize.method="none",
               plot=FALSE,
               verbose=verbose,
               ...);
            weights <- imatrix_v$weights;
            rownames(weights) <- rownames(imatrix);
            colnames(weights) <- colnames(imatrix);
            if (verbose) {
               jamba::printDebug("se_contrast_stats(): ",
                  "   Determined Voom weight matrix.");
            }
         }

         #######################################################
         ## Optionally convert zero (or less than zero) to NA
         if (length(floor_min) == 1 && !is.na(floor_min) && any(imatrix <= floor_min)) {
            if (verbose) {
               jamba::printDebug("se_contrast_stats(): ",
                  c("Applying floor_min:",
                     floor_min,
                     ", replacing with floor_value:",
                     floor_value),
                  sep="");
            }
            imatrix[imatrix <= floor_min] <- floor_value;
         }

         ## Run limma
         rlr_result <- run_limma_replicate(imatrix=imatrix,
            idesign=idesign,
            icontrasts=icontrasts,
            weights=weights,
            robust=robust,
            verbose=verbose,
            adjp_cutoff=adjp_cutoff,
            p_cutoff=p_cutoff,
            fold_cutoff=fold_cutoff,
            mgm_cutoff=mgm_cutoff,
            int_adjp_cutoff=int_adjp_cutoff,
            int_p_cutoff=int_p_cutoff,
            int_fold_cutoff=int_fold_cutoff,
            ave_cutoff=ave_cutoff,
            collapse_by_gene=collapse_by_gene,
            block=block,
            correlation=correlation,
            posthoc_test=posthoc_test,
            posthoc_args=posthoc_args,
            ...);
      }
      ## rlr_result is a list
      ## - statsDF
      ## - statsDFs
      ## - repFits=list(subFit1,subFit2,subFit3)
      return(rlr_result);
   });

   ## Assemble list of statsDF
   stats_df <- lapply(stats_hits_dfs1, function(i){
      i$stats_df;
   });
   ret_list <- list(stats_df=stats_df);

   ## Assemble list of statsDFs
   stats_dfs <- lapply(stats_hits_dfs1, function(i){
      i$stats_dfs;
   });
   ret_list$stats_dfs <- stats_dfs;

   ## list of named lists
   stats_hits <- lapply(stats_dfs, function(iDFs){
      lapply(iDFs, function(iDF){
         iHitCols <- jamba::nameVector(
            jamba::provigrep("^hit[ .]", colnames(iDF)));
         if (verbose) {
            jamba::printDebug("se_contrast_stats(): ",
               "iHitCols:", iHitCols);
         }
         lapply(iHitCols, function(iHitCol){
            iHitRows <- (!is.na(iDF[,iHitCol]) & iDF[,iHitCol] != 0);
            jamba::nameVector(
               iDF[iHitRows,iHitCol],
               rownames(iDF)[iHitRows]);
         });
      });
   });
   if (verbose) {
      jamba::printDebug("se_contrast_stats(): ",
         "ssdim(stats_hits[[1]]):");
      print(jamba::ssdim(stats_hits[[1]]));
   }

   ## array of named lists
   all_cutoffs_list <- lapply(stats_hits[[1]], function(i){
      gsub("[ ]+[^ ]+$", "", names(i))
   })
   all_cutoffs <- unique(unlist(all_cutoffs_list));
   #arrayDim <- c(length(stats_hits[[1]][[1]]),
   arrayDim <- c(length(all_cutoffs),
      length(stats_hits[[1]]),
      length(stats_hits));
   #arrayDimnames <- list(gsub("[ ]+[^ ]+$", "", names(stats_hits[[1]][[1]])),
   arrayDimnames <- list(all_cutoffs,
      names(stats_hits[[1]]),
      names(stats_hits));
   names(arrayDimnames) <- c("Cutoffs",
      "Contrasts",
      "Signal");
   if (verbose) {
      jamba::printDebug("se_contrast_stats(): ",
         "arrayDim:",
         arrayDim);
      jamba::printDebug("se_contrast_stats(): ",
         "arrayDimnames:");
      print(arrayDimnames);
   }

   # assemble hit_array
   # - allow for missing entries
   # first define the array indices
   kji <- jamba::rbindList(lapply(names(stats_hits), function(i){
      jamba::rbindList(lapply(names(stats_hits[[i]]), function(j){
         jamba::rbindList(lapply(names(stats_hits[[i]][[j]]), function(k){
            k1 <- gsub("[ ][^ ]+$", "", k);
            kji <- cbind(match(k1, arrayDimnames[[1]]),
               match(j, arrayDimnames[[2]]),
               match(i, arrayDimnames[[3]]))
            kji
         }))
      }))
   }))
   hit_array <- array(dim=arrayDim,
      dimnames=arrayDimnames);
   # fill data into the array, which somehow converts it into a list
   # - thanks R
   hit_array[kji] <- unlist(recursive=FALSE,
      unlist(recursive=FALSE,
         stats_hits))
   # create the array again using the list data
   hit_array <- array(dim=arrayDim,
      data=hit_array,
      dimnames=arrayDimnames);

   ret_list$hit_array <- hit_array;
   ret_list$hit_list <- stats_hits;

   ## Add design and contrast data used
   ret_list$idesign <- idesign;
   ret_list$icontrasts <- icontrasts;
   if (length(normgroup) > 0 && TRUE %in% do_normgroups) {
      ret_list$normgroup <- normgroup;
   }

   return(ret_list);
}


#' Handle NA values in a numeric matrix
#'
#' Handle NA values in a numeric matrix
#'
#' This function provides reasonable alternatives intended to
#' manage the presence of missing data encoded as `NA` values
#' in a numeric matrix.
#'
#' The alternatives are defined by argument `handle_na`:
#'
#' * `"full"`: Retain `NA` values, except when an entire group is `NA`
#' it is replaced with `na_value`. This option is intended for sparse data
#' where non-NA values are accepted as real measurements
#' for each group, and where a group with all NA values should be
#' retained for statistical contrasts by assigning `na_value`.
#' This method essentially keeps the data as-is, except when groups
#' are otherwise entirely `NA` the value `na_value` is used in
#' order to retain any relevant contrasts that involve this group.
#' * `"full1"`: Similar to `"full"`, retain `NA` values, except
#' when an entire group is `NA`, then replace only one entry
#' with `na_value`. This option is intended to keep data as-is,
#' except to retain groups that are otherwise entirely `NA`.
#' In these groups, only one `na_value` is used in order to
#' prevent the group from contributing toward group variability
#' or dispersion calculations.
#' * `"partial"`: Replace `NA` values with `na_value`, except
#' when an entire group is `NA` the entire group is kept at `NA`.
#' * `"all"`: Replace all `NA` values with `na_value`.
#' * `"none"`: Perform no replacement of `NA` values.
#'
#' @family jamses stats
#'
#' @export
handle_na_values <- function
(x,
 idesign,
 handle_na=c("full1", "full", "partial", "none", "all"),
 na_value=0,
 na_weight=0,
 verbose=FALSE,
 ...)
{
   handle_na <- match.arg(handle_na);
   xNA <- is.na(x);
   groupL <- multienrichjam::im2list(idesign);
   groupV <- jamba::nameVector(
      rep(names(groupL),
         lengths(groupL)),
      unlist(groupL));
   if ("partial" %in% handle_na) {
      ## We replace NA with zero, except when an entire
      ## group is NA, then we leave it as NA
      # x[xNA] <- 0;
      if (verbose) {
         jamba::printDebug("handle_na_values(): ",
            c("Filling singlet NA with ",
               NAvalue,
               ", leaving full group NA as-is."),
            sep="");
      }
      x <- jamba::rowGroupMeans(x[,names(groupV),drop=FALSE],
         groups=rep(names(groupL), lengths(groupL)),
         rowStatsFunc=function(x,...){
            x1 <- x;
            x1[is.na(x)] <- na_value;
            x1[rowMins(is.na(x)*1) == 1,] <- NA;
            x1;
         });
   } else if ("full" %in% handle_na) {
      ## We replace NA with na_value only when an entire group is NA
      ## in order to retain the group in final results.
      ## Otherwise, partial NA values within a group are kept NA.
      if (verbose) {
         jamba::printDebug("handle_na_values(): ",
            c("Leaving singlet NA as-is, filling full-group NA with ",
               na_value),
            sep="");
      }
      x <- jamba::rowGroupMeans(x[,names(groupV),drop=FALSE],
         groups=rep(names(groupL), lengths(groupL)),
         rowStatsFunc=function(x,...){
            x1 <- x;
            x1[is.na(x)] <- NA;
            x1[rowMins(is.na(x)*1) == 1, ] <- na_value;
            x1;
         });
   } else if ("full1" %in% handle_na) {
      ## Same as "full" in that this method replaces NA only
      ## when an entire group contains NA values, however, it
      ## only replaces one entry with na_value in order to
      ## prevent the value from contributing toward estimate of
      ## group variability.
      if (verbose) {
         jamba::printDebug("handle_na_values(): ",
            c("Leaving singlet NA as-is, filling full-group NA with ",
               na_value,
               " once per group."),
            sep="");
      }
      x <- jamba::rowGroupMeans(x[,names(groupV),drop=FALSE],
         groups=rep(names(groupL), lengths(groupL)),
         rowStatsFunc=function(x,...){
            x1 <- x;
            x1[is.na(x)] <- NA;
            x1[rowMins(is.na(x)*1) == 1, 1] <- na_value;
            x1;
         });
   } else if ("all" %in% handle_na) {
      if (verbose) {
         jamba::printDebug("handle_na_values(): ",
            c("Replacing all NA values with ",
               na_value),
            sep="");
      }
      x[xNA] <- na_value;
   }
   ## define weight matrix
   weights <- jamba::noiseFloor(1-(xNA),
      minimum=na_weight);

   return(list(x=x, weights=weights));
}

#' Limma-voom customized for Jam
#'
#' Limma-voom customized for Jam
#'
#' This function is based directly upon `limma::voom()` with a
#' small adjustment to handle the presence of `NA` values, which
#' otherwise causes the `stats::lowess()` output to be clearly
#' incorrect. The correction removes `NA` values during this step,
#' producing a result as expected.
#'
#' @family jamses stats
#'
#' @export
voom_jam <- function
(counts,
 design=NULL,
 lib.size=NULL,
 normalize.method="none",
 block=NULL,
 correlation=NULL,
 weights=NULL,
 span=0.5,
 plot=FALSE,
 save.plot=TRUE,
 verbose=FALSE,
 ...)
{
   out <- list()

   ## Check counts
   if(is(counts,"DGEList")) {
      out$genes <- counts$genes
      out$targets <- counts$samples
      if(is.null(design) && diff(range(as.numeric(counts$sample$group)))>0) {
         design <- model.matrix(~group,
            data=counts$samples)
      }
      if(is.null(lib.size)) {
         lib.size <- with(counts$samples,
            lib.size * norm.factors);
      }
      counts <- counts$counts;
   } else {
      isExpressionSet <- suppressPackageStartupMessages(is(counts,"ExpressionSet"))
      if(isExpressionSet) {
         if(length(Biobase::fData(counts))) {
            out$genes <- Biobase::fData(counts)
         }
         if(length(Biobase::pData(counts))) {
            out$targets <- Biobase::pData(counts)
         }
         counts <- Biobase::exprs(counts)
      } else {
         counts <- as.matrix(counts)
      }
   }

   n <- nrow(counts);
   if (n < 2L) {
      stop("Need at least two rows to fit a mean-variance trend.");
   }

   ## Check design
   if (is.null(design)) {
      design <- matrix(1, ncol(counts), 1);
      rownames(design) <- colnames(counts);
      colnames(design) <- "GrandMean";
   }

   ## Check lib.size
   if (is.null(lib.size)) {
      lib.size <- colSums(counts,
         na.rm=TRUE);
   }
   if (verbose) {
      jamba::printDebug("voom_jam(): ",
         "lib.size:", format(digits=2,
            big.mark=",",
            scientific=FALSE,
            lib.size));
      jamba::printDebug("voom_jam(): ",
         "span:", format(digits=2,
            big.mark=",",
            span));
   }

   ## Fit linear model to log2-counts-per-million
   y <- t(log2(t(counts + 0.5) / (lib.size + 1) * 1e6));
   y <- limma::normalizeBetweenArrays(y,
      method=normalize.method);
   fit <- limma::lmFit(y,
      design,
      block=block,
      correlation=correlation,
      weights=weights,
      ...);
   if (is.null(fit$Amean)) {
      fit$Amean <- rowMeans(y,
         na.rm=TRUE);
   }

   NWithReps <- sum(fit$df.residual > 0L)
   if (NWithReps < 2L) {
      if (NWithReps == 0L) {
         warning("The experimental design has no replication. Setting weights to 1.")
      }
      if (NWithReps == 1L) {
         warning("Only one gene with any replication. Setting weights to 1.")
      }
      out$E <- y;
      out$weights <- y;
      out$weights[] <- 1;
      out$design <- design;
      if (is.null(out$targets)) {
         out$targets <- data.frame(lib.size=lib.size);
      } else {
         out$targets$lib.size <- lib.size;
      }
      return(new("EList", out));
   }

   ## Fit lowess trend to sqrt-standard-deviations by log-count-size
   sx <- fit$Amean + mean(log2(lib.size + 1)) - log2(1e6);
   sy <- sqrt(fit$sigma);
   allzero <- (rowSums(counts, na.rm=TRUE)==0 %in% c(TRUE,NA));
   if (verbose) {
      jamba::printDebug("voom_jam(): ",
         "head(allzero, 10):");
      print(head(allzero, 10));
   }
   if (any(allzero)) {
      sx <- sx[!allzero];
      sy <- sy[!allzero];
   }
   ## 03dec2019: Major change in lowess logic (below)
   ## If there are any NA values in sx or sy, the results
   ## are highly variable, and not correct.
   ## y1 <- rnorm(5000);
   ## y1[sample(1:5000, size=20)] <- NA;
   ## x1 <- seq(from=0.01, to=10, length.out=5000);
   ## plot(x1, y1);
   ## lines(lowess(x1, y1, f=2/3), col="red");
   ## lines(lowess(x1[!is.na(y1)], y1[!is.na(y1)], f=2/3), col="green");
   #l <- stats::lowess(sx,
   #   sy,
   #   f=span);
   s_no_na <- (!is.na(sx) & !is.na(sy));
   l <- stats::lowess(sx[s_no_na],
      sy[s_no_na],
      f=span,
      ...);
   out$sx <- sx;
   out$sy <- sy;
   out$l <- l;
   if (plot) {
      jamba::plotSmoothScatter(x=sx,
         y=sy,
         xlab="log2(count size + 0.5)",
         ylab="sqrt(standard deviation)",
         pch=16,
         cex=0.25);
      title("voom: Mean-variance trend");
      lines(l, col="red");
   }

   ## Make interpolating rule
   f <- tryCatch({
      approxfun(l,
         rule=2,
         ties=list("ordered", mean));
   }, error=function(e){
      jamba::printDebug("Error in voom_jam() approxfun():",
         fgText="red");
      #print(e);
      jamba::printDebug("head(fit$Amean, 20):");
      print(head(fit$Amean, 20));
      jamba::printDebug("head(sx, 20):");
      print(head(sx, 20));
      jamba::printDebug("head(fit$sigma, 20):");
      print(head(fit$sigma, 20));
      jamba::printDebug("head(sy, 20):");
      print(head(sy, 20 ));
      jamba::printDebug("str(l):");
      print(str(l));
      stop("voom_jam() failed.");
   });

   ## Find individual quarter-root fitted counts
   if (fit$rank < ncol(design)) {
      if (verbose) {
         jamba::printDebug("voom_jam(): ",
            "Using subset of fit$rank, length(fit$rank):",
            length(fit$rank));
      }
      j <- fit$pivot[seq_len(fit$rank)];
      fitted.values <- fit$coef[,j,drop=FALSE] %*% t(fit$design[,j,drop=FALSE])
   } else {
      fitted.values <- fit$coef %*% t(fit$design)
   }

   fitted.cpm <- 2^fitted.values;
   fitted.count <- 1e-6 * t(t(fitted.cpm) * (lib.size + 1));
   fitted.logcount <- log2(fitted.count);

   ## Apply trend to individual observations
   w <- 1 / f(fitted.logcount)^4;
   dim(w) <- dim(fitted.logcount);

   ## Output
   rownames(w) <- rownames(y);
   colnames(w) <- colnames(y);
   out$E <- y;
   out$weights <- w;
   out$design <- design;
   out$other <- list(fitted.values=fitted.values);
   if (is.null(out$targets)) {
      out$targets <- data.frame(lib.size=lib.size);
   } else {
      out$targets$lib.size <- lib.size;
   }
   if (save.plot) {
      out$voom.xy <- list(x=sx,
         y=sy,
         xlab="log2( count size + 0.5 )",
         ylab="Sqrt( standard deviation )");
      out$voom.line <- l;
   }
   new("EList", out);
}

#' Run limma contrasts with optional probe replicates
#'
#' Run limma contrasts with optional probe replicates
#'
#' This function is called by `se_contrast_stats()` to perform
#' the comparisons defined as contrasts. The `se_contrast_stats()`
#' function operates on a `SummarizedExperiment` object,
#' this function operates on the `numeric` `matrix` values
#' directly.
#'
#' This function also calls `ebayes2dfs()` which extracts
#' each contrast result as a `data.frame`, whose column names
#' are modified to include the contrast names.
#'
#' This function optionally (not yet ported from previous
#' implementation) detects replicate probes, and performs
#' the internal correlation calculations recommended by
#' `limma user guide` for replicate probes. In that case,
#' it detects each level of probe replication so that
#' each can be properly calculated. For example, Agilent
#' human 4x44 arrays often contain a large number of probes
#' with 8 replicates; a subset of probes with 4 replicates;
#' then the remaining probes (the majority overall) have
#' only one replicate each. In that case, this function
#' splits data into 8-replicate, 4-replicate, and 1-replicate
#' subsets, calculates correlations on the 8-replicate and
#' 4-replicate subsets separately, then runs limma calculations
#' on the three subsets independently, then merges the results
#' into one large table. The end result is that the
#' final table contains one row per unique probe after
#' adjusting for probe replication properly in each scenario.
#' As the Agilent microarray layout is markedly less widely
#' used that in past, the priority to port this methodology
#' is quite low.
#'
#' @return `list` with the following entries:
#'    * "stats_df": `data.frame` with all individual `data.frame` per contrast,
#'    merged together.
#'    * "stats_df": `list` of individual `data.frame` per contrast, each
#'    result is the output from `ebayes2dfs()`.
#'    * "rep_fits": `list` of various intermediate model fits, dependent
#'    upon whether limma, limma-voom, or limma-DEqMS were used.
#'
#' @inheritParams ebayes2dfs
#'
#' @family jamses stats
#'
#' @export
run_limma_replicate <- function
(imatrix,
 idesign,
 icontrasts,
 weights=NULL,
 robust=FALSE,
 adjust.method="BH",
 confint=FALSE,
 trim_colnames=c("t",
    "B",
    "F",
    "sca.t"),
 adjp_cutoff=0.05,
 p_cutoff=NULL,
 fold_cutoff=1.5,
 int_adjp_cutoff=adjp_cutoff,
 int_p_cutoff=p_cutoff,
 int_fold_cutoff=fold_cutoff,
 mgm_cutoff=NULL,
 ave_cutoff=NULL,
 block=NULL,
 collapse_by_gene=FALSE,
 correlation=NULL,
 posthoc_test=c("none",
    "DEqMS"),
 posthoc_args=list(DEqMS=list(
    PSM_counts=NULL,
    fit.method="loess")),
 verbose=FALSE,
 ...)
{
   # validate arguments
   posthoc_test <- match.arg(posthoc_test);

   # prepate ExpressionSet
   imatrixES <- Biobase::ExpressionSet(assayData=imatrix,
      featureData=new("AnnotatedDataFrame",
         data=data.frame(
            probes=rownames(imatrix),
            row.names=rownames(imatrix)
         )
      )
   );

   ## lmFit
   subFit1 <- limma::lmFit(imatrixES,
      design=idesign,
      weights=weights,
      block=block,
      correlation=correlation,
      ...);
   if (length(rownames(subFit1$coefficients)) == 0) {
      rownames(subFit1$coefficients) <- rownames(subFit1$genes);
   }

   ## Add back in the Amean values (probe means)
   if (length(subFit1$Amean) == 0) {
      subFit1$Amean <- jamba::nameVector(
         rowMeans(imatrix, na.rm=TRUE),
         rownames(imatrix));
      if (verbose) {
         jamba::printDebug("run_limma_replicate(): ",
            "Added Amean values since lmFit did not.");
      }
   }

   ## run contrasts on the limma model
   if (verbose) {
      jamba::printDebug("run_limma_replicate(): ",
         "Running contrasts.fit on subFit1");
   }
   subFit2 <- limma::contrasts.fit(subFit1,
      icontrasts);

   ## eBayes adjustment for signal-based noise
   if (verbose) {
      jamba::printDebug("run_limma_replicate(): ",
         "Running eBayes on subFit2");
   }
   subFit3 <- limma::eBayes(subFit2,
      robust=robust);
   #if (nrow(subMatrix) == ndups) {
   #   subFit3 <- subFit3[1,];
   #}

   ## define colnames to rename to include the contrast name
   renameCols <- c("logFC",
      "P.Value",
      "adj.P.Val",
      "sca.p",
      "sca.adj.pval",
      "t",
      "F",
      "B",
      "CI.L",
      "CI.R",
      "AveExpr");

   ## optionally apply post-hoc tests here, driving example is DEqMS
   subFit4 <- NULL;
   if ("DEqMS" %in% posthoc_test) {
      if (!jamba::check_pkg_installed("DEqMS")) {
         stop(paste("The Bioconductor package DEqMS is not installed,",
            "but posthoc_test='DEqMS'"));
      }
      PSM_counts <- posthoc_args$DEqMS$PSM_counts;
      if (length(posthoc_args$DEqMS$fit.method) == 0) {
         posthoc_args$DEqMS$fit.method <- "loess";
      }
      if (length(PSM_counts) == 0) {
         stop(paste("posthoc_args$DEqMS$PSM_counts is empty,",
            "and is required when using posthoc_test='DEqMS'."));
      }
      if (!all(rownames(subFit3$coefficients) %in% names(PSM_counts))) {
         print(head(subFit3$coefficients));
         stop(paste("all rownames(subFit3$coefficients) must be present in",
            "names(PSM_counts) when using posthoc_test='DEqMS'."));
      }
      subFit3$count <- PSM_counts[rownames(subFit3$coefficients)];
      subFit4 <- DEqMS::spectraCounteBayes(subFit3,
         fit.method=posthoc_args$DEqMS$fit.method);
      renameCols <- c(renameCols,
         "sca.p",
         "sca.adj.pval",
         "sca.t",
         "count");
   }

   ## Get summary table for each contrast
   contrastNames <- colnames(subFit3$contrast);

   dimNum <- nrow(subFit3$coefficients);

   ## top table for each contrast
   if (TRUE) {
      stats_dfs <- ebayes2dfs(lmFit3=subFit3,
         lmFit1=subFit1,
         lmFit4=subFit4,
         define_hits=TRUE,
         trim_colnames=trim_colnames,
         adjp_cutoff=adjp_cutoff,
         p_cutoff=p_cutoff,
         fold_cutoff=fold_cutoff,
         int_adjp_cutoff=int_adjp_cutoff,
         int_p_cutoff=int_p_cutoff,
         int_fold_cutoff=int_fold_cutoff,
         mgm_cutoff=mgm_cutoff,
         ave_cutoff=ave_cutoff,
         posthoc_test=posthoc_test,
         verbose=verbose,
         ...);
   } else {
      stats_dfs <- lapply(jamba::nameVector(contrastNames), function(contrastName) {
         tt <- limma::topTable(subFit3,
            coef=contrastName,
            number=Inf,
            adjust.method=adjust.method,
            confint=confint);

         ## Rename columns to include the contrast
         ## Note: this step could rename "probes" or "genes"
         tt <- jamba::renameColumn(tt,
            from=c(renameCols),
            to=c(paste(renameCols, contrastName)));

         if (verbose) {
            jamba::printDebug("run_limma_replicate(): ",
               "contrastName:",
               contrastName,
               " head(topTable):");
            print(head(tt, 5));
         }
         tt;
      } );
   }

   # stats_df <- jamba::mergeAllXY(stats_dfs);
   stats_df_colnames <- unique(
      unlist(lapply(stats_dfs, colnames)));
   stats_df <- jamba::mergeAllXY(stats_dfs);
   stats_df <- stats_df[, stats_df_colnames, drop=FALSE];

   if (length(jamba::tcount(stats_df[,1], minCount=2)) == 0) {
      rownames(stats_df) <- stats_df[,1];
   }
   # subFit4 is included below, but removed by jamba::rmNULL() when empty
   return(
      list(stats_df=stats_df,
         stats_dfs=stats_dfs,
         rep_fits=jamba::rmNULL(list(
            lmFit1=subFit1,
            lmFit2=subFit2,
            lmFit3=subFit3,
            lmFit4=subFit4))
      )
   )
}


#' Convert limma eBayes fit to data.frame with annotated hits
#'
#' Convert limma eBayes fit to data.frame with annotated hits
#'
#' This function is called by `run_limma_replicate()` as
#' an extension to `limma::topTable()`, that differs in that
#' it is performed for each contrast in the input `lmFit3` object.
#'
#' By default the columns include the contrast, so that each `data.frame`
#' is self-described.
#'
#' When `define_hits=TRUE`, then statistical thresholds are applied
#' to define a set of statistical hits. The thresholds available include:
#'
#' 1. `adjp_cutoff` - applied to `"adj.P.Val"` for adjusted P-value.
#' 2. `p_cutoff` - applied to `"P.Value"` for raw, unadjusted P-value.
#' 3. `fold_cutoff` - normal space fold change, applied to `"logFC"`
#'    by using `log2(fold_cutoff)`.
#' 4. `mgm_cutoff` - max group mean, applied to the highest group mean
#'    value involved in each specific contrast.
#' 5. `ave_cutoff` - applied to `"AveExpr"` which represents the mean
#'    value across all sample groups.
#'
#' Note that `mgm_cutoff` requires input `lmFit1` which stores the
#' group mean values used in the limma workflow.
#'
#' Note also there are optional arguments specific to interaction
#' contrasts, which in this context is assumed to be a
#' "fold change of fold changes" style of contrast, for example:
#' `(groupA-groupB)-(groupC-groupD)`. The purpose is distinct interaction
#' thresholds is to enable reasonable data mining, sometimes with
#' somewhat more lenient thresholds for interaction contrasts.
#' For example, one may use `adjp_cutoff=0.01` and `int_adjp_cutoff=0.05`,
#' or `fold_cutoff=2` and `int_fold_cutoff=1.5`.
#'
#' By default, `rename_headers=TRUE` causes colnames to include the
#' contrast, for example renaming colname `"logFC"` to `"logFC contrastA"`.
#' This change helps reinforce the source of the statistical results,
#' and allows the `data.frame` results to be merged together using
#' `base::merge()`.
#'
#' Indeed, `merge_df=TRUE` will cause all `data.frame` results to be
#' merged into one large `data.frame`, using `jamba::mergeAllXY()`.
#'
#' @return `list` with one `data.frame` per contrast defined in
#'    the input `lmFit3` object. When `define_hits=TRUE` there
#'    will be one column per statistical threshold, named `"hit"`
#'    followed by an abbreviation of the statistical thresholds
#'    which were applied.
#'    When `merge_df=TRUE` the returned data will be one
#'    `data.frame` object.
#'
#' @family jamses stats
#'
#' @param lmFit3 object returned by `limma::eBayes()`.
#' @param lmFit1 object returned by `limma::lmFit()`, optional.
#' @param lmFit4 object returned by `posthoc_test="DEqMS"` in
#'    `run_limma_replicate()`.
#' @param define_hits `logical` indicating whether to define hits
#'    using the statistical thresholds.
#' @param adjp_cutoff,p_cutoff,fold_cutoff,mgm_cutoff,ave_cutoff `numeric`
#'    values representing the appropriate statistical threshold,
#'    or `NULL` when a threshold should not be applied.
#' @param int_adjp_cutoff,int_p_cutoff,int_fold_cutoff `numeric`
#'    thresholds to apply only to interaction contrasts.
#' @param confint `logical` passed to `limma::topTable()`, which defines
#'    whether to return confidence intervals for each log2 fold change.
#' @param use_cutoff_colnames `logical` whether to include the
#'    statistical thresholds abbreviated in the `"hit"` colname,
#'    when `define_hits=TRUE`.
#' @param rename_headers `logical` indicating whether to rename
#'    statistical colnames returned by `limma::topTable()` to the
#'    colnames include the contrast name.
#' @param return_fold `logical` whether to return an additional column
#'    with the signed fold change, see `log2fold_to_fold()`.
#' @param merge_df `logical` indicating whether to merge the final
#'    `data.frame` list into one `data.frame`.
#' @param include_ave_expr `logical` indicating whether to retain
#'    the column `"AveExpr"`. This column can be misleading, especially
#'    if the `mgm` (max group mean) threshold is used when determining
#'    statistical hits. This column is mainly useful in reviewing limma
#'    output, since it uses the `"AveExpr"` values to apply its moderated
#'    variance statistic.
#' @param include_group_means `logical` indicating whether to include each
#'    group mean along with the relevant contrast. These values are
#'    helpful, in that they should exactly represent the reported `logFC`
#'    value. Sometimes it is helpful and comforting to see the exact values
#'    used in that calculation.
#' @param collapse_by_gene `logical` indicating whether to apply
#'    `collapse_stats_by_gene` which chooses one "best" exemplar per gene
#'    when there are multiple rows that represent the same gene.
#' @param rename_contrasts `logical` (inactive) which will in future allow
#'    for automated renaming of contrasts.
#' @param sep `character` string used as a delimiter in certain output
#'    colnames.
#' @param int_grep `character` string used to recognize contrasts which
#'    are considered "interaction contrasts". The default pattern recognizes
#'    any contrasts that contain multiple fold changes, recognized by the
#'    presence of more than one hypen `"-"` in the contrast name.
#' @param verbose `logical` indicating whether to print verbose output.
#'
#' @export
ebayes2dfs <- function
(lmFit3=NULL,
 lmFit1=NULL,
 lmFit4=NULL,
 define_hits=TRUE,
 adjp_cutoff=0.05,
 p_cutoff=NULL,
 fold_cutoff=1.5,
 int_adjp_cutoff=adjp_cutoff,
 int_p_cutoff=p_cutoff,
 int_fold_cutoff=fold_cutoff,
 mgm_cutoff=NULL,
 ave_cutoff=NULL,
 confint=FALSE,
 use_cutoff_colnames=TRUE,
 rename_headers=TRUE,
 return_fold=TRUE,
 merge_df=FALSE,
 include_ave_expr=FALSE,
 include_group_means=TRUE,
 transform_means=c("none", "exp2signed", "10^"),
 collapse_by_gene=FALSE,
 rename_contrasts=FALSE,
 sep=" ",
 int_grep="[(].+-.+-.+[)]|-.+-",
 trim_colnames=c("t",
    "B",
    "F",
    "sca.t"),
 posthoc_test=c("none",
    "DEqMS"),
 verbose=FALSE,
 ...)
{
   ## Purpose is to convert the lmFit3 results of eBayes() into a list of data.frames
   ##
   ## Note the cutoffFold is normal space fold change, but converted to log2fold to compare with limma output
   ## Note cutoffMaxGroupmean is in whatever units are sent to limma... usually log2 intensity or log2 counts
   ##
   ## collapseByGene=TRUE will attempt to produce per-gene results, using pre-defined logic to select the best
   ## exemplar(s) among multiple probes for the same gene.
   ##    1. choose statistical hits first
   ##       a. if multiple hits, same direction, choose them all.
   ##       b. if multiple hits, diff direction,
   ##          i.  mark this gene as multi-direction
   ##          ii. choose best P-value, then take hits with same direction
   ##    2. if no statistical hits, choose entry(ies) above maxMean signal cutoff
   ##    3. If no hits, and no entries above signal cutoff, choose all entries
   ##
   ## renameContrasts=TRUE will rename a fully described contrast into a more human-readable
   ## contrast, e.g.
   ## from: dHSA10GR_EtOH-SW13GR_EtOH
   ## to:   dHSA10GR-SW13GR (EtOH)
   ##
   ## intCutoffPVal and intCutoffAdjPVal are intended to allow
   ## applying a different P-value threshold for interaction effects.
   ##
   ## confint is used by topTable(), when FALSE no confidence intervals are
   ## calculated; when TRUE, by default it calculates 0.95 confidence
   ## intervals, reported as logFC upper and lower bounds, CI.L and CI.R,
   ## respectively.
   ##
   ## Optionally transform the AveExpr values, most useful when reporting normal space values which are stored in log space
   transform_means <- match.arg(transform_means);
   posthoc_test <- match.arg(posthoc_test);

   ## cutoffMaxGroupMean requires lmFit1
   if (!define_hits) {
      mgm_cutoff <- NULL;
   }
   if ((length(jamba::rmNA(mgm_cutoff)) > 0 || include_group_means) && length(lmFit1) == 0) {
      #stop("To use cutoffMaxGroupMean, lmFit1 must be supplied, from which the group means are obtained.");
      jamba::printDebug("ebayes2dfs(): ",
         c("lmFit1 is required for mgm_cutoff or include_group_means. Setting ",
            "mgm_cutoff=NULL",
            ", and ",
            "include_group_means=FALSE"),
         sep="");
      include_group_means <- FALSE;
      mgm_cutoff <- NULL;
   }

   ## TODO: Allow some curation of labels here
   contrastNames <- colnames(coef(lmFit3));
   contrastLabels <- jamba::nameVector(contrastNames, contrastNames);
   is_interaction <- grepl(int_grep,
      ignore.case=TRUE,
      contrastNames);

   # define hits is applied only to each relevant contrast
   define_hits <- rep(define_hits,
      length.out=length(contrastNames));
   names(define_hits) <- contrastNames;

   if (any(define_hits)) {
      ## Add cutoff parameters to the colnames, optional
      cutoff_l <- list(
         mgm=mgm_cutoff,
         p=p_cutoff,
         adjp=adjp_cutoff,
         fc=fold_cutoff);
      cutoff_df <- as.data.frame(jamba::rmNULL(cutoff_l));
      if (length(cutoff_df) == 0 || nrow(cutoff_df) == 0) {
         define_hits[!is_interaction] <- FALSE;
         cutoff_df <- NULL;
      }

      int_cutoff_l <- list(
         mgm=mgm_cutoff,
         p=int_p_cutoff,
         adjp=int_adjp_cutoff,
         fc=int_fold_cutoff);
      int_cutoff_df <- as.data.frame(jamba::rmNULL(int_cutoff_l));
      if (length(int_cutoff_df) == 0 || nrow(int_cutoff_df) == 0) {
         define_hits[is_interaction] <- FALSE;
         int_cutoff_df <- NULL;
      }
   }

   # define cutoff strings for contrasts and interaction contrasts
   if (any(define_hits)) {
      if (length(cutoff_df) > 0) {
         cutoff_string_df <- cutoff_df;
         for (i in colnames(cutoff_df)) {
            cutoff_string_df[[i]] <- paste0(i, cutoff_df[[i]]);
         }
         cutoff_string <- jamba::pasteByRow(cutoff_string_df,
            sep=sep);
      }
      if (length(int_cutoff_df) > 0) {
         int_cutoff_string_df <- int_cutoff_df;
         for (i in colnames(int_cutoff_df)) {
            int_cutoff_string_df[[i]] <- paste0(i, int_cutoff_df[[i]]);
         }
         int_cutoff_string <- jamba::pasteByRow(int_cutoff_string_df,
            sep=sep);
      }

      for (i in colnames(cutoff_df)) {
         cutoff_string_df[[i]] <- paste0(i, cutoff_df[[i]]);
         ## if column is "int"
         ## remove if values are all identical to "non-int" cutoff
         if (jamba::igrepHas("^int", i)) {
            j <- gsub("^int", "", i);
            if (j %in% colnames(cutoff_df)) {
               if (all(cutoff_df[[i]] == cutoff_df[[j]])) {
                  cutoff_string_df[,i] <- list(NULL);
               }
            }
         }
      }
      if (verbose) {
         jamba::printDebug("ebayes2dfs(): ",
            "cutoff_df:");
         print(cutoff_df);
         jamba::printDebug("ebayes2dfs(): ",
            "int_cutoff_df:");
         print(int_cutoff_df);
         jamba::printDebug("ebayes2dfs(): ",
            "cutoff_string_df:");
         print(cutoff_string_df);
         jamba::printDebug("ebayes2dfs(): ",
            "cutoff_string:",
            cutoff_string);
         jamba::printDebug("ebayes2dfs(): ",
            "int_cutoff_string_df:");
         print(int_cutoff_string_df);
         jamba::printDebug("ebayes2dfs(): ",
            "int_cutoff_string:",
            int_cutoff_string);
      }
      #return(cutoff_string_df);
   }

   ## assign rownames if not present in lmFit1$coefficients
   if (length(lmFit1) > 0 && length(rownames(lmFit1$coefficients)) == 0) {
      jamba::printDebug("ebayes2dfs(): ",
         c("Note there are no",
            " rownames(lmFit1$coefficients)"),
         sep="",
         fgText=c("darkorange1", "red1"));
      if ("genes" %in% names(lmFit1)) {
         rownames(lmFit1$coefficients) <- rownames(lmFit1$genes);
      }
   }
   if (length(rownames(lmFit3$coefficients)) == 0) {
      jamba::printDebug("ebayes2dfs(): ",
         c("Note there are no",
            " rownames(lmFit3$coefficients)"),
         sep="",
         fgText=c("darkorange1", "red1"));
      if ("genes" %in% names(lmFit3)) {
         rownames(lmFit3$coefficients) <- rownames(lmFit3$genes);
      }
   }

   ## TODO: handle cases with zero residual degrees of freedom,
   ## where we would not have a P-value but still have fold changes.
   ## Examples would be per-patient fold changes.
   ##
   ## lmTopTables is a list:
   ## - named by contrastNames
   ## - containing elements "iTopTableDF"
   ## - if collapseByGene=TRUE
   ##    - "iTopTableByGene"
   ##    - "multiDirProbes"
   lmTopTables <- lapply(jamba::nameVector(contrastNames), function(i){
      retVals <- list();
      iLabel <- contrastLabels[i];
      ## Detect whether the contrast is an interaction effect (2-way or 3-way)
      ## or is a simple pairwise style t-test
      isInteraction <- jamba::igrepHas(int_grep, iLabel);
      if (verbose) {
         if (isInteraction) {
            jamba::printDebug("ebayes2dfs(): ",
               c("Interaction effect detected for contrast:",
                  iLabel),
               sep="",
               fgText=c("darkorange1", "purple"));
         } else {
            jamba::printDebug("ebayes2dfs(): ",
               c("Evaluting contrast:",
                  iLabel),
               sep="");
         }
      }
      #if (rename_contrasts) {
      #   iLabel <- renameTwoFactorComparisons(iLabel);
      #}

      # confirm there is some data available in "genes" by using
      # rownames(coefficients) as a backup plan
      if (!"genes" %in% names(lmFit3)) {
         lmFit3$genes <- rownames(lmFit3$coefficients);
      }

      ## Assemble top table, handling single replicate data in a specific way
      if (!any(lmFit3$df.residual > 0)) {
         jamba::printDebug("ebayes2dfs(): ",
            "No values with df.residual > 0.",
            fgText=c("darkorange1", "red1"));
         if (confint) {
            iTopTable <- data.frame(
               check.names=FALSE,
               stringsAsFactors=FALSE,
               lmFit3$genes,
               logFC=lmFit3$coefficients[,i],
               CI.L=lmFit3$coefficients[,i],
               CI.R=lmFit3$coefficients[,i],
               adj.P.Val=1,
               P.Value=1,
               AveExpr=lmFit3$Amean);
         } else {
            iTopTable <- data.frame(
               check.names=FALSE,
               stringsAsFactors=FALSE,
               lmFit3$genes,
               logFC=lmFit3$coefficients[,i],
               adj.P.Val=1,
               P.Value=1,
               AveExpr=lmFit3$Amean);
         }
      } else {
         if (length(lmFit4) > 0) {
            if ("DEqMS" %in% posthoc_test) {
               if (verbose) {
                  jamba::printDebug("ebayes2dfs(): ",
                     c("DEqMS::outputResult(lmFit4) for contrast: ",
                        iLabel),
                     sep="");
                  print(head(lmFit4$coefficients));
                  print(head(lmFit4$sca.t));
               }
               # determine coefficient index position
               # with only one coefficient the colname gets dropped
               # from the numeric matrix lmFit4$sca.t
               # produced by DEqMS. So DEqMS::outputResult() fails to find the
               # column by name, and must refer to the column by integer number
               coef_col <- match(i,
                  colnames(lmFit4$coefficients))
               iTopTable <- DEqMS::outputResult(lmFit4,
                  coef_col=coef_col);
               if (verbose) {
                  jamba::printDebug("ebayes2dfs(): ",
                     "head(iTopTable): ");
                  print(head(head(iTopTable)));
               }
            } else {
               if (verbose) {
                  jamba::printDebug("ebayes2dfs(): ",
                     c("limma::topTable(lmFit4) for contrast: ",
                        iLabel),
                     sep="");
               }
               iTopTable <- limma::topTable(
                  lmFit4,
                  coef=i,
                  sort.by="none",
                  number=nrow(lmFit4),
                  confint=confint);
            }
         } else {
            if (verbose) {
               jamba::printDebug("ebayes2dfs(): ",
                  c("limma::topTable(lmFit3) for contrast: ",
                     iLabel),
                  sep="");
            }
            iTopTable <- limma::topTable(
               lmFit3,
               coef=i,
               sort.by="none",
               number=nrow(lmFit3),
               confint=confint);
         }
      }
      ## Optionally remove extraneous colnames
      if (length(trim_colnames) > 0) {
         iTopTable <- iTopTable[,setdiff(colnames(iTopTable), trim_colnames), drop=FALSE];
      }

      ## Optionally include group mean and maxgroupmean values
      if (length(lmFit1) > 0) {
         ## Improve maxGroupMean by using the specific groups included in the contrast
         iCoefCols <- names(which(lmFit3$contrasts[,i] != 0));
         iCoefLabs <- paste(iCoefCols,
            sep=sep,
            "mean");
         coef_match <- match(rownames(iTopTable),
            rownames(lmFit1$coefficients));
         gm_m <- jamba::renameColumn(
            lmFit1$coefficients[coef_match, iCoefCols, drop=FALSE],
            from=iCoefCols,
            to=iCoefLabs);
         mgm <- matrixStats::rowMaxs(gm_m,
            na.rm=TRUE);
         if (include_group_means) {
            iTopTable[,colnames(gm_m)] <- gm_m;
         }
         iTopTable[,"mgm"] <- mgm;
      }

      ## Apply statistical thresholds
      if (define_hits[i] || collapse_by_gene) {
         probe_colname <- head(colnames(iTopTable), 1);
         if (collapse_by_gene) {
            if (verbose) {
               jamba::printDebug("ebayes2dfs(): ",
                  "collapse_by_gene.");
            }
            gene_colname <- head(jamba::provigrep(
               c("GeneName", "geneSymbol", "gene"),
               colnames(iTopTable)), 1);
            isGenes <- nameVector(iTopTable[,gene_colname],
               rownames(iTopTable));
         }
         pcol <- "p";
         adjpcol <- "adjp";
         foldcol <- "fc";

         ## iterate each set of thresholds by row in cutoff_df
         ## but only iterate unique set of cutoff thresholds
         if (isInteraction) {
            k_rows <- match(unique(int_cutoff_string),
               int_cutoff_string);
            use_cutoff_string <- int_cutoff_string;
            use_cutoff_df <- int_cutoff_df;
         } else {
            k_rows <- match(unique(cutoff_string),
               cutoff_string);
            use_cutoff_string <- cutoff_string;
            use_cutoff_df <- cutoff_df;
         }
         for (k in k_rows) {
            hit_colname <- paste("hit",
               sep=sep,
               use_cutoff_string[k]);
            mgm_cutoff <- ifelse("mgm" %in% colnames(use_cutoff_df),
               use_cutoff_df[k, "mgm"],
               NA);
            ave_cutoff <- ifelse("ave" %in% colnames(use_cutoff_df),
               use_cutoff_df[k, "ave"],
               NA);
            p_cutoff <- ifelse(pcol %in% colnames(use_cutoff_df),
               use_cutoff_df[k, pcol],
               NA);
            adjp_cutoff <- ifelse(adjpcol %in% colnames(use_cutoff_df),
               use_cutoff_df[k, adjpcol],
               NA);
            fold_cutoff <- ifelse(foldcol %in% colnames(use_cutoff_df),
               use_cutoff_df[k, foldcol],
               NA);

            if (verbose) {
               jamba::printDebug("ebayes2dfs(): ",
                  c("hit_colname:",
                     hit_colname),
                  sep="");
            }

            ## Define colnames to use for each cutoff
            adjp_colname<- "adj.P.Val"
            p_colname <- "P.Value"
            logfc_colname <- "logFC"
            mgm_colname <- "mgm"
            ave_colname <- "AveExpr"
            if ("DEqMS" %in% posthoc_test) {
               adjp_colname <- "sca.adj.pval";
               p_colname <- "sca.P.Value";
            }

            ## Call utility function to get hit flags
            hit_values <- mark_stat_hits(x=iTopTable,
               adjp_cutoff=adjp_cutoff,
               p_cutoff=p_cutoff,
               fold_cutoff=fold_cutoff,
               mgm_cutoff=mgm_cutoff,
               ave_cutoff=ave_cutoff,
               adjp_colname=adjp_colname,
               p_colname=p_colname,
               logfc_colname=logfc_colname,
               mgm_colname=mgm_colname,
               ave_colname=ave_colname);
            iTopTable[[hit_colname]] <- hit_values;
         }
      }

      ## Optionally transform intensities, e.g. exponentiating log2 values
      ## Note that the 2^ conversion now subtracts 1, since the typical
      ## transformation is log2(1+x)
      if (FALSE) {
         if (transformAveExpr %in% c("2^")) {
            #iTopTable[,"AveExpr"] <- 2^iTopTable[,"AveExpr"];
            iTopTable[,"AveExpr"] <- 2^iTopTable[,"AveExpr"] - 1;
            if ("maxGroupMean" %in% colnames(iTopTable)) {
               #iTopTable[,"maxGroupMean"] <- 2^iTopTable[,"maxGroupMean"];
               iTopTable[,"maxGroupMean"] <- 2^iTopTable[,"maxGroupMean"] - 1;
            }
         } else if (transformAveExpr %in% c("10^")) {
            #iTopTable[,"AveExpr"] <- 10^iTopTable[,"AveExpr"];
            iTopTable[,"AveExpr"] <- 10^iTopTable[,"AveExpr"] - 1;
            if ("maxGroupMean" %in% colnames(iTopTable)) {
               #iTopTable[,"maxGroupMean"] <- 10^iTopTable[,"maxGroupMean"];
               iTopTable[,"maxGroupMean"] <- 10^iTopTable[,"maxGroupMean"] - 1;
            }
         }
      }
      ## Optionally convert log2 fold change to normal fold change
      if (return_fold && !"fold" %in% colnames(iTopTable)) {
         iTopTable[,"fold"] <- log2fold_to_fold(iTopTable[,"logFC"]);
      }
      ################################################
      ## Per-gene collapse:
      ##
      ## Note we ultimately define per-gene rows using the first stats hit criteria,
      ## otherwise too many columns will be created.  E.g. changing the stats hit criteria
      ## changes the prioritization of probes to include in the per-gene row, so it affects
      ## the fold change, groupMean, P-Value, and adj-P-Val.
      ## The current solution is to use only the first hit filters, so only one set of these
      ## values is propagated downstream.
      ## However, the "hit" columns can represent multiple stats hit criteria.
      ##
      ## TODO: implement the interaction P-value cutoffs here as well
      if (collapse_by_gene) {
         if (verbose) {
            printDebug("collapseByGene iTopTable:");
            print(head(iTopTable));
            printDebug("geneColname:", geneColname);
         }
         # Note that DEqMS should not ever proceed here since the
         # method is inherently based upon per-gene logic.
         iTopTableByGeneL <- collapseTopTableByGene(iTopTable,
            geneColname=geneColname,
            aveExprColname="AveExpr",
            maxGroupMeanColname="maxGroupMean",
            adjPvalColname="adj.P.Val",
            PvalueColname="P.Value",
            logFCcolname="logFC",
            cutoffAveExpr=ave_cutoff[1],
            cutoffMaxGroupMean=mgm_cutoff[1],
            cutoffAdjPVal=adjp_cutoff[1],
            cutoffPVal=p_cutoff[1],
            cutoffFold=fold_cutoff[1]);
         iTopTableByGene <- iTopTableByGeneL$iTopTableByGene;
         retVals$top_bygene_df <- iTopTableByGene;
         retVals$multidir_probes <- iTopTableByGeneL$multiDirProbes;
      }

      ## Optionally rename column headers to include the contrast name
      gene_colnames <- intersect(colnames(iTopTable),
         c(colnames(lmFit3$genes),
            "gene",
            probe_colname));

      # optionally remove AveExpr
      if (!include_ave_expr) {
         iTopTable <- iTopTable[,jamba::unvigrep("AveExpr", colnames(iTopTable)),drop=FALSE];
      }

      # rename headers to include the contrast name
      if (rename_headers) {
         ## Note we do rename maxGroupMean now,
         ## since we calculate it per each specific contrast
         #
         # rename_from <- jamba::vigrep(c("^AveExpr$|^mgm$|p.value|adj.p.val|sca.adj.pval|^sca.t$|^count$|^logFC$|^fold$"),
         #    colnames(iTopTable));
         rename_from <- setdiff(
            jamba::unvigrep("AveExpr| mean$|^gene|^symbol|^probe",
               colnames(iTopTable)),
            gene_colnames);
         rename_to <- paste(rename_from,
            iLabel,
            sep=sep);
         # remove repeat blank spaces
         rename_to <- gsub("^[ ]+|[ ]+$", "",
            gsub("[ ]+", " ",
               rename_to));
         if (verbose) {
            jamba::printDebug("ebayes2dfs(): ",
               "rename stat columns:");
            print(data.frame(rename_from, rename_to));
         }
         iTopTable <- jamba::renameColumn(iTopTable,
            from=rename_from,
            to=rename_to);
      }

      # re-order colnames
      if (TRUE) {
         neworder_colnames <- jamba::provigrep(c(
            "^probe|^gene|symbol|accession|accnum",
            paste0("^hit", sep),
            paste0("^logFC", sep),
            paste0("^fold", sep),
            "^sca.p.value",
            "^sca.adj.pval",
            "^p.value",
            "^adj.p.val",
            paste0("^mgm", sep),
            " mean$",
            paste0("^AveExpr", sep),
            "."),
            colnames(iTopTable));
         # disable old logic for now, will remove in near future
         if (FALSE) {
            reorder_colnames <- jamba::vgrep("^(t|F|B|logFC|P.Value|adj.P.Val|fold|AveExpr)($| )",
               colnames(iTopTable));
            reorder_match <- match(reorder_colnames, colnames(iTopTable));
            keep_colnames1 <- colnames(iTopTable)[seq_len(min(reorder_match) - 1)];
            keep_colnames2 <- setdiff(colnames(iTopTable),
               c(keep_colnames1, reorder_colnames));
            reorder_sorted <- provigrep(c("logFC",
               "fold",
               "P.Value",
               "adj.P.Val",
               "."),
               reorder_colnames);
            neworder_colnames <- c(keep_colnames1,
               reorder_sorted,
               keep_colnames2);
         }
         iTopTable <- iTopTable[, neworder_colnames, drop=FALSE];
      }
      rownames(iTopTable) <- jamba::makeNames(iTopTable[,1]);

      retVals$top_df <- iTopTable;
      retVals;
   });

   ## Now prepare the data to return
   if (merge_df) {
      lmTopTablesAll <- jamba::mergeAllXY(lapply(lmTopTables, function(i){
         i$top_df;
      }));
      if (collapse_by_gene) {
         lmTopTablesAllG <- jamba::mergeAllXY(lapply(lmTopTables, function(i){
            i$top_bygene_df;
         }));
      }
      ## Re-order columns so the genes, then hits, appear first
      colname_order <- jamba::provigrep(c(gene_colnames, "^hit", "."),
         colnames(lmTopTablesAll));
      lmTopTablesAll <- lmTopTablesAll[,colname_order, drop=FALSE];
      rownames(lmTopTablesAll) <- jamba::makeNames(lmTopTablesAll[,head(gene_colnames, 1)]);
   } else {
      lmTopTablesAll <- lapply(lmTopTables, function(i){
         i$top_df;
      });
      if (collapse_by_gene) {
         lmTopTablesAllG <- lapply(lmTopTables, function(i){
            i$top_bygene_df;
         });
      }
   }

   if (any(define_hits)) {
      if (length(cutoff_df) > 0) {
         attr(lmTopTablesAll, "cutoff_df") <- cutoff_df;
      }
      if (length(int_cutoff_df) > 0) {
         attr(lmTopTablesAll, "int_cutoff_df") <- int_cutoff_df;
      }
   }
   if (collapse_by_gene) {
      if (any(define_hits)) {
         if (length(cutoff_df) > 0) {
            attr(lmTopTablesAllG, "cutoff_df") <- cutoff_df;
         }
         if (length(int_cutoff_df) > 0) {
            attr(lmTopTablesAllG, "int_cutoff_df") <- int_cutoff_df;
         }
      }
      retList <- list(top_df=lmTopTablesAll,
         top_bygene_df=lmTopTablesAllG);
      return(retList);
   }
   return(lmTopTablesAll);
}

