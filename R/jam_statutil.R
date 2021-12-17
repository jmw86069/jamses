
#' Mark statistical hits by threshold cutoffs
#'
#' Mark statistical hits by threshold cutoffs
#'
#' This function is lightweight method of applying one or more
#' statistical thresholds to define "statistical hits". The
#' thresholds are based upon three questions:
#'
#' * Is it "detected"? (above minimum signal)
#' * Is it "changing"? (above minimum fold change)
#' * Is it "significant"? (below defined P-value threshold)
#'
#' The reasoning is roughly described:
#'
#' * If the signal for a measurement is not above a noise threshold,
#' or above a defined level of signature required for adequate
#' confirmation experiments, the other statistical measurements
#' are not relevant.
#' * If the change between two experimental groups is not sufficient
#' for follow-up experiments, or is below a biologically meaningful
#' level of change, the other statistical measures are not relevant.
#' * If the signal is detected, and the change is potentially sufficient
#' for follow-up confirmation experiments, and/or to induce biologically
#' meaningful effects, it must also be statistically robust as
#' defined by the relevant adjusted P-value.
#'
#' These thresholds are dependent upon the experiment itself,
#' and each threshold, if used, must be well-defined and defensible.
#'
#' That is, in order to define a signal threshold, one should evaluate
#' the level of noise below which a measured value is no longer
#' sufficient for follow-up experiments, or no longer reliable
#' based upon the technology being used.
#'
#' In order to impost a minimum fold change threshold, one should
#' have some clear indication of any limitations in follow-up assay
#' techniques, and some indication of the magnitude of change expected
#' for a biologically meaningful response. In some cases, a biologically
#' meaningful change may be defined in other experiments, ideally showing
#' small changes not associated with biologically meaningful effects
#' and changes which are associated with biollogically meaningful effects.
#'
#' A useful technique to review statistical thresholds is a volcano plot,
#' which depicts the relationship of log fold change versus adjusted P-value.
#' The plot can indicate the range of fold changes for which the statistical
#' model found significance. Some technologies or protocols naturally
#' compress the effective fold change, yielding a very narrow volcano plot,
#' while others with high variability may result in a relatively short
#' and wide volcano plot. The range of fold changes with no significant
#' P-value may indicate a reasonable expectation for inherent variability,
#' thus a fold change threshold may be defined above that observed for
#' the majority of non-significant entries.
#'
#' Entries which meet the statistical criteria are marked:
#'
#' * `-1` for entries that meet all criteria, with negative fold change
#' * `0` for entries that do not meet all thresholds
#' * `1` for entries that meet all criteria, with positive fold change
#'
#' "Detected" is defined either by the "max group mean", representing
#' the highest group mean signal intensity, or by "average signal",
#' representing the average group mean signal intensity across all
#' groups. These columns should already be present in the input
#' data `x`.
#'
#' "Changing" is defined by the log2 fold change, which must
#' meet the criteria defined by `fold_cutoff` which is in normal
#' space. For example, `fold_cutoff=1.5` represents 1.5-fold change,
#' and would be applied `abs(log2fc) >= log2(fold_cutoff)`. Most
#' statistical results are reported using log2 fold change, but
#' scientists usually define fold changes in normal space.
#'
#' "Significant" is define using the P-value, and/or adjusted P-value,
#' and requires entries to be at or below the threshold.
#'
#' @return `numeric` vector with length `nrow(x)` and values
#'    defined by argument `assign_value`. The order is identical
#'    to the order of rows in `x` input. The output vector
#'    will be named by `rownames(x)` if rownames exist.
#'
#' @family jamses utilities
#'
#' @param x `data.frame` containing one or more statistical columns
#' @param adjp_cutoff,p_cutoff,fold_cutoff,mgm_cutoff,ave_cutoff `numeric`
#'    value for each cutoff to be enforced, or `NULL` or `NA` to ignore
#'    each threshold. Each argument must have only 1 value assigned to
#'    be enforced.
#' @param adjp_colname,p_colname,logfc_colname,mgm_colname,ave_colname
#'    `character` string for each colname in `x` to be used for the
#'    appropriate statistical threshold.
#' @param assign_value `character` string indicating the value assigned
#'    to hits: `"sign"` uses the sign of the log2 fold change; `"fold"`
#'    uses the normal space fold change, by `log2fold_to_fold()`;
#'    `"logfc"` uses the log2 fold change value. If there is no
#'    matching colname for `logfc_colname` then all hits are assigned `1`.
#'    In all cases, entries which are not hits are assigned `0`.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
mark_stat_hits <- function
(x,
 adjp_cutoff=NULL,
 p_cutoff=NULL,
 fold_cutoff=NULL,
 mgm_cutoff=NULL,
 ave_cutoff=NULL,
 adjp_colname="adj.P.Val",
 p_colname="P.Value",
 logfc_colname="logFC",
 mgm_colname="mgm",
 ave_colname="AveExpr",
 assign_value=c("sign", "fold", "logfc"),
 verbose=FALSE,
 ...)
{
   ## validate assign_value
   assign_value <- match.arg(assign_value);
   ## Check stat colname to make sure at least one can be applied
   stats_applied <- character(0);

   ## Detected rows, using mgm and ave
   if (mgm_colname %in% colnames(x) &&
         length(mgm_cutoff) == 1) {
      is_mgm <- (x[[mgm_colname]] >= jamba::rmNA(naValue=-Inf, mgm_cutoff));
      stats_applied <- c(stats_applied, "mgm");
   } else {
      is_mgm <- rep(TRUE, nrow(x));
   }
   if (ave_colname %in% colnames(x) &&
         length(ave_cutoff) == 1) {
      is_ave <- (x[[ave_colname]] >= jamba::rmNA(naValue=-Inf, ave_cutoff));
      stats_applied <- c(stats_applied, "ave");
   } else {
      is_ave <- rep(TRUE, nrow(x));
   }
   is_detected <- (is_mgm & is_ave) * 1;

   ## significance using P.Value and adj.P.Val
   if (p_colname %in% colnames(x) &&
         length(p_cutoff) == 1) {
      is_p <- (x[[p_colname]] <= jamba::rmNA(naValue=Inf, p_cutoff));
      stats_applied <- c(stats_applied, "p");
   } else {
      is_p <- rep(TRUE, nrow(x));
   }
   if (adjp_colname %in% colnames(x) &&
         length(adjp_cutoff) == 1) {
      is_adjp <- (x[[adjp_colname]] <= jamba::rmNA(naValue=Inf, adjp_cutoff));
      stats_applied <- c(stats_applied, "adjp");
   } else {
      is_adjp <- rep(TRUE, nrow(x));
   }
   is_significant <- (is_p & is_adjp) * 1;

   ## changing above the fold change threshold
   if (logfc_colname %in% colnames(x) &&
         length(fold_cutoff) == 1) {
      is_fc <- (abs(x[[logfc_colname]]) >= log2(jamba::rmNA(naValue=0, fold_cutoff)));
      stats_applied <- c(stats_applied, "fold");
   } else {
      is_fc <- rep(TRUE, nrow(x));
   }
   is_changing <- (is_fc);

   ## direction of change
   if (logfc_colname %in% colnames(x)) {
      if ("sign" %in% assign_value) {
         i_value <- sign(x[[logfc_colname]]);
      } else if ("fold" %in% assign_value) {
         i_value <- log2fold_to_fold(x[[logfc_colname]]);
      } else if ("logfc" %in% assign_value) {
         i_value <- x[[logfc_colname]];
      }
   } else {
      i_value <- rep(1, nrow(x));
   }

   if (verbose) {
      jamba::printDebug("mark_stat_hits(): ",
         "stats applied:",
         stats_applied);
   }
   if (length(stats_applied) == 0) {
      warning("Note: mark_stat_hits() did not apply any stat thresholds.");
   }
   ## Put it all together
   hit_values <- (is_detected & is_significant & is_changing) * i_value;
   if (length(rownames(x)) > 0) {
      names(hit_values) <- rownames(x);
   }
   return(hit_values);

}


#' Convert log2 fold change to signed fold change
#'
#' Convert log2 fold change to signed fold change
#'
#' This function takes log2 fold change values as input, and returns
#' normal space fold change values that retain the positive and negative
#' sign, and the magnitude.
#'
#' For example:
#'
#' * `log2 fold change = 2` becomes `fold change = 4`.
#' * `log2 fold change = -2` becomes `fold change = -4`.
#'
#' This function therefore differs from similar functions that convert
#' log2 fold change into a ratio. Instead, `log2fold_to_fold()`
#' specifically retains the magnitude of negative changes.
#'
#' @param x `numeric` vector
#'
#' @family jamses utilities
#'
#' @return `numeric` vector representing signed fold change values.
#'
#' @examples
#' x <- c(-3, -2, -1, 0, 1, 2, 3);
#' fc <- log2fold_to_fold(x);
#' fc;
#'
#' fold_to_log2fold(fc);
#'
#' @export
log2fold_to_fold <- function
(x,
   ...)
{
   x_sign <- sign(x);
   x_sign <- ifelse(x_sign == 0, 1, x_sign);
   2^(abs(x)) * x_sign;
}

#' Convert normal signed fold change to log2 fold change
#'
#' Convert normal signed fold change to log2 fold change
#'
#' This function takes fold change values as input, and returns
#' log2 fold change values.
#'
#' This function recognizes two forms of input:
#'
#' * ratio, which includes values between 0 and 1, but no negative values;
#' * fold change, as from `log2fold_to_fold()` which includes no values
#' between 0 and 1, but may include negative values.
#'
#' For example, for ratio input:
#'
#' * `ratio = 4` becomes `log2 fold change = 2`.
#' * `ratio = 0.25` becomes `log2 fold change = -2`.
#'
#' For example, for fold change input:
#'
#' * `fold change = 4` becomes `log2 fold change = 2`.
#' * `fold change = -4` becomes `log2 fold change = -2`.
#'
#' @param x `numeric` vector
#'
#' @return `numeric` vector representing log2 fold change values.
#'
#' @family jamses utilities
#'
#' @examples
#' x <- c(-3, -2, -1, 0, 1, 2, 3);
#' fc <- log2fold_to_fold(x);
#' fc;
#'
#' fold_to_log2fold(fc);
#'
#' @export
fold_to_log2fold <- function
(x,
   ...)
{
   is_ratio <- (x > 0 & x < 1);
   ifelse(is_ratio,
      log2(x),
      log2(abs(x)) * sign(x));
}
