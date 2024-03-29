% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jam_secontrasts.R
\name{run_limma_replicate}
\alias{run_limma_replicate}
\title{Run limma contrasts with optional probe replicates}
\usage{
run_limma_replicate(
  imatrix,
  idesign,
  icontrasts,
  weights = NULL,
  robust = FALSE,
  adjust.method = "BH",
  confint = FALSE,
  trim_colnames = c("t", "B", "F", "sca.t"),
  adjp_cutoff = 0.05,
  p_cutoff = NULL,
  fold_cutoff = 1.5,
  int_adjp_cutoff = adjp_cutoff,
  int_p_cutoff = p_cutoff,
  int_fold_cutoff = fold_cutoff,
  mgm_cutoff = NULL,
  ave_cutoff = NULL,
  block = NULL,
  rowData_df = NULL,
  collapse_by_gene = FALSE,
  correlation = NULL,
  posthoc_test = c("none", "DEqMS"),
  posthoc_args = list(DEqMS = list(PSM_counts = NULL, fit.method = "loess")),
  seed = 123,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{confint}{\code{logical} passed to \code{limma::topTable()}, which defines
whether to return confidence intervals for each log2 fold change.}

\item{adjp_cutoff, p_cutoff, fold_cutoff, mgm_cutoff, ave_cutoff}{\code{numeric}
values representing the appropriate statistical threshold,
or \code{NULL} when a threshold should not be applied.}

\item{int_adjp_cutoff, int_p_cutoff, int_fold_cutoff}{\code{numeric}
thresholds to apply only to interaction contrasts.}

\item{rowData_df}{\code{data.frame} representing optional rowData annotation
to be retained in the resulting stat \code{data.frame}. This argument
is usually defined using \code{rowData_colnames} in \code{se_contrast_stats()},
which uses corresponding columns from \code{rowData(se)}.}

\item{collapse_by_gene}{\code{logical} indicating whether to apply
\code{collapse_stats_by_gene} which chooses one "best" exemplar per gene
when there are multiple rows that represent the same gene.}

\item{correlation}{\code{numeric} or \code{NULL} passed to \code{limma::lmFit()}.
Note that when \code{block} is defined (and non-empty), and when
\code{correlation=NULL}, the correlation will be calculated by
calling \code{limma::duplicateCorrelation()}.}

\item{seed}{\code{numeric} value used to define \code{set.seed()} for reproducibility.
To avoid setting seed, use \code{seed=NULL}.}

\item{verbose}{\code{logical} indicating whether to print verbose output.}
}
\value{
\code{list} with the following entries:
\itemize{
\item "stats_df": \code{data.frame} with all individual \code{data.frame} per contrast,
merged together.
\item "stats_df": \code{list} of individual \code{data.frame} per contrast, each
result is the output from \code{ebayes2dfs()}.
\item "rep_fits": \code{list} of various intermediate model fits, dependent
upon whether limma, limma-voom, or limma-DEqMS were used.
}
}
\description{
Run limma contrasts with optional probe replicates
}
\details{
This function is called by \code{se_contrast_stats()} to perform
the comparisons defined as contrasts. The \code{se_contrast_stats()}
function operates on a \code{SummarizedExperiment} object,
this function operates on the \code{numeric} \code{matrix} values
directly.

This function also calls \code{ebayes2dfs()} which extracts
each contrast result as a \code{data.frame}, whose column names
are modified to include the contrast names.

This function optionally (not yet ported from previous
implementation) detects replicate probes, and performs
the internal correlation calculations recommended by
\verb{limma user guide} for replicate probes. In that case,
it detects each level of probe replication so that
each can be properly calculated. For example, Agilent
human 4x44 arrays often contain a large number of probes
with 8 replicates; a subset of probes with 4 replicates;
then the remaining probes (the majority overall) have
only one replicate each. In that case, this function
splits data into 8-replicate, 4-replicate, and 1-replicate
subsets, calculates correlations on the 8-replicate and
4-replicate subsets separately, then runs limma calculations
on the three subsets independently, then merges the results
into one large table. The end result is that the
final table contains one row per unique probe after
adjusting for probe replication properly in each scenario.
As the Agilent microarray layout is markedly less widely
used that in past, the priority to port this methodology
is quite low.
}
\seealso{
Other jamses stats: 
\code{\link{ebayes2dfs}()},
\code{\link{handle_na_values}()},
\code{\link{matrix_normalize}()},
\code{\link{se_contrast_stats}()},
\code{\link{se_normalize}()},
\code{\link{voom_jam}()}
}
\concept{jamses stats}
