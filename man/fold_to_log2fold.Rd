% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jam_statutil.R
\name{fold_to_log2fold}
\alias{fold_to_log2fold}
\title{Convert normal signed fold change to log2 fold change}
\usage{
fold_to_log2fold(x, ...)
}
\arguments{
\item{x}{\code{numeric} vector}
}
\value{
\code{numeric} vector representing log2 fold change values.
}
\description{
Convert normal signed fold change to log2 fold change
}
\details{
This function takes fold change values as input, and returns
log2 fold change values.

This function recognizes two forms of input:
\itemize{
\item ratio, which includes values between 0 and 1, but no negative values;
\item fold change, as from \code{log2fold_to_fold()} which includes no values
between 0 and 1, but may include negative values.
}

For example, for ratio input:
\itemize{
\item \code{ratio = 4} becomes \verb{log2 fold change = 2}.
\item \code{ratio = 0.25} becomes \verb{log2 fold change = -2}.
}

For example, for fold change input:
\itemize{
\item \verb{fold change = 4} becomes \verb{log2 fold change = 2}.
\item \verb{fold change = -4} becomes \verb{log2 fold change = -2}.
}
}
\examples{
x <- c(-3, -2, -1, 0, 1, 2, 3);
fc <- log2fold_to_fold(x);
fc;

fold_to_log2fold(fc);

}
\seealso{
Other jamses utilities: 
\code{\link{log2fold_to_fold}()},
\code{\link{mark_stat_hits}()},
\code{\link{se_collapse_by_column}()},
\code{\link{se_collapse_by_row}()},
\code{\link{se_rbind}()},
\code{\link{se_to_rowcoldata}()},
\code{\link{shortest_unique_abbreviation}()},
\code{\link{shrinkDataFrame}()},
\code{\link{shrink_df}()},
\code{\link{shrink_matrix}()},
\code{\link{strsplitOrdered}()},
\code{\link{sub_split_vector}()},
\code{\link{update_function_params}()},
\code{\link{update_list_elements}()}
}
\concept{jamses utilities}
