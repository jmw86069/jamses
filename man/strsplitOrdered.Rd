% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jam_groupstosedesign.R
\name{strsplitOrdered}
\alias{strsplitOrdered}
\title{Split the elements of an ordered factor vector}
\usage{
strsplitOrdered(
  x,
  split = "_",
  fixed = FALSE,
  perl = FALSE,
  useBytes = FALSE,
  sortFunc = jamba::mixedSort,
  keepOrder = TRUE,
  ...
)
}
\arguments{
\item{x}{character or factor vector.}

\item{split}{character split value sent to \code{base::strsplit()}.}

\item{fixed, perl, useBytes}{additional arguments sent to \code{base::split()}.}

\item{sortFunc}{function used to sort character values when the input
\code{x} is a character vector. The default \code{jamba::mixedSort()} applies
alphanumeric sort.}

\item{keepOrder}{logical indicating whether to keep the order of values
in the input data, for example with character input the values will
be ordered by the first appearance of each term.}

\item{...}{additional arguments are ignored.}
}
\value{
list of factor vectors, where each factor shares the same
global factor levels based upon the input data.
}
\description{
Split the elements of an ordered factor vector
}
\details{
This function performs \code{base::strsplit()} while trying to maintain
the order of factor levels in the output, based upon the order of
factor levels in the input data.
}
\examples{
# first define a vector of sample groups
iGroups <- jamba::nameVector(paste(rep(c("WT", "KO"), each=6),
   rep(c("Control", "Treated"), each=3),
   sep="_"));
iGroups <- factor(iGroups, levels=unique(iGroups));
iGroups;
strsplitOrdered(iGroups, "_");

}
\seealso{
Other jamses utilities: 
\code{\link{fold_to_log2fold}()},
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
\code{\link{sub_split_vector}()},
\code{\link{update_function_params}()},
\code{\link{update_list_elements}()}
}
\concept{jamses utilities}
