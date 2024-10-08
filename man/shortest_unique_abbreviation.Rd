% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jam_abbrevations.R
\name{shortest_unique_abbreviation}
\alias{shortest_unique_abbreviation}
\title{Find the shortest abbrevation to retain unique values}
\usage{
shortest_unique_abbreviation(
  x,
  retain_contig_numbers = TRUE,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{x}{\code{character} vector}

\item{retain_contig_numbers}{\code{logical}, default \code{TRUE}, whether numbers
at the end of an abbreviated string should remain contiguous.
\itemize{
\item When \code{TRUE}, the goal is not to split a numeric value in the middle
of the number.
\item When \code{FALSE} the string will be abbreviated at the first position
of uniqueness.
}}

\item{...}{additional arguments are ignored.}
}
\value{
\code{character} vector named using unique values in \code{x}, and
whose values are the shortest abbreviated substrings which
maintain consistent uniqueness.
}
\description{
Find the shortest abbrevation to retain unique values
}
\details{
This function is intended to abbreviate factor levels used in
statistical contrasts to the smallest substring that uniquely
represents the unique entries provided in \code{x}.

For example, \code{c("one", "two", "three", "four")} would be converted
to \code{c("on", "tw", "th", "fo")}.

The default \code{retain_contig_numbers=TRUE} will attempt to retain
numeric values at the end of a string, to avoid splitting the number
at an intermediate position. This option only applies when the
character substring is not already unique before encountering
the numeric substring.

\if{html}{\out{<div class="sourceCode">}}\preformatted{* For this input:
}\if{html}{\out{</div>}}

\code{c("a", "p6", "p12", "p21")}
the output keeps the contiguous numbers together:
\code{c("a", "p6", "p12", "p21")}
\itemize{
\item For this input:
\code{c("a", "b6", "c12", "d21")}
only the first character is retained,
because it is already unique:
\code{c("a", "b", "c", "d")}
}
}
\section{Todo}{
\itemize{
\item Consider some method to retain contiguous numbers at the end
of a long string, while abbreviating the long string.
\itemize{
\item For this input:
\code{c("adult", "prenatal6", "prenatal12", "prenatal21")}
the ideal output would be:
\code{c("a", "p6", "p12", "p21")}
\item To be fair, I do not know how to describe this logic.
It may required breaking into words by character/non-character
breakpoints, then applying substring to each?
}
}
}

\examples{
x <- c("a", "p6", "p12", "p21");
shortest_unique_abbreviation(x)

shortest_unique_abbreviation(x, retain_contig_numbers=TRUE)

x1 <- c("male", "female");
shortest_unique_abbreviation(x1)

x2 <- c("Control", "Nicotine");
shortest_unique_abbreviation(x2)

x3 <- c("Control", "Nicotine10", "Nicotine12", "Nicotine20");
shortest_unique_abbreviation(x3)

x4 <- c("one", "two", "three", "four");
shortest_unique_abbreviation(x4)

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
\code{\link{shrinkDataFrame}()},
\code{\link{shrink_df}()},
\code{\link{shrink_matrix}()},
\code{\link{strsplitOrdered}()},
\code{\link{sub_split_vector}()},
\code{\link{update_function_params}()},
\code{\link{update_list_elements}()}
}
\concept{jamses utilities}
