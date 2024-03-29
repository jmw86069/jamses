% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jam_sestats_hitarray_to_list.R
\name{hit_array_to_list}
\alias{hit_array_to_list}
\title{Quick conversion of hit_array to hit_list}
\usage{
hit_array_to_list(
  hit_array,
  contrast_names = NULL,
  cutoff_names = NULL,
  assay_names = NULL,
  ...
)
}
\arguments{
\item{hit_array}{\code{array} output from \code{se_contrast_stats()}, list element
\code{"hit_array"}.}

\item{contrast_names}{\code{character} vector of contrasts.}
}
\value{
\code{list} named by \code{contrast_names}, that contains unique statistical
hits by combining entries across the \code{cutoff_names} and \code{assay_names}
for each contrast.
}
\description{
Quick conversion of hit_array to hit_list
}
\details{
This function is mainly useful when there are multiple dimensions
unresolved in a hit_array, in which case this function will combine
hits across the different cutoffs and signals contained
in the \code{hit_array} of \code{sestats} output from \code{se_contrast_stats()}.
}
