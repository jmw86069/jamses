% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jam_save_sestats.R
\name{save_sestats}
\alias{save_sestats}
\title{Save SE contrast stats output}
\usage{
save_sestats(
  sestats,
  file,
  assay_names = NULL,
  contrast_names = NULL,
  cutoff_names = NULL,
  type = c("xlsx"),
  max_nchar_sheetname = 31,
  review_output = TRUE,
  sheet_prefix = NULL,
  use_assay_suffix = TRUE,
  width_factor = 1,
  max_rows = NULL,
  colorSub = NULL,
  rename_contrasts = TRUE,
  se = NULL,
  rowData_colnames = NULL,
  row_type = "gene_name",
  verbose = FALSE,
  ...
)
}
\arguments{
\item{sestats}{\code{list} object output from \code{se_contrast_stats()}}

\item{file}{\code{character} string indicating the filename to save.}

\item{assay_names}{\code{character} string indicating which assay names
to save, stored in \code{dimnames(sestats$hit_array)$Signal}.
When \code{NULL} then all assay names are saved.}

\item{contrast_names}{\code{character} string indicating which contrasts
to save, stored in \code{dimnames(sestats$hit_array)$Contrasts}.
When \code{NULL} then all contrasts are saved.}

\item{type}{\code{character} string indicating the type of file to save.
\itemize{
\item \code{"xlsx"}: saves an Excel xlsx file using \code{jamba::writeOpenxlsx()}.
}}

\item{max_nchar_sheetname}{\code{integer} number of characters allowed in
MS Excel worksheet names, currently 31 characters.}

\item{review_output}{\code{logical} indicating whether a summary of output
should be returned as a \code{data.frame} without exporting data. This
summary will indicate all worksheets to be saved, in addition
to the sheetName for each worksheet.}

\item{sheet_prefix}{\code{character} string with optional character prefix
to use when creating worksheet names.}

\item{use_assay_suffix}{\code{logical} indicating whether to include
\code{assay_names} as suffix when forming sheet names, when there is more
than one unique assay name to be saved. This step will
attempt to abbreviate \code{assay_names} by taking up to 4 characters
from each word in the assay name, where each word is
delimited by \code{"[-.:_ ]+"}.
Otherwise, sheet names are forced to be unique by taking a substring
of the contrast name of up to \code{max_nchar_sheetname}, passing any
duplicate strings to \code{jamba::makeNames()} with suffix \code{"_v"}
followed by an integer number.}

\item{width_factor}{\code{numeric} used to adjust relative column widths
in the output Excel worksheets.}

\item{colorSub}{\code{character} vector of colors, optional, used to define
categorical background colors for text string fields in Excel.
The \code{names(colorSub)} are matched to character strings to assign
colors.}

\item{verbose}{\code{logical} indicating whether to print verbose output.}

\item{...}{additional arguments are passed to \code{jamba::writeOpenxlsx()}}
}
\description{
Save SE contrast stats output
}
\details{
This function is intended as a convenient method to export
a series of statistical tables into organized, formatted Excel
worksheets.
}
