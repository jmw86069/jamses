% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jam_plot_sedesign.R
\name{plot_sedesign}
\alias{plot_sedesign}
\title{Plot sedesign object contrasts}
\usage{
plot_sedesign(
  sedesign,
  se = NULL,
  factor_names = NULL,
  factor_sep = "_",
  contrast_sep = "-",
  axis1 = NULL,
  axis2 = NULL,
  axis3 = NULL,
  axis4 = NULL,
  which_contrasts = NULL,
  contrast_style = c("comp", "contrast", "none"),
  contrast_labels = NULL,
  label_position = 0.9,
  twoway_label_position = 0.5,
  contrast_position = NULL,
  sestats = NULL,
  assay_names = NULL,
  cutoff_names = NULL,
  label_cex = 1,
  flip_twoway = FALSE,
  colorset = NULL,
  twoway_lwd = 5,
  extend_ex = 0.5,
  extend_angle = 10,
  bump_factor = 1,
  group_buffer = 0.02,
  group_border = "grey65",
  group_fill = "grey95",
  replicate_color = "grey40",
  replicate_cex = 0.8,
  do_plot = TRUE,
  plot_margins = c(0.1, 0.1, 0.1, 0.1),
  plot_type = c("base", "grid"),
  verbose = FALSE,
  debug = FALSE,
  ...
)
}
\arguments{
\item{sedesign}{\code{SEDesign} object as returned by \code{groups_to_sedesign()}.}

\item{se}{\code{SummarizedExperiment} (optional) and not yet used by this
function. In future this object may be used to assign factor level
order to factor values.}

\item{factor_names}{\code{character} vector equal to the number of delimited
values in each group name, recognized in the group names of the
design matrix of \code{sedesign} as in \code{colnames(sedesign@design)}.}

\item{factor_sep}{\code{character} string separator between factor values
in each group name, typically \code{factor_sep="_"}.}

\item{contrast_sep}{\code{character} string separator between group names
in each contrast name, typically \code{contrast_sep="-"}.}

\item{axis1, axis2, axis3, axis4}{\code{character} vectors which define the
factors to represent on each axis, with axes defined in order
\code{1=bottom}, \code{2=left}, \code{3=top}, \code{4=right}. All factors in
\code{factor_names} must be represented.}

\item{which_contrasts}{\code{numeric} index of contrasts defined in \code{sedesign},
or \code{character} vector of values present in \code{contrast_names(sedesign)}.
When a two-way contrast is defined, its component one-way
contrasts are also included.}

\item{contrast_style}{\code{character} string deciding how to format
the contrast:
\itemize{
\item \code{"comp"}: calls contrast2comp()
\item \code{"contrast"}: uses the contrast as-is
\item \code{"none"}: hides the contrast label, appending \code{contrast_labels}
when provided
}}

\item{contrast_labels}{\code{character} vector of labels named by contrast}

\item{label_position, twoway_label_position}{\code{numeric} value between 0 and 1,
which define the default position of each contrast label for one-way
and two-way contrasts, respectively.
This adjustment can be helpful to minimize overlapping labels
on the resulting figure.
\itemize{
\item \code{0} places the label toward the beginning of the arrow, which also
applies right/top justification of text at the start of the arrow.
\item \code{1} places the label at the end of the arrow, which also
applies left/bottom justification of text at the end of the arrow.
}}

\item{contrast_position}{\code{numeric} vector named by contrast, whose
values position the contrast label as described for \code{label_position}
and \code{twoway_label_position}. Any values not defined by this argument
will use the default values in \code{label_position} or
\code{twoway_label_position}.}

\item{sestats}{\code{list} object that contains element \code{"hit_array"} as
produced by \code{se_contrast_stats()}, with statistical hits for
each contrast, after applying statistical cutoffs.
When supplied, statistical hits are included in each contrast label.
The three relevant optional values used to specify specific hits:
\itemize{
\item \code{"assay_name"}: this argument defines the values from
\code{SummarizedExperiment::assays()} that were used in the contrasts.
\item \code{"cutoff_name"}: this argument defines a specific cutoff to
use, otherwise hits from any applied cutoff are included.
\item \code{"contrast_name"}: this value uses argument \code{which_contrasts}
}}

\item{assay_names, cutoff_names}{\code{character} values used with \code{sestats}
to define the statistical hits to use when \code{sestats} is supplied.}

\item{do_plot}{\code{logical} indicating whether to render the plot,
or when \code{do_plot=FALSE} only the underlying data is returned.}

\item{plot_type}{\code{character} string (experimental) to define one of
multiple plot output types:
\itemize{
\item \code{"base"} uses base R graphics.
\item \code{"grid"} uses R grid graphics, specifically the \code{vwline} package
for variable width lines.
}}

\item{verbose}{\code{logical} indicating whether to print verbose output.}

\item{...}{additional arguments are ignored.}

\item{contrast}{\code{character} contrast name}
}
\description{
Plot contrasts from sedesign object (in development), showing
one-way contrasts as block arrow, and two-way contrasts as
two block arrows connected in proper order.
}
\details{
TODO:
\itemize{
\item Mostly done: Confirm functionality with different combinations
of axis values.
\item Confirm functionality with only one factor on one axis.
\item Adjust drawing order:
\itemize{
\item Group contrasts into sets of two-way contrasts which share
any one one-way contrast with each other. This group should
be drawn together as a set, to minimize weird effects of overlaps.
\item Do not draw a one-way contrast independently when it is already
being rendered as part of a two-way contrast.
}
\item Implement method to assign colors to contrasts.
\itemize{
\item Simplest option: Allow \code{color_sub} whose names match values in
\code{contrast_names(sedesign)}.
\item Next potential option: Use \code{color_sub} to match each group name,
then define either solid color from \code{colorjam::blend_colors()}, or
using a color gradient for each one-way block arrow.
\item Two-way connectors use the first contrast end color, and
second contrast start color as a gradient.
\item If colors are not defined per group, call \code{design2colors()}?
}
\item DONE. Confirm/implement method to display fewer factors on axes
than are present in the underlying group labels.
\item Improve location of axis labels - currently uses \code{jamba::groupedAxis()}
however they appear too distant from the figure itself.
\item Determine method to "recognize" factor order from \code{colData(se)}.
\itemize{
\item Simplest option is to use argument \code{factor_names} to match
\code{colnames(colData(se))} (when supplied) and use that to define
factor order.
\item One option is to update \code{group_to_sedesign()} so it stores the
factor design data as a \code{data.frame} with proper factor level order.
This update could also benefit \code{platjam::design2colors()} by
informing the necessary \code{colData()} colnames to use.
\item A simpler option is to update \code{sedesign} to include
\code{colnames(colData(se))} as a \code{character} vector, without having to store
the full \code{data.frame}. It would requiring passing both the
\code{sedesign} and \code{se} objects together.
\item Another option is to require \code{colnames(colData(se))} to match
the order in the group names defined in \code{colnames(sedesign@design)}.
}
\item DONE. implement block arrows functions for improved quality output.
\item consider \code{grid} graphics (package \code{vwline}) or \code{ggplot2} output.
\item implement sensible method to display a subset of one-way or
two-way contrasts. For example, two-way contrasts are easier to see
when showing only a subset, perhaps only along one common axis.
\item Consider implementing gradient colors for block arrows.
\itemize{
\item This enhancement requires changing block arrow from one polygon
to a list of polygons, so each smaller polygon has its own
color from the color gradient.
}
\item Two-way contrasts:
\itemize{
\item Handle two-way contrasts for which the one-way contrasts may
not also be defined.
\item Change drawing order so the one-way block arrow label is not
overdrawn by the two-way connector.
\itemize{
\item This step probably requires grouping one-way and two-way
contrasts so that for each two-way contrast,
each one-way contrast is drawn, then the two-way connector,
then the one-way labels.
\item Probably need helper function \code{draw_twoway_contrast()}
which calls \code{draw_oneway_contrast()}, \code{draw_twoway_connector()},
and \code{draw_oneway_label()}. The one-way steps can be "skipped".
\item To be "fancy", when a one-way contrast would be rendered
multiple times, the rendering should be "skipped" and rendered
only the last time, so the label would always be rendered
after the incoming two-way connector, and so the one-way
contrast (and its label) would only need to be rendered
once overall.
}
}
}
}
\examples{
isamples_1 <- paste0(
   rep(c("DMSO", "Etop", "DMSO", "Etop"), each=6),
   "_",
   rep(c("NF", "Flag"), each=12),
   "_",
   rep(c("WT", "KO", "WT", "KO", "WT", "D955N", "WT", "D955N"), each=3),
   "_",
   LETTERS[1:3])
# simple data.frame with group information
idf <- data.frame(jamba::rbindList(strsplit(isamples_1, "_")))[,1:3]
rownames(idf) <- isamples_1;
# convert to sedesign
sedesign_1 <- groups_to_sedesign(idf)

# plot the contrasts
plot_sedesign(sedesign_1)

# re-order the factors along each axis
plot_sedesign(sedesign_1, axis1=1, axis2=3, axis3=2)

# flip the group ordering for two-way contrasts
# (These are mathematically equivalent, but shown in flipped orientation)
plot_sedesign(sedesign_1, axis1=1, axis2=3, axis3=2, flip_twoway=TRUE)

# plot only the two-way contrasts
is_twoway <- grepl("[(]", contrast_names(sedesign_1))
plot_sedesign(sedesign_1, which_contrasts=which(is_twoway),
   axis1=1, axis2=3, axis3=2)

group_names <- paste0(
   rep(c("DMSO", "Etop"), each=4),
   "_",
   rep(c("NF", "Flag"), each=2),
   "_",
   rep(c("WT", "KO"), 4))
sedesign <- groups_to_sedesign(group_names)
plot_sedesign(sedesign)

# plot only the two-way contrasts
is_twoway <- grepl("[(]", contrast_names(sedesign))
plot_sedesign(sedesign, which_contrasts=which(is_twoway),
   axis1=1, axis2=3, axis3=2)

}