% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jam_hmgroupcolumnlabels.R
\name{heatmap_column_group_labels}
\alias{heatmap_column_group_labels}
\title{Add Heatmap column group labels}
\usage{
heatmap_column_group_labels(
  hm_group_list,
  hm_drawn = NULL,
  se = NULL,
  add_group_label = TRUE,
  add_group_line = TRUE,
  add_group_box = FALSE,
  group_line_buffer = grid::unit(1, "mm"),
  group_line_lwd = 2,
  group_line_requires_label = TRUE,
  group_box_lwd = 2,
  group_box_outer = TRUE,
  font_cex = 1,
  hm_title_base = NULL,
  hm_body_base = NULL,
  hm_title_base_default = "centered\\nexpression_column_title_",
  hm_body_base_default = "centered\\nexpression_heatmap_body_1_",
  y_offset_lines = 0,
  endlines = 1,
  use_gridtext = TRUE,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{hm_group_list}{\code{character} or \code{list} with one of the
following types of content:
\itemize{
\item \code{character} colnames in the column data of \code{se} argument, specifically
\code{colnames(colData(se))} to define groupings.
\item \code{list} of \code{integer} vectors, where list names become labels for
each group.
}}

\item{hm_drawn}{\code{ComplexHeatmap::HeatmapList} object which is returned
by the function \code{ComplexHeatmap::draw()}. This object is a pointer
to the grid graphical elements required.}

\item{se}{\code{SummarizedExperiment} object, required only when \code{hm_group_list}
is supplied in the form of colnames of the
\code{SummarizedExperiment::colData(se)}.}

\item{add_group_label}{\code{logical} indicating whether to draw the group
\code{character} label above each group defined in \code{hm_group_list}.}

\item{add_group_line}{\code{logical} indicating whether to draw the group
line for each group defined in \code{hm_group_list}. The line is intended
to appear below the text label when \code{add_group_label=TRUE}.}

\item{add_group_box}{\code{logical} indicating whether to draw a box around
the heatmap region defined by each group in \code{hm_group_list}.}

\item{group_line_buffer}{\code{grid::unit} object indicating the whitespace
buffer region betwen adjacent lines when \code{add_group_line=TRUE}.
The default enforces 1-mm white space between lines so they
do not touch.}

\item{group_line_lwd}{\code{numeric} line width when \code{add_group_line=TRUE}.}

\item{group_line_requires_label}{\code{logical} whether to require
the associated group label to contain non-whitespace visible
text.
\itemize{
\item \code{group_line_requires_label=TRUE} (default) requires group label
to have visible characters, which means the presence of an
empty labels will cause the group line not to be drawn.
\item \code{group_line_requires_label=FALSE} does not require a group label,
therefore all group lines are drawn.
}}

\item{group_box_lwd}{\code{numeric} line width when \code{add_group_box=TRUE}.}

\item{group_box_outer}{\code{logical} indicating whether to draw the group
box as an outer border, which means the border will be drawn only
on the very outer edge of the heatmap body, and will not overlap
the heatmap contents. This option is only active when
\code{add_group_box=TRUE}.}

\item{font_cex}{\code{numeric} adjustment for font sizes overall.}

\item{hm_title_base}{\code{character} string used to search for matching
heatmap column title grid layout regions. It is derived from the
\code{Heatmap} argument \code{name}, usually followed by \verb{"_column_title_". When }NULL\verb{the default values are defined by}detect_heatmap_components()\code{element}column_title_components`.}

\item{hm_body_base}{\code{character} string used to search for matching
heatmap body grid layout regions.  It is derived from the
\code{Heatmap} argument \code{name}, usually followed by \verb{"_body_1_". When }NULL\verb{the default values are defined by}detect_heatmap_components()\code{element}heatmap_body_components`.
Note: Currently the group box only includes the first row_split,
although it is intended to include the entire set of heatmap
rows in future iterations.}

\item{hm_title_base_default, hm_body_base_default}{\code{character} string
used as reference for hm_title_base,hm_body_base.}

\item{y_offset_lines}{\code{numeric} adjustment used to shift the text
label and underline by this many lines of character height.
It is mainly used internally for iterative calls, when \code{hm_group_list}
includes multiple layers of groups. The text and underline are
intended to fill 2 character height lines, although this argument
can be used to make manual adjustments.}

\item{endlines}{\code{integer} number of blank lines to append to the end
of each group label, which has the effect of shifting the group
label upwards slightly.}

\item{use_gridtext}{\code{logical} (default TRUE) whether to render text using
\code{gridtext::richtext_grob()} in order to enable markdown and
limited HTML-based formatting.}

\item{verbose}{\code{logical} indicating whether to print verbose output.}

\item{...}{additional arguments are ignored.}
}
\description{
Add Heatmap column group labels, specifically for ComplexHeatmap output.
}
\details{
This function is currently experimental, and is intended only for
a specific scenario, to augment a \code{ComplexHeatmap::Heatmap} object,
as produced by \code{heatmap_se()}, that used the argument \code{column_split}
to sub-divide the heatmap columns into subgroups.
This function draws labels above the heatmap to describe group
factor values associated with column splits.

When there are multiple layers of grouping, this function will
also draw multiple layers. For example, when columns are split
by \code{column_split=c("Treatment", "Time")}, it will produce a heatmap
where each column slice has one combination of Treatment and Time.
This function can be used to add a layer of labels by \code{"Time"},
and a layer of labels by \code{"Treatment"}. The labels are shown by default
with a broad underline to indicate contiguous column slices that
contain the same \code{"Treatment"} or \code{"Time"} values.

The output can provide a cleaner visualization than the alternative of
displaying colorized annotation boxes at the top of the heatmap.
\subsection{Input Formats}{

There are two strategies for defining the column group data to display:
\subsection{Input using colnames from colData(se)}{

This option is recommended as the "easiest" method. It requires:
\itemize{
\item \code{hm_group_list} provides a \code{character} vector of \code{colnames(colData(se))}.
Note that the columns are applied bottom-to-top, so it is sometimes
helpful to supply columns in reverse order for \code{column_split}.
\item \code{se} must be provided, since it supplies \code{colData(se)}
}

Each value in \code{hm_group_list} is applied from bottom-to-top, to define
a row of labels. By default, \code{add_group_line=TRUE}, so labels are also
underlined to indicate samples included in each group.
}

\subsection{Input using a named list}{

This option is recommended when labeling a subset of column slices,
or when the column slice groups need to be customized in some way.
}

}

\subsection{Additional requirements}{
\itemize{
\item The Heatmap \code{column_title} should (usually) be empty, so that no text
labels are drawn which may overlap the labels drawn by this function.
Use \code{column_title=" "} (with a whitespace character) to prevent
\code{ComplexHeatmap::Heatmap()} from using its internal default text label.
\item There usually needs to be whitespace above the heatmap, which
can be accomplished when drawing a Heatmap object \code{hm} like this:
\code{ComplexHeatmap::draw(hm, column_title="\\n\\n\\n\\n")}
\item When also displaying a heatmap title as defined by \code{heatmap_se()},
the adjustment can be done like this:
}

\if{html}{\out{<div class="sourceCode r">}}\preformatted{ComplexHeatmap::draw(hm,
   column_title=paste0(attr(hm, "hm_title"),
      "\\n\\n\\n\\n"))
}\if{html}{\out{</div>}}
\itemize{
\item The arguments \code{hm_title_base} and \code{hm_body_base} should match the
heatmap name, which is defined in \code{heatmap_se()} with \code{data_type},
usually prefixed \code{"centered\\n"} when the data is centered by that
function. For example, when \code{data_type="abundance"}, the corresponding
argument value should be
\code{hm_title_base="centered\\nabundance_column_title_"}.
\item use \code{detect_heatmap_components()} to help identify the appropriate
\code{grid} elements available to be used by this function.
}

When the heatmap contains \code{"column_title"} elements defined in
\code{ComplexHeatmap::list_components()}, and there is no element
\code{"global_column_title"}, the \code{y_offset_lines} is adjusted down so that
the position is inside the column_title region, typically below the
column_title when using trailing whitespace (see argument
\code{hm_title_buffer} in \code{heatmap_se()}). However, when element
\code{"global_column_title"} is present, the \code{y_offset_lines} are not adjusted,
so that the position is above the column_title region. In this case, to
position the labels below the column_title, you can adjust
\code{y_offset_lines} down manually like this: \code{y_offset_lines=-9}.
}

\subsection{Todo}{
\itemize{
\item Automate determining the column_title grid layout name.
\itemize{
\item Specifically when only one \code{heatmapname_column_title_1} is present,
but there are multiple column groups, it should take the x-axis
coordinate values (left and right boundaries) from the heatmap body
instead of the column_title region.
}
\item When \code{add_group_box=TRUE}, and \code{row_split} indicates multiple row
groups, it should calculate the y-axis coordinate values (top and bottom
boundaries) using the full set of row groups.
\item Enable blank annotations, either by passing a subset \code{se}, or by
annotations with no associated label.
}
}
}
