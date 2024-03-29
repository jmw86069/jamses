% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jam_filter_contrast_names.R
\name{filter_contrast_names}
\alias{filter_contrast_names}
\title{Filter contrast names}
\usage{
filter_contrast_names(
  contrast_names = NULL,
  sedesign = NULL,
  factor_controls = NULL,
  apply_to_sedesign = FALSE,
  must_be_control = TRUE,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{contrast_names}{\code{character} vector of contrast names}

\item{sedesign}{\code{SEDesign} object}

\item{factor_controls}{optional \code{character} vector with one value per
experiment factor column used in the contrast names. Typically the
first observed value in each column is taken as the overall control,
using only the baseline control group.}

\item{apply_to_sedesign}{\code{logical} indicating when \code{sedesign} is supplied,
whether to return an \code{SEDesign} object with update \code{contrasts()}
consistent with the contrast_names generated by this function.}

\item{must_be_control}{\code{logical} indicating whether the control factor
level must also be the control in factor comparisons.
When \code{must-be_control=TRUE} (default) the factor comparison must
be in the form (Test-Control).
When \code{must_be_control=FALSE} it permits the control factor level
to be present as the test or control: (Control-Test) or (Test-Control).}

\item{verbose}{\code{logical} indicating whether to print verbose output.}

\item{...}{additional arguments are ignored.}
}
\description{
Filter contrast names
}
\details{
Utility function to help filter a large set of contrasts to
a smaller set of comparisons versus specific factor level
controls.

This function was motivated by a design with 5 treatments and 3 genotypes,
the "all-by-all" pairwise strategy produced 45 oneway contrasts, and
75 contrasts overall.

After filtering for comparisons versus the first level per
experimental factor, it produced 22 oneway contrasts, and 30 overall.

This function is still being tested to determine effective ways
to pare down an enormous number of potential contrasts to the
minimal set of "useful contrasts". The examples show a few use cases.
}
\examples{
group_names <- paste0(
   rep(c("UL3", "dH1A", "dH1B"), each=5), "_",
   c("Veh", "DEX", "PMA", "SF", "Ins"))
sedesign <- groups_to_sedesign(group_names)
plot_sedesign(sedesign,
   which_contrasts=jamba::unvigrep("[(]", contrast_names(sedesign)),
   arrow_ex=0.2, twoway_lwd=1, contrast_style="none")
title(main="45 total oneway and twoway contrasts")

new_contrast_names <- filter_contrast_names(contrast_names(sedesign))
new_contrast_names

# apply to sedesign directly
sedesign2 <- filter_contrast_names(sedesign=sedesign, apply_to_sedesign=TRUE)
plot_sedesign(sedesign2,
   which_contrasts=jamba::unvigrep("[(]", contrast_names(sedesign2)),
   arrow_ex=0.2, twoway_lwd=1, contrast_style="none")
title(main="22 oneway contrasts (vs UL3 / vs Veh)")

# subset for contrasts involving certain control values
# SF is the control for factors Ins,SF, so we filter for comparisons
# where SF is the control, using `must_be_control=TRUE` (default)
# but the contrast must involve SF or Ins directly
contrast_names_sf <- filter_contrast_names(
   contrast_names=jamba::vigrep("SF|Ins",
   contrast_names(sedesign)))
data.frame(jamba::rbindList(strsplit(contrast2comp(contrast_names_sf), ":")))

# must_be_control=FALSE allows any orientation involving "UL3" or "SF"
# it therefore also permits "SF-Veh"
contrast_names_sf <- filter_contrast_names(
   contrast_names=jamba::vigrep("SF|Ins", contrast_names(sedesign)),
   must_be_control=FALSE,
   factor_controls=c(Genotype="dH1A", Treatment="SF"))
data.frame(jamba::rbindList(strsplit(contrast2comp(contrast_names_sf), ":")))

# use must_be_control=c(FALSE, TRUE) to enforce only on the second factor
# 1. do not require UL3 to be the control in the comparison
# 2. do require SF to be the control in the comparison
contrast_names_sf <- filter_contrast_names(
   contrast_names=jamba::vigrep("SF|Ins", contrast_names(sedesign)),
   must_be_control=c(FALSE, TRUE),
   factor_controls=c(Genotype="dH1A", Treatment="SF"))
data.frame(jamba::rbindList(strsplit(contrast2comp(contrast_names_sf), ":")))

# for the example below, use this form
contrast_names_sf <- filter_contrast_names(jamba::vigrep("SF|Ins", contrast_names(sedesign)))
data.frame(jamba::rbindList(strsplit(contrast2comp(contrast_names_sf), ":")))

# Veh is the control for all other groups, omit Ins in these comparisons
contrast_names_veh <- filter_contrast_names(jamba::unvigrep("Ins", contrast_names(sedesign)))
data.frame(jamba::rbindList(strsplit(contrast2comp(contrast_names_veh), ":")))

use_contrasts <- unique(c(contrast_names_veh, contrast_names_sf))
sedesign3 <- sedesign;
contrast_names(sedesign3) <- use_contrasts

# show contrasts after filtering
jamba::rbindList(strsplit(contrast2comp(contrast_names(sedesign3)), ":"))

# show all one-way contrasts
plot_sedesign(sedesign3,
   which_contrasts=jamba::unvigrep("[(]", contrast_names(sedesign3)),
   arrow_ex=0.5, twoway_lwd=1, contrast_style="none")
title(main="22 oneway contrasts (filtered)")

# show only the two-way contrasts
plot_sedesign(sedesign3,
   which_contrasts=jamba::vigrep("[(]", contrast_names(sedesign3)),
   arrow_ex=0.5, twoway_lwd=1, contrast_style="none")
title(main="8 twoway contrasts (filtered)")

# same two-way contrasts, showing flipped orientation
plot_sedesign(sedesign3,
   flip_twoway=TRUE,
   which_contrasts=jamba::vigrep("[(]", contrast_names(sedesign3)),
   arrow_ex=0.5, twoway_lwd=1, contrast_style="none")
title(main="8 twoway contrasts (flipped)")

}
\seealso{
Other jam experiment design: 
\code{\link{check_sedesign}()},
\code{\link{contrast_colors_by_group}()},
\code{\link{contrasts_to_factors}()},
\code{\link{groups_to_sedesign}()},
\code{\link{plot_sedesign}()},
\code{\link{sedesign_to_factors}()}
}
\concept{jam experiment design}
