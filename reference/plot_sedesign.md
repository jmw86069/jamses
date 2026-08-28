# Plot sedesign object contrasts

Plot contrasts from sedesign object (in development), showing one-way
contrasts as block arrow, and two-way contrasts as two block arrows
connected in proper order.

## Usage

``` r
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
  oneway_position = 0.9,
  twoway_position = 0.5,
  contrast_position = NULL,
  contrast_depths = NULL,
  sestats = NULL,
  sestats_style = c("label", "number", "simple label"),
  assay_names = NULL,
  cutoff_names = NULL,
  label_cex = 1,
  arrow_ex = NULL,
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
```

## Arguments

- sedesign:

  `SEDesign` object as returned by
  [`groups_to_sedesign()`](https://jmw86069.github.io/jamses/reference/groups_to_sedesign.md).

- se:

  `SummarizedExperiment` (optional) and not yet used by this function.
  In future this object may be used to assign factor level order to
  factor values.

- factor_names:

  `character` vector equal to the number of delimited values in each
  group name, recognized in the group names of the design matrix of
  `sedesign` as in `colnames(sedesign@design)`.

- factor_sep:

  `character` string separator between factor values in each group name,
  typically `factor_sep="_"`.

- contrast_sep:

  `character` string separator between group names in each contrast
  name, typically `contrast_sep="-"`.

- axis1, axis2, axis3, axis4:

  `character` vectors which define the factors to represent on each
  axis, with axes defined in order `1=bottom`, `2=left`, `3=top`,
  `4=right`. All factors in `factor_names` must be represented.

- which_contrasts:

  one of the following:

  - `numeric` index of contrasts defined in `sedesign`, or

  - `character` vector of values present in `contrast_names(sedesign)`.

  When a two-way contrast is defined, its component one-way contrasts
  are also included.

- contrast_style:

  `character` string deciding how to format the contrast:

  - `"comp"`: calls contrast2comp()

  - `"contrast"`: uses the contrast as-is

  - `"none"`: hides the contrast label, appending `contrast_labels` when
    provided

- contrast_labels:

  `character` vector of labels named by contrast

- oneway_position, twoway_position:

  `numeric` value between 0 and 1, which define the default position of
  each contrast label for one-way and two-way contrasts, respectively.
  These values are overridden by optional argument `contrast_position`
  when supplied.

  - `0` places the label toward the beginning of the arrow, which also
    applies right/top justification of text at the start of the arrow.

  - `1` places the label at the end of the arrow, which also applies
    left/bottom justification of text at the end of the arrow.

- contrast_position:

  `numeric` vector named by contrast, whose values position each
  contrast label given. Default values are defined by `oneway_position`
  and `twoway_position`, except where defined by `contrast_position`.

- contrast_depths:

  `numeric` with one or more values used to display only contrasts of
  the given depth: 1=oneway, 2=twoway. When `NULL` all contrasts are
  displayed.

- sestats:

  `SEStats` or `list` object with element `'hit_array'` produced by
  [`se_contrast_stats()`](https://jmw86069.github.io/jamses/reference/se_contrast_stats.md),
  with statistical hits for each contrast, after applying statistical
  cutoffs. When supplied, statistical hits are included in each contrast
  label. The three relevant optional values used to specify specific
  hits:

  - `'assay_name'`: this argument defines the values from
    [`SummarizedExperiment::assays()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
    that were used in the contrasts.

  - `'cutoff_name'`: this argument defines a specific cutoff to use,
    otherwise hits from any applied cutoff are included.

  - `'contrast_name'`: this value uses argument `which_contrasts`

- sestats_style:

  `character` string indicating how to present the number of hits for
  each contrast.

  - `"label"`: uses the full label: number hits (number up, number
    down)"

  - `"number"`: uses only the number of hits "number"

  - `"simple label"`: uses a simple label: "number hits" without the
    number of hits up and down.

- assay_names, cutoff_names:

  `character` values used with `sestats` to define the statistical hits
  to use when `sestats` is supplied.

- label_cex:

  `numeric` expansion factor to adjust contrast label font sizes.

- arrow_ex:

  `numeric` (default NULL) to adjust arrow width and arrow head size
  together. When `NULL` it is adjusted starting at `arrow_ex=1` and
  reduced proportional to the number of contrast bumps (parallel
  contrasts that would otherwise overlap). When provided as a numeric
  value, it is used without adjustment

- flip_twoway:

  `logical` indicating whether to flip the orientation of two-way
  contrasts, for example `"(A-B)-(C-D)"` would be flipped to equivalent
  form `"(A-C)-(B-D)"`, which will alter the orientation of the two-way
  contrast connection. Note that the individual one-way contrasts will
  be added if they did not already exist in the data.

- colorset:

  `character` vector of colors used for one-way contrasts. When the
  vector contains names, they are assigned to
  `contrast_names(sedesign)`, and any missing colors are assigned using
  [`colorjam::group2colors()`](https://jmw86069.github.io/colorjam/reference/group2colors.html).
  When the vector does not contain names, it is recycled to the number
  of contrast names.

- twoway_lwd:

  `numeric` line width for two-way contrasts, passed to
  [`draw_twoway_contrast()`](https://jmw86069.github.io/jamses/reference/draw_twoway_contrast.md).

- extend_ex:

  `numeric` expansion factor to define control points for two-way
  contrasts, beyond each one-way contrast by this fraction of the group
  width in the diagram.

- extend_angle:

  `numeric` angle in degrees to define control points for two-way
  contrasts, using this angle from the end of each one-way contrast
  toward the other one-way contrast in the set.

  - When the control point crosses the midpoint between the two
    contrasts, half the angle is used to re-define control points.

  - When the control point crosses the other contrast in the set, the
    first control point is retained, and the second control point uses
    the opposite angle so the resulting bezier curve from the first
    contrast "loops around" the far side of the second contrast, then
    connects from the opposite side.

- bump_factor:

  `numeric` factor applied to the relative amount of "bump" used to
  adjust contrasts which would otherwise overlap on the same x- or
  y-axis intercept. It can optionally accept two values, applied to the
  y-axis, then x-axis bump. Values range from 0 (no bump) to 1.5
  (expands to full width of each group), with default 1 covering roughly
  80% the size of each group box. This value is also adjusted by
  `group_buffer`.

- group_buffer:

  `numeric` value (default 0.02) indicating the relative buffer in
  between each group square, as a fraction of total width. The range is
  restricted to minimum 0 (no buffer) to 0.5 (groups are drawn as a
  point) with recommended values between 0 and 0.1,

- group_border, group_fill:

  `character` color used for the border and fill colors, respectively,
  to draw a square for each experimental group defined in `sedesign`.
  These values can be supplied as a named vector, whose names match the
  group names defined in `sedesign`, and they will be applied to each
  group. Any missing groups will used recycled values.

- replicate_color:

  `character` string with R color, used for the label in each group for
  the number of replicates as defined in `sedesign`.

- replicate_cex:

  `numeric` expansion factor used to adjust the font size for the
  replicate label in each group defined by `sedesign`.

- do_plot:

  `logical` indicating whether to render the plot, or when
  `do_plot=FALSE` only the underlying data is returned.

- plot_margins:

  `numeric` value applied to `par("mar")` around the plot, to define
  minimal whitespace around the plot.

- plot_type:

  `character` string (experimental) to define one of multiple plot
  output types:

  - `"base"` uses base R graphics.

  - `"grid"` (not yet implemented) uses R grid graphics. This option is
    expected to enable more methods to reduce overlapping labels, and
    potentially labels with markdown markup.

- verbose:

  `logical` indicating whether to print verbose output.

- debug:

  `logical` indicating whether to print very detailed debug output.

- ...:

  additional arguments are ignored.

- contrast:

  `character` contrast name

## Value

invisible `list` of `data.frame` representing individual contrasts to be
rendered. Mainly useful for reviewing the data used to produce the
figure.

## Details

TODO:

- Mostly done: Confirm functionality with different combinations of axis
  values.

- Confirm functionality with only one factor on one axis.

- Adjust drawing order:

  - Group contrasts into sets of two-way contrasts which share any one
    one-way contrast with each other. This group should be drawn
    together as a set, to minimize weird effects of overlaps.

  - Do not draw a one-way contrast independently when it is already
    being rendered as part of a two-way contrast.

- Implement method to assign colors to contrasts.

  - Simplest option: Allow `color_sub` whose names match values in
    `contrast_names(sedesign)`.

  - Next potential option: Use `color_sub` to match each group name,
    then define either solid color from
    [`colorjam::blend_colors()`](https://jmw86069.github.io/colorjam/reference/blend_colors.html),
    or using a color gradient for each one-way block arrow.

  - Two-way connectors use the first contrast end color, and second
    contrast start color as a gradient.

  - If colors are not defined per group, call `design2colors()`?

- DONE. Confirm/implement method to display fewer factors on axes than
  are present in the underlying group labels.

- Improve location of axis labels - currently uses
  [`jamba::groupedAxis()`](https://jmw86069.github.io/jamba/reference/groupedAxis.html)
  however they appear too distant from the figure itself.

- Determine method to "recognize" factor order from `colData(se)`.

  - Simplest option is to use argument `factor_names` to match
    `colnames(colData(se))` (when supplied) and use that to define
    factor order.

  - One option is to update `group_to_sedesign()` so it stores the
    factor design data as a `data.frame` with proper factor level order.
    This update could also benefit `platjam::design2colors()` by
    informing the necessary
    [`colData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
    colnames to use.

  - A simpler option is to update `sedesign` to include
    `colnames(colData(se))` as a `character` vector, without having to
    store the full `data.frame`. It would requiring passing both the
    `sedesign` and `se` objects together.

  - Another option is to require `colnames(colData(se))` to match the
    order in the group names defined in `colnames(sedesign@design)`.

- DONE. implement block arrows functions for improved quality output.

- consider `grid` graphics (package `vwline`) or `ggplot2` output.

- implement sensible method to display a subset of one-way or two-way
  contrasts. For example, two-way contrasts are easier to see when
  showing only a subset, perhaps only along one common axis.

- Consider implementing gradient colors for block arrows.

  - This enhancement requires changing block arrow from one polygon to a
    list of polygons, so each smaller polygon has its own color from the
    color gradient.

- Two-way contrasts:

  - Handle two-way contrasts for which the one-way contrasts may not
    also be defined.

  - Change drawing order so the one-way block arrow label is not
    overdrawn by the two-way connector.

    - This step probably requires grouping one-way and two-way contrasts
      so that for each two-way contrast, each one-way contrast is drawn,
      then the two-way connector, then the one-way labels.

    - Probably need helper function
      [`draw_twoway_contrast()`](https://jmw86069.github.io/jamses/reference/draw_twoway_contrast.md)
      which calls
      [`draw_oneway_contrast()`](https://jmw86069.github.io/jamses/reference/draw_oneway_contrast.md),
      `draw_twoway_connector()`, and `draw_oneway_label()`. The one-way
      steps can be "skipped".

    - To be "fancy", when a one-way contrast would be rendered multiple
      times, the rendering should be "skipped" and rendered only the
      last time, so the label would always be rendered after the
      incoming two-way connector, and so the one-way contrast (and its
      label) would only need to be rendered once overall.

## See also

Other jam experiment design:
[`SEDesign()`](https://jmw86069.github.io/jamses/reference/SEDesign.md),
[`[,SEDesign-method`](https://jmw86069.github.io/jamses/reference/sub-SEDesign-method.md),
[`check_sedesign()`](https://jmw86069.github.io/jamses/reference/check_sedesign.md),
[`contrast2comp()`](https://jmw86069.github.io/jamses/reference/contrast2comp.md),
[`contrast_colors_by_group()`](https://jmw86069.github.io/jamses/reference/contrast_colors_by_group.md),
[`contrast_names_to_sedesign()`](https://jmw86069.github.io/jamses/reference/contrast_names_to_sedesign.md),
[`contrastnames()`](https://jmw86069.github.io/jamses/reference/contrastnames.md),
[`contrasts()`](https://jmw86069.github.io/jamses/reference/contrasts.md),
[`contrasts<-()`](https://jmw86069.github.io/jamses/reference/contrasts-set.md),
[`contrasts_to_factors()`](https://jmw86069.github.io/jamses/reference/contrasts_to_factors.md),
[`contrasts_to_venn_setlists()`](https://jmw86069.github.io/jamses/reference/contrasts_to_venn_setlists.md),
[`design,SEDesign-method`](https://jmw86069.github.io/jamses/reference/design-SEDesign-method.md),
`design<-,SEDesign,matrix-method`,
[`draw_oneway_contrast()`](https://jmw86069.github.io/jamses/reference/draw_oneway_contrast.md),
[`draw_twoway_contrast()`](https://jmw86069.github.io/jamses/reference/draw_twoway_contrast.md),
[`factors()`](https://jmw86069.github.io/jamses/reference/factors.md),
[`filter_contrast_names()`](https://jmw86069.github.io/jamses/reference/filter_contrast_names.md),
[`groups()`](https://jmw86069.github.io/jamses/reference/groups.md),
[`groups_to_sedesign()`](https://jmw86069.github.io/jamses/reference/groups_to_sedesign.md),
[`plot.SEDesign()`](https://jmw86069.github.io/jamses/reference/plot.SEDesign.md),
[`print,SEDesign-method`](https://jmw86069.github.io/jamses/reference/print-SEDesign-method.md),
[`samples()`](https://jmw86069.github.io/jamses/reference/samples.md),
[`sedesign_to_factors()`](https://jmw86069.github.io/jamses/reference/sedesign_to_factors.md),
[`sort_contrasts()`](https://jmw86069.github.io/jamses/reference/sort_contrasts.md),
[`validate_sedesign()`](https://jmw86069.github.io/jamses/reference/validate_sedesign.md)

## Examples

``` r
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

```
