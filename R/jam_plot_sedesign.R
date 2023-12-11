
#' Plot sedesign object contrasts
#'
#' Plot contrasts from sedesign object (in development), showing
#' one-way contrasts as block arrow, and two-way contrasts as
#' two block arrows connected in proper order.
#'
#' TODO:
#' * Mostly done: Confirm functionality with different combinations
#' of axis values.
#' * Confirm functionality with only one factor on one axis.
#' * Adjust drawing order:
#'
#'    * Group contrasts into sets of two-way contrasts which share
#'    any one one-way contrast with each other. This group should
#'    be drawn together as a set, to minimize weird effects of overlaps.
#'    * Do not draw a one-way contrast independently when it is already
#'    being rendered as part of a two-way contrast.
#'
#' * Implement method to assign colors to contrasts.
#'
#'    * Simplest option: Allow `color_sub` whose names match values in
#'    `contrast_names(sedesign)`.
#'    * Next potential option: Use `color_sub` to match each group name,
#'    then define either solid color from `colorjam::blend_colors()`, or
#'    using a color gradient for each one-way block arrow.
#'    * Two-way connectors use the first contrast end color, and
#'    second contrast start color as a gradient.
#'    * If colors are not defined per group, call `design2colors()`?
#'
#' * DONE. Confirm/implement method to display fewer factors on axes
#' than are present in the underlying group labels.
#' * Improve location of axis labels - currently uses `jamba::groupedAxis()`
#' however they appear too distant from the figure itself.
#' * Determine method to "recognize" factor order from `colData(se)`.
#'
#'    * Simplest option is to use argument `factor_names` to match
#'    `colnames(colData(se))` (when supplied) and use that to define
#'    factor order.
#'    * One option is to update `group_to_sedesign()` so it stores the
#'    factor design data as a `data.frame` with proper factor level order.
#'    This update could also benefit `platjam::design2colors()` by
#'    informing the necessary `colData()` colnames to use.
#'    * A simpler option is to update `sedesign` to include
#'    `colnames(colData(se))` as a `character` vector, without having to store
#'    the full `data.frame`. It would requiring passing both the
#'    `sedesign` and `se` objects together.
#'    * Another option is to require `colnames(colData(se))` to match
#'    the order in the group names defined in `colnames(sedesign@design)`.
#'
#' * DONE. implement block arrows functions for improved quality output.
#' * consider `grid` graphics (package `vwline`) or `ggplot2` output.
#' * implement sensible method to display a subset of one-way or
#' two-way contrasts. For example, two-way contrasts are easier to see
#' when showing only a subset, perhaps only along one common axis.
#' * Consider implementing gradient colors for block arrows.
#'
#'    * This enhancement requires changing block arrow from one polygon
#'    to a list of polygons, so each smaller polygon has its own
#'    color from the color gradient.
#' * Two-way contrasts:
#'
#'    * Handle two-way contrasts for which the one-way contrasts may
#'    not also be defined.
#'    * Change drawing order so the one-way block arrow label is not
#'    overdrawn by the two-way connector.
#'
#'       * This step probably requires grouping one-way and two-way
#'       contrasts so that for each two-way contrast,
#'       each one-way contrast is drawn, then the two-way connector,
#'       then the one-way labels.
#'       * Probably need helper function `draw_twoway_contrast()`
#'       which calls `draw_oneway_contrast()`, `draw_twoway_connector()`,
#'       and `draw_oneway_label()`. The one-way steps can be "skipped".
#'       * To be "fancy", when a one-way contrast would be rendered
#'       multiple times, the rendering should be "skipped" and rendered
#'       only the last time, so the label would always be rendered
#'       after the incoming two-way connector, and so the one-way
#'       contrast (and its label) would only need to be rendered
#'       once overall.
#'
#' @family jam experiment design
#'
#' @returns invisible `list` of `data.frame` representing individual
#'    contrasts to be rendered. Mainly useful for reviewing the
#'    data used to produce the figure.
#'
#' @param sedesign `SEDesign` object as returned by `groups_to_sedesign()`.
#' @param se `SummarizedExperiment` (optional) and not yet used by this
#'    function. In future this object may be used to assign factor level
#'    order to factor values.
#' @param factor_names `character` vector equal to the number of delimited
#'    values in each group name, recognized in the group names of the
#'    design matrix of `sedesign` as in `colnames(sedesign@design)`.
#' @param factor_sep `character` string separator between factor values
#'    in each group name, typically `factor_sep="_"`.
#' @param contrast_sep `character` string separator between group names
#'    in each contrast name, typically `contrast_sep="-"`.
#' @param axis1,axis2,axis3,axis4 `character` vectors which define the
#'    factors to represent on each axis, with axes defined in order
#'    `1=bottom`, `2=left`, `3=top`, `4=right`. All factors in
#'    `factor_names` must be represented.
#' @param which_contrasts one of the following:
#'    * `numeric` index of contrasts defined in `sedesign`, or
#'    * `character` vector of values present in `contrast_names(sedesign)`.
#'
#'    When a two-way contrast is defined, its component one-way
#'    contrasts are also included.
#' @param contrast_style `character` string indicating how to format
#'    contrast labels displayed on each one-way contrast:
#'    * `"comp"`: the abbreviated 'comp' from `contrast2comp()`
#'    * `"contrast"`: the full contrast name
#'    * `"none"`: hides the contrast name
#' @param contrast_labels `character` vector of optional custom labels
#'    to append to each contrast name. The `names(contrast_labels)`
#'    are expected to match contrast names obtained
#'    from `contrast_names(sedesign)`.
#'    Note: When `sestats` is supplied, `contrast_labels` are populated
#'    with the number of statistical hits for each contrast.
#' @param oneway_position,twoway_position `numeric` value between 0 and 1,
#'    which define the default position of each contrast label for one-way
#'    and two-way contrasts, respectively. These values are overridden
#'    by optional argument `contrast_position` when supplied.
#'    * `0` places the label toward the beginning of the arrow, which also
#'    applies right/top justification of text at the start of the arrow.
#'    * `1` places the label at the end of the arrow, which also
#'    applies left/bottom justification of text at the end of the arrow.
#' @param contrast_position `numeric` vector named by contrast, whose
#'    values position each contrast label given. Default values are
#'    defined by `oneway_position` and `twoway_position`, except
#'    where defined by `contrast_position`.
#' @param contrast_depths `numeric` with one or more values used to
#'    display only contrasts of the given depth: 1=oneway, 2=twoway.
#'    When `NULL` all contrasts are displayed.
#' @param sestats `list` object that contains element `"hit_array"` as
#'    produced by `se_contrast_stats()`, with statistical hits for
#'    each contrast, after applying statistical cutoffs.
#'    When supplied, statistical hits are included in each contrast label.
#'    The three relevant optional values used to specify specific hits:
#'    * `"assay_name"`: this argument defines the values from
#'    `SummarizedExperiment::assays()` that were used in the contrasts.
#'    * `"cutoff_name"`: this argument defines a specific cutoff to
#'    use, otherwise hits from any applied cutoff are included.
#'    * `"contrast_name"`: this value uses argument `which_contrasts`
#' @param sestats_style `character` string indicating how to present the
#'    number of hits for each contrast.
#'    * `"label"`: uses the full label: number hits (number up, number down)"
#'    * `"number"`: uses only the number of hits "number"
#'    * `"simple label"`: uses a simple label: "number hits" without
#'    the number of hits up and down.
#' @param assay_names,cutoff_names `character` values used with `sestats`
#'    to define the statistical hits to use when `sestats` is supplied.
#' @param label_cex `numeric` expansion factor to adjust contrast label
#'    font sizes.
#' @param arrow_ex `numeric` (default NULL) to adjust arrow width and
#'    arrow head size together. When `NULL` it is adjusted starting at
#'    `arrow_ex=1` and reduced proportional to the number of contrast
#'    bumps (parallel contrasts that would otherwise overlap). When
#'    provided as a numeric value, it is used without adjustment
#' @param flip_twoway `logical` indicating whether to flip the orientation
#'    of two-way contrasts, for example `"(A-B)-(C-D)"` would be flipped
#'    to equivalent form `"(A-C)-(B-D)"`, which will alter the orientation
#'    of the two-way contrast connection. Note that the individual one-way
#'    contrasts will be added if they did not already exist in the data.
#' @param colorset `character` vector of colors used for one-way contrasts.
#'    When the vector contains names, they are assigned to
#'    `contrast_names(sedesign)`, and any missing colors are assigned
#'    using `colorjam::group2colors()`.
#'    When the vector does not contain names, it is recycled to
#'    the number of contrast names.
#' @param twoway_lwd `numeric` line width for two-way contrasts,
#'    passed to `draw_twoway_contrast()`.
#' @param extend_ex `numeric` expansion factor to define control points
#'    for two-way contrasts, beyond each one-way contrast by this fraction
#'    of the group width in the diagram.
#' @param extend_angle `numeric` angle in degrees to define control points
#'    for two-way contrasts, using this angle from the end of each one-way
#'    contrast toward the other one-way contrast in the set.
#'    * When the control point crosses the midpoint between the two contrasts,
#'    half the angle is used to re-define control points.
#'    * When the control point crosses the other contrast in the set, the
#'    first control point is retained, and the second control point
#'    uses the opposite angle so the resulting bezier curve
#'    from the first contrast "loops around" the far side of the second
#'    contrast, then connects from the opposite side.
#' @param bump_factor `numeric` factor applied to the relative
#'    amount of "bump" used to adjust contrasts which would otherwise
#'    overlap on the same x- or y-axis intercept. It can optionally accept
#'    two values, applied to the y-axis, then x-axis bump. Values
#'    range from 0 (no bump) to 1.5 (expands to full width of each group),
#'    with default 1 covering roughly 80% the size of each group box.
#'    This value is also adjusted by `group_buffer`.
#' @param group_buffer `numeric` value (default 0.02) indicating the
#'    relative buffer in between each group square, as a fraction of
#'    total width. The range is restricted to minimum 0 (no buffer) to
#'    0.5 (groups are drawn as a point) with recommended values between 0
#'    and 0.1,
#' @param group_border,group_fill `character` color used for the border
#'    and fill colors, respectively, to draw a square
#'    for each experimental group defined in `sedesign`.
#'    These values can be supplied as a named vector, whose names match
#'    the group names defined in `sedesign`, and they will be applied
#'    to each group. Any missing groups will used recycled values.
#' @param replicate_color `character` string with R color, used for
#'    the label in each group for the number of replicates as defined
#'    in `sedesign`.
#' @param replicate_cex `numeric` expansion factor used to adjust the font
#'    size for the replicate label in each group defined by `sedesign`.
#' @param do_plot `logical` indicating whether to render the plot,
#'    or when `do_plot=FALSE` only the underlying data is returned.
#' @param plot_margins `numeric` value applied to `par("mar")` around the
#'    plot, to define minimal whitespace around the plot.
#' @param plot_type `character` string (experimental) to define one of
#'    multiple plot output types:
#'    * `"base"` uses base R graphics.
#'    * `"grid"` (not yet implemented) uses R grid graphics. This option
#'    is expected to enable more methods to reduce overlapping labels,
#'    and potentially labels with markdown markup.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param debug `logical` indicating whether to print very detailed
#'    debug output.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' isamples_1 <- paste0(
#'    rep(c("DMSO", "Etop", "DMSO", "Etop"), each=6),
#'    "_",
#'    rep(c("NF", "Flag"), each=12),
#'    "_",
#'    rep(c("WT", "KO", "WT", "KO", "WT", "D955N", "WT", "D955N"), each=3),
#'    "_",
#'    LETTERS[1:3])
#' # simple data.frame with group information
#' idf <- data.frame(jamba::rbindList(strsplit(isamples_1, "_")))[,1:3]
#' rownames(idf) <- isamples_1;
#' # convert to sedesign
#' sedesign_1 <- groups_to_sedesign(idf)
#'
#' # plot the contrasts
#' plot_sedesign(sedesign_1)
#'
#' # re-order the factors along each axis
#' plot_sedesign(sedesign_1, axis1=1, axis2=3, axis3=2)
#'
#' # flip the group ordering for two-way contrasts
#' # (These are mathematically equivalent, but shown in flipped orientation)
#' plot_sedesign(sedesign_1, axis1=1, axis2=3, axis3=2, flip_twoway=TRUE)
#'
#' # plot only the two-way contrasts
#' is_twoway <- grepl("[(]", contrast_names(sedesign_1))
#' plot_sedesign(sedesign_1, which_contrasts=which(is_twoway),
#'    axis1=1, axis2=3, axis3=2)
#'
#' group_names <- paste0(
#'    rep(c("DMSO", "Etop"), each=4),
#'    "_",
#'    rep(c("NF", "Flag"), each=2),
#'    "_",
#'    rep(c("WT", "KO"), 4))
#' sedesign <- groups_to_sedesign(group_names)
#' plot_sedesign(sedesign)
#'
#' # plot only the two-way contrasts
#' is_twoway <- grepl("[(]", contrast_names(sedesign))
#' plot_sedesign(sedesign, which_contrasts=which(is_twoway),
#'    axis1=1, axis2=3, axis3=2)
#'
#' @export
plot_sedesign <- function
(sedesign,
 se=NULL,
 factor_names=NULL,
 factor_sep="_",
 contrast_sep="-",
 axis1=NULL,
 axis2=NULL,
 axis3=NULL,
 axis4=NULL,
 which_contrasts=NULL,
 contrast_style=c("comp",
    "contrast",
    "none"),
 contrast_labels=NULL,
 oneway_position=0.9,
 twoway_position=0.5,
 contrast_position=NULL,
 contrast_depths=NULL,
 sestats=NULL,
 sestats_style=c("label",
    "number",
    "simple label"),
 assay_names=NULL,
 cutoff_names=NULL,
 label_cex=1,
 arrow_ex=NULL,
 flip_twoway=FALSE,
 colorset=NULL,
 twoway_lwd=5,
 extend_ex=0.5,
 extend_angle=10,
 bump_factor=1,
 group_buffer=0.02,
 group_border="grey65",
 group_fill="grey95",
 replicate_color="grey40",
 replicate_cex=0.8,
 do_plot=TRUE,
 plot_margins=c(0.1, 0.1, 0.1, 0.1),
 plot_type=c("base",
    "grid"),
 verbose=FALSE,
 debug=FALSE,
 ...)
{
   plot_type <- match.arg(plot_type);
   contrast_style <- match.arg(contrast_style);
   sestats_style <- match.arg(sestats_style);

   # convert group names to factor summary table
   # TODO: call sedesign_to_factors()
   group_names <- groups(sedesign);
   factors_df <- data.frame(check.names=FALSE,
      stringsAsFactors=FALSE,
      jamba::rbindList(strsplit(group_names, factor_sep)))

   # determine group sizes
   im2list_internal <- function
   (x,
    empty=0,
    ...)
   {
      # the reciprocal of list2im()
      x_rows <- rownames(x);
      x_cols <- colnames(x);
      l <- lapply(jamba::nameVector(x_cols), function(i){
         i_empty <- as(empty, class(x[,i]));
         has_value <- (!x[,i] %in% i_empty);
         x_rows[has_value];
      });
      return(l);
   }
   groups_list <- im2list_internal(design(sedesign))
   groups_n <- lengths(groups_list);
   if (verbose) {
      jamba::printDebug("plot_sedesign(): ",
         "groups_n:");
      print(groups_n);
   }

   # factor_names
   if (length(factor_names) == 0) {
      # assign temporary factor_names
      factor_names <- jamba::makeNames(
         rep("factor", ncol(factors_df)),
         suffix="");
   } else {
      # assign given factor_names
      factor_names <- jamba::makeNames(
         rep(factor_names, length.out=ncol(factors_df)),
         startN=2,
         renameFirst=FALSE)
   }
   colnames(factors_df) <- factor_names;
   rownames(factors_df) <- jamba::pasteByRow(factors_df);

   # attempt to recognize colnames(colData(se))
   if (length(se) > 0 && inherits(se, what="SummarizedExperiment")) {
      colData_colnames <- colnames(SummarizedExperiment::colData(se));
      use_colData_colnames <- intersect(factor_names, colData_colnames);
      # if any colnames match, use factor levels from colData(se)
      if (length(use_colData_colnames) > 0) {
         for (icol in use_colData_colnames) {
            if (is.factor(SummarizedExperiment::colData(se)[[icol]])) {
               # use defined factor level order,
               # but allow additional values if present
               use_levels <- unique(c(
                  levels(SummarizedExperiment::colData(se)[[icol]]),
                  factors_df[[icol]]));
            } else {
               # use values in the order they appear in group names,
               # include additional values from colData(se) if present
               use_levels <- unique(c(
                  factors_df[[icol]],
                  unique(SummarizedExperiment::colData(se)[[icol]])));
            }
            factors_df[[icol]] <- factor(factors_df[[icol]],
               levels=use_levels)
         }
      }
   }

   # Convert remaining non-factor columns to factor
   # by coercing to character, then assigning levels in the order they appear
   for (icol in factor_names) {
      if (!is.factor(factors_df[[icol]])) {
         use_levels <- unique(as.character(factors_df[[icol]]));
         factors_df[[icol]] <- factor(as.character(factors_df[[icol]]),
            levels=use_levels)
      }
   }

   # add group sizes to factors_df
   group_match <- match(rownames(factors_df), names(groups_n));
   factors_df$n <- jamba::rmNA(groups_n[group_match],
      naValue=0)

   # now all columns in factors_df should be factors with levels
   if (verbose) {
      jamba::printDebug("plot_sedesign(): ",
         "factors_df:");
      print(factors_df);
   }

   # define axis factor values
   # custom function to apply logic
   handle_axis_values <- function
   (axisN, axis_values, factor_names, assign_default=FALSE)
   {
      #
      if (length(axisN) > 0) {
         if (is.numeric(axisN)) {
            axisN <- factor_names[axisN]
         } else {
            axisN <- intersect(axisN, factor_names)
         }
      }
      if (TRUE %in% assign_default && length(axisN) == 0) {
         axisN <- head(setdiff(factor_names, axis_values), 1)
      }
      axisN <- setdiff(axisN, c(axis_values, NA));
      return(axisN)
   }
   axis_values <- character(0);
   # if all are empty, define some suitable defaults
   assign_default <- FALSE;
   if (length(c(axis1, axis2, axis3, axis3)) == 0) {
      assign_default <- TRUE;
   }
   axis1 <- handle_axis_values(axis1,
      axis_values,
      factor_names,
      assign_default)
   axis_values <- c(axis_values, axis1);
   axis2 <- handle_axis_values(axis2,
      axis_values,
      factor_names,
      assign_default)
   axis_values <- c(axis_values, axis2);
   axis3 <- handle_axis_values(axis3,
      axis_values,
      factor_names,
      assign_default)
   axis_values <- c(axis_values, axis3);
   axis4 <- handle_axis_values(axis4,
      axis_values,
      factor_names,
      assign_default)
   axis_values <- c(axis_values, axis4);

   axis_values <- c(axis1, axis2, axis3, axis4)
   remaining_names <- setdiff(factor_names,
      axis_values)
   # add convenient label to align with axis labels
   factors_df$label <- jamba::pasteByRowOrdered(
      factors_df[, axis_values, drop=FALSE])
   # sort factors which helps axis1,axis3 ordering, not axis2,axis4
   factors_df <- jamba::mixedSortDF(factors_df,
      byCols=c("label", factor_names))


   if (verbose) {
      jamba::printDebug("plot_sedesign(): ",
         "axis_values:",
         axis_values);
      if (length(remaining_names) > 0) {
         jamba::printDebug("plot_sedesign(): ",
            "remaining_names:",
            remaining_names);
      }
   }

   # assign group names with factor orders
   # group_labels <- jamba::pasteByRowOrdered(
   #    factors_df[,c(axis_values, remaining_names), drop=FALSE])
   # rownames(factors_df) <- as.character(group_labels);
   if (verbose) {
      # jamba::printDebug("plot_sedesign(): ",
      #    "group_labels:",
      #    group_labels);
      jamba::printDebug("plot_sedesign(): ",
         "factors_df:");
      print(factors_df);
   }

   axis1values <- jamba::pasteByRowOrdered(factors_df[, axis1, drop=FALSE])
   axis2values <- jamba::pasteByRowOrdered(factors_df[, axis2, drop=FALSE])
   axis3values <- jamba::pasteByRowOrdered(factors_df[, axis3, drop=FALSE])
   axis4values <- jamba::pasteByRowOrdered(factors_df[, axis4, drop=FALSE])
   remaining_values <- jamba::pasteByRowOrdered(
      factors_df[, remaining_names, drop=FALSE])
   # optionally print verbose summary
   if (verbose) {
      jamba::printDebug("plot_sedesign(): ",
         "axis1values:");
      print(axis1values);
      jamba::printDebug("plot_sedesign(): ",
         "axis2values:");
      print(axis2values);
      jamba::printDebug("plot_sedesign(): ",
         "axis3values:");
      print(axis3values);
      jamba::printDebug("plot_sedesign(): ",
         "axis4values:");
      print(axis4values);
      jamba::printDebug("plot_sedesign(): ",
         "remaining_values:");
      print(remaining_values);
   }

   # bottom/top axis
   axis13df <- jamba::mixedSortDF(
      unique(as.data.frame(jamba::rmNULL(nullValue=NA,
         list(axis1=axis1values, axis3=axis3values)))))
   axis13df[,"x_coord"] <- seq_len(nrow(axis13df));
   # left/right axis
   axis24df <- jamba::mixedSortDF(
      unique(as.data.frame(jamba::rmNULL(nullValue=NA,
         list(axis2=axis2values, axis4=axis4values)))))
   axis24df[,"y_coord"] <- seq_len(nrow(axis24df));

   axis_df <- data.frame(check.names=FALSE,
      stringsAsFactors=FALSE,
      row.names=NULL,
      axis13df,
      axis24df[rep(axis24df$y_coord, each=nrow(axis13df)), , drop=FALSE])
   # if (length(remaining_values) > 0) {
   #    axis_df$remaining_names <- remaining_values;
   # }
   use_axis_colnames <- intersect(c(paste0("axis", 1:4), "remaining_values"),
      colnames(axis_df));
   axis_df$label <- jamba::pasteByRow(
      axis_df[, c("axis1", "axis2", "axis3", "axis4"), drop=FALSE]);

   # merge group counts into axis_df
   axis_match <- match(axis_df$label, factors_df$label);
   axis_df$n <- factors_df$n[axis_match];

   # add group name to axis_df for convenience
   fmatch <- match(as.character(axis_df$label),
      as.character(factors_df$label))
   axis_df$group <- rownames(factors_df)[fmatch];
   if (names(group_fill) == 0 ||
         !all(axis_df$group %in% names(group_fill))) {
      missing_groups <- setdiff(axis_df$group, names(group_fill))
      group_fill[missing_groups] <- rep(group_fill,
         length.out=length(missing_groups));
      group_fill <- rep(group_fill, length.out=length(axis_df$group));
      names(group_fill) <- axis_df$group;
   }
   if (names(group_border) == 0 ||
         !all(axis_df$group %in% names(group_border))) {
      missing_groups <- setdiff(axis_df$group, names(group_border))
      group_border[missing_groups] <- rep(group_border,
         length.out=length(missing_groups));
      group_border <- rep(group_border, length.out=length(axis_df$group));
      names(group_border) <- axis_df$group;
   }
   # jamba::printDebug("factors_df:");print(factors_df);# debug
   # jamba::printDebug("axis_df:");print(axis_df);# debug

   # layout is defined here
   if (verbose) {
      jamba::printDebug("plot_sedesign(): ",
         "axis_df:");
      print(axis_df);
   }

   # define contrasts
   contrast_names <- unique(contrast_names(sedesign))
   contrast_numbers <- jamba::nameVector(seq_along(contrast_names),
      contrast_names);
   if (length(which_contrasts) > 0) {
      if (is.logical(which_contrasts)) {
         contrast_names <- contrast_names[which_contrasts]
      } else if (is.numeric(which_contrasts)) {
         contrast_names <- jamba::rmNA(contrast_names[which_contrasts])
         contrast_numbers <- jamba::rmNA(contrast_numbers[which_contrasts])
         # jamba::printDebug("contrast_numbers: ");print(contrast_numbers);# debug
         if (length(contrast_names) == 0) {
            stop("Numeric which_contrasts did not include any contrast_names.")
         }
      } else {
         contrast_names <- intersect(as.character(which_contrasts),
            contrast_names)
         if (length(contrast_names) == 0) {
            stop("Character which_contrasts did not match any contrast_names.")
         }
      }
   }

   # define contrast_labels when sestats is supplied
   if (length(sestats) > 0) {
      if ("list" %in% class(sestats) && "hit_array" %in% names(sestats)) {
         hit_list <- hit_array_to_list(sestats,
            assay_names=assay_names,
            cutoff_names=cutoff_names)
         contrast_labels_hits <- NULL;
         if (sestats_style %in% c("label", "simple label")) {
            contrast_labels_hits <- format_hits(hits=hit_list,
               style="text");
            if (sestats_style %in% "simple label") {
               # remove directional part of the label
               contrast_labels_hits <- jamba::nameVector(
                  gsub(" [(].+", "",
                     contrast_labels_hits),
                  names(contrast_labels_hits));
            }
         } else if (sestats_style %in% "number") {
            # use only the number of hits
            contrast_labels_hits <- jamba::nameVector(
               jamba::formatInt(lengths(hit_list)),
               names(hit_list));
         }
         if (length(contrast_labels) == 0) {
            contrast_labels <- contrast_labels_hits;
         } else {
            # combine custom and contrast labels into a list
            all_labels <- jamba::nameVector(unique(c(
               names(contrast_labels),
               names(contrast_labels_hits))));
            contrast_labels <- lapply(all_labels, function(iname){
               c(contrast_labels[[iname]],
                  contrast_labels_hits[[iname]])
            })
         }
      }
   }

   # handle contrast_position
   if (length(oneway_position) == 0) {
      oneway_position <- 0.5;
   }
   if (length(twoway_position) == 0) {
      twoway_position <- 0.5;
   }
   is_twoway <- grepl("[(]", contrast_names)
   if (length(contrast_position) == 0) {
      if (verbose) {
         jamba::printDebug("plot_sedesign(): ",
            "Applying one-way and two-way label position defaults.")
      }
      contrast_position <- jamba::nameVector(
         ifelse(is_twoway,
            twoway_position,
            oneway_position),
         contrast_names)
   } else if (length(names(contrast_position)) == 0) {
      if (verbose) {
         jamba::printDebug("plot_sedesign(): ",
            "Expanding ",
            "contrast_position",
            " as provided to all contrasts.")
      }
      contrast_position <- jamba::nameVector(
         rep(contrast_position,
            length.out=length(contrast_names)),
         contrast_names)
   } else {
      if (verbose) {
         jamba::printDebug("plot_sedesign(): ",
            "Applying ",
            "contrast_position",
            " to matching contrast_names.")
      }
      # define 0.5 for all contrasts
      use_contrast_position <- jamba::nameVector(
         ifelse(is_twoway,
            twoway_position,
            oneway_position),
         contrast_names)
      # substitute custom values where defined
      use_names <- setdiff(names(use_contrast_position),
         names(contrast_position))
      if (length(use_names) > 0) {
         contrast_position[use_names] <- use_contrast_position;
      }
   }
   if (debug) {
      jamba::printDebug("resolved contrast_position:");print(contrast_position);# debug
   }

   # optionally "flip" two-way contrasts
   is_twoway <- grepl("[(]", contrast_names)
   if (TRUE %in% flip_twoway && TRUE %in% is_twoway) {
      twoway_split <- strsplit(
         gsub("[()]", "", contrast_names[is_twoway]),
         "-")
      twoway_flipped <- jamba::cPaste(lapply(twoway_split, function(icon){
         paste0(c("(", ""), icon[c(1, 3, 2, 4)], c("", ")"))
      }), sep="-")
      names(twoway_flipped) <- contrast_names[is_twoway];
      contrast_names[is_twoway] <- twoway_flipped;
      # flip the contrast name format in names(contrast_labels)
      if (length(contrast_labels) > 0) {
         labeled_flip <- names(contrast_labels) %in% names(twoway_flipped)
         if (any(labeled_flip)) {
            flip_match <- match(names(contrast_labels)[labeled_flip],
               names(twoway_flipped))
            names(contrast_labels)[labeled_flip] <- twoway_flipped[flip_match];
         }
      }
      # flip the contrast name format in names(contrast_position)
      if (length(contrast_position) > 0) {
         labeled_flip <- names(contrast_position) %in% names(twoway_flipped)
         # jamba::printDebug("pre-flipped contrast_position:");print(contrast_position);# debug
         if (any(labeled_flip)) {
            flip_match <- match(names(contrast_position)[labeled_flip],
               names(twoway_flipped))
            names(contrast_position)[labeled_flip] <- twoway_flipped[flip_match];
            # jamba::printDebug("flipped contrast_position:");print(contrast_position);# debug
         }
      }
   }

   # expand two-way contrasts to one-way contrasts to ensure each are present
   is_twoway <- grepl("[(]", contrast_names)
   if (any(is_twoway)) {
      oneway_contrast_names <- unlist(
         strsplit(
            gsub("^[(]|[)]$", "", contrast_names[is_twoway]),
            "[)]-[(]"))
      # ensure one-way contrasts appear before two-way contrasts
      new_oneway_contrast_names <- setdiff(oneway_contrast_names,
         contrast_names[!is_twoway]);
      if (length(new_oneway_contrast_names) > 0) {
         contrast_names <- unique(c(contrast_names[!is_twoway],
            new_oneway_contrast_names,
            contrast_names[is_twoway]))
      }
      new_oneway_contrast_position <- setdiff(oneway_contrast_names,
         names(contrast_position));
      if (length(new_oneway_contrast_position) > 0) {
         contrast_position[new_oneway_contrast_position] <- jamba::nameVector(
            rep(oneway_position,
               length.out=length(new_oneway_contrast_position)),
            new_oneway_contrast_position)
      }
   }

   split_contrast_names <- function(contrast_names, factors_df) {
      is_twoway <- grepl("[(]", contrast_names);
      split_names <- as.list(rep(NA, length(contrast_names)))
      if (any(is_twoway)) {
         twoway_list <- strsplit(contrast_names[is_twoway], "[)]-[(]")
         twoway_split <- rep(seq_along(twoway_list), lengths(twoway_list));
         oneway_contrasts <- gsub("^[(]|[)]$", "", unlist(twoway_list))
         oneway_split_contrasts <- split_contrast_names(oneway_contrasts,
            factors_df)
         twoway_split_list <- split(oneway_split_contrasts, twoway_split)
         split_names[is_twoway] <- twoway_split_list;
      }
      split_contrasts <- strsplit(contrast_names[!is_twoway], "-")
      # apply labels
      split_contrasts <- lapply(split_contrasts, function(igroups){
         jamba::nameVector(
            factors_df[match(igroups, rownames(factors_df)), "label"],
            igroups);
      })
      split_names[!is_twoway] <- split_contrasts;
      names(split_names) <- contrast_names;
      return(split_names)
   }
   contrast_group_list <- split_contrast_names(contrast_names, factors_df)
   contrast_group_dfs <- lapply(names(contrast_group_list), function(iname){
      i <- contrast_group_list[[iname]];
      if (is.list(i)) {
         # two-way contrasts
         idf <- jamba::rbindList(lapply(i, function(j){
            jdf <- data.frame(check.names=FALSE,
               stringsAsFactors=FALSE,
               contrast=iname,
               axis_df[match(j, axis_df$label), , drop=FALSE])
            if (length(unique(jdf$x_coord)) > 1) {
               if (length(unique(jdf$y_coord)) > 1) {
                  jdf$angle <- jamba::rad2deg(atan2(y=diff(jdf$y_coord),
                     x=diff(jdf$x_coord)))
                  islope <- diff(jdf$y_coord) / diff(jdf$x_coord)
                  jdf$intercept <- (islope * (-jdf$x_coord[1])) + jdf$y_coord;
               } else {
                  if (jdf$x_coord[1] > jdf$x_coord[2]) {
                     jdf$angle <- 180
                  } else {
                     jdf$angle <- 0
                  }
                  jdf$intercept <- jdf$y_coord;
               }
            } else if (length(unique(jdf$y_coord)) > 1) {
               if (jdf$y_coord[1] > jdf$y_coord[2]) {
                  jdf$angle <- 270;
               } else {
                  jdf$angle <- 90;
               }
               jdf$intercept <- jdf$x_coord;
            } else {
               jdf$angle <- 0;
               jdf$intercept <- jdf$y_coord;
            }
            jdf$depth <- 2;
            jdf$oneway_contrast <- jamba::cPaste(j, sep=contrast_sep);
            # Manhattan distance
            jdf$distance <- abs(diff(range(jdf$x_coord))) +
               abs(diff(range(jdf$y_coord)));
            jdf
         }))
      } else {
         idf <- data.frame(check.names=FALSE,
            stringsAsFactors=FALSE,
            contrast=iname,
            axis_df[match(i, axis_df$label), , drop=FALSE])
         if (length(unique(idf$x_coord)) > 1) {
            if (length(unique(idf$y_coord)) > 1) {
               idf$angle <- jamba::rad2deg(atan2(y=diff(idf$y_coord),
                  x=diff(idf$x_coord)))
               islope <- diff(idf$y_coord) / diff(idf$x_coord)
               idf$intercept <- (islope * (-idf$x_coord[1])) + idf$y_coord;
               # idf$distance <- sqrt(diff(idf$x_coord)^2 + diff(idf$y_coord)^2)
            } else {
               if (idf$x_coord[1] > idf$x_coord[2]) {
                  idf$angle <- 180
               } else {
                  idf$angle <- 0
               }
               idf$intercept <- idf$y_coord;
               # idf$distance <- abs(diff(idf$x_coord))
            }
         } else if (length(unique(idf$y_coord)) > 1) {
            if (idf$y_coord[1] > idf$y_coord[2]) {
               idf$angle <- 270;
            } else {
               idf$angle <- 90;
            }
            idf$intercept <- idf$x_coord;
            # idf$distance <- abs(diff(idf$y_coord))
         } else {
            idf$angle <- 0;
            idf$intercept <- idf$y_coord;
         }
         idf$depth <- 1;
         idf$oneway_contrast <- jamba::cPaste(i, sep=contrast_sep);
         # Manhattan distance (x + y)
         idf$distance <- abs(diff(idf$x_coord)) + abs(diff(idf$y_coord))
         idf;
      }
   })
   contrast_group_df <- jamba::rbindList(contrast_group_dfs)

   # sort contrasts by increasing distance,
   # so shorter contrasts are not bumped as much
   contrast_group_df <- jamba::mixedSortDF(contrast_group_df,
      byCols="distance");


   # determine which contrast should be labeled,
   # assuming a contrast may be present multiple times, and
   # only needs a label the last time it is rendered.
   oneway_seq <- seq(from=1, to=nrow(contrast_group_df), by=2);
   oneway_render <- rep(rev(!duplicated(
      rev(contrast_group_df$oneway_contrast[oneway_seq]))), each=2);
   contrast_group_df$render_contrast <- oneway_render;
   # jamba::printDebug("contrast_group_df:");print(contrast_group_df);# debug
   # if (nrow(contrast_group_df) == 48) {contrast_group_df <- contrast_group_df[-17:-24,]}
   # return(contrast_group_df);

   # apply "bump"
   contrast_group_df$bump_set <- jamba::pasteByRow(
      contrast_group_df[,c("angle", "intercept"), drop=FALSE]);
   contrast_group_df$bump <- 0;
   for (ibump in unique(contrast_group_df$bump_set)) {
      ibump_rows <- which(contrast_group_df$bump_set %in% ibump &
            !grepl("[(]", contrast_group_df$contrast))
      if (length(ibump_rows) > 1) {
         ibump_values <- seq_len(length(ibump_rows)/2)
         ibump_values <- ibump_values - mean(ibump_values)
         contrast_group_df[ibump_rows, "bump"] <- rep(ibump_values, each=2);
      }
   }

   contrast_group_split <- split(contrast_group_df,
      factor(contrast_group_df$contrast,
         levels=unique(contrast_group_df$contrast)))
   xlim <- range(contrast_group_df$x_coord, na.rm=TRUE)
   ylim <- range(contrast_group_df$y_coord, na.rm=TRUE)
   xlim <- range(axis_df$x_coord, na.rm=TRUE)
   ylim <- range(axis_df$y_coord, na.rm=TRUE)

   # define buffer between groups
   group_buffer <- jamba::noiseFloor(group_buffer,
      minimum=0,
      ceiling=0.5);
   use_group_buffer <- 0.5 - group_buffer;

   # adjust bump_factor proportional to group_buffer
   if (length(bump_factor) == 0) {
      bump_factor <- 1;
   }
   bump_factor <- rep(bump_factor, length.out=2);
   use_bump_factor <- jamba::noiseFloor(bump_factor,
      minimum=0,
      ceiling=1.5) * (use_group_buffer * 2)
   max_x_bump <- (1.5 / use_bump_factor[2]) * (2 + max(subset(contrast_group_df,
      !angle %in% c(0, 180))$bump));
   max_y_bump <- (1.5 / use_bump_factor[1]) * (2 + max(subset(contrast_group_df,
      angle %in% c(0, 180))$bump));
   # define arrow_ex
   if (length(arrow_ex) == 0 || !is.numeric(arrow_ex)) {
      max_x_bumps <- jamba::noiseFloor(
         2 * (0.5 + max(subset(contrast_group_df, !angle %in% c(0, 180))$bump)),
         minimum=1,
         ceiling=20);
      max_y_bumps <- jamba::noiseFloor(
         2 * (0.5 + max(subset(contrast_group_df, angle %in% c(0, 180))$bump)),
         minimum=1,
         ceiling=20);
      use_arrow_ex <- 1 / sqrt(jamba::rmNA(naValue=1,
         max(c(max_x_bumps, max_y_bumps), na.rm=TRUE)))
   } else {
      use_arrow_ex <- arrow_ex;
   }
   # return(invisible(contrast_group_list));
   if (debug) {
      jamba::printDebug("contrast_group_df:");print(contrast_group_df);# debug
   }

   # custom function to handle contrast labeling
   #' @param contrast `character` contrast name
   #' @param contrast_labels `character` vector of labels named by contrast
   #' @param contrast_style `character` string deciding how to format
   #'    the contrast:
   #'    * `"comp"`: calls contrast2comp()
   #'    * `"contrast"`: uses the contrast as-is
   #'    * `"none"`: hides the contrast label, appending `contrast_labels`
   #'    when provided
   handle_contrast_label <- function
   (contrast,
    contrast_labels,
    contrast_style=c("comp", "contrast", "none"))
   {
      contrast_style <- match.arg(contrast_style);
      # initial label to use
      if ("comp" %in% contrast_style) {
         use_label1 <- contrast2comp(contrast)
      } else if ("contrast" %in% contrast_style) {
         use_label1 <- contrast
      } else {
         use_label1 <- NULL;
      }
      if (length(contrast_labels) > 0 &&
            contrast %in% names(contrast_labels)) {
         use_label1 <- jamba::cPaste(c(
            use_label1,
            contrast_labels[[contrast]]),
            sep="\n")
      }
      if (length(use_label1) == 0) {
         use_label1 <- ""
      }
      return(use_label1)
   }


   if ("base" %in% plot_type && TRUE %in% do_plot) {
      if (length(plot_margins) > 0) {
         plot_margins <- rep(plot_margins, length.out=4)
         opar <- par(mar=plot_margins)
         on.exit(par(opar), add=TRUE)
      }
      jamba::nullPlot(xlim=xlim + c(-1, 1),
         ylim=ylim + c(-1, 1),
         doBoxes=FALSE,
         asp=1)

      # jamba::printDebug("axis_df:");print(axis_df);# debug
      if (any(!is.na(axis_df$axis1))) {
         jamba::groupedAxis(
            side=1,
            las=1,
            pos=0.45,
            x=unique(axis_df[,c("axis1", "x_coord")])$axis1,
            nudge=use_group_buffer,
            group_style="grouped")
      }
      if (any(!is.na(axis_df$axis2))) {
         jamba::groupedAxis(
            side=2,
            pos=0.45,
            x=unique(axis_df[,c("axis2", "y_coord")])$axis2,
            nudge=use_group_buffer,
            group_style="grouped")
      }
      if (any(!is.na(axis_df$axis3))) {
         if (debug) {
            jamba::printDebug("axis3:");# debug
            print(unique(axis_df[,c("axis3", "x_coord")]));# debug
         }
         jamba::groupedAxis(
            side=3,
            las=1,
            pos=max(axis_df$y_coord, na.rm=TRUE) + 0.55,
            x=unique(axis_df[,c("axis3", "x_coord")])$axis3,
            nudge=use_group_buffer,
            group_style="grouped")
      }
      if (any(!is.na(axis_df$axis4))) {
         jamba::groupedAxis(
            side=4,
            pos=max(axis_df$x_coord, na.rm=TRUE) + 0.55,
            x=unique(axis_df[,c("axis4", "y_coord")])$axis4,
            nudge=use_group_buffer,
            group_style="grouped")
      }

      is_twoway <- grepl("[(]", contrast_names);
      if (length(colorset) == 0) {
         colorset <- colorjam::rainbowJam(
            n=sum(!is_twoway),
            Crange=c(60, 90),
            Lrange=c(44, 80))
         names(colorset) <- contrast_names[!is_twoway];
      } else {
         if (length(names(colorset)) > 0) {
            missing_names <- setdiff(contrast_names[!is_twoway],
               names(colorset))
            if (length(missing_names) > 0) {
               new_colorset <- colorjam::rainbowJam(
                  n=length(missing_names),
                  Crange=c(60, 90),
                  Lrange=c(44, 80))
               names(new_colorset) <- missing_names;
               colorset[names(missing_names)] <- new_colorset;
            }
         } else {
            colorset <- rep(colorset,
               length.out=sum(!is_twoway))
            names(colorset) <- contrast_names[!is_twoway];
         }
      }
      # jamba::printDebug(colorset, sep=",\n");# debug

      # iterate each group to draw a visible square around each group
      if (verbose > 1) {
         jamba::printDebug("axis_df:");print(axis_df);
         print(data.frame(xleft=axis_df$x_coord - use_group_buffer,
            xright=axis_df$x_coord + use_group_buffer,
            ybottom=axis_df$y_coord - use_group_buffer,
            ytop=axis_df$y_coord + use_group_buffer))
      }
      rect(xleft=axis_df$x_coord - use_group_buffer,
         xright=axis_df$x_coord + use_group_buffer,
         ybottom=axis_df$y_coord - use_group_buffer,
         ytop=axis_df$y_coord + use_group_buffer,
         col=ifelse(is.na(axis_df$n),
            "transparent",
            group_fill[axis_df$group]),
         border=ifelse(is.na(axis_df$n),
            "transparent",
            group_border[axis_df$group]))
      # optionally add number of replicates as labels
      text(
         x=axis_df$x_coord - use_group_buffer + 0.04,
         y=axis_df$y_coord + use_group_buffer - 0.04,
         adj=c(0, 1),
         labels=ifelse(is.na(axis_df$n), "",
            paste0("n=", axis_df$n)),
         col=replicate_color,
         cex=replicate_cex)

      # generate summary of pairwise contrasts
      contrast_summary_df <- jamba::rbindList(
         lapply(seq_along(contrast_group_split), function(i){
            idf <- contrast_group_split[[i]];
            idf$color <- colorset[as.character(idf$oneway_contrast)];
            if (nrow(idf) == 2) {
               data.frame(from=idf$label[1],
                  to=idf$label[2],
                  contrast=idf$contrast[1],
                  full_contrast=idf$contrast[1],
                  color=colorset[idf$contrast[1]],
                  depth="one-way",
                  contrast_depth=1)
            } else if (nrow(idf) == 4) {
               idf$group <- strsplit(gsub("[()]", "", idf$contrast[1]),
                  contrast_sep)[[1]];
               colorset_key <- jamba::cPaste(
                  list(idf$group[c(1, 3)],
                     idf$group[c(2, 4)]),
                  sep="-");
               use_colors <- jamba::rmNA(colorset[colorset_key],
                  naValue="grey");
               names(use_colors) <- colorset_key;
               data.frame(from=idf$label[c(1, 3)],
                  to=idf$label[c(2, 4)],
                  contrast=paste0(idf$group[c(1, 3)], "-",
                     idf$group[c(2, 4)]),
                  full_contrast=idf$contrast[1],
                  color=use_colors,
                  depth="two-way",
                  contrast_depth=2)
            } else {
               data.frame(from="a",
                  to="a",
                  contrast="a",
                  full_contrast="a",
                  color="a",
                  depth="a",
                  contrast_depth=0)[0,]
            }
         }));

      # iterate each contrast
      use_contrasts <- seq_along(contrast_group_split)
      # jamba::printDebug("contrast_summary_df:");print(contrast_summary_df);# debug
      # jamba::printDebug("contrast_group_split:");print(contrast_group_split);# debug
      if (debug) {
         jamba::printDebug("contrast_group_split cgs_df:");
         # print(contrast_group_split);# debug
         cgs_df <- jamba::rbindList(contrast_group_split);
         rownames(cgs_df) <- NULL;
         cgs_df$num <- as.numeric(factor(cgs_df$contrast,
            levels=unique(cgs_df$contrast)))
         print(cgs_df);# debug
      }
      for (i in seq_along(contrast_group_split)[use_contrasts]) {
         idf <- contrast_group_split[[i]];
         if (length(contrast_depths) > 0 && is.numeric(contrast_depths)) {
            if (!contrast_depths %in% idf$depth) {
               if (verbose) {
                  jamba::printDebug("plot_sedesign(): ",
                     "Skipping contrast with depth=", head(idf$depth, 1));
               }
               next;
            }
         }
         if (idf$angle[1] %in% c(0, 180)) {
            idf$y_coord <- idf$y_coord + idf$bump / max_y_bump;
         } else if (idf$angle[1] %in% c(90, 270)) {
            idf$x_coord <- idf$x_coord + idf$bump / max_x_bump;
         }
         if (nrow(idf) == 2) {
            if (TRUE %in% idf$render_contrast) {
               # assemble custom contrast label as relevant
               # jamba::printDebug("idf oneway contrast:");print(idf);# debug
               use_label <- handle_contrast_label(
                  contrast=idf$contrast[1],
                  contrast_labels=contrast_labels,
                  contrast_style=contrast_style)
               # draw the contrast
               draw_oneway_contrast(
                  x=idf$x_coord[1],
                  x1=idf$x_coord[2],
                  y=idf$y_coord[1],
                  y1=idf$y_coord[2],
                  plot_type="base",
                  color=colorset[idf$contrast[1]],
                  label=use_label,
                  label_cex=label_cex,
                  oneway_position=contrast_position[idf$contrast[1]],
                  arrow_ex=use_arrow_ex,
                  verbose=verbose,
                  ...);
            }
         } else if (nrow(idf) == 4) {
            # match with each pairwise contrast
            # jamba::printDebug("contrast_summary_df:");print(contrast_summary_df);# debug
            summary_df1 <- subset(contrast_summary_df,
               full_contrast %in% idf$contrast[1])
            # jamba::printDebug("summary_df1:");print(summary_df1);# debug
            summary_df <- subset(contrast_summary_df,
               contrast %in% summary_df1$contrast &
               depth %in% "one-way")
            if (nrow(summary_df) == 0) {
               summary_df <- subset(contrast_summary_df,
                  contrast %in% summary_df1$contrast &
                     depth %in% "two-way")
            }
            # jamba::printDebug("summary_df:");print(summary_df);# debug
            summary_df <- summary_df[match(summary_df1$contrast,
               summary_df$contrast), , drop=FALSE];
            colorset_twoway <- summary_df$color;
            if (any("grey" %in% colorset_twoway) && TRUE %in% debug) {
               jamba::printDebug("contrast_summary_df:");print(contrast_summary_df);# debug
               jamba::printDebug("idf:");print(idf);# debug
               jamba::printDebug("summary_df1:");print(summary_df1);# debug
               jamba::printDebug("summary_df:");print(summary_df);# debug
            }
            contrast_dfs <- lapply(jamba::nameVector(as.character(summary_df1$contrast)), function(jname){
               jdf <- contrast_group_split[[jname]];
               if (debug && length(jdf) == 0) {
                  jamba::printDebug("contrast_group_split:");print(contrast_group_split);# debug
                  jamba::printDebug("summary_df1$contrast:");print(summary_df1$contrast);# debug
                  jamba::printDebug("contrast_group_split (jdf):");print(jdf);# debug
               }
               jmatch <- match(jdf$oneway_contrast, idf$oneway_contrast);
               jdf[, "render_contrast"] <- idf[jmatch, "render_contrast"];
               if (jdf$angle[1] %in% c(0, 180)) {
                  jdf$y_coord <- jdf$y_coord + jdf$bump / max_y_bump;
               } else if (idf$angle[1] %in% c(90, 270)) {
                  jdf$x_coord <- jdf$x_coord + jdf$bump / max_x_bump;
               }
               jdf;
            })
            # call draw_twoway_contrast()

            # assemble custom contrast label as relevant
            icontrast1 <- head(contrast_dfs[[1]]$contrast, 1);
            use_label1 <- handle_contrast_label(
               contrast=icontrast1,
               contrast_labels=contrast_labels,
               contrast_style=contrast_style)
            icontrast2 <- head(contrast_dfs[[2]]$contrast, 1);
            use_label2 <- handle_contrast_label(
               contrast=icontrast2,
               contrast_labels=contrast_labels,
               contrast_style=contrast_style)
            if (debug) {
               jamba::printDebug("icontrast1:", icontrast1,
                  ", use_label1:", use_label1,
                  ", icontrast2:", icontrast2,
                  ", use_label2:", use_label2)# debug
            }

            # draw the contrast
            if (verbose) {
               jamba::printDebug("plot_sedesign(): ",
                  "draw_twoway_contrast()");
            }
            if (debug) {
               jamba::printDebug("contrast_dfs:");print(contrast_dfs);# debug
            }
            use_twoway_label <- NULL;
            if (length(contrast_labels) > 0 &&
                  idf$contrast[1] %in% names(contrast_labels)) {
               use_twoway_label <- contrast_labels[idf$contrast[1]];
            }
            # jamba::printDebug("idf$contrast[1]:");print(idf$contrast[1]);# debug
            # jamba::printDebug("use_twoway_label:", use_twoway_label);# debug
            # jamba::printDebug("contrast_dfs[[1]]$contrast[1]:");print(contrast_dfs[[1]]$contrast[1]);# debug
            # jamba::printDebug("contrast_position:");print(contrast_position);# debug
            # jamba::printDebug("contrast_dfs[[1]]:");print(contrast_dfs[[1]]);# debug
            # jamba::printDebug("contrast_dfs[[2]]:");print(contrast_dfs[[2]]);# debug
            # print(list(
            #    x0=c(contrast_dfs[[1]]$x_coord[1],
            #       contrast_dfs[[2]]$x_coord[1]),
            #    x1=c(contrast_dfs[[1]]$x_coord[2],
            #       contrast_dfs[[2]]$x_coord[2]),
            #    y0=c(contrast_dfs[[1]]$y_coord[1],
            #       contrast_dfs[[2]]$y_coord[1]),
            #    y1=c(contrast_dfs[[1]]$y_coord[2],
            #       contrast_dfs[[2]]$y_coord[2]),
            #    color=colorset_twoway))
            use_position <- contrast_position[c(
               contrast_dfs[[1]]$contrast[1],
               contrast_dfs[[2]]$contrast[1])];
            use_twoway_position <- contrast_position[idf$contrast[1]];
            draw_twoway_contrast(
               x0=c(contrast_dfs[[1]]$x_coord[1],
                  contrast_dfs[[2]]$x_coord[1]),
               x1=c(contrast_dfs[[1]]$x_coord[2],
                  contrast_dfs[[2]]$x_coord[2]),
               y0=c(contrast_dfs[[1]]$y_coord[1],
                  contrast_dfs[[2]]$y_coord[1]),
               y1=c(contrast_dfs[[1]]$y_coord[2],
                  contrast_dfs[[2]]$y_coord[2]),
               color=colorset_twoway,
               # border=border,
               draw_oneway=TRUE,
               twoway_label=use_twoway_label,
               extend_ex=extend_ex,
               arrow_ex=use_arrow_ex,
               extend_angle=extend_angle,
               label=c(use_label1,
                  use_label2),
               label_cex=label_cex,
               oneway_position=use_position,
               twoway_position=use_twoway_position,
               twoway_lwd=twoway_lwd,
               ...)
            # jamba::printDebug("contrast_position[idf$contrast[1]]:", contrast_position[idf$contrast[1]]);# debug
            # jamba::printDebug("contrast_position:");print(contrast_position);# debug
            # jamba::printDebug("idf$contrast[1]:");print(idf$contrast[1]);
            # jamba::printDebug("twoway_position:");print(contrast_position[idf$contrast[1]]);# debug
         }
      }
      attr(contrast_group_split, "max_x_bump") <- max_x_bump;
      attr(contrast_group_split, "max_y_bump") <- max_y_bump;
      return(invisible(contrast_group_split))
   } else if ("grid" %in% plot_type) {
   }

}

#' Make block arrow polygon coordinates for line segments
#'
#' Make block arrow polygon coordinates for line segments
#'
#' This function defines a block arrow defined by line segments.
#'
#' The block arrow is defined with a fixed arrow head size,
#' in order to preserve the aspect ratio for the arrow head.
#' The arrow stem is extended to the length of the line segment,
#' unless the arrow head itself is larger than the line segment
#' in which case only the arrow head is shown.
#'
#' The overall arrow size can be adjusted with `arrow_ex`,
#' which adjusts the arrow stem width and the arrow head
#' proportionally.
#' The size of the arrow head relative to the arrow stem width
#' can be adjusted with `head_ex`.
#'
#' When any line segments have zero width, the final point in the line
#' segment is shifted `0.1` to create a small horizontal
#' line segment with length `0.1`.
#'
#' TODO: Prepare `list` of polygons for gradient color fill.
#'
#' * Optionally return a `list` of polygons usable to render gradient
#' color fill along the arrow stem, with final color on the arrow head.
#' The list should include color and border for each polygon, so the border
#' is used to fill the tiny "gap" between adjacent polygons, to prevent
#' visual artifacts.
#' The final element in the list should be the full arrow, with `NA` color,
#' and proper border to be drawn atop the series of gradient colors.
#'
#' Note: Unfortunately even `grid::grid.rect()` used with
#' `grid::linearGradient()` is unable to fill this purpose,
#' since it is not available for `quartz()` and `windows()` devices,
#' and only works properly with a subset of graphics devices such as
#' PDF, SVG, PNG, etc. That said, it would otherwise be the preferred
#' approach, along with the `vwline` variable grid line package.
#'
#' @param x `numeric` vector with the start value for each line segment,
#'    or when `x1=NULL` then `x` must contain two values per line segment.
#' @param y `numeric` vector with the start value for each line segment,
#'    or when `y1=NULL` then `y` must contain two values per line segment.
#' @param x1 `numeric` with the end value for each line segment.
#' @param y1 `numeric` with the end value for each line segment.
#' @param reference `character` string indicating which point is the
#'    reference when the arrow head is longer than the line segment.
#'    * `"last"` (default) fixes the last point in the line to the
#'    point of the block arrow head, so the arrow base (stem) will
#'    extend past the first point in the line.
#'    * `"first"` fixes the first point in the line to the end
#'    of the arrow stem, so the arrow head will extend beyond the
#'    final point in the line.
#' @param data_format `character` string indicating the data format to
#'    return, based upon the type of plot used in subsequent steps:
#'    * `"grid"` (default) returns a `list` with `x`, `y`, and `id`,
#'    usable with `grid::grid.polygon(x, y, id)` to plot multiple
#'    polygons in vectorized form.
#'    * `"base"` returns a `list` with `x` and `y` values, where each
#'    vector contains `NA` values to separate each polygon, usable
#'    with `graphics::polygon(x, y)` for vectorized plotting.
#'    * `"list"` returns a `list` where `x` is a `list` separated by
#'    each polygon, and `y` is a `list` separated by each polygon.
#'    One would iterate the coordinate lists to plot the polygons.
#' @param arrow_ex `numeric` expansion factor applied to the arrow width,
#'    and by default to the arrow head, relative to the arrow width.
#' @param head_ex `numeric` expansion factor applied to the arrow head,
#'    by default applied relative to the arrow width.
#' @param arrow_w `numeric` width for the arrow stem, actually the half-width
#'    since the width is applied to each side perpendicular to the line.
#' @param head_w `numeric` head width, added to `arrow_w` for each side
#'    of the arrow.
#' @param head_l `numeric` head length, a fixed distance from the end of
#'    the line described between x,y points.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' arrow_coords <- matrix(ncol=4, byrow=TRUE, c(
#'    0, 0, 1, 0,
#'    0, 1, 2, 1,
#'    0, 2, 1, 2,
#'    0, 3, 2, 3,
#'    0, 0, -1, 0,
#'    0, 0, 0, 1,
#'    0, 0, 0, -1,
#'    0, 2, 1, 1,
#'    0, -2, 1, -2,
#'    0, -3, 0, -3,
#'    0, 0, -2, -2,
#'    1, -1, 1.1, -1.1))
#' x <- arrow_coords[,1];
#' y <- arrow_coords[,2];
#' x1 <- arrow_coords[,3];
#' y1 <- arrow_coords[,4];
#'
#' colorset <- colorjam::rainbowJam(nrow(arrow_coords), Crange=c(60, 90), Lrange=c(46, 80));
#'
#' jamba::nullPlot(xlim=c(-3, 3), ylim=c(-3, 3), asp=1, doBoxes=FALSE);
#' axis(1, las=2); axis(2);
#' box();
#' abline(h=-3:3, v=-3:3, col="grey90", lty=2)
#'
#' k <- seq_len(nrow(arrow_coords));
#' arrowxy <- make_block_arrow_polygon(
#'    x=x[k], y=y[k],
#'    x1=x1[k], y1=y1[k],
#'    head_ex=rep(c(1, 3, 1), c(1, 1, 9)),
#'    arrow_ex=rep(c(1, 2, 1), c(2, 1, 8)),
#'    data_format="base")
#' polygon(x=arrowxy$x, y=arrowxy$y,
#'    col=colorset[k])
#'
#' for (i in k) {
#'    irad <- atan2(y=(y1 - y)[i], x=(x1 - x)[i]);
#'    iangle <- ((jamba::rad2deg(irad) + 90) %% 180 - 90) %% 360;
#'    jamba::shadowText(x=(x + x1)[i] / 2,
#'       y=(y + y1)[i] / 2,
#'       col="white",
#'       srt=attr(arrowxy, "text_angle")[i],
#'       label=paste0("label ", i))
#' }
#'
#' @export
make_block_arrow_polygon <- function
(x,
 y,
 x1=NULL,
 y1=NULL,
 reference="last",
 data_format=c("grid",
    "base",
    "list"),
 arrow_ex=1,
 head_ex=1,
 arrow_w=0.1 * arrow_ex,
 head_w=arrow_w * 0.75 * head_ex,
 head_l=(head_w + arrow_w) * 1.5,
 color=NA,
 border=NA,
 gradient_n=15,
 verbose=FALSE,
 ...)
{
   #
   data_format <- match.arg(data_format);

   if (!all(reference %in% c("first", "last"))) {
      stop("reference must only contain values 'first' and 'last'.")
   }

   if (length(x1) == 0) {
      is_odd <- (seq_along(x) %% 2 == 1)
      x1 <- x[!is_odd];
      x <- x[is_odd];
   }
   if (length(y1) == 0) {
      is_odd <- (seq_along(y) %% 2 == 1)
      y1 <- y[!is_odd];
      y <- y[is_odd];
   }

   # check for duplicate start/end points
   dupe_points <- (x == x1 & y == y1)
   if (any(dupe_points)) {
      x1[dupe_points] <- x1[dupe_points] + 0.1
   }

   # line length
   line_len <- sqrt((x - x1)^2 + (y - y1)^2);
   # angle
   xy_angle <- atan2(y=(y1 - y), x=(x1 - x))
   text_angle <- ((jamba::rad2deg(xy_angle) + 90) %% 180 - 90) %% 360

   arrow_w <- rep(arrow_w, length.out=length(x));
   head_w <- rep(head_w, length.out=length(x));
   head_l <- rep(head_l, length.out=length(x));
   reference <- rep(reference, length.out=length(x));

   if (verbose) {
      jamba::printDebug("make_block_arrow_polygon(): ",
         "input data:");
      print(data.frame(x, x1, y, y1, line_len, xy_angle,
         arrow_w, head_w, head_l));
   }

   # vectorized format
   arrow_w_use <- rep(arrow_w, each=7)
   head_w_use <- rep(head_w, each=7)
   head_l_use <- rep(head_l, each=7)
   line_len_use <- rep(line_len, each=7)

   # standard block arrow with sufficient line length
   arrow_x_list <- lapply(seq_along(x), function(i){
      if (line_len[i] < head_l[i]) {
         iarrow_x <- c(-head_l[i],
            -head_l[i],
            -head_l[i],
            0,
            -head_l[i],
            -head_l[i],
            -head_l[i])
      } else {
         iarrow_x <- c(0,
            line_len[i] - head_l[i],
            line_len[i] - head_l[i],
            line_len[i],
            line_len[i] - head_l[i],
            line_len[i] - head_l[i],
            0)
      }
   })
   arrow_y_list <- lapply(seq_along(x), function(i){
      if (line_len[i] < head_l[i]) {
         iarrow_y <- c(0,
            arrow_w[i],
            arrow_w[i] + head_w[i],
            0,
            -arrow_w[i] - head_w[i],
            -arrow_w[i],
            0)
      } else {
         iarrow_y <- c(arrow_w[i],
            arrow_w[i],
            arrow_w[i] + head_w[i],
            0,
            -arrow_w[i] - head_w[i],
            -arrow_w[i],
            -arrow_w[i])
      }
   })
   if (verbose > 1) {
      jamba::printDebug("make_block_arrow_polygon(): ",
         "arrow_x_list:");
      print(arrow_x_list);
      jamba::printDebug("make_block_arrow_polygon(): ",
         "arrow_y_list:");
      print(arrow_y_list);
   }
   arrow_x <- unlist(arrow_x_list);
   arrow_y <- unlist(arrow_y_list);

   # Previous single-entry block arrow logic
   if (FALSE) {
      if (line_len < head_l) {
         arrow_x <- c(-head_l, -head_l, -head_l,
            0,
            -head_l, -head_l, -head_l)
         arrow_y <- c(0, arrow_w,
            arrow_w + head_w, 0, -arrow_w - head_w,
            -arrow_w, 0)
      } else if (FALSE) {
         arrow_x <- c(0,
            line_len - head_l,
            line_len - head_l,
            line_len,
            line_len - head_l,
            line_len - head_l,
            0)
         arrow_y <- c(arrow_w,
            arrow_w,
            arrow_w + head_w,
            0,
            -arrow_w - head_w,
            -arrow_w, -arrow_w)
      }
   }

   # rotate matrix coordinates
   co <- cos(xy_angle)
   si <- sin(xy_angle)
   use_x <- (rep(co, each=7) * arrow_x -
         rep(si, each=7) * arrow_y);
   use_y <- (rep(si, each=7) * arrow_x +
         rep(co, each=7) * arrow_y);

   index7 <- seq_along(x) * 7;
   index1 <- index7 - 6;
   index4 <- index7 - 3;
   start_x <- (use_x[index1] + use_x[index7]) / 2
   start_y <- (use_y[index1] + use_y[index7]) / 2
   end_x <- use_x[index4]
   end_y <- use_y[index4]
   ## previous single-arrow logic
   # start_x <- mean(use_x[c(1, 7)])
   # start_y <- mean(use_y[c(1, 7)])
   # end_x <- use_x[4]
   # end_y <- use_y[4]

   # adjust each polygon so the final point is the reference
   # new_x <- (use_x + rep(x1 - end_x, each=7))
   # new_y <- use_y + (rep(y1 - end_y, each=7))
   # adjust points relative to the reference
   new_x <- ifelse(rep(reference, each=7) %in% "last",
      (use_x + rep(x1 - end_x, each=7)),
      (use_x + rep(x - start_x, each=7)))
   new_y <- ifelse(rep(reference, each=7) %in% "last",
      (use_y + rep(y1 - end_y, each=7)),
      (use_y + rep(y - start_y, each=7)))
   if (verbose > 1) {
      print(data.frame(start_x, start_y,
         end_x,
         end_y=round(digits=3, end_y),
         use_x, use_y, new_x, new_y))
   }

   ## optionally convert to base R polygon() format
   if ("grid" %in% data_format) {
      id <- rep(seq_along(x), each=7)
      retlist <- list(
         x=new_x,
         y=new_y,
         id=id);
      attr(retlist, "text_angle") <- text_angle;
      return(retlist);
   }

   new_x_list <- split(new_x,
      rep(seq_along(x), each=7))
   new_y_list <- split(new_y,
      rep(seq_along(x), each=7))
   if ("list" %in% data_format) {
      retlist <- list(
         x=new_x_list,
         y=new_y_list);
      attr(retlist, "text_angle") <- text_angle;
      return(retlist);
   }

   if ("base" %in% data_format) {
      new_x <- head(unlist(lapply(new_x_list, function(ix){
         c(ix, NA)
      })), -1)
      new_y <- head(unlist(lapply(new_y_list, function(iy){
         c(iy, NA)
      })), -1)
   }

   retlist <- list(x=new_x,
      y=new_y);
   attr(retlist, "text_angle") <- text_angle;
   return(retlist)
}

#' convert point-slope to axis intercept
#'
#' @return `numeric` intercept on the y-axis, except with slope is
#'    infinite in which case the value returned is the x-axis intercept.
#'
#' @param pt `numeric` matrix with two columns, or coerced to matrix
#'    with two columns using `byrow=TRUE`.
#' @param slope `numeric` slope.
#' @param angle `numeric` value, currently ignored.
#' @param ... additional arguments are ignored.
#'
#' @export
point_slope_intercept <- function
(pt,
 slope,
 angle,
 ...)
{
   #
   if (length(dim(pt)) == 0) {
      pt <- matrix(ncol=2,
         byrow=TRUE,
         pt)
   }
   slope <- rep(slope, length.out=nrow(pt));
   intercept <- ifelse(is.infinite(slope),
      pt[,1],
      (slope * (-pt[,1])) + pt[,2]);
   return(intercept)
}

#' Determine which side one point is to another, given a slope or angle
#'
#' Determine which side one point is to another, given a slope or angle
#'
#' The result describes the position of the first line relative
#' to the second line, assuming both lines are parallel with identical slope.
#' For example "left" indicates that line 1 is on the left side of line 2.
#' Note that the `angle` is a more accurate measure of directionality,
#' otherwise `slope` is always assumed to desribe an angle moving to the right.
#'
#' When both lines are exactly overlapping, the result may be unstable,
#' however the result tends to favor "right" by default.
#'
#' @param pt1 `numeric` matrix of 2 columns, with x and y coordinates.
#' @param pt2 `numeric` matrix of 2 columns, with x and y coordinates.
#' @param slope `numeric` slope for each point in pt1 and pt2.
#' @param do_plot `logical` indicating whether to plot the result.
#' @param ... additional arguments are ignored.
#'
#' @return `character` vector equal to the number of points, `nrow(pt1)`:
#'    * `"right"` indicates `pt1` is on the right side of `pt2`
#'    * `"left"` indicates `pt1` is on the left side of `pt2`
#'
#' @examples
#' pt1 <- matrix(ncol=2, c(1, 1))
#' pt2 <- matrix(ncol=2, c(2, 2))
#' point_handedness(pt1, pt2, angle=0, do_plot=TRUE)
#' point_handedness(pt1=c(1, 1), pt2=c(2, 2), angle=0, do_plot=TRUE)
#'
#' pt1 <- matrix(ncol=2, c(1, 1))
#' pt2 <- matrix(ncol=2, c(0, 1))
#'
#' point_handedness(pt1, pt2, angle=45, do_plot=TRUE)
#' point_handedness(pt2, pt1, angle=45, do_plot=TRUE)
#'
#' point_handedness(rbind(pt2, pt2-1), rbind(pt1, pt1-1),
#'    angle=c(45, 45+180), do_plot=TRUE)
#'
#' point_handedness(pt1, pt2, angle=45, do_plot=TRUE)
#' title(main="angle = 45,\n(slope = 1)")
#' point_handedness(pt1, pt2, angle=45 + 180, do_plot=TRUE)
#' title(main="angle = 225,\n(slope = 1)")
#'
#' point_handedness(pt1, pt2, slope=Inf, do_plot=TRUE)
#' title(main="slope = Inf")
#' point_handedness(pt1, pt2, angle=90, do_plot=TRUE)
#' title(main="angle = 90")
#' point_handedness(pt2, pt1, slope=-Inf, do_plot=TRUE)
#' title(main="slope = -Inf")
#' point_handedness(pt2, pt1, angle = 270, do_plot=TRUE)
#' title(main="angle = 270")
#'
#' point_handedness(pt1, pt2, slope=0, do_plot=TRUE)
#' point_handedness(pt2, pt1, slope=0, do_plot=TRUE)
#'
#' pt1 <- matrix(ncol=2, c(-2, 5))
#' pt2 <- matrix(ncol=2, c(2, 3))
#' point_handedness(pt1, pt2, slope=1, do_plot=TRUE)
#' point_handedness(pt1, pt2, slope=-1, do_plot=TRUE)
#' point_handedness(pt1, pt2, slope=-1/3, do_plot=TRUE)
#' point_handedness(pt1, pt2, slope=-1/3, do_plot=TRUE)
#'
#' point_handedness(pt1, pt2, slope=Inf, do_plot=TRUE)
#' point_handedness(pt1, pt2, slope=-Inf, do_plot=TRUE)
#'
#' @export
point_handedness <- function
(pt1,
 pt2,
 slope=NULL,
 angle=NULL,
 do_plot=FALSE,
 verbose=FALSE,
 ...)
{
   # coerce pt1 and pt2 to two-column matrix form with byrow=TRUE
   if (length(dim(pt1)) == 0) {
      pt1 <- matrix(pt1, ncol=2, byrow=TRUE)
   }
   if (length(dim(pt2)) == 0) {
      pt2 <- matrix(pt2, ncol=2, byrow=TRUE)
   }
   if (length(slope) == 0) {
      if (length(angle) == 0) {
         stop("Either slope or angle must be supplied.")
      }
      angle <- rep(angle, length.out=nrow(pt1));
      angle <- angle %% 360;
      slope <- sin(jamba::deg2rad(angle)) / cos(jamba::deg2rad(angle))
      slope <- ifelse(slope > 1e8, Inf * sign(slope), slope)
   } else if (length(angle) == 0) {
      slope <- rep(slope, length.out=nrow(pt1));
      angle <- jamba::rad2deg(atan2(y=slope, x=1))
   }
   angle <- angle %% 360;

   pt1_intercept <- point_slope_intercept(pt1, slope=slope)
   pt2_intercept <- point_slope_intercept(pt2, slope=slope)
   handedness <- ifelse(is.infinite(slope),
      ifelse(angle <= 90,
         ifelse(pt1_intercept >= pt2_intercept, "right", "left"),
         ifelse(pt1_intercept < pt2_intercept, "right", "left")),
      ifelse(slope == 0,
         ifelse(pt2[,2] >= pt1[,2], "right", "left"),
         ifelse((angle > 0 & angle < 90) | (angle > 270),
            ifelse(pt1_intercept <= pt2_intercept, "right", "left"),
            ifelse(pt1_intercept >= pt2_intercept, "right", "left"))
      )
   )
   if (verbose) {
      jamba::printDebug("point_handedness(): ",
         "angle: ", angle);
      jamba::printDebug("point_handedness(): ",
         "pt1_intercept: ", pt1_intercept);
      jamba::printDebug("point_handedness(): ",
         "pt2_intercept: ", pt2_intercept);
      jamba::printDebug("point_handedness(): ",
         "handedness: ", handedness);
   }

   # optionally make a visual plot
   if (TRUE %in% do_plot) {
      xrange <- range(c(pt1[,1], pt2[,1]), na.rm=TRUE)
      yrange <- range(c(pt1[,2], pt2[,2]), na.rm=TRUE)
      xydiff <- max(abs(c(diff(xrange), diff(yrange))))
      xlim <- xrange + c(-1, 1) * xydiff * 4;
      ylim <- yrange + c(-1, 1) * xydiff * 4;
      jamba::nullPlot(doBoxes=FALSE,
         asp=1,
         xlim=xlim,
         ylim=ylim)
      box(); axis(1, las=2); axis(2);
      is_inf <- is.infinite(slope);
      if (any(is_inf)) {
         abline(v=pt1_intercept[is_inf], col="grey")
         abline(v=pt2_intercept[is_inf], col="grey")
      }
      if (any(!is_inf)) {
         abline(a=pt1_intercept[!is_inf], b=slope[!is_inf], col="grey")
         abline(a=pt2_intercept[!is_inf], b=slope[!is_inf], col="grey")
      }
      # draw arrow atop the sloped lines indicating direction
      arrows(
         x0=pt1[,1],
         y0=pt1[,2],
         x1=pt1[,1] + cos(jamba::deg2rad(angle)) * xydiff * 3,
         y1=pt1[,2] + sin(jamba::deg2rad(angle)) * xydiff * 3,
         col="grey20")
      arrows(
         x0=pt2[,1],
         y0=pt2[,2],
         x1=pt2[,1] + cos(jamba::deg2rad(angle)) * xydiff * 3,
         y1=pt2[,2] + sin(jamba::deg2rad(angle)) * xydiff * 3,
         col="grey20")
      print(pt1)
      if (nrow(pt1) == 1) {
         prefix <- "";
      } else {
         prefix <- paste0(
            jamba::colNum2excelName(seq_len(nrow(pt1))), ":")
      }
      print(prefix);
      points(pt1, pch=21, col="grey", bg="white", cex=3);
      text(pt1, paste0(prefix, "1"))
      pt1a <- pt2 - (pt2 - pt1) / 7;
      pt2a <- (pt2 - pt1) / 7 + pt1;
      arrows(x0=pt1a[,1], x1=pt2a[,1],
         y0=pt1a[,2], y1=pt2a[,2],
         col="grey20")
      points(pt2, pch=21, col="grey", bg="white", cex=3);
      text(pt2, paste0(prefix, "2"))
      jamba::drawLabels(preset="bottomright",
         drawBox=FALSE,
         labelCex=2,
         txt=jamba::cPaste(
            paste0(prefix, '"', handedness, '"'), "\n"))
   }
   return(handedness)
}
