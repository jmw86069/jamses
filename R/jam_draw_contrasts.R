
#' Draw one-way contrast using block arrows
#'
#' Draw one-way contrast using block arrows
#'
#' @param x0,x1,y0,y1 `numeric` values with the start and end coordinates,
#'    for the x and y axes, respectively.
#' @param color,border `character` R colors used to define color fill
#'    and border, respectively, for each block arrow,
#' @param plot_type `character` string indicating the type of plot output:
#'    * `"base"`: base R graphics
#'    * `"grid"`: grid graphics (not yet implemented)
#' @param label `character` vector or `list` with optional label to
#'    display atop each block arrow.
#'    For base R graphics, the label is drawn using
#'    `jamba::shadowText()` to render an outline around the text.
#'    * When `label` is a `character` vector, it is converted to a `list`
#'    in two ways depending upon the number of block arrows (`length(x0)`):
#'       * `length == 1`:  `label` is converted to `list` with length == 1.
#'       * `length > 1`: `label` is converted to `list` using `as.list`,
#'       then expanded to `length(x0)`.
#'
#'    * When `label` is passed as a `list`, or after `label` is converted
#'    to a `list`:
#'
#'       * Each block arrow label uses one concatenated string after
#'       calling `jamba::cPaste(..., sep=label_sep)` which joins values
#'       by default using newline `"\n"` between each value.
#' @param label_sep `character` string used as separator, passed to
#'    `jamba::cPaste()`, so that each block arrow may contain a vector
#'    which is concatenated using `label_sep` between each value.
#'    By default `label_sep="\n"` which prints each value on a new line.
#' @param na.rm `logical` passed to `jamba::cPaste()` to define how to
#'    display NA labels:
#'    * `na.rm=FALSE`: `"NA"`
#'    * `na.rm=TRUE`: `""`.
#' @param label_color `character` color used for the `label`.
#' @param label_cex `numeric` label font expansion factor, used to adjust
#'    the font size of the text label.
#' @param label_font `numeric` indicating the font face, defined as:
#'    * 1 = normal font
#'    * 2 = bold font
#'    * 3 = italic font
#'    * 4 = bold, italic font
#' @param do_plot `logical` indicating whether to draw the block arrow.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param debug `logical` indicating whether to print additional debug info.
#' @param ... additional arguments are passed to `make_block_arrow_polygon()`,
#'    including `arrow_ex` the arrow size expansion factor, and
#'    `head_ex` the arrow head expansion factor, which is adjusted
#'    relative to the arrow stem width.
#'
#' @examples
#' plot(NULL, xlim=c(0, 5), ylim=c(0, 4), asp=1, xlab="", ylab="")
#' draw_oneway_contrast(1, 4, 1, 1, label="contrast label")
#' draw_oneway_contrast(1, 4, 2, 2, head_ex=2, label_cex=1, label="contrast label")
#' draw_oneway_contrast(1, 4, 3, 3, arrow_ex=2, label_cex=2, label="contrast label")
#' draw_oneway_contrast(3, 0, 4, 1, arrow_ex=2, label_cex=2, label="contrast label")
#'
#' @export
draw_oneway_contrast <- function
(x0,
 x1,
 y0,
 y1,
 color="peachpuff",
 border="black",
 plot_type=c("base",
    "grid"),
 label=NULL,
 label_sep="\n",
 na.rm=TRUE,
 label_color="white",
 label_cex=1,
 label_font=2,
 oneway_position=0.5,
 do_plot=TRUE,
 verbose=FALSE,
 debug=FALSE,
 ...)
{
   # validate arguments
   plot_type <- match.arg(plot_type);

   # ensure vectorized values are proper lengths
   n <- max(c(length(x0), length(x1), length(y0), length(y1)));
   if (n == 0) {
      # consider adding a warning
      return(NULL)
   }
   x0 <- rep(x0, length.out=n)
   x1 <- rep(x1, length.out=n)
   y0 <- rep(y0, length.out=n)
   y1 <- rep(y1, length.out=n)

   # handle label argument
   if (length(label) == 0) {
      label <- rep(list(""), length.out=n);
   }

   # TODO: recognize color as list type,
   # so each block arrow can use a gradient color fill.
   # border should still be vector with length n, one border per block arrow.
   if (n == 1) {
      if (is.atomic(label)) {
         label <- list(label)
      }
   } else if (is.atomic(label)) {
      label <- as.list(label)
   }

   # define polygon coordinates
   arrow_coords <- make_block_arrow_polygon(
      x=x0,
      x1=x1,
      y=y0,
      y1=y1,
      data_format=plot_type,
      # color=color,
      # border=border,
      ...)

   # render polygon(s)
   if ("base" %in% plot_type) {
      if (TRUE %in% do_plot) {
         graphics::polygon(x=arrow_coords$x,
            y=arrow_coords$y,
            col=color,
            border=border)
      }
   }

   # optional text label
   if (length(label) > 0) {
      # ensure vectorized values are proper lengths
      label <- rep(label, length.out=length(x0))
      label_cex <- rep(label_cex, length.out=length(x0))
      label_color <- rep(label_color, length.out=length(x0))
      label_font <- rep(label_font, length.out=length(x0))
      oneway_position <- rep(oneway_position, length.out=length(x0))
      oneway_position <- jamba::noiseFloor(oneway_position,
         minimum=0,
         ceiling=1);

      # base R graphics
      if ("base" %in% plot_type) {
         arrow_angle <- jamba::rad2deg(atan2(y=(y1 - y0), x=(x1 - x0))) %% 360;
         # round to integer values
         label_angle <- round(((arrow_angle + 90) %% 180 - 90) %% 360);
         # label_angle <- arrow_angle;
         angle_flip <- (arrow_angle != label_angle)
         if (verbose) {
            jamba::printDebug("draw_oneway_contrast(): ",
               "x0: ", x0,
               ", x: ", (x0 + x1) / 2,
               ", x1: ", x1,
               ", y0: ", y0,
               ", y: ", (y0 + y1) / 2,
               ", y1: ", y1,
               ", arrow_angle: ", arrow_angle,
               ", label_angle: ", label_angle)
         }
         use_label <- jamba::cPaste(label,
            sep=label_sep,
            na.rm=TRUE)
         which_labels <- which(!is.na(label) & nchar(label) > 0);
         if (length(which_labels) > 0) {
            # k_list <- split(which_labels, label_angle[which_labels])
            label_split <- paste0(label_angle[which_labels], "_",
               oneway_position[which_labels], "_",
               as.character(angle_flip[which_labels]))
            # print(data.frame(label_split, label_angle[which_labels], oneway_position[which_labels], which_labels));# debug
            k_list <- split(which_labels, label_split)
            # TODO: use jamba::drawLabels() with rotation and shadowText
            # or gridtext::richtext_grob() with rounded corners and styling
            for (k in k_list) {
               use_oneway_position <- 1 - head(oneway_position[k], 1)
               use_label_angle <- head(label_angle[k], 1);
               use_arrow_angle <- head(arrow_angle[k], 1);
               use_angle_flip <- head(angle_flip[k], 1);
               use_adj <- c(
                  ifelse(use_angle_flip,
                     round(use_oneway_position * 2) / 2,
                     1 - (round(use_oneway_position * 2) / 2)),
                  0.5);
               jamba::shadowText(
                  x=((x0 * ((1-oneway_position)*0.9 + 0.1) + x1 * (1 - ((1-oneway_position)*0.9 + 0.1))) )[k],
                  y=((y0 * ((1-oneway_position)*0.9 + 0.1) + y1 * (1 - ((1-oneway_position)*0.9 + 0.1))) )[k],
                  col=label_color[k],
                  srt=use_label_angle,
                  cex=label_cex[k],
                  alphaOutline=getOption("jam.alphaOutline", 0.95),
                  r=getOption("jam.shadow.r", 0.1),
                  n=getOption("jam.shadow.n", 8),
                  adj=use_adj,
                  font=label_font[k],
                  jamba::cPaste(label[k], sep="\n"))
               if (debug) {
                  jamba::printDebug("use_adj:", use_adj,
                     ", use_oneway_position:", use_oneway_position,
                     ", use_label_angle:", use_label_angle,
                     ", use_arrow_angle:", use_arrow_angle,
                     ", use_angle_flip:", use_angle_flip,
                     ", label[k]:", label[k]);# debug
               }
            }
         }
      } else if ("grid" %in% plot_type) {
         # TODO: consider using shadowtext::grid.shadowtext()
      }
   }
   return(invisible(arrow_coords));
}

#' Draw two one-way contrasts using block arrows, showing a two-way connector
#'
#' Draw two one-way contrasts using block arrows, showing a two-way connector
#'
#' This function essentially draws two one-way contrasts
#' in the form: `(group1-group2)-(group3-group4)`
#'
#' This two-way contrast is represented by two one-way contrasts:
#' 1. group1-group2
#' 2. group3-group4
#'
#' This function renders two individual one-way contrasts, then
#' draws a connector from the end of group2, to the beginning of group3.
#'
#' ## TODO:
#'
#' * Change the order of drawing:
#'
#'    1. Draw the two-way connector border.
#'    2. Draw the one-way contrasts.
#'    3. Draw the two-way connector fill.
#'    4. These steps would ensure the connector line does not
#'    overlap the one-way contrasts, and would allow the connector to
#'    connect directly to the contrast at the ends of the block arrows.
#'
#'
#' @inheritParams draw_oneway_contrast
#' @param extend_ex `numeric` which defines the amount to extend each
#'    contrast arrow when defining control points for a bezier spline
#'    from the end of one contrast, to the beginning of the next contrast.
#'    The value `extend_ex` extends the contrast by 20% of the contrast
#'    line length.
#' @param extend_angle `numeric` value in degrees (values from 0 to 360)
#'    which defines the relative angle from the first contrast, toward
#'    the second contrast, when defining control points for the bezier
#'    spline. The default value `extend_angle=15` angles the extended line
#'    15 degrees toward the next contrast, starting at the end of the
#'    first contrast. The same angle is used, rotated 180 degrees, for
#'    the second control point at the beginning of the second contrast.
#' @param twoway_label `character` or `list` handled the same as `label`.
#'    This label is placed atop the two-way connector, at an angle aligned
#'    with the middle of the connector curved line.
#' @param twoway_label_color `character` color used for the `label`.
#' @param twoway_label_cex `numeric` label font expansion factor, used to adjust
#'    the font size of the text label.
#' @param twoway_label_font `numeric` indicating the font face, defined as:
#'    * 1 = normal font
#'    * 2 = bold font
#'    * 3 = italic font
#'    * 4 = bold, italic font
#'
#'
#' @examples
#' plot(NULL, xlim=c(0, 5), ylim=c(0, 4), asp=1, xlab="", ylab="")
#' draw_twoway_contrast(c(1, 1), c(4, 4), c(1, 2), c(1, 2),
#'    label=c("contrast one", "contrast two"))
#'
#' plot(NULL, xlim=c(0, 5), ylim=c(0, 6), asp=1, xlab="", ylab="")
#' draw_twoway_contrast(c(1, 2, 1, 2, 2, 2), c(4, 5, 4, 5, 5, 5),
#'    c(1, 1.3, 2, 3, 4, 5), c(1, 1.3, 2, 3, 4, 5),
#'    label=c("contrast one", "contrast two"),
#'    extend_ex=0.5, extend_angle=10, color=colorjam::rainbowJam(6))
#'
#' plot(NULL, xlim=c(0, 5), ylim=c(0, 6), asp=1, xlab="", ylab="")
#' draw_twoway_contrast(c(1, 2, 1, 2, 2, 2), c(4, 5, 4, 5, 5, 5),
#'    c(1, 1.3, 2, 2.2, 4, 5), c(1, 1.3, 2, 2.2, 4, 5),
#'    label=c("contrast one", "contrast two"),
#'    extend_ex=0.3, extend_angle=20, color=colorjam::rainbowJam(6))
#'
#' plot(NULL, xlim=c(-2, 8), ylim=c(1, 6), asp=1, xlab="", ylab="")
#' draw_twoway_contrast(
#'    x0=c(1, 2, 1, 2, 1, 1, 0, -2, 4, 3, 2, 4),
#'    x1=c(4, 5, 4, 5, 6, 6, 0, -2, 7, 6, 2, 4),
#'    y0=c(1, 1.3, 2, 2.2, 4, 5, 5, 5, 3.2, 3, 6, 6),
#'    y1=c(1, 1.3, 2, 2.2, 4, 5, 3, 3, 3.2, 3, 3.5, 3.5),
#'    label=c("contrast one", "contrast two"), verbose=TRUE,
#'    extend_ex=0.5, extend_angle=10,
#'    color=colorjam::rainbowJam(12, Crange=c(90, 120)))
#'
#' draw_twoway_contrast(x0=c(4, 3), x1=c(4, 3), y0=c(2, 2), y1=c(1, 1),
#'    label=c("contrast one", "contrast two"),
#'    color=c("gold", "dodgerblue"), extend_angle=-20)
#' @export
draw_twoway_contrast <- function
(x0,
 x1,
 y0,
 y1,
 color="peachpuff",
 border="black",
 extend_ex=0.3,
 extend_angle=10,
 plot_type=c("base",
    "grid"),
 label=NULL,
 label_sep="\n",
 na.rm=TRUE,
 label_color="white",
 label_cex=1,
 label_font=2,
 oneway_position=0.5,
 twoway_label=NULL,
 twoway_label_color=label_color,
 twoway_label_cex=label_cex,
 twoway_label_font=label_font,
 twoway_position=0.5,
 twoway_lwd=5,
 contingency=c("loop",
    "scrunch"),
 draw_oneway=TRUE,
 drawing_order=c(
    "two-two-one",
    "one-two",
    "two-one-two"),
 do_plot=TRUE,
 verbose=FALSE,
 ...)
{
   # validate arguments
   plot_type <- match.arg(plot_type);
   contingency <- match.arg(contingency);
   drawing_order <- match.arg(drawing_order);

   # ensure vectorized values are proper lengths
   n <- max(c(length(x0), length(x1), length(y0), length(y1)));
   if (n == 0) {
      # consider adding a warning
      return(NULL)
   }
   x0 <- rep(x0, length.out=n)
   x1 <- rep(x1, length.out=n)
   y0 <- rep(y0, length.out=n)
   y1 <- rep(y1, length.out=n)
   color <- rep(color, length.out=n);
   colorset_twoway <- split(color, rep(seq_len(n / 2), each=2))
   border <- rep(border, length.out=n);

   # first draw the two one-way contrasts
   if (TRUE %in% draw_oneway && "one-two" %in% drawing_order) {
      doc <- draw_oneway_contrast(
         x0=x0,
         x1=x1,
         y0=y0,
         y1=y1,
         color=color,
         border=border,
         plot_type=plot_type,
         label=label,
         label_sep=label_sep,
         na.rm=na.rm,
         label_color=label_color,
         label_cex=label_cex,
         label_font=label_font,
         do_plot=do_plot,
         verbose=verbose,
         ...)
   }

   # index values for first and second contrast
   k1 <- seq(from=1, to=length(x0), by=2)
   k2 <- k1 + 1;

   # determine the two-way connector
   contrast_angles <- jamba::rad2deg(atan2(y=(y1 - y0),
      x=(x1 - x0)));
   # flip each second contrast
   contrast_angles[k2] <- (contrast_angles[k2] + 180) %% 360

   # determine handedness of the two contrasts
   hand <- point_handedness(
      pt1=cbind(c(x0[k1]), c(y0[k1])),
      pt2=cbind(c(x0[k2]), c(y0[k2])),
      angle=contrast_angles[k1])
   if (verbose) {
      jamba::printDebug("draw_twoway_contrast(): ",
         "hand: ", hand);
   }
   contrast_lengths <- sqrt(
      (x1 - x0)^2 + (y1 - y0)^2)

   # determine appropriate angle
   angle_diff <- -(extend_angle);
   angle_diff <- ifelse(rep(hand, each=2) %in% "right",
      extend_angle,
      -extend_angle)

   new_contrast_angles <- contrast_angles + angle_diff;
   if (verbose) {
      jamba::printDebug("draw_twoway_contrast(): ",
         "new_contrast_angles: ");
      print(data.frame(contrast_angles, hand=rep(hand, each=2),
         angle_diff, new_contrast_angles));
   }

   # make bezier control points
   make_bezier_lines_df <- function
   (x0, x1, y0, y1, k1, k2, contrast_lengths, extend_ex, xadj, yadj)
   {
      # Note: using length=1 instead of contrast_lengths[k1]
      # so the extension has consistent length.
      lines_df <- data.frame(
         x=c(
            # x0=x0[k1],
            `x0.5`=(x0[k1] +
                  x1[k1])/2,
            `x0.75`=(x0[k1] +
                  x1[k1]*3)/4,
            x1=x1[k1],
            `x1.2`=(x1[k1] +
                  xadj[k1] * 1 * extend_ex),
            `x1.5`=(x1[k1] +
                  x0[k2])/2,
            `x1.8`=(x0[k2] +
                  xadj[k2] * 1 * extend_ex),
            x2=x0[k2],
            `x2.25`=(x0[k2]*3 +
                  x1[k2])/4,
            `x2.5`=(x0[k2] +
                  x1[k2])/2),
         # x3=x1[k2]),
         y=c(
            # y0=y0[k1],
            `y0.5`=(y0[k1] +
                  y0[k1])/2,
            `y0.75`=(y0[k1] +
                  y0[k1]*3)/4,
            y1=y1[k1],
            `x1.2`=(y1[k1] +
                  yadj[k1] * 1 * extend_ex),
            `y1.5`=(y1[k1] +
                  y0[k2])/2,
            `y1.8`=(y0[k2] +
                  yadj[k2] * 1 * extend_ex),
            y2=y0[k2],
            `y2.25`=(y0[k2]*3 +
                  y1[k2])/4,
            `y2.5`=(y0[k2] +
                  y1[k2])/2)
      );
      return(lines_df);
   }

   # optionally mirror the angle when lines are too close?
   # mirror the first contrast? might be useful when lines are too close
   # new_contrast_angles[k1] <- (360 - new_contrast_angles[k1]) %% 360;
   # mirror the second contrast? might be useful when lines are too close
   # new_contrast_angles[k2] <- (360 - new_contrast_angles[k2]) %% 360;
   # new_contrast_angles <- (360 - new_contrast_angles) %% 360;

   # jamba::printDebug("new_contrast_angles:", new_contrast_angles); # debug
   xadj <- cos(jamba::deg2rad(new_contrast_angles))
   yadj <- sin(jamba::deg2rad(new_contrast_angles))

   # define bezier control points for typical cases
   lines_df <- make_bezier_lines_df(x0, x1, y0, y1,
      k1, k2,
      contrast_lengths, extend_ex,
      xadj, yadj);

   # check for weird cases:
   # - when control point x1.2 crosses either the x1.5 or x2.5
   # - test whether "handedness" changes relative to those points
   x12_rows <- jamba::igrep("^x1.2", rownames(lines_df));
   x15_rows <- jamba::igrep("^x1.5", rownames(lines_df));
   x25_rows <- jamba::igrep("^x2.5", rownames(lines_df));
   new_hand_25 <- point_handedness(
      pt1=as.matrix(lines_df[x12_rows, , drop=FALSE]),
      pt2=as.matrix(lines_df[x25_rows, , drop=FALSE]),
      angle=contrast_angles[k1])
   new_hand_15 <- point_handedness(
      pt1=as.matrix(lines_df[x12_rows, , drop=FALSE]),
      pt2=as.matrix(lines_df[x15_rows, , drop=FALSE]),
      angle=contrast_angles[k1])
   # if it crosses x25, we should loop
   # if it crosses x15 but not x25, maybe we scrunch
   changed_hands_25 <- (hand != new_hand_25);
   changed_hands_15 <- (hand != new_hand_15);
   if (any(c(changed_hands_25, changed_hands_15))) {
      # mirror the second angle?
      if (any(changed_hands_25)) {
         k1_changed <- k1[changed_hands_25];
         k2_changed <- k2[changed_hands_25];
         new_contrast_angles[k2_changed] <- (360 -
               new_contrast_angles[k2_changed]) %% 360;
      }
      if (any(changed_hands_15 & !changed_hands_25)) {
         k1_changed <- k1[changed_hands_15 & !changed_hands_25];
         k2_changed <- k2[changed_hands_15 & !changed_hands_25];
         new_contrast_angles[k1_changed] <- (contrast_angles[k1_changed] +
               angle_diff / 3);
         new_contrast_angles[k2_changed] <- (contrast_angles[k2_changed] +
               angle_diff / 3);
      }
      xadj <- cos(jamba::deg2rad(new_contrast_angles))
      yadj <- sin(jamba::deg2rad(new_contrast_angles))
      # define bezier control points for typical cases
      lines_df <- make_bezier_lines_df(x0, x1, y0, y1,
         k1, k2,
         contrast_lengths, extend_ex,
         xadj, yadj);
   }

   # optionally adjust x1.5,y1.5 to be midpoint between x1.2-x1.8 and y1.2-y1.8
   x12_rows <- jamba::igrep("^x1.2", rownames(lines_df));
   x15_rows <- jamba::igrep("^x1.5", rownames(lines_df));
   x18_rows <- jamba::igrep("^x1.8", rownames(lines_df));
   lines_df[x15_rows, "x"] <- (lines_df[x12_rows, "x"] + lines_df[x18_rows, "x"]) / 2;
   lines_df[x15_rows, "y"] <- (lines_df[x12_rows, "y"] + lines_df[x18_rows, "y"]) / 2;

   # bezier::bezier() is not vectorized, so iterate each set of control points
   k_multiple <- (length(x0) / 2)
   # iterate to calculate each bezier curve
   bezier_curves <- lapply(seq_len(k_multiple), function(k){
      k_seq <- seq(from=k, length.out=9, by=k_multiple)
      use_lines_df <- lines_df[k_seq, , drop=FALSE]

      # optionally display the control points with labels for debugging
      if (FALSE) {
         points(use_lines_df,
            pch=21,
            cex=4,
            col="grey",
            bg="white")
         text(use_lines_df, rownames(lines_df), pch=21, cex=0.8)
      }

      # calculate bezier curve
      bezier_points <- bezier::bezier(
         t=seq(from=0, to=4, length.out=100),
         p=as.matrix(use_lines_df), deg=2)
      bezier_points
   })

   # iterate to draw two-way border lines for each connector
   for (k in seq_along(bezier_curves)) {
      bezier_points <- bezier_curves[[k]];

      # render wider line to display a border
      if ("one-two" %in% drawing_order) {
         ik <- 25:76
      } else {
         ik <- 24:76;
      }
      bezier_line_use <- bezier_points[ik, , drop=FALSE];
      lines(bezier_line_use,
         lwd=twoway_lwd + 2,
         col=border[k * 2 - 1]);
   }

   # optionally draw one-way contrasts here
   if (TRUE %in% draw_oneway && "two-one-two" %in% drawing_order) {
      doc <- draw_oneway_contrast(
         x0=x0,
         x1=x1,
         y0=y0,
         y1=y1,
         color=color,
         border=border,
         plot_type=plot_type,
         label=label,
         label_sep=label_sep,
         na.rm=na.rm,
         label_color=label_color,
         label_cex=label_cex,
         label_font=label_font,
         oneway_position=oneway_position,
         do_plot=do_plot,
         verbose=verbose,
         ...)
   }

   # iterate to draw two-way border lines for each connector
   for (k in seq_along(bezier_curves)) {
      bezier_points <- bezier_curves[[k]];

      # render wider line to display a border
      if ("one-two" %in% drawing_order) {
         ik <- 23:77;
      } else {
         ik <- 23:77;
      }
      # render thinner line with colored segments
      bezier_use <- bezier_points[ik, , drop=FALSE];
      graphics::segments(
         x0=head(bezier_use, -1)[,1],
         x1=tail(bezier_use, -1)[,1],
         y0=head(bezier_use, -1)[,2],
         y1=tail(bezier_use, -1)[,2],
         col=jamba::getColorRamp(
            # c("dodgerblue", "firebrick2"),
            colorset_twoway[[k]],
            n=nrow(bezier_use) - 1),
         lwd=twoway_lwd)
   }

   # optionally draw one-way contrasts here
   if (TRUE %in% draw_oneway && "two-two-one" %in% drawing_order) {
      doc <- draw_oneway_contrast(
         x0=x0,
         x1=x1,
         y0=y0,
         y1=y1,
         color=color,
         border=border,
         plot_type=plot_type,
         label=NULL,
         label_sep=label_sep,
         na.rm=na.rm,
         label_color=label_color,
         label_cex=label_cex,
         label_font=label_font,
         oneway_position=oneway_position,
         do_plot=do_plot,
         verbose=verbose,
         ...)
      # optional "flourish" to re-draw the start/end
      # iterate to draw two-way border lines for each connector
      for (k in seq_along(bezier_curves)) {
         bezier_points <- bezier_curves[[k]];

         # render wider line to display a border
         ik <- 23:77;
         # then take first 5, last 5 points
         # render thinner line with colored segments
         bezier_use <- bezier_points[ik, , drop=FALSE];
         bezier_use_start <- head(bezier_use, 5);
         bezier_use_end <- tail(bezier_use, 5);
         bezier_colors <- jamba::getColorRamp(
            colorset_twoway[[k]],
            n=nrow(bezier_use) - 1);
         graphics::segments(
            x0=head(bezier_use_start, -1)[,1],
            x1=tail(bezier_use_start, -1)[,1],
            y0=head(bezier_use_start, -1)[,2],
            y1=tail(bezier_use_start, -1)[,2],
            col=head(bezier_colors, 5),
            lwd=twoway_lwd)
         graphics::segments(
            x0=head(bezier_use_end, -1)[,1],
            x1=tail(bezier_use_end, -1)[,1],
            y0=head(bezier_use_end, -1)[,2],
            y1=tail(bezier_use_end, -1)[,2],
            col=tail(bezier_colors, 5),
            lwd=twoway_lwd)
      }

      # last step is to draw labels so they are not covered by the bezier curve
      # oneway labels
      doc <- draw_oneway_contrast(
         x0=x0,
         x1=x1,
         y0=y0,
         y1=y1,
         color=NA,
         border=NA,
         plot_type=plot_type,
         label=label,
         label_sep=label_sep,
         na.rm=na.rm,
         label_color=label_color,
         label_cex=label_cex,
         label_font=label_font,
         oneway_position=oneway_position,
         do_plot=do_plot,
         verbose=verbose,
         ...)
      # twoway labels
      for (k in seq_along(bezier_curves)) {
         bezier_points <- bezier_curves[[k]];

         # define where to place the label
         kn <- floor(twoway_position * 50) + 25;
         # grab points in the middle
         klen <- 7;
         kn_start <- min(c(max(c(kn - klen, 15)), 85 - klen * 2))
         ik <- seq(kn_start, length.out=15);
         kn <- floor(mean(ik));
         # ik <- 43:57;

         bezier_use <- bezier_points[ik, , drop=FALSE];
         # optionally render thinner line with colored segments
         # graphics::segments(
         #    x0=head(bezier_use, -1)[,1],
         #    x1=tail(bezier_use, -1)[,1],
         #    y0=head(bezier_use, -1)[,2],
         #    y1=tail(bezier_use, -1)[,2],
         #    col="grey25",
         #    lwd=twoway_lwd * 2)
         vector_angle <- jamba::rad2deg(atan2(
            y=(tail(bezier_use, 1)[, 2] - head(bezier_use, 1)[, 2]),
            x=(tail(bezier_use, 1)[, 1] - head(bezier_use, 1)[, 1]))) %% 360;
         # round to integer values
         twoway_label_angle <- round(
            ((vector_angle + 90) %% 180 - 90) %% 360);
         use_angle_flip <- (twoway_label_angle != vector_angle)
         use_adj <- c(
            ifelse(use_angle_flip,
               1 - (round(twoway_position * 8) / 8),
               round(twoway_position * 8) / 8),
            0.5);
         jamba::shadowText(
            x=bezier_points[kn, 1],
            y=bezier_points[kn, 2],
            col=twoway_label_color,
            srt=twoway_label_angle,
            cex=twoway_label_cex,
            alphaOutline=getOption("jam.alphaOutline", 0.95),
            r=getOption("jam.shadow.r", 0.1),
            n=getOption("jam.shadow.n", 8),
            # adj=c(0.5, 0.5),
            adj=use_adj,
            font=twoway_label_font,
            jamba::cPaste(twoway_label, sep="\n"))

      }


   }


   invisible(NULL)
}
