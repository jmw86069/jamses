
#' Detect ComplexHeatmap Heatmap grid layout components
#'
#' Detect ComplexHeatmap Heatmap grid layout components
#'
#' This function is a wrapper around `ComplexHeatmap::list_components()`
#' which also creates a `list` of component names based upon common
#' naming conventions used in the `ComplexHeatmap` package.
#'
#' @return `list` of `character` vectors, where `list` names represent
#'    different features of the heatmap, and each `character` vector
#'    includes the grid layout component name stem suitable for use
#'    in `heatmap_column_group_labels()`.
#'
#' @param ... additional arguments are ignored.
#'
#' @export
detect_heatmap_components <- function
(...)
{
   # column title components
   hm_components <- ComplexHeatmap::list_components();
   column_title_components <- unique(gsub("_[0-9]+$", "_",
      jamba::vigrep("column_title_[1-9]", hm_components)))

   # row title components
   row_title_components <- unique(gsub("_[0-9]+$", "_",
      jamba::vigrep("row_title_[1-9]", hm_components)))

   # heatmap body components
   heatmap_body_components <- unique(gsub("_[0-9]+$", "_",
      jamba::vigrep("body_[1-9]", hm_components)))

   # annotation components
   annotation_components <- unique(gsub("_[0-9]+$", "_",
      jamba::vigrep("annotation_.+_[1-9]", hm_components)))

   # assemble into a list
   component_sets <- list();
   if (length(column_title_components) > 0) {
      component_sets$column_title_components <- column_title_components
   }
   if (length(row_title_components) > 0) {
      component_sets$row_title_components <- row_title_components
   }
   if (length(heatmap_body_components) > 0) {
      component_sets$heatmap_body_components <- heatmap_body_components
   }
   if (length(annotation_components) > 0) {
      component_sets$annotation_components <- annotation_components
   }
   return(component_sets)
}


#' Add Heatmap column group labels
#'
#' Add Heatmap column group labels, specifically for ComplexHeatmap output.
#'
#' This function is currently experimental, and is intended only for
#' a specific scenario, to augment a `ComplexHeatmap::Heatmap` object,
#' as produced by `heatmap_se()`, that used the argument `column_split`
#' to sub-divide the heatmap columns into subgroups.
#' This function draws labels above the heatmap to describe group
#' factor values associated with column splits.
#'
#' When there are multiple layers of grouping, this function will
#' also draw multiple layers. For example, when columns are split
#' by `column_split=c("Treatment", "Time")`, it will produce a heatmap
#' where each column slice has one combination of Treatment and Time.
#' This function can be used to add a layer of labels by `"Time"`,
#' and a layer of labels by `"Treatment"`. The labels are shown by default
#' with a broad underline to indicate contiguous column slices that
#' contain the same `"Treatment"` or `"Time"` values.
#'
#' The output can provide a cleaner visualization than the alternative of
#' displaying colorized annotation boxes at the top of the heatmap.
#'
#' ## Input Formats
#'
#' There are two strategies for defining the column group data to display:
#'
#' ### Input using colnames from colData(se)
#'
#' This option is recommended as the "easiest" method. It requires:
#'
#' * `hm_group_list` provides a `character` vector of `colnames(colData(se))`.
#' Note that the columns are applied bottom-to-top, so it is sometimes
#' helpful to supply columns in reverse order for `column_split`.
#' * `se` must be provided, since it supplies `colData(se)`
#'
#' Each value in `hm_group_list` is applied from bottom-to-top, to define
#' a row of labels. By default, `add_group_line=TRUE`, so labels are also
#' underlined to indicate samples included in each group.
#'
#' ### Input using a named list
#'
#' This option is recommended when labeling a subset of column slices,
#' or when the column slice groups need to be customized in some way.
#'
#' ## Additional requirements
#'
#' * The Heatmap `column_title` should (usually) be empty, so that no text
#' labels are drawn which may overlap the labels drawn by this function.
#' Use `column_title=" "` (with a whitespace character) to prevent
#' `ComplexHeatmap::Heatmap()` from using its internal default text label.
#' * There usually needs to be whitespace above the heatmap, which
#' can be accomplished when drawing a Heatmap object `hm` like this:
#' `ComplexHeatmap::draw(hm, column_title="\n\n\n\n")`
#' * When also displaying a heatmap title as defined by `heatmap_se()`,
#' the adjustment can be done like this:
#' ```r
#' ComplexHeatmap::draw(hm,
#'    column_title=paste0(attr(hm, "hm_title"),
#'       "\n\n\n\n"))
#' ```
#'
#' * The arguments `hm_title_base` and `hm_body_base` should match the
#' heatmap name, which is defined in `heatmap_se()` with `data_type`,
#' usually prefixed `"centered\n"` when the data is centered by that
#' function. For example, when `data_type="abundance"`, the corresponding
#' argument value should be
#' `hm_title_base="centered\nabundance_column_title_"`.
#' * use `detect_heatmap_components()` to help identify the appropriate
#' `grid` elements available to be used by this function.
#'
#' When the heatmap contains `"column_title"` elements defined in
#' `ComplexHeatmap::list_components()`, and there is no element
#' `"global_column_title"`, the `y_offset_lines` is adjusted down so that
#' the position is inside the column_title region, typically below the
#' column_title when using trailing whitespace (see argument
#' `hm_title_buffer` in `heatmap_se()`). However, when element
#' `"global_column_title"` is present, the `y_offset_lines` are not adjusted,
#' so that the position is above the column_title region. In this case, to
#' position the labels below the column_title, you can adjust
#' `y_offset_lines` down manually like this: `y_offset_lines=-9`.
#'
#' ## Todo
#'
#' * Automate determining the column_title grid layout name.
#'
#'    * Specifically when only one `heatmapname_column_title_1` is present,
#'    but there are multiple column groups, it should take the x-axis
#'    coordinate values (left and right boundaries) from the heatmap body
#'    instead of the column_title region.
#'
#' * When `add_group_box=TRUE`, and `row_split` indicates multiple row
#' groups, it should calculate the y-axis coordinate values (top and bottom
#' boundaries) using the full set of row groups.
#' * Enable blank annotations, either by passing a subset `se`, or by
#' annotations with no associated label.
#'
#'
#' @param hm_group_list `character` or `list` with one of the
#'    following types of content:
#'    * `character` colnames in the column data of `se` argument, specifically
#'    `colnames(colData(se))` to define groupings.
#'    * `list` of `integer` vectors, where list names become labels for
#'    each group.
#' @param hm_drawn `ComplexHeatmap::HeatmapList` object which is returned
#'    by the function `ComplexHeatmap::draw()`. This object is a pointer
#'    to the grid graphical elements required.
#' @param se `SummarizedExperiment` object, required only when `hm_group_list`
#'    is supplied in the form of colnames of the
#'    `SummarizedExperiment::colData(se)`.
#' @param add_group_label `logical` indicating whether to draw the group
#'    `character` label above each group defined in `hm_group_list`.
#' @param add_group_line `logical` indicating whether to draw the group
#'    line for each group defined in `hm_group_list`. The line is intended
#'    to appear below the text label when `add_group_label=TRUE`.
#' @param add_group_box `logical` indicating whether to draw a box around
#'    the heatmap region defined by each group in `hm_group_list`.
#' @param group_line_buffer `grid::unit` object indicating the whitespace
#'    buffer region betwen adjacent lines when `add_group_line=TRUE`.
#'    The default enforces 1-mm white space between lines so they
#'    do not touch.
#' @param group_line_lwd `numeric` line width when `add_group_line=TRUE`.
#' @param group_line_requires_label `logical` whether to require
#'    the associated group label to contain non-whitespace visible
#'    text.
#'    * `group_line_requires_label=TRUE` (default) requires group label
#'    to have visible characters, which means the presence of an
#'    empty labels will cause the group line not to be drawn.
#'    * `group_line_requires_label=FALSE` does not require a group label,
#'    therefore all group lines are drawn.
#' @param group_box_lwd `numeric` line width when `add_group_box=TRUE`.
#' @param group_box_outer `logical` indicating whether to draw the group
#'    box as an outer border, which means the border will be drawn only
#'    on the very outer edge of the heatmap body, and will not overlap
#'    the heatmap contents. This option is only active when
#'    `add_group_box=TRUE`.
#' @param font_cex `numeric` adjustment for font sizes overall.
#' @param hm_title_base `character` string used to search for matching
#'    heatmap column title grid layout regions. It is derived from the
#'    `Heatmap` argument `name`, usually followed by `"_column_title_".
#'    When `NULL` the default values are defined by
#'    `detect_heatmap_components()` element `column_title_components`.
#' @param hm_body_base `character` string used to search for matching
#'    heatmap body grid layout regions.  It is derived from the
#'    `Heatmap` argument `name`, usually followed by `"_body_1_".
#'    When `NULL` the default values are defined by
#'    `detect_heatmap_components()` element `heatmap_body_components`.
#'    Note: Currently the group box only includes the first row_split,
#'    although it is intended to include the entire set of heatmap
#'    rows in future iterations.
#' @param hm_title_base_default,hm_body_base_default `character` string
#'    used as reference for hm_title_base,hm_body_base.
#' @param y_offset_lines `numeric` adjustment used to shift the text
#'    label and underline by this many lines of character height.
#'    It is mainly used internally for iterative calls, when `hm_group_list`
#'    includes multiple layers of groups. The text and underline are
#'    intended to fill 2 character height lines, although this argument
#'    can be used to make manual adjustments.
#' @param endlines `integer` number of blank lines to append to the end
#'    of each group label, which has the effect of shifting the group
#'    label upwards slightly.
#' @param use_gridtext `logical` (default TRUE) whether to render text using
#'    `gridtext::richtext_grob()` in order to enable markdown and
#'    limited HTML-based formatting.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @export
heatmap_column_group_labels <- function
(hm_group_list,
 hm_drawn=NULL,
 se=NULL,
 add_group_label=TRUE,
 add_group_line=TRUE,
 add_group_box=FALSE,
 group_line_buffer=grid::unit(1, "mm"),
 group_line_lwd=2,
 group_line_requires_label=TRUE,
 group_box_lwd=2,
 group_box_outer=TRUE,
 font_cex=1,
 hm_title_base=NULL,
 hm_body_base=NULL,
 hm_title_base_default="centered\nexpression_column_title_",
 hm_body_base_default="centered\nexpression_heatmap_body_1_",
 y_offset_lines=0,
 endlines=1,
 use_gridtext=TRUE,
 verbose=FALSE,
 ...)
{
   #
   # tip: find column_title names with
   # vigrep("column_title_[1-9]", list_components())
   #
   # tip: find heatmap_body names with
   # vigrep("body_[1-9]", list_components())
   #
   if (length(group_line_buffer) == 0) {
      group_line_buffer <- grid::unit(0, "mm")
   }
   if (length(group_box_outer) == 0) {
      group_box_outer <- FALSE
   }
   if (is.logical(group_box_outer)) {
      group_box_outer <- group_box_outer * 1;
   }

   # try to detect heatmap components when possible
   dhc <- detect_heatmap_components();
   if (length(hm_title_base) == 0 || length(hm_body_base) == 0) {
      if (length(hm_title_base) == 0) {
         if (verbose) {
            jamba::printDebug("heatmap_column_group_labels(): ",
               "Assigning hm_title_base using detect_heatmap_components(): '",
               dhc$column_title_components, "'")
         }
         hm_title_base <- dhc$column_title_components;
      }
      if (length(hm_body_base) == 0) {
         if (verbose) {
            jamba::printDebug("heatmap_column_group_labels(): ",
               "Assigning hm_body_base using detect_heatmap_components(): '",
               dhc$heatmap_body_components, "'")
         }
         hm_body_base <- dhc$heatmap_body_components;
         # 0.0.64.900 - warning about multiple entries in grep() used later
         hm_body_base <- head(hm_body_base, 1);
      }
   }

   # determine whether there is also a "global_column_title"
   # - if no global title, assume the hm_title is applied to column_title
   hlc <- ComplexHeatmap::list_components();
   if (!jamba::igrepHas("global_column_title", hlc)) {
      y_offset_lines <- y_offset_lines - 4;
   }

   # small helper function
   bbv2seq <- function(bbv, ...){
      bbvl <- lapply(seq_along(bbv$breakPoints), function(i){
         j1 <- 1;
         if (i > 1) {
            j1 <- bbv$breakPoints[[i - 1]] + 1
         }
         seq(j1, bbv$breakPoints[[i]])
      })
      names(bbvl) <- bbv$useLabels;
      bbvl;
   }

   # hm_group_list can be a list of integer vectors, named by label
   # hm_group_list can be a character vector of colnames(colData(se))
   # - hm_drawn and se must also be supplied
   if ("character" %in% class(hm_group_list)) {
      if (length(se) == 0) {
         stop(paste0("When hm_group_list is character colnames(colData(se)), ",
            "se must be supplied."))
      }
      if (length(hm_drawn) == 0) {
         stop(paste0("When hm_group_list is character colnames(colData(se)), ",
            "hm_drawn must be supplied."))
      }
      icolnames_list <- jamba::heatmap_column_order(hm_drawn);
      icolnames_list_num <- rep(seq_along(icolnames_list),
         lengths(icolnames_list))
      icolnames <- unlist(unname(icolnames_list));
      # iterate through each colData column
      # hm_group_list <- c("Treatment", "IP")
      for (icol in hm_group_list) {
         # get column values
         gv1 <- SummarizedExperiment::colData(
            se[,icolnames])[,c(icol)]
         # convert values to character at this step
         bbv1 <- jamba::breaksByVector(as.character(gv1))
         bbvl1 <- bbv2seq(bbv1)
         bbvl1num <- lapply(bbvl1, function(ibbvl1){
            unique(icolnames_list_num[ibbvl1])
         })
         # call this function
         if (verbose) {
            jamba::printDebug("heatmap_column_group_labels(): ",
               "Calling function for column: ", icol);
         }
         heatmap_column_group_labels(
            hm_group_list=bbvl1num,
            add_group_label=add_group_label,
            add_group_line=add_group_line,
            group_line_requires_label=group_line_requires_label,
            add_group_box=add_group_box,
            group_line_buffer=group_line_buffer,
            group_line_lwd=group_line_lwd,
            group_box_lwd=group_box_lwd,
            font_cex=font_cex,
            y_offset_lines=y_offset_lines,
            hm_title_base=hm_title_base,
            hm_body_base=hm_body_base,
            ...)
         y_offset_lines <- y_offset_lines + 2;
         add_group_box <- FALSE;
      }
      return(invisible(NULL))
   }

   # define grid layout components
   all_hm_components <- ComplexHeatmap::list_components();
   # printDebug("hm_group_list: ");print(hm_group_list);
   # use_hm_components <- paste0(hm_title_base, seq_along(hm_group_list))
   use_hm_components <- paste0(hm_title_base,
      sort(unique(unlist(hm_group_list))))
   # use_hm_body_components <- paste0(hm_body_base, seq_along(hm_group_list))
   use_hm_body_components <- paste0(hm_body_base,
      sort(unique(unlist(hm_group_list))))

   if (any(grepl("_1_$", hm_body_base))) {
      hm_body_test <- gsub("_1_$", "_2_", hm_body_base)
      use_hm_test_components <- paste0(hm_body_test, seq_along(hm_group_list))
      if (all(use_hm_test_components %in% all_hm_components)) {
         hm_body_base_t <- use_hm_body_components;
         hm_body_pattern <- gsub("_1_", "_[0-9]+_", hm_body_base);
         hm_body_base_b1 <- tail(unique(gsub("_[0-9]+$", "_",
            setdiff(
               jamba::vigrep(hm_body_pattern, all_hm_components),
               use_hm_body_components))), 1)
         hm_body_base_b <- paste0(hm_body_base_b1, seq_along(hm_group_list))
      } else {
         hm_body_base_t <- use_hm_body_components;
         hm_body_base_b <- use_hm_body_components;
      }
   } else {
      hm_body_base_t <- use_hm_body_components;
      hm_body_base_b <- use_hm_body_components;
   }

   if (all(use_hm_components %in% all_hm_components)) {
      hm_title_base_lr <- use_hm_components;
      hm_title_base_tb <- use_hm_components;
   } else {
      use_hm_components <- rep(head(use_hm_components, 1),
         length.out=length(use_hm_components))
      hm_title_base_lr <- use_hm_body_components;
      hm_title_base_tb <- use_hm_components;
   }
   if (verbose) {
      jamba::printDebug("heatmap_column_group_labels(): ",
         "use_hm_components %in% all_hm_components: ",
         use_hm_components %in% all_hm_components);
      jamba::printDebug("heatmap_column_group_labels(): ",
         "use_hm_body_components %in% all_hm_components: ",
         use_hm_body_components %in% all_hm_components);
   }

   for (inum in seq_along(hm_group_list)) {
      i <- names(hm_group_list)[inum];
      if (endlines > 0) {
         i <- paste0(
            i,
            paste0(rep("\n", endlines), collapse=""))
      }
      hm_num <- hm_group_list[[inum]];
      if (verbose) {
         jamba::printDebug("heatmap_column_group_labels(): ",
            "hm_num: ", hm_num);
      }

      # left
      grid::seekViewport(hm_title_base_lr[[head(hm_num, 1)]])
      if (verbose) {
         jamba::printDebug("heatmap_column_group_labels(): ",
            "seekViewport() left-right: '",
            gsub("\n", "\\\\n", hm_title_base_lr[[inum]]), "'");
      }
      loc1_l <- grid::deviceLoc(
         x=grid::unit(0, "npc"),# + group_line_buffer,
         y=grid::unit(1, "npc") + grid::unit(y_offset_lines, "lines"))
      # right
      if (head(hm_num, 1) != tail(hm_num, 1)) {
         grid::seekViewport(hm_title_base_lr[[tail(hm_num, 1)]])
      }
      loc1_r <- grid::deviceLoc(
         x=grid::unit(1, "npc"),# - group_line_buffer,
         y=grid::unit(1, "npc") + grid::unit(y_offset_lines, "lines"))
      # top/bottom label boundaries
      if (!hm_title_base_tb[[inum]] == hm_title_base_lr[[inum]]) {
         grid::seekViewport(hm_title_base_tb[[inum]])
         if (verbose) {
            jamba::printDebug("heatmap_column_group_labels(): ",
               "seekViewport() top-bottom '",
               gsub("\n", "\\\\n", hm_title_base_tb[[inum]]), "'");
         }
         # top label
         loc2_t <- grid::deviceLoc(
            x=grid::unit(0, "npc"),
            y=grid::unit(1, "npc") + grid::unit(y_offset_lines, "lines"))
         # bottom label
         loc2_b <- grid::deviceLoc(
            x=grid::unit(1, "npc"),
            y=grid::unit(1, "npc") + grid::unit(y_offset_lines, "lines"))
      } else {
         loc2_t <- loc1_l;
         loc2_b <- loc1_r;
      }

      # grid::seekViewport(paste0(hm_title_base, head(hm_num, 1)))
      # loc1 <- grid::deviceLoc(
      #    x=grid::unit(0, "npc") + group_line_buffer,
      #    y=grid::unit(1, "npc") + grid::unit(y_offset_lines, "lines"))
      # grid::seekViewport(paste0(hm_title_base, tail(hm_num, 1)))
      # loc2 <- grid::deviceLoc(
      #    x=grid::unit(1, "npc") - group_line_buffer,
      #    y=grid::unit(1, "npc") + grid::unit(y_offset_lines, "lines"))

      grid::seekViewport("global")
      # label above each group only if label exists and has non-whitespace
      if (add_group_label &&
            !all(is.na(i)) &&
            length(i) > 0 &&
            jamba::igrepHas("[^ \t\n]", i)) {
         if (TRUE %in% use_gridtext) {
            gtg <- gridtext::richtext_grob(text=i,
               gp=grid::gpar(
                  fontsize=10 * font_cex,
                  col="black",
                  family="Arial"),
               hjust=0.5,
               vjust=-0.8,
               # halign=0.5,
               # valign=0,
               x=(loc1_l$x + loc1_r$x)*0.5,
               y=(loc2_t$y + loc2_b$y)*0.5)
            grid::grid.draw(gtg);
         } else {
            grid::grid.text(i,
               just="bottom",
               gp=grid::gpar(
                  fontsize=10 * font_cex,
                  col="red",
                  family="Arial"),
               x=(loc1_l$x + loc1_r$x)*0.5,
               y=(loc2_t$y + loc2_b$y)*0.5)
         }
      }
      # line below group label
      if (TRUE %in% add_group_line) {
         if (FALSE %in% group_line_requires_label ||
               (!any(is.na(i)) &&
                  length(i) > 0 &&
                  jamba::igrepHas("[^ \t\n]", i))) {
            if (verbose) {
               jamba::printDebug("heatmap_column_group_labels(): ",
                  "Drawing group line for i: '", i, "'");
            }
            grid::grid.lines(
               gp=grid::gpar(
                  lwd=group_line_lwd,
                  col="black"),
               default.units="native",
               x=grid::unit.c(
                  loc1_l$x + group_line_buffer,
                  loc1_r$x - group_line_buffer),
               y=grid::unit.c(loc2_t$y, loc2_b$y) + grid::unit(2, "mm"))
         } else {
            if (verbose) {
               jamba::printDebug("heatmap_column_group_labels(): ",
                  "Hiding group line for i: '", i, "'");
            }
         }
      }

      # box around heatmap
      if (TRUE %in% add_group_box) {
         # top box
         grid::seekViewport(hm_body_base_t[[inum]])
         if (verbose) {
            jamba::printDebug("heatmap_column_group_labels(): ",
               "seekViewport() box top '",
               gsub("\n", "\\\\n", hm_body_base_t[[inum]]), "'");
         }
         loc3_t <- grid::deviceLoc(
            x=grid::unit(0, "npc"),
            y=grid::unit(1, "npc"))
         if (!hm_body_base_t[[inum]] == hm_body_base_b[[inum]]) {
            grid::seekViewport(hm_body_base_b[[inum]])
            if (verbose) {
               jamba::printDebug("heatmap_column_group_labels(): ",
                  "seekViewport() box bottom '",
                  gsub("\n", "\\\\n", hm_body_base_b[[inum]]), "'");
            }
         }
         loc3_b <- grid::deviceLoc(
            x=grid::unit(0, "npc"),
            y=grid::unit(0, "npc"))

         # grid::seekViewport(paste0(hm_body_base, head(hm_num, 1)))
         # loc1 <- grid::deviceLoc(
         #    x=grid::unit(0, "npc"),
         #    y=grid::unit(0, "npc"))
         # grid::seekViewport(paste0(hm_body_base, tail(hm_num, 1)))
         # loc2 <- grid::deviceLoc(
         #    x=grid::unit(1, "npc"),
         #    y=grid::unit(1, "npc"))

         grid::seekViewport("global")
         grid::grid.rect(
            gp=grid::gpar(fill=NA,
               col="black",
               linejoin="mitre",
               lwd=group_box_lwd),
            width=loc1_r$x - loc1_l$x +
               grid::unit(group_box_lwd/2 * group_box_outer, "points"),
            height=loc3_t$y - loc3_b$y +
               grid::unit(group_box_lwd/2 * group_box_outer, "points"),
            x=(loc1_l$x + loc1_r$x)*0.5,
            y=(loc3_t$y + loc3_b$y)*0.5)
         # width=loc1$x - loc1$x +
         #    grid::unit(group_box_lwd/2 * group_box_outer, "points"),
         # height=loc2$y - loc1$y +
         #    grid::unit(group_box_lwd/2 * group_box_outer, "points"),
         # x=(loc1$x + loc2$x)*0.5,
         # y=(loc1$y + loc2$y)*0.5)
      }
   }

}
