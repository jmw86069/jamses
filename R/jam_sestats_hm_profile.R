
#' Profile Line Plot for Heatmap of SummarizedExperiment data
#'
#' Convert Heatmap of SummarizedExperiment data into Profile Line Plot,
#' with row and column groups and facets as relevant.
#'
#' It uses the data within the heatmap to produce the plot,
#' therefore all data centering is used as defined in `heatmap_se()`.
#' This plot is intended to produce a line plot equivalent to
#' the same data shown in heatmap form. Sometimes the linear plot
#' is more visually intuitive, and certainly is more effective
#' at conveying magnitudes than color gradient.
#'
#' Row groups are detected by calling `jamba::heatmap_row_order()`.
#' If there are no row groups, all rows are placed into the same group.
#'
#' Column groups are also detected by calling `jamba::heatmap_column_order()`.
#' The argument `summarize_column_split=TRUE` will optionally calculate
#' group mean per column group, but the default is `FALSE`.
#'
#' @param hm `ComplexHeatmap::Heatmap` object, intended to use output
#'    from `heatmap_se()`.
#'    It uses the data within the heatmap to produce the plot,
#'    therefore all data centering is used as defined in `heatmap_se()`.
#' @param display_profiles `character` indicating which profiles to include:
#'    * 'rows': include a profile line for each row
#'    * 'group': include a profile line for the row group.
#' @param summarize_column_split `logical` default TRUE, whether to summarize
#'    values in `column_split` by taking the mean across replicates (columns)
#'    within each column grouping.
#'    * When `summarize_column_split=FALSE` it will keep replicates as
#'    individual points, which may appear more consistent with the heatmap.
#' @param rowData_colnames `character` vector, currently ignored.
#'    In future it may be used to apply different row groupings.
#' @param profile_color `character` string, default is partially transparent,
#'    light grey. This color is used for the line profile color.
#' @param include_count `logical` default TRUE, whether to include the
#'    row count in the label for each row group.
#'    * It uses `sep="\n"` by default between the row group label,
#'    and the row count label. So to make the label all appear on one line,
#'    use `sep=" "`.
#' @param row_type `character` string, default `"Gene"` used when
#'    `include_count=TRUE`, so the label becomes `"11 Genes"`
#'    The label should be singular, and the 's' is added when there
#'    is more than 1 row in the group.
#' @param strip.position `character` string, default 'right', places the
#'    ggplot2 facet strip label on the right side.
#'    The other common option is 'top', which is better for long labels,
#'    but requires a taller graphics device to accomodate the number of
#'    row groups.
#' @param base_size `numeric` base default font size for the ggplot2 theme,
#'    default 10 is relatively small, suitable for analysis.
#'    Use a larger number for figures, but adjust the graphics device size
#'    accordingly.
#' @param panel.background.fill `character` color used for the ggplot2
#'    panel background color, default is very light grey (lighter than
#'    ggplot2 default). The grey color makes it slightly easier to recognize
#'    each plot panel.
#' @param ... additional arguments are passed to internal functions,
#'    specifically `colorjam::theme_jam()` which defines the ggplot2 theme.
#'    The theme can be replaced by adding `ggplot2::theme_bw()` or something
#'    else for example.
#'
#' @returns `ggplot` object suitable for plotting or further customization.
#'    * When `return_type='list'` it returns a `list` with 'hmtall'
#'    which is a `data.frame` containing data used for the ggplot.
#'    If display_profiles includes 'group' the `list` will contain 'hmtall4'
#'    which is a `data.frame` containing the row group mean values.
#'    * When return_type='ggplot' (default) the `output@data`
#'    can also be used to inspect the underlying data.
#'
#' @family jamses heatmaps
#'
#' @examples
#' se <- make_se_test()
#' hm1 <- heatmap_se(se,
#'    rowData_colnames="Class",
#'    row_split=9,
#'    column_split="group",
#'    show_left_annotation_name="top",
#'    left_annotation_name_rot=150,
#'    sample_color_list=list(group=c(groupA="red")))
#' hm1drawn <- ComplexHeatmap::draw(hm1)
#' heatmap_profile_plot(hm1drawn, summarize_column_split=FALSE);
#'
#' se <- jamses::make_se_test(nrow=200, ngroups=4, nreps=8)
#' hm <- jamses::heatmap_se(se,
#'    rowData_colnames="Class",
#'    apply_hm_column_title=TRUE,
#'    column_split="group",
#'    cluster_row_slices=TRUE,
#'    controlSamples=head(colnames(se), 8),
#'    row_split=12)
#' ComplexHeatmap::draw(hm)
#'
#' heatmap_profile_plot(hm, strip.position="top")
#'
#' heatmap_profile_plot(hm, summarize_column_split=FALSE, strip.position="top")
#'
#' heatmap_profile_plot(hm, summarize_column_split=FALSE)
#'
#' @export
heatmap_profile_plot <- function
(hm,
 display_profiles=c("rows", "group"),
 summarize_column_split=TRUE,
 se=NULL,
 rowData_colnames=NULL,
 profile_color="#88888844",
 mean_color="#000000FF",
 include_count=TRUE,
 row_type="Gene",
 strip.position="right",
 base_size=10,
 panel.background.fill="grey97",
 facet_ncol=1,
 return_type=c("ggplot",
    "list"),
 ...)
{
   # validate arguments
   display_profiles <- match.arg(display_profiles,
      several.ok=TRUE);
   return_type <- match.arg(return_type);
   # verify ggplot2 is available
   if ("ggplot" %in% return_type &&
         !requireNamespace("ggplot2", quietly=TRUE)) {
      stop(paste0("heatmap_profile_plot() requires ggplot2 when ",
         "return_type='ggplot'."));
   }

   # custom pivot function
   pivotmatrix <- function
   (x,
      rows_to="row",
      names_to="name",
      values_to="value",
      ...)
   {
      idf <- data.frame(check.names=FALSE,
         tidyr::pivot_longer(
            data.frame(row=rownames(x), x),
            names_to=names_to,
            values_to=values_to,
            cols=colnames(x)))
      colnames(idf)[1] <- rows_to;
      idf;
   }

   # custom function
   # convert heatmap row/column order to named factor vector
   column_order_to_vector <- function
   (hro,
      add_sizes=FALSE,
      type="item",
      sep="\n",
      ...)
   {
      if (all(lengths(hro) == 1)) {
         names(hro) <- NULL;
         hro <- list(Group=unlist(hro))
      }
      if (length(names(hro)) == 0) {
         names(hro) <- as.character(seq_along(hro));
      }
      hrov <- jamba::nameVector(
         factor(rep(names(hro), lengths(hro)),
            levels=names(hro)),
         unlist(hro));
      if (isTRUE(add_sizes)) {
         hrovt <- table(hrov);
         hrovtl <- jamba::nameVector(
            paste0(names(hrovt),
               sep, "(",
               jamba::formatInt(hrovt)," ",
               type,
               ifelse(hrovt > 1, "s", ""), ")"),
            names(hrovt))
         levels(hrov) <- hrovtl[levels(hrov)];
      }
      return(hrov)
   }

   # accept HeatmapList
   if (inherits(hm, "HeatmapList")) {
      hml <- hm;
      hm <- hml@ht_list[[1]];
   }
   # row order and groupings
   hro <- jamba::heatmap_row_order(hm)
   hrov <- column_order_to_vector(hro,
      add_sizes=include_count,
      type=row_type,
      ...)

   # column order and groupings
   hco <- jamba::heatmap_column_order(hm)
   hcov <- column_order_to_vector(hco)

   # make tall data.frame
   hmtall <- pivotmatrix(hm@matrix,
      "Gene",
      "Sample",
      "Centered")

   # add row groups
   hmtall$row_group <- hrov[hmtall$Gene]
   hmtall$Sample <- factor(hmtall$Sample,
      levels=names(hcov))
   hmtall$Gene <- factor(hmtall$Gene,
      levels=names(hrov))

   if (length(hco) > 1 && length(hco) < ncol(hm@matrix)) {
      # column groups
      hmtall$Group <- hcov[hmtall$Sample]
   }
   # Optionally use column groups
   if (isTRUE(summarize_column_split) &&
         length(hco) < ncol(hm@matrix)) {

      library(dplyr)
      hmtall2 <- hmtall %>%
         group_by(Group, Gene) %>%
         mutate(GroupMean=mean(Centered, na.rm=TRUE)) %>%
         ungroup()
      hmtall <- subset(hmtall2, !duplicated(paste0(Group, Gene)))
      hmtall$Sample <- hmtall$Group
   }

   ret_vals <- list();

   # row group mean
   if ("group" %in% display_profiles) {
      hmtall3 <- hmtall |>
         dplyr::group_by(row_group, Sample) |>
         dplyr::mutate(ProfileMean=mean(Centered, na.rm=TRUE)) |>
         dplyr::ungroup()
      hmtall4 <- subset(hmtall3, !duplicated(paste0(row_group, Sample)))
      hmtall4$Centered <- hmtall4$ProfileMean;
      hmtall4$ProfileType <- "Mean";
      hmtall$ProfileType <- row_type;
      hmtall4 <- hmtall4[, colnames(hmtall), drop=FALSE];
      hmtall4$Gene <- paste0(hmtall4$Gene, "_mean");
      hmtall <- rbind(hmtall, hmtall4)
      # hmtall4$Sample <- hmtall$Group
   }

   if (!"rows" %in% display_profiles) {
      hmtall <- subset(hmtall, ProfileType %in% "Mean")
   }

   # return list of data objects
   if ("list" %in% return_type) {
      ret_vals$hmtall <- hmtall;
      if ("group" %in% display_profiles) {
         ret_vals$hmtall4 <- hmtall4;
      }
      return(ret_vals);
   }

   # return ggplot2
   gg <- ggplot2::ggplot(hmtall,
      ggplot2::aes(
         x=Sample,
         y=Centered,
         group=Gene)) +
      ggplot2::geom_hline(yintercept=0,
         color="firebrick4") +
      colorjam::theme_jam(
         base_size=base_size,
         panel.background=ggplot2::element_rect(
            fill=panel.background.fill,
            colour=NA),
         ...)
   if (length(display_profiles) == 1) {
      if ("rows" %in% display_profiles) {
         gg <- gg +
            ggplot2::geom_line(color=profile_color)
      } else {
         gg <- gg +
            ggplot2::geom_line(color=mean_color)
      }
   } else {
      gg <- gg +
         ggplot2::geom_line(ggplot2::aes(color=ProfileType)) +
         ggplot2::scale_color_manual(values=jamba::nameVector(
            c(profile_color, mean_color),
            c(row_type, "Mean"))
         )
   }
   # if (length(hco) > 1 && isFALSE(summarize_column_split)) {
   #    gg <- gg +
   #       ggplot2::facet_grid(row_group ~ Group, scales="free_x")
   # } else {
   gg <- gg +
      ggplot2::facet_wrap(~row_group,
         strip.position=strip.position,
         ncol=facet_ncol)
   # }
   gg
}
