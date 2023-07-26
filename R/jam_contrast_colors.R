
#' Define contrast colors by group colors
#'
#' Define contrast colors by group colors
#'
#' This function takes colors assigned to groups, and applies them
#' to contrasts based upon the groups represented in each contrast.
#' Group colors are blended into one output color per contrast,
#' then the color is adjusted for minimum chroma (C) and maximum
#' luminance (L) to ensure the colors have acceptable vibrance.
#'
#' When `sample_color_list` is supplied, which defines colors per group,
#' these colors are used as appropriate.
#'
#' When `sample_color_list` is not supplied, categorical colors are
#' defined for all observed group names using `colorjam::group2colors()`.
#'
#' @returns `character` vector of R colors, with names defined
#'    by `contrast_names(sedesign)`.
#'
#' @family jam experiment design
#'
#' @param sedesign `SEDesign` object
#' @param sample_color_list `list` of colors, with one list element with
#'    name defined by `group_name`, default `"Group"`, so the colors
#'    use `sample_color_list[[group_name]]`. This element should be
#'    a `character` vector of colors, whose names are group names
#'    used in the `contrast_names(sedesign)`.
#'    The `sample_color_list` is produced by `platjam::design2colors()`.
#' @param group_name `character` string for the element in `sample_color_list`
#'    that contains group colors.
#' @param contrast_sep `character` string used as delimiter between
#'    group names in each contrast name, default is `"-"`.
#' @param C_min `numeric` value to impose a minimum chroma (C) value for
#'    the resulting contrast color.
#' @param L_max `numeric` value to impose a maximum luminance (C) value for
#'    the resulting contrast color.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to `colorjam::group2colors()`
#'    for any groups for which no colors are assigned.
#'
#' @export
contrast_colors_by_group <- function
(sedesign,
 sample_color_list=NULL,
 group_name="Group",
 contrast_sep="-",
 C_min=80,
 L_max=70,
 verbose=FALSE,
 ...)
{
   #
   use_contrast_names <- contrast_names(sedesign);
   group_names <- strsplit(
      gsub("[()]", "", use_contrast_names),
      contrast_sep);
   names(group_names) <- use_contrast_names;
   all_group_names <- unique(unlist(group_names))

   group_colors <- sample_color_list[[group_name]];
   if (length(sample_color_list) == 0 ||
         !all(all_group_names %in% names(group_colors))) {
      missing_group_colors <- setdiff(all_group_names,
         names(group_colors));
      if (verbose) {
         jamba::printDebug("contrast_colors_by_group(): ",
            c("Using ",
            "colorjam::group2colors()",
            " to assign colors to ",
            jamba::formatInt(length(missing_group_colors)),
            " missing group colors."),
            sep="")
      }
      new_group_colors <- colorjam::group2colors(missing_group_colors,
         ...);
      group_colors[names(new_group_colors)] <- new_group_colors;
   }

   contrast_colors <- lapply(group_names, function(igroups) {
      icolors <- group_colors[igroups];
      blended <- c(blended=colorjam::blend_colors(icolors));
      new_blended <- farver::raise_channel(colour=blended,
         value=C_min,
         channel="c",
         space="hcl")
      new_blended2 <- farver::cap_channel(colour=new_blended,
         value=L_max,
         channel="l",
         space="hcl")
      unname(new_blended2)
   })
   # remove any empty entries and convert to vector
   use_contrast_colors <- unlist(jamba::rmNULL(contrast_colors))
   use_contrast_colors
}
