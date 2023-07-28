
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
#'    Alternatively, a `character` vector of colors, named by contrast.
#' @param group_name `character` string for the element in `sample_color_list`
#'    that contains group colors.
#' @param contrast_sep `character` string used as delimiter between
#'    group names in each contrast name, default is `"-"`.
#' @param C_min,C_max `numeric` value to impose a minimum and maximum
#'    chroma (C) value for the resulting contrast color.
#' @param L_min,L_max `numeric` value to impose a minimum and maximum
#'    luminance (C) value for the resulting contrast color.
#' @param default_color_style `character` string indicating how to assign
#'    colors to groups without a color assignment.
#'    * `"categorical"`: calls `colorjam::group2colors()` which calls
#'    `colorjam::rainbowJam()` for rainbow categorical colors.
#'    * any single R color: this one R color is applied to all other groups.
#'    * any non-R color: assumed to be the name of a color gradient or set,
#'    and is resolved by calling `jamba::getColorRamp()`. It recognizes
#'    RColorBrewer and viridis color sets for example.
#'    * vector of multiple colors: also resolved by calling
#'    `jamba::getColorRamp()` which creates a gradient across the colors.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are passed to `colorjam::group2colors()`,
#'    or `jamba::getColorRamp()` as needed.
#'
#' @export
contrast_colors_by_group <- function
(sedesign,
 sample_color_list=NULL,
 group_name="Group",
 contrast_sep="-",
 C_min=80,
 C_max=100,
 L_min=35,
 L_max=70,
 default_color_style=c("categorical",
    "navy"),
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

   # handle empty sample_color_list
   if (length(sample_color_list) == 0) {
      sample_color_list <- "navy";
   }
   if (is.atomic(sample_color_list)) {
      if (length(names(sample_color_list)) == 0) {
         group_colors <- rep(sample_color_list,
            length.out=length(all_group_names));
         names(group_colors) <- all_group_names;
      } else {
         group_colors <- sample_color_list;
      }
   } else {
      group_colors <- sample_color_list[[group_name]];
   }
   # for any
   if (!all(all_group_names %in% names(group_colors))) {
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
      if ("categorical" %in% default_color_style) {
         new_group_colors <- colorjam::group2colors(missing_group_colors,
            ...);
      } else {
         if (length(default_color_style) == 1 &&
               jamba::isColor(default_color_style)) {
            # assign single color to all contrasts
               new_group_colors <- jamba::nameVector(
                  rep(default_color_style,
                     length.out=length(missing_group_colors)),
                  missing_group_colors);
         } else {
            # if multiple colors, or non-color, use getColorRamp
            new_group_colors <- jamba::nameVector(
               jamba::getColorRamp(default_color_style,
                  gradientN=length(missing_group_colors),
                  n=length(missing_group_colors)),
               missing_group_colors);
         }
      }
      group_colors[names(new_group_colors)] <- new_group_colors;
   }

   contrast_colors <- lapply(group_names, function(igroups) {
      icolors <- group_colors[igroups];
      blended <- c(blended=colorjam::blend_colors(icolors));
      new_blended <- farver::raise_channel(colour=blended,
         value=C_min,
         channel="c",
         space="hcl")
      new_blended <- farver::raise_channel(colour=blended,
         value=L_min,
         channel="l",
         space="hcl")
      new_blended2 <- farver::cap_channel(colour=new_blended,
         value=L_max,
         channel="l",
         space="hcl")
      new_blended2 <- farver::cap_channel(colour=new_blended,
         value=C_max,
         channel="c",
         space="hcl")
      unname(new_blended2)
   })
   # remove any empty entries and convert to vector
   use_contrast_colors <- unlist(jamba::rmNULL(contrast_colors))
   use_contrast_colors
}
