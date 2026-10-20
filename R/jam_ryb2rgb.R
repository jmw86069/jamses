
#' Convert red-yellow-blue color matrix to RGB, internal
#' 
#' Convert red-yellow-blue color matrix to RGB, and interpret
#' the RGB components of a hex color as RRYYBB. Then convert
#' the colors defined using a red-yellow-blue color definition into
#' a standard red-green-blue (RGB) color.
#' 
#' This function is currently only internal, intended to
#' help combine three-factor experiments into RGB style
#' colors, except using red-yellow-blue to improve
#' the range of colors.
#' 
#' @returns `numeric` red,green,blue matrix when input is
#'    `numeric`, or `character` hex color when input is `character`.
#' 
#' @param use_ryb `numeric` matrix with columns representing
#'    red, yellow, blue. It assumes max value 1 in each column,
#'    and does not validate `colnames()`. It assumes RYB order.
#'    * Alternatively, `character` hex colors, which are assumed
#'    to represent RRYYBB (red, yellow, blue) color definitions.
#' @param ... additional arguments are ignored.
#' 
#' @examples
#' group_names <- paste0(
#'    rep(c("DMSO", "Etop"), each=2), "_",
#'    rep(c("NF", "Flag"), each=4), "_",
#'    c(rep(c("WT", "KO"), 2), rep(c("WT", "D955N"), 2)))
#' group_fill <- (c(
#'    "gray50",
#'    "slateblue",
#'    "salmon3",
#'    "orchid4",
#'    "lightgoldenrod3",
#'    "aquamarine3",
#'    "orange2",
#'    "lightgoldenrod1"))
#' names(group_fill) <- group_names;
#' 
#' use_ryb <- matrix(ncol=3, data=c(
#'    0,0,1,1,0,0,1,1,
#'    0,0,0,0,1,1,1,1,
#'    0,1,0,1,0,1,0,1
#' ) /4 + 0.5)
#' rownames(use_ryb) <- names(group_fill)
#' use_rgb <- ryb2rgb(use_ryb)
#' use_rgb <- ryb2rgb(use_ryb, do_plot=TRUE)
#' use_fill_rgb <- jamba::rgb2col(use_ryb)
#' use_fill_ryb <- jamba::rgb2col(use_rgb)
#' jamba::showColors(list(
#'    use_fill_rgb=use_fill_rgb,
#'    use_fill_ryb=use_fill_ryb,
#'    group_fill=group_fill
#' ))
#' 
#' ryb2rgb(c("gold", "green", "blue"), do_plot=TRUE)
#' 
#' @keywords internal
#' @noRd
ryb2rgb <- function(
   use_ryb,
   use_model=c("hsl", "hcl"),
   do_plot=FALSE,
   ...
) {
   #
   use_model <- match.arg(use_model)
   #
   if (is.numeric(use_ryb) && ncol(use_ryb) >= 3) {
      use_fill <- jamba::rgb2col(use_ryb)
   } else if (is.character(use_ryb)) {
      use_fill <- use_ryb;
   }
   use_hcl <- jamba::col2hcl(
      use_fill
   )

   # print("use_hcl:");print(use_hcl);# debug
   use_hsl <- jamba::col2hsl(
      use_fill
   )
   use_h_ryb <- colorjam::hw2h(
      h=use_hcl["H", ],
      preset="ryb"
   )
   use_hcl_ryb <- use_hcl;
   use_hcl_ryb["H", ] <- use_h_ryb;
   use_hsl_ryb <- use_hsl;
   hsl_lwt <- 0.8;
   if ("hcl" %in% use_model) {
      # print("use_hcl_ryb:");print(use_hcl_ryb);# debug
      use_fill_ryb <- jamba::hcl2col(use_hcl_ryb)
   } else {
      use_hsl_ryb["H", ] <- (use_hcl_ryb["H", ] +
         colorjam::hcl_to_hsl_hue(use_hcl_ryb["H", ])) / 2;
      use_hsl_ryb["H", ] <- colorjam::hcl_to_hsl_hue(use_hcl_ryb["H", ]);
      new_s <- ifelse(
         use_hsl["S", ] > 10,
         (pmax(use_hsl["S", ], use_hcl["C", ]) + 25) / 1.25,
         pmax(use_hsl["S", ], use_hcl["C", ])
      )
      use_hsl_ryb["S", ] <- new_s;
      use_hsl_ryb["L", ] <- (use_hsl["L", ]*hsl_lwt + use_hcl["L", ]*(1 - hsl_lwt));
      # print("use_hsl_ryb:");print(use_hsl_ryb);# debug
      use_fill_ryb <- jamba::hsl2col(use_hsl_ryb)
   }

   # optional plot
   if (isTRUE(do_plot)) {
         jamba::showColors(list(
            input=use_fill,
            output=use_fill_ryb
         ))
   }
   if (is.character(use_ryb)) {
      return(use_fill_ryb)
   }
   use_rgb_ryb <- t(col2rgb(use_fill_ryb))
   return(use_rgb_ryb)
}
