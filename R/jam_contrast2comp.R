
#' Convert contrast to short-form comp, convert comp to contrast
#'
#' Convert contrast to short-form comp, and convert comp back to original contrast.
#'
#' These functions are intended to reduce the number of characters
#' required to represent a statistical contrast.
#' `contrast2comp()` converts long to short form, and
#' `comp2contrast()` converts short to long form.
#'
#' 1. `"contrast"`: the fully-defined contrast
#' 2. `"comp"`: equivalent abbreviated form, a short comparison
#'
#' Note that one goal is to reduce characters in Excel worksheet names,
#' currently limited to 31 characters. Also note, the `":"` delimiter is
#' not permitted in Excel sheet names, thus `save_sestats()` uses
#' semicolon `";"`. This limitation may warrant using a different default
#' delimiter between factors, such as comma `","`, or pipe `"|"`,
#' or forward-slash `"/"`.
#'
#' ## Assumptions
#'
#' The key assumption is that an experimental group name is a
#' `character` string composed of its factor levels, with a delimiter
#' between factors. For example:
#'
#' * `CellA_Treated` - is interpreted as `"CellA"` and `"Treated"`
#' * `CellA_Control` - is interpreted as `"CellA"` and `"Control"`
#'
#' A contrast therefore:
#'
#' * `CellA_Treated-CellA_Control` can be re-written
#' * `CellA:Treated-Control`
#'
#' Factors must be in identical order for all groups, and there
#' must be no empty factor levels.
#' Do not use: `"CellA_Treated_Time0"`, `"CellA_Time0"`.
#'
#' Finally, the overall assumption is that contrasts are composed
#' of reasonable comparisons between factor levels, with no
#' more factors being compared than the depth of contrast. For example,
#' a one-way contrast can compare one factor,
#' a two-way contrast can compared two factors, and so on.
#' In most cases where the assumptions above are broken, the
#' output should be the same as the input, with no change.
#'
#' When using `groups_to_sedesign()`, the output contrasts should all
#' meet these requirements, therefore all contrasts can be
#' converted to `"comp"` form for plot labels, and converted
#' back to `"contrasts"` as needed.
#'
#' Delimiters can be customized, however they must all be
#' single-character values, avoiding `()[]` which are reserved.
#' For example, sometimes factors are separated by `"."` such as
#' in the contrast: `"A.B-C.B"`. In this case, use:
#' `contrast2comp("A.B-C.B", contrast_factor_delim=".")`.
#' The corresponding conversion back to contrast would be:
#' `comp2contrast("A-C:B", contrast_factor_delim=".")`
#'
#' ## Design goals for conversion to short form comp
#'
#' 1. `"comp"` should be interchangeable with `"contrast"`
#'
#'    * use `contrast2comp()` and `comp2contrast()`
#'
#' 2. when a contrast cannot be abbreviated, comp will use contrast
#'
#'    * see examples
#'    * when more factors are being compared than the contrast order,
#'    the function will leave the contrast as-is
#'    * Consider `"CellA_Treated-CellB_Control"`.
#'    Both `"CellA-CellB"` and `"Treated-Control"` are compared
#'    in a one-way contrast, therefore it cannot be abbreviated.
#'
#' 3. `"comp"` shall not create any whitespace
#'
#'    * factors will be delimited with `":"`
#'    * factor levels will be delimited with `"-"`
#'    * other potential delimiters `"*"`, `"+"`
#'    already have meaning in formula context.
#'
#' 4. `"comp"` shall not use parentheses `"()"`, where possible
#'
#'    * the goal is to reduce characters
#'    * parentheses are not necessary for balanced contrasts
#'    * unbalanced contrasts (see point 2) will retain the original syntax
#'
#' 5. the order of factors should be maintained in `"comp"`
#'
#'    * goal is to reproduce the original correct group name in contrast form
#'    * the original group name is necessary for the design matrix
#'
#' ## Worked examples
#'
#' 1. One-way contrast
#'
#'    * contrast: `CellA_Treated-CellA_Control`
#'    * comment: `CellA` is unchanged, `Treated-Control` is changed
#'    * comp: `CellA:Treated-Control`
#'
#' 2. Two-way contrast
#'
#'    * contrast: `(CellA_Treated-CellA_Control)-(CellB_Treated-CellB_Control)`
#'    * comment: `CellA-CellB` is changed, `Treated-Control` is changed
#'    * comp: `CellA-CellB:Treated-Control`
#'    * note: when converting comp `CellA-CellB:Treated-Control` back to
#'    contrast, two forms are mathematically equivalent:
#'
#'    ```
#'    # form 1
#'    (CellA_Treated-CellA_Control)-(CellB_Treated-CellB_Control)
#'    # form 2
#'    (CellA_Treated-CellB_Treated)-(CellA_Control-CellB_Control)
#'    # both are equivalent
#'    CellA_Treated - CellB_Treated - CellA_Control + CellB_Control
#'    ```
#'    * These two forms can be controlled in `comp2contrast()` with
#'    argument `factor_order`.
#'
#' 3. Three-way contrast (it happens rarely, but does happen)
#'
#'    * contrast:
#'    ```
#'    (CellA_Treated_Mut-CellA_Control_Mut)-(CellB_Treated_Mut-CellB_Control_Mut) -
#'    (CellA_Treated_WT-CellA_Control_WT)-(CellB_Treated_WT-CellB_Control_WT)
#'    ```
#'    * comment: `CellA-CellB`, `Treated-Control`, `Mut-WT` are changed
#'    * comp: `CellA-CellB:Treated-Control:Mut-WT`
#'
#' 4. One-way contrast with additional unchanged factors
#'
#'    * contrast: `CellA_Treated_WT-CellA_Control_WT`
#'    * comment: `CellA`, `WT` are unchanged, `Treated-Control` is changed
#'    * comp: `CellA:Treated-Control:WT`
#'
#' 5. Unbalanced one-way contrast
#'
#'    * contrast: `CellA_Treated-CellB_Control`
#'    * comment: `CellA-CellB` and `Treated-Control` are changed
#'    * comp: `CellA_Treated-CellB_Control`
#'
#' 6. Mis-directed two-way contrast
#'
#'    * contrast: `(CellA_Treated-CellA_Control)-(CellB_Control-CellB_Treated)`
#'    * comment: `CellA-CellB` are changed, `Treated-Control/Control-Treated` are changed
#'    * Note: The `Treated-Control` and `Control-Treated` do not agree in
#'    direction. The output is partially abbreviated, and maintains the
#'    original direction to prevent loss of information.
#'    * comp: `(CellA:Treated-Control)-(CellB:Control-Treatment)`
#'
#' @examples
#' contrast_names <- c(
#'    "CellA_Treated-CellA_Control",
#'    "CellB_Treated-CellB_Control",
#'    "CellB_Treated-CellA_Control",
#'    "(CellA_Treated-CellA_Control)-(CellB_Treated-CellB_Control)",
#'    "(CellB_Treated-CellB_Control)-(CellA_Treated-CellA_Control)",
#'    "(CellA_Treated-CellB_Treated)-(CellA_Control-CellB_Control)"
#' );
#' contrast2comp(contrast_names)
#'
#' contrast2comp(contrast_names, comp_factor_delim=";")
#'
#' comps <- contrast2comp(contrast_names)
#' data.frame(contrast_names,
#'    nchar_contrasts=nchar(contrast_names),
#'    comps,
#'    nchar_comps=nchar(comps))
#'
#' # compare conversion back to contrast
#' data.frame(contrast_names,
#'    comps=comps,
#'    contrast_again=comp2contrast(comps),
#'    changed=contrast_names != comp2contrast(comps))
#'
#' # factors can be ordered by contrast
#' contrasts2 <- comp2contrast(comps,
#'    factor_order=list(1:2, 1:2, 1:2,
#'       2:1, 2:1, 1:2))
#' # compare conversion back to contrast
#' data.frame(contrast_names,
#'    comps=comps,
#'    contrasts2,
#'    changed=contrast_names != contrasts2)
#'
#' # note change in direction for two-way contrasts
#' # Treated-Control and Control-Treated
#' contrast_diff <- "(CellA_Treated-CellA_Control)-(CellB_Control-CellB_Treated)";
#' comp_diff <- contrast2comp(contrast_diff)
#' # partially abbreviated comp
#' comp_diff
#' # it is converted back to original form
#' comp2contrast(comp_diff)
#'
#' data.frame(contrast_diff,
#'    nchar_contrasts=nchar(contrast_diff),
#'    comp_diff,
#'    nchar_comps=nchar(comp_diff))
#'
#' # evaluate the rare three-way contrast
#' contrast_names_3way <- c(
#'    contrast_names[4],
#'    gsub("([a-zA-Z])([-)])", "\\1_Mut\\2", contrast_names[4]),
#'    gsub("([a-zA-Z])([-)])", "\\1_WT\\2", contrast_names[4]),
#'    paste0("(",
#'    gsub("([a-zA-Z])([-)])", "\\1_Mut\\2", contrast_names[4]),
#'    ")-(",
#'    gsub("([a-zA-Z])([-)])", "\\1_WT\\2", contrast_names[4]),
#'    ")"))
#' contrast_names_3way <- c(
#'    paste0("(CellA_Treated-CellA_Control)-",
#'       "(CellB_Treated-CellB_Control)"),
#'    paste0("(CellA_Treated_Mut-CellA_Control_Mut)-",
#'       "(CellB_Treated_Mut-CellB_Control_Mut)"),
#'    paste0("(CellA_Treated_WT-CellA_Control_WT)-",
#'       "(CellB_Treated_WT-CellB_Control_WT)"),
#'    paste0("((CellA_Treated_Mut-CellB_Treated_Mut)-",
#'       "(CellA_Control_Mut-CellB_Control_Mut))-",
#'       "((CellA_Treated_WT-CellB_Treated_WT)-",
#'       "(CellA_Control_WT-CellB_Control_WT))"),
#'    paste0("((CellA_Treated_Mut-CellA_Control_Mut)-",
#'       "(CellB_Treated_Mut-CellB_Control_Mut))-",
#'       "((CellA_Treated_WT-CellA_Control_WT)-",
#'       "(CellB_Treated_WT-CellB_Control_WT))"))
#' comp_3way <- contrast2comp(contrast_names_3way);
#' data.frame(contrast_names_3way,
#'    nchar_contrasts=nchar(contrast_names_3way),
#'    comp_3way,
#'    nchar_comps=nchar(comp_3way));
#'
#' # compare to input
#' contrasts2_3way <- comp2contrast(comp_3way);
#' # mathematically correct contrasts but in different order from input
#' data.frame(contrast_names_3way,
#'    contrasts2_3way,
#'    changed=contrast_names_3way != contrasts2_3way);
#'
#' # custom factor order produces the same contrasts as input
#' contrasts2_3way_v2 <- comp2contrast(comp_3way,
#'    factor_order=list(c(2,1,3), c(2,1,3), c(2,1,3),
#'       c(1,2,3), c(2,1,3)));
#' data.frame(contrast_names_3way,
#'    contrasts2_3way_v2,
#'    changed=contrast_names_3way != contrasts2_3way_v2);
#'
#' @param contrast_names `character` vector of statistical contrasts
#' @param contrast_delim `character` string delimiter between groups,
#'    typically `"-"` to indicate subtraction of group means.
#' @param contrast_factor_delim `character` string delimiter between
#'    design factors in a contrast.
#' @param comp_factor_delim `character` string delimiter between
#'    design factors in a comp.
#' @param factor_order `integer`, `list` of `integer` vectors, or `NULL`.
#'    When supplied as `integer` vector, it is converted to a `list`
#'    and expanded to `length()` of the input. The `integer` values
#'    are used by `comp2contrast()` to force the order of factor
#'    comparisons for two-way and higher order contrasts.
#' @param add_attr `logical` indicating whether to add attributes to
#'    the output, containing the input values provided.
#' @param verbose `logical` indicating whether to print verbose output,
#'    or for much more verbose output use `verbose=2`.
#' @param ... additional arguments are ignored.
#'
#' @export
contrast2comp <- function
(contrast_names,
 contrast_delim="-",
 contrast_factor_delim="_",
 comp_factor_delim=":",
 add_attr=FALSE,
 verbose=FALSE,
 ...)
{
   # validate delim all have 1 character
   check_delim <- function(x){
      (!is.list(x) && length(x) == 1 && !is.na(x) && nchar(x) == 1 &&
            !x %in% c("(", ")", "[", "]"))
   }
   if (!all(check_delim(contrast_delim) &&
         check_delim(contrast_factor_delim) &&
         check_delim(comp_factor_delim))) {
      stop("All delim arguments must be single character strings.")
   }
   # check for two-way contrasts
   # contrast <- gsub("[(]([^()]+)-([^()]+)[)]",
   #    "<\\1=\\2>",
   #    contrast_names)
   # contrast
   # contrast2 <- gsub("[(]([^()]+)-([^()]+)[)]",
   #    "<\\1=\\2>",
   #    contrast)
   # contrast2

   split_pattern <- paste0("[", contrast_delim, "()", "]+");
   contrast_delim_pattern <- paste0("[", contrast_delim, "]");
   contrast_factor_delim_pattern <- paste0("[", contrast_factor_delim, "]")

   contrast_names_split <- strsplit(contrast_names,
      paste0("[", contrast_delim, "()]+"));
   if (verbose > 1) {jamba::printDebug("contrast_names_split:");print(contrast_names_split);}
   comps <- unname(sapply(contrast_names_split, function(i){
      if (length(i) <= 1) {
         return(i)
      }
      im <- jamba::rbindList(strsplit(setdiff(i, ""), contrast_factor_delim_pattern));
      ksplit <- rep(head(LETTERS, nrow(im)/2), each=2);
      # code below keeps contrast order intact
      if (verbose) {jamba::printDebug("im:");print(im);}
      imdf_split <- split(
         data.frame(
            stringsAsFactors=FALSE,
            check.names=FALSE,
            im),
         ksplit)
      imnew <- jamba::rbindList(
         lapply(imdf_split, function(idf){
            do.call(cbind, lapply(seq_len(ncol(idf)), function(j){
               jamba::cPasteU(idf[,j],
                  sep=contrast_delim)
            }))
         }))
      if (verbose) {jamba::printDebug("imnew:");print(imnew);}
      # check for more factors being compared than expected order contrasts
      num_comp_columns <- sum(sapply(seq_len(ncol(imnew)), function(icol){
         any(grepl(contrast_delim_pattern, imnew[,icol]))
      }));
      if (num_comp_columns > log2(nrow(im))) {
         return(stack_contrasts(
            jamba::pasteByRow(im,
               sep=contrast_factor_delim),
            contrast_delim=contrast_delim))
      }

      diff_cols <- which(sapply(seq_len(ncol(im)), function(j){
         (any(grepl(contrast_delim_pattern, imnew[,j])) &&
            length(unique(imnew[,j])) > 1)
      }));
      if (length(diff_cols) > 0) {
         diff_split <- jamba::pasteByRow(imnew[,diff_cols,drop=FALSE]);
         diff_split <- factor(diff_split, levels=unique(diff_split));
      } else {
         diff_split <- rep(1, nrow(imnew));
      }
      imnew_split <- split(
         data.frame(
            stringsAsFactors=FALSE,
            check.names=FALSE,
            imnew),
         diff_split);
      if (verbose > 1) {jamba::printDebug("imnew_split:");print(imnew_split);}
      imcontrasts <- unique(jamba::pasteByRow(sep=comp_factor_delim,
         jamba::rbindList(lapply(imnew_split, function(imnew1) {
            do.call(cbind, lapply(seq_len(ncol(imnew1)), function(j){
               if (!any(grepl(contrast_delim_pattern, imnew1[,j]))) {
                  jout <- jamba::cPasteU(imnew1[,j],
                     sep=contrast_delim)
               } else {
                  jout <- imnew1[,j]
               }
               jout;
            }))
         }))
      ))
      imcontrasts <- stack_contrasts(imcontrasts,
         contrast_delim=contrast_delim);
      imcontrasts
   }))
   if (add_attr) {
      attr(comps, "contrast_names") <- contrast_names;
   }
   names(comps) <- names(contrast_names);
   comps
}

#' Convert short-form comparison to statistical contrast
#'
#' @rdname contrast2comp
#'
#' @export
comp2contrast <- function
(comps,
 contrast_delim="-",
 contrast_factor_delim="_",
 comp_factor_delim=":",
 factor_order=NULL,
 add_attr=FALSE,
 verbose=FALSE,
 ...)
{
   # validate delim all have 1 character
   check_delim <- function(x){
      (!is.list(x) && length(x) == 1 && !is.na(x) && nchar(x) == 1 &&
            !x %in% c("(", ")", "[", "]"))
   }
   if (!all(check_delim(contrast_delim) &&
         check_delim(contrast_factor_delim) &&
         check_delim(comp_factor_delim))) {
      stop("All delim arguments must be single character strings.")
   }

   # factor_order
   if (length(factor_order) > 0) {
      if ("list" %in% class(factor_order)) {
         factor_order <- rep(factor_order, length(comps));
      } else {
         factor_order <- rep(list(factor_order), length(comps));
      }
   }

   contrast_delim_pattern <- paste0("[", contrast_delim, "]");
   comp_factor_delim_pattern <- paste0("[", comp_factor_delim, "]");

   comps1 <- strsplit(comps, comp_factor_delim_pattern)
   contrast_names <- sapply(seq_along(comps1), function(i){
      # if input contains parentheses, try to re-build piecewise
      if (grepl("[()]", comps[[i]])) {
         comppar <- comps[[i]];
         comppar1 <- strsplit(
            gsub("^[()]+|[()]+$", "", comppar),
            paste0("[)]+[", contrast_delim, "][(]+"))
         imcontrasts <- sapply(comppar1, function(icomppar1){
            stack_contrasts(comp2contrast(icomppar1),
               contrast_delim=contrast_delim)
         });
         return(imcontrasts);
      }
      icomp <- comps1[[i]];
      if (verbose) {
         jamba::printDebug("comp2contrast(): ",
            "comp: ", icomp);
      }
      im <- matrix(nrow=1, icomp)
      if (verbose > 1) {
         jamba::printDebug("comp2contrast(): ",
            "           initial im:");
         print(im);
      }
      jseq <- seq_len(ncol(im));
      if (length(factor_order) > 0 && all(jseq %in% factor_order[[i]])) {
         jseq <- intersect(factor_order[[i]], jseq);
      }
      for (j in jseq) {
         if (any(grepl(contrast_delim_pattern, im[,j]))) {
            iv <- unlist(strsplit(unique(im[,j]), contrast_delim_pattern));
            irow <- rep(seq_len(nrow(im)), length(iv));
            im <- im[irow, , drop=FALSE];
            im[,j] <- rep(iv, each=nrow(im)/length(iv));
            if (verbose > 1) {
               jamba::printDebug("comp2contrast(): ",
                  "      intermediate im:");
               print(im);
            }
         }
      }
      # imcontrasts <- stack_contrasts(im, sep=contrast_factor_delim)
      imgroups <- jamba::pasteByRow(im,
         sep=contrast_factor_delim);
      imcontrasts <- stack_contrasts(imgroups,
         contrast_delim=contrast_delim)

      imseq <- seq(from=1, to=length(imgroups), by=2);
      imcontrasts <- paste0(imgroups[imseq],
         contrast_delim,
         imgroups[imseq + 1])
      imcseq <- rev(seq_len(log2(length(imcontrasts))));
      for (k in imcseq) {
         kn <- (2 ^ k) / 2;
         ksplit <- rep(head(LETTERS, kn), each=2);
         imcontrasts <- sapply(split(imcontrasts, ksplit), function(kcontrast){
            paste0(paste0("(", kcontrast, ")"),
               collapse=contrast_delim)
         })
      }
      imcontrasts
   })
   if (add_attr) {
      attr(contrast_names, "comps") <- comps;
   }
   names(contrast_names) <- names(comps);
   contrast_names
}

#' Stack a vector of contrasts into a full string
#'
#' Stack a vector of contrasts into a full string
#'
#' This is an internal function called by `contrast2comp()` and
#' `comp2contrast()`. It takes a single vector of groups, with
#' 2^(integer) number of entries, and combines them into a
#' series of pairwise contrasts.
#'
#' For example:
#' * `A, B` becomes `A-B`
#' * `A, B, C, D` becomes `(A-B)-(C-D)`
#' * `A, B, C, D, E, F, G, H` becomes `((A-B)-(C-D))-((E-F)-(G-H))`
#'
#' @return `character` vector with fully written contrasts.
#'
#' @param imgroups `character` vector of group names.
#' @param contrast_delim `character` string used as a delimiter between
#'    group names to indicate a statistical contrast.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' stack_contrasts(head(LETTERS, 2))
#' stack_contrasts(head(LETTERS, 4))
#' stack_contrasts(head(LETTERS, 8))
#' contrast2comp(stack_contrasts(head(LETTERS, 16)))
#'
#' # more typical case
#' stacked <- stack_contrasts(c("KO_Treated", "KO_Control",
#'    "WT_Treated", "WT_Control"))
#' stacked
#'
#' # this contrast can be rewritten as a "comp"
#' stacked_comp <- contrast2comp(stacked)
#' stacked_comp
#'
#' # of course the "comp" can be converted back to full contrast name
#' comp2contrast(stacked_comp)
#'
#' # note when the contrast cannot be stacked, an error is thrown
#' tryCatch({
#'    stack_contrasts(head(LETTERS, 5))
#' }, error=function(e){
#'    print(e)
#' })
#'
#' @export
stack_contrasts <- function
(imgroups,
 contrast_delim="-",
 ...)
{
   if ("list" %in% class(imgroups)) {
      return(sapply(imgroups, function(i){
         stack_contrasts(imgroups=i,
            contrast_delim=contrast_delim,
            ...)
      }))
   }
   contrast_delim_pattern <- paste0("[", contrast_delim, "]");

   if (length(imgroups) <= 1) {
      return(imgroups)
   }
   if ((log2(length(imgroups)) %% 1) != 0) {
      print(imgroups);
      stop("length(imgroups) must be a value in 2^(int)");
   }
   imseq <- seq(from=1, to=length(imgroups), by=2);
   imgroups <- ifelse(grepl(contrast_delim_pattern, imgroups),
      paste0("(", imgroups, ")"),
      imgroups);
   imcontrasts <- paste0(imgroups[imseq],
      contrast_delim,
      imgroups[imseq + 1])
   imcseq <- rev(seq_len(log2(length(imcontrasts))));
   for (k in imcseq) {
      kn <- (2 ^ k) / 2;
      ksplit <- rep(head(LETTERS, kn), each=2);
      imcontrasts <- sapply(split(imcontrasts, ksplit), function(kcontrast){
         paste0(paste0("(", kcontrast, ")"),
            collapse=contrast_delim)
      })
   }
   imcontrasts
}


#' Convert contrast to short-form comparison using object names
#'
#' @rdname contrast2comp
#'
#' @export
names_contrast2comp <- function
(contrast_names,
 contrast_delim="-",
 contrast_factor_delim="_",
 comp_factor_delim=":",
 add_attr=FALSE,
 verbose=FALSE,
 ...)
{
   #
   xnames <- names(contrast_names);
   xcomp <- contrast2comp(contrast_names=xnames,
      contrast_delim=contrast_delim,
      contrast_factor_delim=contrast_factor_delim,
      comp_factor_delim=comp_factor_delim,
      add_attr=add_attr,
      verbose=verbose,
      ...)
   names(contrast_names) <- xcomp;
   return(contrast_names);
}


#' Convert short-form comparison to statistical contrast using object names
#'
#' @rdname contrast2comp
#'
#' @export
names_comp2contrast <- function
(comps,
 contrast_delim="-",
 contrast_factor_delim="_",
 comp_factor_delim=":",
 factor_order=NULL,
 add_attr=FALSE,
 verbose=FALSE,
 ...)
{
   #
   xcomp <- names(comps);
   xnames <- comp2contrast(comps=xcomp,
      contrast_delim=contrast_delim,
      contrast_factor_delim=contrast_factor_delim,
      comp_factor_delim=comp_factor_delim,
      factor_order=factor_order,
      add_attr=add_attr,
      verbose=verbose,
      ...)
   names(comps) <- xnames;
   return(comps);
}
