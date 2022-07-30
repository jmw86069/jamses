
#' Convert contrast to short-form comp, convert comp to contrast (DEV)
#'
#' Convert contrast to short-form comp, convert comp to contrast (DEV)
#'
#' This method is developmental, intended to provide a vectorized
#' approach to improve speed.
#'
#' @examples
#' contrast_names1 <- c("(A_c-B_c)-(A_d-B_d)",
#'    "A_c-B_c",
#'    "GroupA-GroupB",
#'    "((A_c_J-B_c_J)-(A_d_J-B_d_J))-((A_c_W-B_c_W)-(A_d_W-B_d_W))")
#'
#' contrast2comp_dev(contrast_names1)
#'
#' # Note: this function fails when contrasts are not balanced
#' contrast2comp_dev(c("(A_c-B_c)-(A_d-C_d)"), verbose=TRUE)
#' contrast2comp(c("(A_c-B_c)-(A_d-C_d)"), verbose=TRUE)
#'
#' contrast2comp_dev(c(contrast_names1[1:2], "(A_c-B_c)-(A_d-C_d)"), verbose=TRUE)
#'
#' @export
contrast2comp_dev <- function
(contrast_names,
   verbose=FALSE,
   ...)
{
   # validate each contrast has a unique name
   names_contrast_names <- names(contrast_names);
   if (length(names(contrast_names)) > 0 && length(jamba::tcount(names(contrast_names), 2)) > 0) {
      names(contrast_names) <- jamba::makeNames(names(contrast_names));
   }
   if (length(names(contrast_names)) == 0) {
      names(contrast_names) <- jamba::makeNames(contrast_names);
   }

   # custom function to expand factor
   # takes factor with repeated values
   # makes values in which_x into unique not-repeated values
   expand_factor_x <- function(x, which_x) {
      x <- as.character(x);
      x[which_x] <- jamba::makeNames(x)[which_x]
      x <- factor(x,
         levels=unique(x));
   }

   # split into list of groups
   ilist <- strsplit(contrast_names, "[()-]+")
   if (verbose) {
      jamba::printDebug("contrast2comp_dev(): ",
         "ilist:");
      print(ilist);
   }

   # convert to data.frame
   idf <- data.frame(check.names=FALSE,
      C=rep(names(ilist), lengths(ilist))[nchar(unlist(ilist)) > 0],
      B=rep(names(ilist), lengths(ilist))[nchar(unlist(ilist)) > 0],
      A=jamba::rbindList(strsplit(unlist(ilist), "_")))
   idf$C <- factor(idf$C,
      levels=unique(idf$C))
   idf$B <- factor(idf$B,
      levels=unique(idf$B))

   # assign split to each unique pairwise contrast
   idf_A_colnames <- jamba::nameVector(colnames(idf));
   A_colnames <- tail(idf_A_colnames, -2);
   idf_B_sizes <- table(idf$B)
   idf_B_factors <- rep(names(idf_B_sizes), idf_B_sizes / 2)
   idf_B_split <- rep(jamba::colNum2excelName(seq_along(idf_B_factors)),
      each=2)
   idf_B_split <- factor(idf_B_split, levels=unique(idf_B_split))

   # review verbose output
   idf$split1 <- as.numeric(idf_B_split);
   idf$idf_B_split <- idf_B_split;
   if (verbose) {
      jamba::printDebug("contrast2comp_dev(): ",
         "idf:");
      print(idf);
   }

   # combine pairwise contrasts
   idf_new_m <- do.call(cbind, lapply(idf_A_colnames, function(idf_A_colname){
      jamba::cPasteU(split(idf[[idf_A_colname]], idf_B_split), sep="-")
   }))
   idf_new <- data.frame(check.names=FALSE, idf_new_m);

   # review output for proper balance
   # only one factor column should produce comparison
   idf_new_check <- matrix(ncol=length(A_colnames),
      grepl("-", idf_new_m[, A_colnames, drop=FALSE]))
   if (verbose > 1) {
      jamba::printDebug("contrast2comp_dev(): ",
         "idf_new_m:");
      print(idf_new_m);
      jamba::printDebug("contrast2comp_dev(): ",
         "idf_new_m[, A_colnames, drop=FALSE]:");
      print(idf_new_m[, A_colnames, drop=FALSE]);
      jamba::printDebug("contrast2comp_dev(): ",
         "idf_A_colnames: ", idf_A_colnames);
      jamba::printDebug("contrast2comp_dev(): ",
         "A_colnames: ", A_colnames);
      jamba::printDebug("contrast2comp_dev(): ",
         "idf_new_check:");
      print(idf_new_check);
      print(rowSums(idf_new_check))
   }
   if (any(rowSums(idf_new_check) > 1)) {
      idf_B_split_imbalanced <- levels(idf_B_split)[rowSums(idf_new_check) > 1];
      if (verbose > 1) {
         jamba::printDebug("contrast2comp_dev(): ",
            "idf_B_split_imbalanced: ", idf_B_split_imbalanced);
      }
      which_imbalanced <- idf_B_split %in% idf_B_split_imbalanced;
      new_imbalanced <- jamba::pasteByRow(idf[which_imbalanced, A_colnames, drop=FALSE], sep="_");
      if (verbose > 1) {
         jamba::printDebug("contrast2comp_dev(): ",
            "new_imbalanced: ", new_imbalanced);
      }
      idf[which_imbalanced, A_colnames] <- "";
      idf[which_imbalanced, head(A_colnames, 1)] <- new_imbalanced;
      # idf_B_split <- expand_factor_x(idf_B_split, which_imbalanced);
      # review verbose output
      idf$split1 <- as.numeric(idf_B_split);
      idf$idf_B_split <- idf_B_split;
      if (verbose) {
         jamba::printDebug("contrast2comp_dev(): ",
            "idf:");
         print(idf);
      }
      # combine pairwise contrasts
      idf_new_m <- do.call(cbind, lapply(idf_A_colnames, function(idf_A_colname){
         jamba::cPasteU(split(idf[[idf_A_colname]], idf_B_split), sep="-")
      }))
      idf_new <- data.frame(check.names=FALSE, idf_new_m);
      # return(invisible(idf_new_check));
   }

   # prepare to split into original contrast_names
   idf_new$B <- factor(idf_new$B, levels=unique(idf_new$B))
   idf_new_B_split <- idf_new$B;
   idf_new$split2 <- as.numeric(idf_new_B_split);

   # review verbose output
   if (verbose) {
      jamba::printDebug("contrast2comp_dev(): ",
         "idf_new:");
      print(idf_new);
   }

   # combine remaining factor comparisons
   # TODO: verify what happens when factors are unbalanced
   idf_newer <- data.frame(check.names=FALSE,
      do.call(cbind, lapply(jamba::nameVector(idf_A_colnames), function(idf_A_colname){
         jamba::cPasteU(split(idf_new[[idf_A_colname]], idf_new$B), sep="-")
      })))

   # review verbose output
   if (verbose) {
      jamba::printDebug("contrast2comp_dev(): ",
         "idf_newer:");
      print(idf_newer);
   }

   # check for imbalance at this step
   idf_newer_imbalanced <- sapply(tail(idf_A_colnames, -2), function(idf_A_colname){
      (grepl("-.+-", idf_newer[[idf_A_colname]]))
   });
   if (any(idf_newer_imbalanced)) {
      if (verbose) {
         jamba::printDebug("contrast2comp_dev(): ",
            "two-way imbalance:");
         print(idf_newer_imbalanced);
      }
      which_imbalanced <- Reduce("|", data.frame(idf_newer_imbalanced));
      B_tofix <- idf_newer$B[which_imbalanced];
      if (verbose > 1) {
         jamba::printDebug("contrast2comp_dev(): ",
            "which_imbalanced: ", which_imbalanced);
         jamba::printDebug("contrast2comp_dev(): ",
            "B_tofix: ", B_tofix);
      }
      idf_new_B_tofix <- subset(idf_new, B %in% B_tofix)$B;
      idf_new$C <- "";
      idf_new[idf_new$B %in% B_tofix, "C"] <- gsub("^.+_(v[0-9]+$)", "\\1",
         jamba::makeNames(idf_new_B_tofix))
      idf_new$newB <- jamba::pasteByRowOrdered(idf_new[,c("B", "C"), drop=FALSE]);
      idf_new <- jamba::renameColumn(idf_new,
         from=c("B", "newB"),
         to=c("newB", "B"));
      if (verbose > 1) {
         jamba::printDebug("contrast2comp_dev(): ",
            "idf_new:");
         print(idf_new);
      }
      # re-calculate idf_newer
      idf_newer <- data.frame(check.names=FALSE,
         do.call(cbind, lapply(jamba::nameVector(c(idf_A_colnames, "newB")), function(idf_A_colname){
            jamba::cPasteU(split(idf_new[[idf_A_colname]], idf_new$B), sep="-")
         })))
      # review verbose output
      if (verbose) {
         jamba::printDebug("contrast2comp_dev(): ",
            "idf_newer:");
         print(idf_newer);
      }
   }

   # combine into comps by row
   idf_newer$comp <- jamba::pasteByRow(idf_newer[, tail(idf_A_colnames, -2), drop=FALSE],
      sep=":");

   # review verbose output
   if (verbose) {
      jamba::printDebug("contrast2comp_dev(): ",
         "idf_newer:");
      print(idf_newer);
   }

   # if any are duplicated, we have imbalanced contrast
   if (any(duplicated(idf_newer$newB))) {
      idf_newer$newB <- factor(idf_newer$newB,
         levels=unique(idf_newer$newB));
      dupeB <- idf_newer$newB[duplicated(idf_newer$newB)]
      dupeB_rows <- idf_newer$newB %in% dupeB;
      idf_newer$comp[dupeB_rows] <- paste0("(", idf_newer$comp[dupeB_rows], ")");
      idf_contrasts <- jamba::cPaste(split(idf_newer$comp, idf_newer$newB),
         sep="-")
   } else {
      idf_contrasts <- idf_newer$comp;
   }
   names(idf_contrasts) <- names_contrast_names;
   return(idf_contrasts)
}
