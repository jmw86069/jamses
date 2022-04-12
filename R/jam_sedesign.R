#
# contrast design object

#' Check SEDesign object
#'
#' Check whether a SEDesign object is valid
#'
#' This function checks whether an `SEDesign` object is valid:
#'
#' * if `samples` is provided, and if `design` is provided,
#'   `samples` must match `rownames(design)`.
#' * if `design` is provided, and if `contrasts` is provided,
#'   `colnames(design)` must match `rownames(contrasts)`.
#' * if `contrasts` is provided, `design` must also be provided.
#'
#' Note that `samples` can be a subset of `rownames(design)`,
#' in which case the `design` will also be subset accordingly.
#'
#' Similarly, `colnames(design)` can be a subset of
#' `rownames(contrasts)`, which would force `contrasts`
#' to be subset accordingly.
#'
#' Typically the order of `samples` should match the order of
#' `rownames(design)` but this is not required. Downstream methods
#' should confirm this order.
#'
#' Typically the order of `colnames(design)` should match the order of
#' `rownames(contrast)` but this is not required. Downstream methods
#' should confirm this order.
#'
#' @param object `SEDesign` object
#'
#' @family jam experiment design
#'
#' @export
check_sedesign <- function
(object)
{
   errors <- character();
   if (length(object@samples) > 0 && !all(is.na(object@samples))) {
      if (length(object@design) > 0) {
         if (!all(object@samples) %in% rownames(object@design)) {
            msg <- paste0("Failed: all(samples %in% rownames(design))");
            errors <- c(errors, msg);
         }
      }
   }
   if (length(object@design) > 0) {
      if (length(object@contrasts) > 0) {
         if (!all(colnames(object@design) %in% rownames(object@contrasts))) {
            msg <- paste0("Failed: colnames(design) %in% rownames(contrasts)");
            errors <- c(errors, msg);
         }
      }
   }
   if (length(object@contrasts) > 0) {
      if (length(object@design) == 0) {
         msg <- paste0("Error: contrasts is provided without design, design is required");
         errors <- c(errors, msg);
      }
   }
   if (length(errors) == 0) {
      TRUE
   } else {
      errors
   }
}


setClass("SEDesign",
   slots=c(
      design="matrix",
      contrasts="matrix",
      samples="character"
   ),
   prototype=prototype(
      design=NULL,
      contrasts=NULL,
      samples=character(0)
   ),
   validity=check_sedesign
);



#' Validate SEDesign object contents
#'
#' Validate SEDesign object contents
#'
#' This function validates and enforces constraints on `SEDesign`
#' objects:
#'
#' * `samples` must match `rownames(design)`
#' * `colnames(design)` must match `rownames(contrasts)`
#'
#' If `samples` does not exist, and `rownames(design)` does exist,
#' then `samples` will be defined as `rownames(design)`.
#'
#' If `design` and `samples` are provided, but `rownames(design)`
#' is empty, it must be the same length as `samples`.
#' then `rownames(design)` will be defined as `samples`.
#'
#' @return `SEDesign` object after validation updates have been applied.
#'
#' @param object `SEDesign` object
#' @param min_reps `integer` indicating the minimum required replicate
#'    samples per design group to be used during analysis. Any design
#'    groups with fewer replicates will be removed from the design matrix,
#'    and subsequently will be removed from the contrasts matrix.
#' @param samples,groups,contrasts `character` vectors used to subset
#'    the samples, groups, or contrasts.
#' @param verbose `logical` indicating whether to print verbose output.
#' @param ... additional arguments are ignored.
#'
#' @examples
#' factors2 <- rep(c("one", "two", "three", "four"), each=3)
#' factors2 <- factor(factors2,
#'    levels=unique(factors2))
#' names(factors2) <- paste0("sample", seq_along(factors2))
#' factors2
#'
#' mm2 <- model.matrix(~0 + factors2)
#' rownames(mm2) <- names(factors2)
#' colnames(mm2) <- levels(factors2);
#' mm2
#'
#' icontrastnames <- c("two-one",
#'    "four-three",
#'    "(four-three)-(two-one)");
#' icon <- c(-1, 1, 0, 0,
#'    0, 0, -1, 1,
#'    1, -1, -1, 1)
#' icontrasts2 <- matrix(icon,
#'    ncol=3,
#'    dimnames=list(levels(factors2),
#'       icontrastnames))
#' icontrasts2
#'
#' condes2 <- new("SEDesign",
#'    design=mm2,
#'    contrasts=icontrasts2)
#' condes2
#'
#' # now subset samples
#' validate_sedesign(condes2,
#'    samples=paste0("sample", 2:12))
#'
#' # now subset enough samples to remove one group
#' validate_sedesign(condes2,
#'    samples=paste0("sample", 4:12))
#'
#' validate_sedesign(condes2, groups=c("one", "two"))
#'
#' condes2[, c("one", "two")]
#'
#' @export
validate_sedesign <- function
(object,
 min_reps=1,
 samples=NULL,
 groups=NULL,
 contrasts=NULL,
 verbose=TRUE,
 ...)
{
   newmsg <- character();

   # convert NA to NULL
   if (length(object@samples) > 0 && all(is.na(object@samples))) {
      object@samples <- NULL;
   }
   if (length(object@design) > 0 && all(is.na(object@design))) {
      object@design <- NULL;
   }
   if (length(object@contrasts) > 0 && all(is.na(object@contrasts))) {
      object@contrasts <- NULL;
   }

   # check samples
   #
   # fill missing samples with rownames(design) if possible
   if (length(object@samples) == 0 || all(is.na(object@samples))) {
      # no samples provided
      if (length(object@design) > 0) {
         if (length(rownames(object@design)) > 0) {
            # rownames(design) used to populate samples
            newmsg <- c(newmsg,
               paste0("assigned: samples <- rownames(design)"));
            object@samples <- rownames(object@design);
         } else {
            # no samples, no rownames(design)
            # therefore only integer subsetting is permitted
         }
      } else {
         # no samples, no design, what even is this object
         #return(object)
      }
   }
   #
   # check samples integrity
   if (length(object@samples) > 0) {
      # optionally subset by samples provided
      if (length(samples) > 0) {
         if (is.numeric(samples)) {
            if (max(samples) > length(object@samples)) {
               stop("samples is greater than length(object@samples)");
            }
            if (!all(samples == round(samples))) {
               stop("samples must contain only integer values when supplied as numeric");
            }
            samples <- object@samples[round(samples)];
         }
         if (!all(samples %in% object@samples)) {
            stop("samples supplied must be present in object@samples")
         }
         # subset by this method in order to retain names
         newmsg <- c(newmsg,
            paste0("subset, re-order: object@samples <- samples"));
         object@samples <- object@samples[match(samples, object@samples)];
      }

      # check design
      if (length(object@design) > 0) {
         if (length(rownames(object@design)) == 0) {
            if (nrow(object@design) == length(object@samples)) {
               # assign rownames(design) using samples
               newmsg <- c(newmsg,
                  paste0("assigned: rownames(design) <- samples"));
               rownames(object@design) <- object@samples;
            } else {
               stop("rownames(design) is empty, and nrow(design) must equal length(samples)");
            }
         } else if (length(object@samples) != nrow(object@design) ||
            !all(object@samples == rownames(object@design))) {
            if (all(object@samples %in% rownames(object@design))) {
               # re-order design rows using samples
               newmsg <- c(newmsg,
                  paste0("re-ordered rows: design <- design[samples, , drop=FALSE]"));
               object@design <- object@design[object@samples, , drop=FALSE];
            } else {
               stop("all values in samples must be present in rownames(design)");
            }
         } else {
            # all samples == rownames(design)
         }
      } else {
         # no design provided
      }
   } else if (FALSE) {
      # no samples provided
      if (length(object@design) > 0) {
         if (length(rownames(object@design)) > 0) {
            # rownames(design) used to populate samples
            newmsg <- c(newmsg,
               paste0("assigned: samples <- rownames(design)"));
            object@samples <- rownames(object@design);
         } else {
            # no samples, no rownames(design)
            # therefore only integer subsetting is permitted
         }
      } else {
         # no samples, no design, what even is this object
         return(object)
      }
   }

   # check design
   if (length(object@design) > 0) {
      # design is provided

      # check design groups have at least one sample
      design_reps <- colSums(abs(object@design) > 0, na.rm=TRUE);
      if (any(design_reps < min_reps)) {
         # remove some groups that are not represented by min_reps samples
         # (this step will trigger a filtering step with contrasts later)
         design_group_drop <- colnames(object@design)[design_reps < min_reps];
         jamba::printDebug("Dropped design groups: ",
            design_group_drop);
         newmsg <- c(newmsg,
            paste0("dropped design groups: ",
               jamba::cPaste(design_group_drop,
                  sep=", ")));
         object@design <- object@design[, design_reps >= min_reps, drop=FALSE];
      } else {
         # all groups are represented
      }

      # optionally subset by group
      if (length(groups) > 0) {
         if (is.numeric(groups)) {
            if (max(groups) > ncol(object@design)) {
               stop("groups is greater than ncol(object@design)");
            }
            if (!all(groups == round(groups))) {
               stop("groups must contain only integer values when supplied as numeric");
            }
            groups <- colnames(object@design)[round(groups)];
         }
         if (!all(groups %in% colnames(object@design))) {
            stop("groups must be present in colnames(object@design)")
         }
         if (!all(colnames(object@design) %in% groups) ||
               !all(colnames(object@design) == groups)) {
            # re-order object@contrasts
            newmsg <- c(newmsg,
               paste0("subset: design <- design[, groups, drop=FALSE]"));
            object@design <- object@design[, groups, drop=FALSE];
         } else {
            # groups == colnames(object@design)
            # no further action is necessary
         }
      }

      # check design and contrasts
      if (length(object@contrasts) > 0) {
         # contrasts is provided
         if (!all(colnames(object@design) %in% rownames(object@contrasts))) {
            # design groups are not represented in contrasts
            # (we are choosing not to subset design groups by contrast groups)
            stop("colnames(design) must be defined in rownames(contrasts)")
         } else {
            if (!all(colnames(object@design) == rownames(object@contrasts))) {
               if (!length(colnames(object@design)) == length(rownames(object@contrasts))) {
                  # design groups are a subset of contrast groups
                  contrast_group_drop <- setdiff(rownames(object@contrasts),
                     colnames(object@design));
                  jamba::printDebug("Dropped contrast groups: ",
                     contrast_group_drop);
                  if (length(contrast_group_drop) > 0) {
                     newmsg <- c(newmsg,
                        paste0("dropped contrast groups: ",
                           jamba::cPaste(contrast_group_drop,
                              sep=", ")));
                     contrast_drop <- Reduce("|", lapply(contrast_group_drop, function(i){
                        !object@contrasts[i,] %in% c(0, NA)
                     }));
                     if (any(contrast_drop)) {
                        jamba::printDebug("Dropped contrasts: ",
                           colnames(object@contrasts)[contrast_drop]);
                        newmsg <- c(newmsg,
                           paste0("dropped contrasts: ",
                              jamba::cPaste(colnames(object@contrasts)[contrast_drop],
                                 sep=", ")));
                     }
                  } else {
                     contrast_drop <- rep(FALSE, ncol(object@contrasts));
                     newmsg <- c(newmsg,
                        paste0("re-ordered: contrasts <- contrasts[colnames(design), , drop=FALSE]"));
                  }
                  object@contrasts <- object@contrasts[colnames(object@design), !contrast_drop, drop=FALSE];
               } else {
                  # re-order object@contrasts
                  newmsg <- c(newmsg,
                     paste0("re-ordered: contrasts <- contrasts[colnames(design), , drop=FALSE]"));
                  object@contrasts <- object@contrasts[colnames(object@design), , drop=FALSE];
               }
            } else {
               # all design groups match contrast groups
            }
         }
      } else {
         # no contrasts
      }
   } else {
      # no design
      if (length(object@contrasts) > 0) {
         stop("design must be present when SEDesign contains contrasts.");
      }
   }

   # optional steps regarding counts per contrast

   # print messages
   if (length(newmsg) > 0) {
      for (i in newmsg) {
         jamba::printDebug(i);
      }
   }

   #
   return(object)
}

#' @export
setMethod("[",
   signature=c(x="SEDesign",
      i="ANY",
      j="ANY"),
   definition=function(x, i=NULL, j=NULL, ...) {
      if (missing(i)) {
         i <- NULL;
      }
      if (missing(j)) {
         j <- NULL;
      }
      validate_sedesign(x,
         samples=i,
         groups=j)
   }
)


setGeneric("samples", function(object) {standardGeneric("samples")})

#' @export
setMethod("samples",
   signature=c(object="SEDesign"),
   definition=function(object) {
      rownames(object@design);
   }
)


setGeneric("samples<-", function(object, value) {standardGeneric("samples<-")})

#' @export
setMethod("samples<-",
   signature=c(object="SEDesign",
      value="character"),
   definition=function(object, value) {
      if (length(object@samples) > 0) {
         if (length(value) != length(object@samples)) {
            stop("length(value) must equal length(object@samples)")
         }
         # update object@samples
         object@samples <- value;
         if (length(object@design) > 0) {
            rownames(object@design) <- value;
         }
      } else if (length(object@design) == 0) {
         stop("cannot assign samples when length(object@samples) == 0 and length(object@design) == 0")
      } else if (length(value) != nrow(object@design)) {
         stop("length(value) must equal nrow(object@design)")
      } else {
         # update rownames(object@design)
         rownames(object@design) <- value;
         object@samples <- value;
      }
      object;
   }
)

setGeneric("groups", function(object) {standardGeneric("groups")})

#' @export
setMethod("groups",
   signature=c(object="SEDesign"),
   definition=function(object) {
      colnames(object@design);
   }
)

setGeneric("groups<-", function(object, value) {standardGeneric("groups<-")})

#' @export
setMethod("groups<-",
   signature=c(object="SEDesign",
      value="character"),
   definition=function(object, value) {
      if (length(value) != ncol(object@design)) {
         stop("length(value) must equal ncol(object@design)")
      }
      colnames(object@design) <- value;
      if (length(object@contrasts) > 0) {
         rownames(object@contrasts) <- value;
      }
      object;
   }
)

setGeneric("contrastnames", function(object) {standardGeneric("contrastnames")})

#' @export
setMethod("contrastnames",
   signature=c(object="SEDesign"),
   definition=function(object) {
      colnames(object@contrasts);
   }
)

setGeneric("contrastnames<-", function(object, value) {standardGeneric("contrastnames<-")})

#' @export
setMethod("contrastnames<-",
   signature=c(object="SEDesign",
      value="character"),
   definition=function(object, value) {
      if (length(value) != ncol(object@contrasts)) {
         stop("length(value) must equal ncol(object@contrasts)")
      }
      colnames(object@contrasts) <- value;
      object;
   }
)

setGeneric("design", function(object) {standardGeneric("design")})

#' @export
setMethod("design",
   signature=c(object="SEDesign"),
   definition=function(object) {
      object@design;
   }
)

#setGeneric("design<-", function(object) {standardGeneric("design<-")})

#' @import BiocGenerics
#' @export
setMethod("design<-",
   signature=c(object="SEDesign",
      value="matrix"),
   definition=function(object, value) {
      object@design <- value;
      validate_sedesign(object);
   }
)


#' @export
setMethod("contrasts",
   signature=c(x="SEDesign"),
   definition=function(x) {
      x@contrasts;
   }
)

#' @export
setMethod("contrasts<-",
   signature=c(x="SEDesign",
      how.many="ANY",
      value="ANY"),
   definition=function(x, how.many, value) {
      if (!missing(value)) {
         x@contrasts <- value;
      } else if (!missing(how.many)) {
         x@contrasts <- how.many;
      }
      validate_sedesign(x);
   }
)
