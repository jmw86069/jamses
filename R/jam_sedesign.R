#
# contrast design object

#' @import methods
NULL

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
#' @returns `logical` TRUE if valid, or `character` vector of errors.
#'
#' @export
check_sedesign <- function(object) {
   errors <- character()
   if (
      length(object@samples) > 0 &&
         !all(is.na(object@samples)) &&
         length(object@design) > 0 &&
         !all(object@samples %in% rownames(object@design))
   ) {
      msg <- paste0("Failed: all(samples %in% rownames(design))")
      errors <- c(errors, msg)
   }
   if (
      length(object@design) > 0 &&
         length(object@contrasts) > 0 &&
         !all(colnames(object@design) %in% rownames(object@contrasts))
   ) {
      msg <- paste0("Failed: colnames(design) %in% rownames(contrasts)")
      errors <- c(errors, msg)
   }
   if (
      length(object@contrasts) > 0 &&
         length(object@design) == 0
   ) {
      msg <- paste0(
         "Error: contrasts is provided without design, design is required"
      )
      errors <- c(errors, msg)
   }
   if (length(errors) == 0) {
      TRUE
   } else {
      errors
   }
}


# S7 class type used for the design/contrasts matrix properties
.class_matrix <- S7::new_S3_class("matrix");

#' Build internal per-group factor decomposition data.frame
#'
#' For internal use. Splits `group_names` (`colnames(design)`, i.e.
#' `groups(object)`) on `"_"` into one column per underlying experimental
#' factor. Column names default to `"factor1"`, `"factor2"`, etc., unless
#' `factor_labels` is supplied with a matching length, in which case it
#' is used instead (this is how user-customized `factors()` labels are
#' applied).
#'
#' @noRd
.sedesign_build_design_df <- function
(group_names,
 factor_labels=character(0))
{
   if (length(group_names) == 0) {
      return(data.frame());
   }
   design_matrix <- jamba::rbindList(strsplit(group_names, "_", fixed=TRUE));
   design_matrix <- matrix(design_matrix,
      nrow=length(group_names));
   if (length(factor_labels) == ncol(design_matrix)) {
      colnames(design_matrix) <- factor_labels;
   } else {
      colnames(design_matrix) <- paste0("factor", seq_len(ncol(design_matrix)));
   }
   design_df <- data.frame(check.names=FALSE,
      stringsAsFactors=FALSE,
      design_matrix);
   rownames(design_df) <- group_names;
   return(design_df);
}

#' Build (or refresh) the cached contrasts_to_factors() data.frame
#'
#' For internal use. This step can be slow, so its result is cached
#' on the `SEDesign` object in the `contrasts_df` property, and is
#' only recomputed when `factors`, `design`, or `contrasts` change.
#'
#' @noRd
.sedesign_build_contrasts_df <- function
(object)
{
   if (length(object@contrasts) == 0 || ncol(object@contrasts) == 0) {
      return(data.frame());
   }
   tryCatch({
      contrasts_to_factors(object);
   }, error=function(e){
      data.frame();
   });
}

#' Refresh internal design_df and contrasts_df caches
#'
#' `design_df` is rebuilt from `colnames(design(object))` (i.e.
#' `groups(object)`), reusing existing `factors()` labels when their
#' count still matches the number of underlying factor columns; the
#' `factors` property is then synced to whatever labels were actually
#' used (custom, or auto-generated `"factor1"`, `"factor2"`, etc.).
#' `contrasts_df` depends only on `design`/`contrasts`, not `factors`.
#'
#' @noRd
.sedesign_refresh_caches <- function(object) {
   group_names <- colnames(object@design)
   if (is.null(group_names)) {
      group_names <- character(0)
   }
   design_df <- .sedesign_build_design_df(group_names, object@factors)
   object@design_df <- design_df
   object@factors <- colnames(design_df)
   object@contrasts_df <- .sedesign_build_contrasts_df(object)
   return(object)
}


#' SEDesign: experiment design and contrasts object (S7 class)
#'
#' `SEDesign` enforces the relationship between individual samples,
#' design groups, and groups involved in statistical contrasts.
#'
#' In addition to the four constructor arguments below, `SEDesign` objects
#' carry two internal-use properties, accessible via `@`:
#'
#' * `design_df`: `data.frame` with one row per design group
#'    (`groups(object)`, i.e. `colnames(design)`), derived by splitting
#'    each group name on `"_"` into one column per underlying
#'    experimental factor. Column names default to `"factor1"`,
#'    `"factor2"`, etc., and can be customized via `factors()<-`.
#' * `contrasts_df`: `data.frame`, cached result of `contrasts_to_factors()`,
#'    refreshed whenever `design` or `contrasts` are updated.
#'
#' @param design `matrix` with rownames as samples, colnames as
#'    design groups, containing 0/1 sample-to-group association values.
#' @param contrasts `matrix` with rownames matching `colnames(design)`,
#'    and colnames as contrast names, containing contrast coefficients.
#' @param samples `character` vector of sample identifiers, typically
#'    equal to `rownames(design)`.
#' @param factors `character` vector of labels for the underlying
#'    experimental factors, equivalent to `colnames(design_df)`, i.e.
#'    the columns produced by splitting each design group name
#'    (`groups(object)`) on `"_"`. When empty, or when its length does
#'    not match the number of underlying factors, it is auto-populated
#'    with `"factor1"`, `"factor2"`, etc. Editing `factors()` only
#'    renames these labels (used primarily by `print()`); it has no
#'    effect on `design`, `contrasts`, or `groups()`.
#'
#' @family jam experiment design
#'
#' @export
SEDesign <- S7::new_class("SEDesign",
   package=NULL,
   properties=list(
      design=.class_matrix,
      contrasts=.class_matrix,
      samples=S7::class_character,
      factors=S7::class_character,
      design_df=S7::class_data.frame,
      contrasts_df=S7::class_data.frame
   ),
   constructor=function(
      design=matrix(nrow=0, ncol=0),
      contrasts=matrix(nrow=0, ncol=0),
      samples=character(0),
      factors=character(0)
   ) {
      if (length(samples) == 0) {
         samples <- rownames(design)
      }
      object <- S7::new_object(S7::S7_object(),
         design=design,
         contrasts=contrasts,
         samples=samples,
         factors=factors,
         design_df=data.frame(),
         contrasts_df=data.frame());
      .sedesign_refresh_caches(object);
   },
   validator=function(self) {
      result <- check_sedesign(self);
      if (isTRUE(result)) {
         NULL
      } else {
         result
      }
   }
);

# register SEDesign as an S4 class, so existing S4 setGeneric/setMethod
# machinery (used for design(), contrasts(), "[", etc.) continues to work.
S7::S4_register(SEDesign);


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
#' @family jam experiment design
#'
#' @param object `SEDesign` object
#' @param min_reps `integer` indicating the minimum required replicate
#'    samples per design group to be used during analysis. Any design
#'    groups with fewer replicates will be removed from the design matrix,
#'    and subsequently will be removed from the contrasts matrix.
#' @param samples,groups,contrasts `character` vectors used to subset
#'    the samples, groups, or contrasts.
#' @param verbose `logical` indicating whether to print verbose output,
#'    where `verbose=TRUE` will print messages at the end of operations,
#'    and `verbose=2` will also print messages during operations.
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
#' condes2 <- SEDesign(
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
#' condes2[, c("one", "two"), ]
#'
#' @export
validate_sedesign <- function
(object,
 min_reps=1,
 samples=NULL,
 groups=NULL,
 contrasts=NULL,
 verbose=FALSE,
 ...)
{
   # S7 validates on every single `@<-` assignment (unlike S4, which only
   # validates at construction). .validate_sedesign_core() below performs
   # a sequence of property mutations that can be transiently "invalid"
   # partway through (e.g. design and samples briefly out of sync), so
   # the whole operation is wrapped in valid_eventually() to suppress
   # validation until the final, fully-reconciled object is produced.
   S7::valid_eventually(object, function(object) {
      .validate_sedesign_core(object,
         min_reps=min_reps,
         samples=samples,
         groups=groups,
         contrasts=contrasts,
         verbose=verbose,
         ...)
   })
}

#' @noRd
.validate_sedesign_core <- function(
   object,
   min_reps = 1,
   samples = NULL,
   groups = NULL,
   contrasts = NULL,
   verbose = FALSE,
   ...
) {
   newmsg <- character()

   # convert NA to empty
   if (length(object@samples) > 0 && all(is.na(object@samples))) {
      object@samples <- character(0)
   }
   if (length(object@design) > 0 && all(is.na(object@design))) {
      object@design <- matrix(nrow = 0, ncol = 0)
   }
   if (length(object@contrasts) > 0 && all(is.na(object@contrasts))) {
      object@contrasts <- matrix(nrow = 0, ncol = 0)
   }

   # check samples
   #
   # fill missing samples with rownames(design) if possible
   if (length(object@samples) == 0 || all(is.na(object@samples))) {
      # no samples provided
      if (length(object@design) > 0) {
         if (length(rownames(object@design)) > 0) {
            # rownames(design) used to populate samples
            newmsg <- c(newmsg, paste0("assigned: samples <- rownames(design)"))
            object@samples <- rownames(object@design)
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
               cli::cli_abort(paste0(
                  "{.var samples} is greater than {.code length(object@samples)}"
               ))
            }
            if (!all(samples == round(samples))) {
               cli::cli_abort(paste0(
                  "{.var samples} must contain only integer values when ",
                  "supplied as {.cls numeric}"
               ))
            }
            samples <- object@samples[round(samples)]
         }
         if (!all(samples %in% object@samples)) {
            cli::cli_abort(paste0(
               "{.var samples} supplied must be present in {.code object@samples}"
            ))
         }
         # subset by this method in order to retain names
         newmsg <- c(
            newmsg,
            paste0("subset, re-order: object@samples <- samples")
         )
         object@samples <- object@samples[match(samples, object@samples)]
      }

      # check design
      if (length(object@design) > 0) {
         if (length(rownames(object@design)) == 0) {
            if (nrow(object@design) == length(object@samples)) {
               # assign rownames(design) using samples
               newmsg <- c(
                  newmsg,
                  paste0("assigned: rownames(design) <- samples")
               )
               rownames(object@design) <- object@samples
            } else {
               cli::cli_abort(paste0(
                  "{.code rownames(design)} is empty, and ",
                  "{.code nrow(design)} must equal {.code length(samples)}"
               ))
            }
         } else if (
            length(object@samples) != nrow(object@design) ||
               !all(object@samples == rownames(object@design))
         ) {
            # confirm all samples are present in rownames
            if (all(object@samples %in% rownames(object@design))) {
               # re-order design rows using samples
               newmsg <- c(
                  newmsg,
                  paste0(
                     "re-ordered rows: design <- design[samples, , drop=FALSE]"
                  )
               )
               object@design <- object@design[object@samples, , drop = FALSE]
            } else {
               cli::cli_abort(paste0(
                  "all values in {.var samples} must be present in {.code rownames(design)}"
               ))
            }
         } else {
            # all samples == rownames(design)
         }
      } else {
         # no design provided
      }
   } else if (sqrt(5) == 2) {
      # no samples provided
      if (length(object@design) > 0) {
         if (length(rownames(object@design)) > 0) {
            # rownames(design) used to populate samples
            newmsg <- c(newmsg, paste0("assigned: samples <- rownames(design)"))
            object@samples <- rownames(object@design)
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
      design_reps <- colSums(abs(object@design) > 0, na.rm = TRUE)
      if (any(design_reps < min_reps)) {
         # remove some groups that are not represented by min_reps samples
         # (this step will trigger a filtering step with contrasts later)
         design_group_drop <- colnames(object@design)[design_reps < min_reps]
         if (verbose > 1) {
            jamba::printDebug("Dropped design groups: ", design_group_drop)
         }
         newmsg <- c(
            newmsg,
            paste0(
               "dropped design groups: ",
               jamba::cPaste(design_group_drop, sep = ", ")
            )
         )
         object@design <- object@design[, design_reps >= min_reps, drop = FALSE]
      } else {
         # all groups are represented
      }

      # optionally subset by group
      if (length(groups) > 0) {
         if (is.numeric(groups)) {
            if (max(groups) > ncol(object@design)) {
               cli::cli_abort(paste0(
                  "{.var groups} is greater than {.code ncol(object@design)}"
               ))
            }
            if (!all(groups == round(groups))) {
               cli::cli_abort(paste0(
                  "{.var groups} must contain only integer values when supplied as {.cls numeric}"
               ))
            }
            groups <- colnames(object@design)[round(groups)]
         }
         if (!all(groups %in% colnames(object@design))) {
            cli::cli_abort(paste0(
               "{.var groups} must be present in {.code colnames(object@design)}"
            ))
         }
         if (
            !all(colnames(object@design) %in% groups) ||
               !all(colnames(object@design) == groups)
         ) {
            # re-order object@contrasts
            newmsg <- c(
               newmsg,
               paste0("subset: design <- design[, groups, drop=FALSE]")
            )
            object@design <- object@design[, groups, drop = FALSE]
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
            if (
               ncol(object@design) != nrow(object@contrasts) ||
                  !all(colnames(object@design) == rownames(object@contrasts))
            ) {
               if (
                  length(colnames(object@design)) !=
                     length(rownames(object@contrasts))
               ) {
                  # design groups are a subset of contrast groups
                  contrast_group_drop <- setdiff(
                     rownames(object@contrasts),
                     colnames(object@design)
                  )
                  if (verbose > 1) {
                     jamba::printDebug(
                        "Dropped contrast groups: ",
                        contrast_group_drop
                     )
                  }
                  if (length(contrast_group_drop) > 0) {
                     newmsg <- c(
                        newmsg,
                        paste0(
                           "dropped contrast groups: ",
                           jamba::cPaste(contrast_group_drop, sep = ", ")
                        )
                     )
                     contrast_drop <- Reduce(
                        "|",
                        lapply(contrast_group_drop, function(i) {
                           !object@contrasts[i, ] %in% c(0, NA)
                        })
                     )
                     if (any(contrast_drop)) {
                        if (verbose > 1) {
                           jamba::printDebug(
                              "Dropped contrasts: ",
                              colnames(object@contrasts)[contrast_drop]
                           )
                        }
                        newmsg <- c(
                           newmsg,
                           paste0(
                              "dropped contrasts: ",
                              jamba::cPaste(
                                 colnames(object@contrasts)[contrast_drop],
                                 sep = ", "
                              )
                           )
                        )
                     }
                  } else {
                     contrast_drop <- rep(FALSE, ncol(object@contrasts))
                     newmsg <- c(
                        newmsg,
                        paste0(
                           "re-ordered: contrasts <- contrasts[colnames(design), , drop=FALSE]"
                        )
                     )
                  }
                  object@contrasts <- object@contrasts[
                     colnames(object@design),
                     !contrast_drop,
                     drop = FALSE
                  ]
               } else {
                  # re-order object@contrasts
                  newmsg <- c(
                     newmsg,
                     paste0(
                        "re-ordered: contrasts <- contrasts[colnames(design), , drop=FALSE]"
                     )
                  )
                  object@contrasts <- object@contrasts[
                     colnames(object@design),
                     ,
                     drop = FALSE
                  ]
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
         stop("design must be present when SEDesign contains contrasts.")
      }
   }

   # optional steps regarding counts per contrast

   # # subset by contrasts, which does not trigger subset of design
   if (length(contrasts) > 0 && length(object@contrasts) > 0) {
      #
      if (!all(contrasts %in% colnames(object@contrasts))) {
         missing_contrasts <- setdiff(
            colnames(object@contrasts),
            contrasts
         )
         missing_contrasts <- missing_contrasts;
         cli::cli_abort(c(
            "Not all {.var contrasts} found in {.var colnames(object@contrasts)}.",
            " Missing contrasts include: {missing_contrasts}"
         ))
      }
      object@contrasts <- object@contrasts[, contrasts, drop = FALSE]
   }

   # print messages
   if (length(newmsg) > 0 && verbose) {
      for (i in newmsg) {
         jamba::printDebug(i)
      }
   }

   # keep design_df, contrasts_df caches synchronized with any
   # changes made to design/contrasts above.
   object <- .sedesign_refresh_caches(object)

   #
   return(object)
}

#' Subset a SEDesign object by samples and/or design groups
#'
#' @param x `SEDesign` object
#' @param i,j optional `character` or `integer` vectors used to subset
#'    samples (`i`) and/or design groups (`j`). See `validate_sedesign()`.
#' @param ... additional arguments are ignored
#'
#' @family jam experiment design
#'
#' @export
S7::method(`[`, SEDesign) <- function(x, i=NULL, j=NULL, k=NULL, ...) {
   if (missing(i)) {
      i <- NULL;
   }
   if (missing(j)) {
      j <- NULL;
   }
   validate_sedesign(x,
      samples=i,
      groups=j,
      contrasts=k)
}


#' Sample identifiers for SEDesign objects
#'
#' `samples()` returns sample identifiers, equivalent to
#' `rownames(design(object))`. `samples()<-` renames samples, and
#' updates `rownames(design(object))` to match. It accepts a `character`
#' vector `value` with length equal to `length(samples(object))`.
#'
#' @param object `SEDesign` object
#' @param ... additional arguments are ignored
#'
#' @family jam experiment design
#'
#' @export
samples <- S7::new_generic("samples", "object")

#' @export
S7::method(samples, SEDesign) <- function(object) {
   rownames(object@design);
}


#' @rdname samples
#' @param value `character` vector, length equal to `length(samples(object))`
#'
#' @export
`samples<-` <- S7::new_generic("samples<-", "object")

#' @export
S7::method(`samples<-`, SEDesign) <- function(object, value) {
   value <- as.character(value);
   # renaming samples and rownames(design) together can be transiently
   # inconsistent between the two assignments, so both are performed
   # inside valid_eventually() and only validated once, at the end.
   S7::valid_eventually(object, function(object) {
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
   })
}

#' Design group names for SEDesign objects
#'
#' `groups()` returns the design group names, equivalent to
#' `colnames(design(object))`. `groups()<-` renames design groups,
#' updating `colnames(design)` and `rownames(contrasts)` together.
#'
#' In principle, design groups are derived identifiers and should not
#' need to be renamed directly (prefer re-running `groups_to_sedesign()`
#' with updated inputs). `groups()<-` is provided for convenience, but
#' use it with care: `design_df` (and therefore `factors()` labels
#' derived from group names) is rebuilt from the new group names, which
#' resets any customized `factors()<-` labels if the number of
#' underlying factors changes.
#'
#' @param object `SEDesign` object
#' @param ... additional arguments are ignored
#'
#' @family jam experiment design
#'
#' @export
groups <- S7::new_generic("groups", "object")

#' @export
S7::method(groups, SEDesign) <- function(object) {
   colnames(object@design);
}

#' @rdname groups
#' @param value `character` vector, length equal to `ncol(design(object))`
#'
#' @export
`groups<-` <- S7::new_generic("groups<-", "object")

#' @export
S7::method(`groups<-`, SEDesign) <- function(object, value) {
   value <- as.character(value);
   if (anyDuplicated(value)) {
      stop("groups must not contain duplicated values.")
   }
   # colnames(design) and rownames(contrasts) are transiently out of
   # sync between these two assignments, so both are performed inside
   # valid_eventually() and only validated once, at the end.
   S7::valid_eventually(object, function(object) {
      if (length(object@design) > 0) {
         if (length(value) != ncol(object@design)) {
            stop("length(value) must equal ncol(object@design)")
         }
         colnames(object@design) <- value;
         if (length(object@contrasts) > 0) {
            rownames(object@contrasts) <- value;
         }
      }
      .sedesign_refresh_caches(object);
   })
}

#' Experimental factor labels for SEDesign objects
#'
#' `factors()` returns the labels of the underlying experimental
#' factors, equivalent to `colnames(x@design_df)`, i.e. the columns
#' produced by splitting each design group name (`groups(object)`) on
#' `"_"`. `factors()<-` renames these labels; it has no effect on
#' `design`, `contrasts`, `groups()`, or `contrasts_df` -- it is used
#' only for internal bookkeeping and for the `print()` summary.
#'
#' @param object `SEDesign` object
#' @param ... additional arguments are ignored
#'
#' @family jam experiment design
#'
#' @export
factors <- S7::new_generic("factors", "object")

#' @export
S7::method(factors, SEDesign) <- function(object) {
   colnames(object@design_df);
}

#' @rdname factors
#' @param value `character` vector, length equal to `ncol(x@design_df)`
#'    (the number of underlying experimental factors)
#'
#' @export
`factors<-` <- S7::new_generic("factors<-", "object")

#' @export
S7::method(`factors<-`, SEDesign) <- function(object, value) {
   value <- as.character(value)
   if (anyDuplicated(value)) {
      stop("factors must not contain duplicated values.")
   }
   if (ncol(object@design_df) > 0) {
      if (length(value) != ncol(object@design_df)) {
         cli::cli_abort(paste0(
            "{.val length(value)} must equal ",
            "{.val ncol(object@design_df)},",
            " the number of underlying experimental factors."
         ))
      }
      colnames(object@design_df) <- value
   }
   if (ncol(object@contrasts_df) > 0) {
      colnames(object@contrasts_df) <- value
   }
   object@factors <- value
   object
}

#' Contrast names for SEDesign objects
#'
#' `contrastnames()` and `contrast_names()` are equivalent, returning
#' `colnames(contrasts(object))`. `contrastnames()<-` renames contrast
#' columns directly, accepting a `character` vector `value` with length
#' equal to `ncol(contrasts(object))`. `contrast_names()<-` instead
#' rebuilds the contrasts matrix using `limma::makeContrasts()`, accepting
#' a `character` vector `value` of contrast names (must not contain
#' duplicated values).
#'
#' @param object `SEDesign` object
#' @param ... additional arguments are ignored
#' @param value `character` vector: contrast column names (for
#'    `contrastnames<-`) or contrast names to rebuild via
#'    `limma::makeContrasts()` (for `contrast_names<-`)
#'
#' @family jam experiment design
#'
#' @export
contrastnames <- S7::new_generic("contrastnames", "object")

#' @export
S7::method(contrastnames, SEDesign) <- function(object) {
   colnames(object@contrasts);
}

#' @rdname contrastnames
#' @export
`contrastnames<-` <- S7::new_generic("contrastnames<-", "object")

#' @export
S7::method(`contrastnames<-`, SEDesign) <- function(object, value) {
   value <- as.character(value);
   if (length(value) != ncol(object@contrasts)) {
      stop("length(value) must equal ncol(object@contrasts)")
   }
   colnames(object@contrasts) <- value;
   object <- .sedesign_refresh_caches(object);
   object;
}

#' Design matrix accessor for SEDesign objects
#'
#' @param object `SEDesign` object
#' @param ... additional arguments are ignored
#'
#' @family jam experiment design
#'
#' @import BiocGenerics
#' @export
setMethod("design",
   signature=c(object="SEDesign"),
   definition=function(object, ...) {
      object@design;
   }
)

#' Design matrix setter for SEDesign objects
#'
#' @param object `SEDesign` object
#' @param value `matrix` whose `colnames` must match `groups(object)`
#'    exactly (when `groups(object)` is already defined); use
#'    `groups(object) <- ...` to rename design groups instead.
#' @param ... additional arguments are ignored
#'
#' @family jam experiment design
#'
#' @importFrom BiocGenerics design
#' @export
setMethod("design<-",
   signature=c(object="SEDesign",
      value="matrix"),
   definition=function(object, ..., value) {
      current_groups <- groups(object);
      value_colnames <- colnames(value);
      if (length(current_groups) > 0 && length(value_colnames) > 0 &&
            !identical(value_colnames, current_groups)) {
         stop(paste0("colnames(value) must match groups(object) exactly. ",
            "Use groups(object) <- ... to rename design groups."));
      }
      # replacing design can transiently break the samples/design
      # alignment (fixed up by .validate_sedesign_core()), so both
      # steps happen inside valid_eventually() and are only validated
      # once, at the end.
      S7::valid_eventually(object, function(object) {
         object@design <- value;
         .validate_sedesign_core(object);
      })
   }
)

#' Contrast matrix accessor for SEDesign objects
#'
#' @param x `SEDesign` object
#' @param contrasts,sparse arguments retained for compatibility with the
#'    base `stats::contrasts()` generic signature; ignored for `SEDesign`.
#'
#' @family jam experiment design
#'
#' @export
# S4 generic to dispatch S3 classes
setGeneric("contrasts")

#' @export
setMethod("contrasts",
   signature=c(x="SEDesign"),
   definition=function(x, contrasts=TRUE, sparse=FALSE) {
      x@contrasts;
   }
)

#' Contrast matrix setter for SEDesign objects
#'
#' @param x `SEDesign` object
#' @param how.many,value arguments retained for compatibility with the
#'    base `stats::contrasts<-()` generic signature; `value` (or
#'    `how.many` when `value` is not supplied) must be a `matrix` whose
#'    `rownames` match `groups(x)` exactly.
#'
#' @family jam experiment design
#'
#' @export
# S4 generic to dispatch S3 classes
setGeneric("contrasts<-")

#' @export
setMethod("contrasts<-",
   signature=c(x="SEDesign",
      how.many="ANY",
      value="ANY"),
   definition=function(x, how.many, value) {
      new_contrasts <- if (!missing(value)) value else how.many;
      current_groups <- groups(x);
      new_rownames <- rownames(new_contrasts);
      if (length(current_groups) > 0 && length(new_rownames) > 0 &&
            !identical(new_rownames, current_groups)) {
         stop(paste0("rownames(value) must match groups(x) exactly. ",
            "Use groups(x) <- ... to rename design groups."));
      }
      S7::valid_eventually(x, function(x) {
         x@contrasts <- new_contrasts;
         .validate_sedesign_core(x);
      })
   }
)


#' @rdname contrastnames
#' @export
contrast_names <- S7::new_generic("contrast_names", "object")

#' @export
S7::method(contrast_names, SEDesign) <- function(object) {
   colnames(object@contrasts);
}

#' @rdname contrastnames
#' @export
`contrast_names<-` <- S7::new_generic("contrast_names<-", "object")

#' @export
S7::method(`contrast_names<-`, SEDesign) <- function(object, value) {
   if (anyDuplicated(value)) {
      stop("contrast_names cannot be duplicated.")
   }
   object <- tryCatch({
      S7::valid_eventually(object, function(object) {
         object@contrasts <- limma::makeContrasts(
            contrasts=value,
            levels=object@design)
         .validate_sedesign_core(object);
      })
   }, error=function(e){
      cli::cli_alert_danger(paste(
         "{.field contrast_names} were not compatible with",
         "{.field design(sedesign)}.",
         "No subsetting was performed."))
      cli::cli_alert_danger(
         cli::format_error(e))
      object
   })
   object
}

#' Print / show method for SEDesign objects
#'
#' S7 objects bypass the classic S4 `show` generic for console
#' auto-printing; a plain S3 `print.SEDesign` method is used instead,
#' since S7 objects retain "SEDesign" in their S3 `class()` vector.
#' 
#' The summary includes:
#' 
#' * the number of samples, groups, and contrasts
#' * Each factor name, with factor levels.
#' 
#' Factor levels are printed in order, so that the first entry
#' is prioritized as the control in subsequent contrasts.
#' 
#' For example: 'level1', 'level2', 'level3' should always produce
#' contrasts 'level2-level1', 'level3-level1', 'level3-level2'.
#' It should not produce a contrast 'level2-level3'.
#'
#' @param x `SEDesign` object
#' @param ... additional arguments are ignored
#'
#' @family jam experiment design
#' @returns This function is called for its output to console.
#'    It invisibly returns the input object.
#' @export
S7::method(print, SEDesign) <- function(x, ...) {
   cat(sprintf("<SEDesign> %d samples, %d groups, %d contrasts\n",
      length(samples(x)),
      length(groups(x)),
      ncol(x@contrasts)));
   design_df <- x@design_df;
   if (ncol(design_df) > 0) {
      cat("factors:\n");
      for (icol in colnames(design_df)) {
         ilevels <- unique(design_df[[icol]]);
         cat(sprintf("  - %s: %s\n",
            icol,
            jamba::cPaste(ilevels, sep=", ")));
      }
   }
   invisible(x);
}


#' Get dimnames from SEDesign
#'
#' Extract the dimension names (samples, groups, contrast_names) from an
#' SEDesign object, similar to how `dimnames()` works on arrays or matrices.
#'
#' @param x An `SEDesign` object
#'
#' @details
#' Returns a `list` with three named elements:
#' - `samples`: Sample names (from `samples(x)`)
#' - `groups`: Group names (from `groups(x)`)
#' - `contrasts`: Contrast names (from `contrast_names(x)`)
#'
#' This provides a unified way to access all dimension names from a SEDesign
#' object, consistent with base R's `dimnames()` function for arrays and
#' matrices.
#'
#' @returns `list` with named elements: `$samples`, `$groups`, `$contrasts`
#'
#' @examples
#' if (FALSE) {
#'   dimnames(sedesign)
#' }
#'
#' @export
S7::method(dimnames, SEDesign) <- function(x) {
  list(
    samples = samples(x),
    groups = groups(x),
    contrasts = contrast_names(x)
  )
}

#' Get dimensions for SEDesign
#'
#' Summarize the dimensions using samples, groups, contrast_names from an
#' `SEDesign` object.
#'
#' @param x An `SEDesign` object
#'
#' @details
#' Returns a `list` with three named elements:
#' - `samples`: Sample names (from `samples(x)`)
#' - `groups`: Group names (from `groups(x)`)
#' - `contrasts`: Contrast names (from `contrast_names(x)`)
#'
#' This method provides a unified way to access all dimensions for an
#' `SEDesign` object.
#'
#' @returns `integer` vector with lengths of `samples`, `groups`, `contrasts`
#'
#' @examples
#' if (FALSE) {
#'   dim(sedesign)
#' }
#'
#' @export
S7::method(dim, SEDesign) <- function(x) {
   lengths(dimnames(x))
}
