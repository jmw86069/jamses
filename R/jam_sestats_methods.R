# #' SEStats and SEDesign Accessor Methods
# #'
# #' Accessor methods for SEStats and SEDesign objects.
# #'
# #' @family SEStats objects
# #' @family jam experiment design
# #'
# NULL

# Create S7 generics if they don't already exist
# These allow methods to be defined for samples, groups, factors, contrast_names

# if (!exists("samples") || !inherits(samples, "S7_generic")) {
#   samples <- S7::new_generic("samples", "object")
# }

# if (!exists("groups") || !inherits(groups, "S7_generic")) {
#   groups <- S7::new_generic("groups", "object")
# }

# if (!exists("factors") || !inherits(factors, "S7_generic")) {
#   factors <- S7::new_generic("factors", "object")
# }

# if (!exists("contrast_names") || !inherits(contrast_names, "S7_generic")) {
#   contrast_names <- S7::new_generic("contrast_names", "object")
# }


#' Get contrast names from SEStats
#'
#' Extract contrast names from the SEStats object's sedesign.
#'
#' This is a read-only accessor that retrieves the contrast names from
#' the SEDesign object stored in the SEStats object. Values cannot be
#' modified through this accessor, as the design is fixed when results
#' are created.
#'
#' @param x An `SEStats` object
#'
#' @return Character vector of contrast names
#'
#' @family SEStats objects
#'
#' @examples
#' if (FALSE) {
#'   contrast_names(sestats_object)
#' }
#'
#' @export
S7::method(contrast_names, SEStats) <- function(object) {
  sedesign <- object@sedesign
  if (is.null(sedesign)) {
    return(character(0))
  }
  # Use the existing contrast_names generic on the SEDesign object
  tryCatch({
    contrast_names(sedesign)
  }, error = function(e) {
    # Fallback if sedesign doesn't have contrast_names defined
    character(0)
  })
}

#' Get factor names from SEStats
#'
#' Extract factor names from the SEStats object's sedesign.
#'
#' This is a read-only accessor that retrieves the factor names from
#' the SEDesign object stored in the SEStats object. Values cannot be
#' modified through this accessor, as the design is fixed when results
#' are created.
#'
#' @param x An `SEStats` object
#'
#' @return Character vector of factor names
#'
#' @family SEStats objects
#'
#' @examples
#' if (FALSE) {
#'   factors(sestats_object)
#' }
#'
#' @export
S7::method(factors, SEStats) <- function(object) {
  sedesign <- object@sedesign
  if (is.null(sedesign)) {
    return(character(0))
  }
  # Use the existing factors generic on the SEDesign object
  tryCatch({
    factors(sedesign)
  }, error = function(e) {
    # Fallback if sedesign doesn't have factors defined
    character(0)
  })
}

#' Get samples from SEStats
#'
#' Extract sample names from the SEStats object's sedesign.
#'
#' This is a read-only accessor that retrieves the sample names from
#' the SEDesign object stored in the SEStats object.
#'
#' @param x An `SEStats` object
#'
#' @return Character vector of sample names
#'
#' @family SEStats objects
#'
#' @export
S7::method(samples, SEStats) <- function(object) {
  sedesign <- object@sedesign
  if (is.null(sedesign)) {
    return(character(0))
  }
  tryCatch({
    samples(sedesign)
  }, error = function(e) {
    character(0)
  })
}

#' Get groups from SEStats
#'
#' Extract group names from the SEStats object's sedesign.
#'
#' This is a read-only accessor that retrieves the group names from
#' the SEDesign object stored in the SEStats object.
#'
#' @param x An `SEStats` object
#'
#' @return Character vector of group names
#'
#' @family SEStats objects
#'
#' @export
S7::method(groups, SEStats) <- function(object) {
  sedesign <- object@sedesign
  if (is.null(sedesign)) {
    return(character(0))
  }
  tryCatch({
    groups(sedesign)
  }, error = function(e) {
    character(0)
  })
}

#' Get dimnames from SEStats
#'
#' Extract the dimension names from an `SEStats`
#' object, how `dimnames()` works on 'hit_array' data.
#'
#' @param x An `SEStats` object
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
#' @returns `list` with named elements: `'Cutoffs'`, `'Contrasts'`, `'Signal'`
#'
#' @examples
#' if (FALSE) {
#'   dimnames(sestats)
#' }
#'
#' @export
S7::method(dimnames, SEStats) <- function(x) {
   if (length(x@hit_array) > 0) {
      dimnames(x@hit_array)
   } else {
      list()
   }
#   list(
#     samples = samples(x),
#     groups = groups(x),
#     contrasts = contrast_names(x)
#   )
}
