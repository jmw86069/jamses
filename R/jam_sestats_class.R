#' SEStats S7 Object
#'
#' S7 object for storing statistical results from comparative analysis.
#'
#' The `SEStats` object is an S7 class designed to store comprehensive
#' statistical results from differential analysis on SummarizedExperiment data.
#' It maintains relationships between the experimental design, statistical
#' results at different cutoff thresholds, and a 3-dimensional hit array
#' that summarizes significant findings across multiple dimensions.
#'
#' @slot sedesign An `SEDesign` object containing the experimental design
#'   and contrast information.
#' @slot stats_dfs A nested `list` structure organizing data.frames of
#'   statistical results. The structure is `stats_dfs[[signal]][[contrast]]`
#'   where each data.frame contains detailed statistical metrics for a
#'   specific signal type (assay) and contrast comparison.
#' @slot stats_objects An optional `list` for storing raw statistical
#'   objects or intermediate result objects used to generate `stats_dfs`.
#'   Currently unused but provided for future extensibility.
#' @slot hit_array A 3-dimensional `array` with dimnames:
#'   - Dimension 1 (Cutoffs): Statistical threshold variations
#'   - Dimension 2 (Contrasts): Comparison between sample groups
#'   - Dimension 3 (Signal): Different assay names/normalization methods
#'   Cell values are named numeric vectors indicating hit status.
#' @slot metadata An optional `list` containing analysis parameters
#'   and metadata associated with the statistical analysis.
#'
#' @details
#' The `hit_array` is the central data structure, a 3-dimensional array where:
#' - **Cutoffs** are derived from column names in `stats_dfs` data.frames
#'   (columns beginning with "hit ")
#' - **Contrasts** are the names within each stats_dfs element
#' - **Signal** values are the top-level names in stats_dfs
#'
#' The dimensions should satisfy these relationships:
#' - All Signal values in hit_array should appear in `names(stats_dfs)`
#' - All Contrast values should appear in `names(stats_dfs[[signal]])`
#'   for at least one signal
#' - Cutoff names are extracted from data.frame column names
#'
#' @family SEStats objects
#'
#' @docType class
#' @export

# Define S3 class wrappers for non-S7 types
.class_list <- S7::new_S3_class("list")

SEStats <- S7::new_class(
   package = NULL,
   name = "SEStats",
   properties = list(
      sedesign = S7::new_property(
         class = S7::new_union(S7::new_S3_class("NULL"), S7::new_S3_class("SEDesign")),
         default = NULL
      ),
      stats_dfs = S7::new_property(
         class = .class_list,
         default = list()
      ),
      stats_objects = S7::new_property(
         class = .class_list,
         default = list()
      ),
      hit_array = S7::new_property(
         class = S7::new_union(S7::new_S3_class("NULL"), S7::new_S3_class("array")),
         default = NULL
      ),
      metadata = S7::new_property(
         class = .class_list,
         default = list()
      )
   ),
   constructor = function(sedesign = NULL,
                           stats_dfs = list(),
                           stats_objects = list(),
                           hit_array = NULL,
                           metadata = list()) {
      S7::new_object(
      S7::S7_object(),
      sedesign = sedesign,
      stats_dfs = stats_dfs,
      stats_objects = stats_objects,
      hit_array = hit_array,
      metadata = metadata
      )
   }
)

#' Print SEStats Object
#'
#' Display a summary of an SEStats object including design information
#' and hit_array structure.
#'
#' @param x An `SEStats` object to print.
#' @param ... Additional arguments passed to other methods (unused).
#'
#' @details
#' The print method displays:
#' - Number of samples, groups, and contrasts from the SEDesign
#' - Dimensions of the hit_array with their names
#' - Structure of stats_dfs
#' - Number of metadata items
#'
#' @family SEStats objects
#'
#' @export
S7::method(print, SEStats) <- function(x, ...) {
  cat("<SEStats> object summary:\n")
  
  # Design information
  sedesign <- x@sedesign
  if (!is.null(sedesign) && length(sedesign) > 0) {
    # Try to extract design info - handle both S4 and other objects
    n_samples <- tryCatch({
      if (is.S4(sedesign)) {
        length(sedesign@samples)
      } else {
        length(sedesign$samples)
      }
    }, error = function(e) 0)
    
    n_groups <- tryCatch({
      if (is.S4(sedesign)) {
        length(sedesign@groups)
      } else {
        length(sedesign$groups)
      }
    }, error = function(e) 0)
    
    n_contrasts <- tryCatch({
      if (is.S4(sedesign)) {
        ncol(sedesign@contrasts)
      } else {
        ncol(sedesign$contrasts)
      }
    }, error = function(e) 0)
    
    if (n_samples > 0 || n_groups > 0 || n_contrasts > 0) {
      cat(sprintf("- sedesign: %d samples, %d groups, %d contrasts\n",
        n_samples, n_groups, n_contrasts))
    } else {
      cat("- sedesign: (present but empty)\n")
    }
  } else {
    cat("- sedesign: (empty or NULL)\n")
  }
  
  # Hit array dimensions
  hit_array <- x@hit_array
  if (!is.null(hit_array) && length(hit_array) > 0) {
    dimnames_list <- dimnames(hit_array)
    dims <- dim(hit_array)
    
    # Dimension names
    dim_names <- names(dimnames_list)
    if (is.null(dim_names) || all(is.na(dim_names))) {
      dim_names <- c("Dim1", "Dim2", "Dim3")[seq_along(dims)]
    }
    
    # Print each dimension
    for (i in seq_along(dims)) {
      dim_name <- dim_names[i]
      dim_size <- dims[i]
      dim_values <- dimnames_list[[i]]
      
      if (length(dim_values) > 0) {
        values_str <- jamba::cPaste(dim_values, sep=", ", maxItems=5)
        if (length(dim_values) > 5) {
          values_str <- paste0(values_str, ", ... (", length(dim_values), " total)")
        }
      } else {
        values_str <- "(unnamed)"
      }
      
      cat(sprintf("- hit_array %s (%d): %s\n",
        dim_name, dim_size, values_str))
    }
  } else {
    cat("- hit_array: (empty or NULL)\n")
  }
  
  # Stats_dfs structure
  stats_dfs <- x@stats_dfs
  if (length(stats_dfs) > 0) {
    total_dfs <- sum(lengths(stats_dfs))
    signal_names <- names(stats_dfs)
    signals_str <- jamba::cPaste(signal_names, sep=", ", maxItems=5)
    if (length(signal_names) > 5) {
      signals_str <- paste0(signals_str, ", ... (", length(signal_names), " total)")
    }
    cat(sprintf("- stats_dfs: %d signal types, %d total data.frames\n",
      length(stats_dfs), total_dfs))
    cat(sprintf("  signals: %s\n", signals_str))
  } else {
    cat("- stats_dfs: (empty)\n")
  }
  
  # Metadata
  metadata <- x@metadata
  if (length(metadata) > 0) {
    cat(sprintf("- metadata: %d items\n", length(metadata)))
  } else {
    cat("- metadata: (empty)\n")
  }
  
  invisible(x)
}
