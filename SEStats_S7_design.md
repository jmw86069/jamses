
# SEStats design ideas

* Implement as an S7 object.
Properties should include:

   * sedesign: `SEDesign`
   * stats_dfs: `list` which contains a `list` of `data.frame` objects.
   * stats_objects: `list` which is empty for now, but in
   future may contain the specific result object or objects
   depending upon the methodology used.
   It would be the data used to create the 'stats_dfs'
   `data.frame` objects, but would be optional.
   * hit_array: `array` with three dimensions
   * metadata: `list` with optional additional data. For example
   it would hold useful parameters from the analysis, most
   values would be optional.

* Some features of the object.

   * hit_array is a three-dimensional array, with Cutoffs,
   Contrasts, Signal.
   The Signal values should all be present in names(stats_dfs).
   They do not need to be in the same order.
   The Contrast values should all be present in names(stats_dfs[[1]])
   as the next layer of names. They do not need to be in the same
   order, and not every contrast needs to appear in each stats_dfs
   list.
   * The Cutoffs are derived from specific colnames in the `data.frame`
   objects within the stats_dfs, for example stats_dfs[[1]][[1]] would
   contain a `data.frame`.
   The cutoffs are stored in a column "hit cutoffname contrastname".
   You can determine the cutoff used in that `data.frame` by looking
   for columns that begin "hit ", then remove the trailing contrastname,
   with the assumption that the contrastname will never contain
   whitespace.
   In this example "hit cutoffname" is the cutoff name.
   Each `data.frame` may have one or more cutoff columns,
   but most often it only has one column.
   * You do not need to create hit_array, but it would be helpful
   to be able to confirm that the dimnames of hit_array are valid,
   by ensuring each value is also present in the corresponding stats_dfs data.

