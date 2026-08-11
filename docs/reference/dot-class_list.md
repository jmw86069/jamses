# SEStats S7 Object

S7 object for storing statistical results from comparative analysis.

## Usage

``` r
.class_list
```

## Details

The `SEStats` object is an S7 class designed to store comprehensive
statistical results from differential analysis on SummarizedExperiment
data. It maintains relationships between the experimental design,
statistical results at different cutoff thresholds, and a 3-dimensional
hit array that summarizes significant findings across multiple
dimensions.

The `hit_array` is the central data structure, a 3-dimensional array
where:

- **Cutoffs** are derived from column names in `stats_dfs` data.frames
  (columns beginning with "hit ")

- **Contrasts** are the names within each stats_dfs element

- **Signal** values are the top-level names in stats_dfs

The dimensions should satisfy these relationships:

- All Signal values in hit_array should appear in `names(stats_dfs)`

- All Contrast values should appear in `names(stats_dfs[[signal]])` for
  at least one signal

- Cutoff names are extracted from data.frame column names

## Slots

- `sedesign`:

  An `SEDesign` object containing the experimental design and contrast
  information.

- `stats_dfs`:

  A nested `list` structure organizing data.frames of statistical
  results. The structure is `stats_dfs[[signal]][[contrast]]` where each
  data.frame contains detailed statistical metrics for a specific signal
  type (assay) and contrast comparison.

- `stats_objects`:

  An optional `list` for storing raw statistical objects or intermediate
  result objects used to generate `stats_dfs`. Currently unused but
  provided for future extensibility.

- `hit_array`:

  A 3-dimensional `array` with dimnames:

  - Dimension 1 (Cutoffs): Statistical threshold variations

  - Dimension 2 (Contrasts): Comparison between sample groups

  - Dimension 3 (Signal): Different assay names/normalization methods
    Cell values are named numeric vectors indicating hit status.

- `metadata`:

  An optional `list` containing analysis parameters and metadata
  associated with the statistical analysis.

## See also

Other SEStats objects:
[`contrast_names,SEStats-method`](contrast_names-SEStats-method.md),
[`factors,SEStats-method`](factors-SEStats-method.md),
[`groups,SEStats-method`](groups-SEStats-method.md),
[`print,SEStats-method`](print-SEStats-method.md),
[`samples,SEStats-method`](samples-SEStats-method.md)
