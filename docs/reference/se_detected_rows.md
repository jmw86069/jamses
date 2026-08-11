# SummarizedExperiment heuristics to define detected rows

SummarizedExperiment heuristics to define detected rows

## Usage

``` r
se_detected_rows(
  se,
  assay_name = 1,
  group_colnames,
  normgroup_colname = NULL,
  detect_mincounts = 0,
  detect_totalreps = 1,
  detect_minreps = 2,
  detect_minpct = 0.65,
  detect_mingroups = 1,
  isamples = colnames(se),
  verbose = FALSE,
  ...
)
```

## Arguments

- se:

  `SummarizedExperiment` object

- assay_name:

  `character` or `integer` index, referring to the entry in `assays(se)`
  to use when determining valid measurements.

- group_colnames:

  `character` vector of colnames in `colData(se)` which defines sample
  grouping.

- normgroup_colname:

  `character` string with optional colname in `colData(se)` to use for
  normgroups.

- detect_mincounts:

  `numeric` value at or above which a measurement is considered "valid".

- detect_totalreps:

  `numeric` minimum total number of replicates that must contain "valid"
  measurements.

- detect_minreps:

  `numeric` minimum replicates which must contain "valid" measurements
  in any given sample group.

- detect_minpct:

  `numeric` minimum fraction of available replicates in a sample group
  that must contain "valid" measurements.

- detect_mingroups:

  `numeric` minimum number of sample groups that are considered "valid"
  based upon other criteria.

- isamples:

  `character` optional vector of `colnames(se)` to use during this
  analysis. This vector is useful for example, when excluding outlier
  samples that were defined by other methods.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are ignored.

## Value

`list` with the following elements:

- `detected_rows` is a `character` vector of detected `rownames(se)`

- `detected_normgroup` is a `list` of `logical` vectors for each
  normgroup, where the vectors encode whether a row is detected within
  each normgroup.

- `detected_df` is a `data.frame` with summary information for each
  normgroup.

## Details

This function is intended to help apply common logical rules to define
valid, "detected" rows for downstream analysis.

The rules:

- minimum value at or above which a measurement is "valid"

- minimum total replicates with "valid" measurement, across all sample
  columns

- minimum replicates with "valid" measurement required in any sample
  group

- minimum percent replicates with "valid" measurement required in any
  sample group

- minimum sample groups with "valid" criteria above required

### Example

Consider an experiment with 7 groups, and n=3 replicates, which contains
21 total samples.

Assume one row of data that contains 6 "valid" measurements.

- If these 6 "valid" measurements are found in only 2 groups, both
  groups contain n=3 "valid" measurements. This row may have sufficient
  data for a statistical comparison across these two groups.

- However, if the 6 "valid" measurements are also found across 6
  different groups, it may not be suitable for statistical testing.

### Use of normgroups

Detection can be carried out within `"normgroups"`, which are
independent subsets of sample columns. In most cases this method is not
necessary, but is intended when the detected rows should be
independently calculated for two or more subsets of sample columns.

A specific example might be an experiment that measures treatment
effects in two very different tissue types, like lung and muscle. The
detected genes in lung may well not be the same as detected genes in
lung. And in fact, statistical comparisons may not be intended to
compare muscle and lung directly. (That judgement is left to the
analyst.) One may define a column in `colData(se)` that represents
tissue type, with values `"muscle"`, and `"lung"`, then define this
column with argument `normgroup_colname`. The detection will be done
within each independent normgroup, returned as a `list` named
`"detected_normgroup"`. The detected rows are also combined into
`"detected_rows"` which returns rows detected across all normgroups.

## See also

Other jamses SE utilities: [`geomx_to_se()`](geomx_to_se.md),
[`make_se_test()`](make_se_test.md),
[`se_collapse_by_column()`](se_collapse_by_column.md),
[`se_collapse_by_row()`](se_collapse_by_row.md),
[`se_normalize()`](se_normalize.md), [`se_rbind()`](se_rbind.md),
[`se_to_assay_data()`](se_to_assay_data.md),
[`se_to_assay_names()`](se_to_assay_names.md),
[`se_to_rowcoldata()`](se_to_rowcoldata.md)
