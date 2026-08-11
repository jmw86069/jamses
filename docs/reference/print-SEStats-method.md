# Print SEStats Object

Display a summary of an SEStats object including design information and
hit_array structure.

## Usage

``` r
## S7 method for class <SEStats>
print(x, ...)
```

## Arguments

- x:

  An `SEStats` object to print.

- ...:

  Additional arguments passed to other methods (unused).

## Details

The print method displays:

- Number of samples, groups, and contrasts from the SEDesign

- Dimensions of the hit_array with their names

- Structure of stats_dfs

- Number of metadata items

## See also

Other SEStats objects: [`.class_list`](dot-class_list.md),
[`contrast_names,SEStats-method`](contrast_names-SEStats-method.md),
[`factors,SEStats-method`](factors-SEStats-method.md),
[`groups,SEStats-method`](groups-SEStats-method.md),
[`samples,SEStats-method`](samples-SEStats-method.md)
