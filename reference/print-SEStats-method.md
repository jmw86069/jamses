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

Other SEStats objects:
[`SEStats()`](https://jmw86069.github.io/jamses/reference/SEStats.md),
[`contrast_names,SEStats-method`](https://jmw86069.github.io/jamses/reference/contrast_names-SEStats-method.md),
[`factors,SEStats-method`](https://jmw86069.github.io/jamses/reference/factors-SEStats-method.md),
[`groups,SEStats-method`](https://jmw86069.github.io/jamses/reference/groups-SEStats-method.md),
[`samples,SEStats-method`](https://jmw86069.github.io/jamses/reference/samples-SEStats-method.md)
