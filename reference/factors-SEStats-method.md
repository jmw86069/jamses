# Get factor names from SEStats

Extract factor names from the SEStats object's sedesign.

## Usage

``` r
## S7 method for class <SEStats>
factors(object)
```

## Arguments

- x:

  An `SEStats` object

## Value

Character vector of factor names

## Details

This is a read-only accessor that retrieves the factor names from the
SEDesign object stored in the SEStats object. Values cannot be modified
through this accessor, as the design is fixed when results are created.

## See also

Other SEStats objects:
[`SEStats()`](https://jmw86069.github.io/jamses/reference/SEStats.md),
[`contrast_names,SEStats-method`](https://jmw86069.github.io/jamses/reference/contrast_names-SEStats-method.md),
[`groups,SEStats-method`](https://jmw86069.github.io/jamses/reference/groups-SEStats-method.md),
[`print,SEStats-method`](https://jmw86069.github.io/jamses/reference/print-SEStats-method.md),
[`samples,SEStats-method`](https://jmw86069.github.io/jamses/reference/samples-SEStats-method.md)

## Examples

``` r
if (FALSE) {
  factors(sestats_object)
}
```
