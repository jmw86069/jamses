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

Other SEStats objects: [`.class_list`](dot-class_list.md),
[`contrast_names,SEStats-method`](contrast_names-SEStats-method.md),
[`groups,SEStats-method`](groups-SEStats-method.md),
[`print,SEStats-method`](print-SEStats-method.md),
[`samples,SEStats-method`](samples-SEStats-method.md)

## Examples

``` r
if (FALSE) {
  factors(sestats_object)
}
```
