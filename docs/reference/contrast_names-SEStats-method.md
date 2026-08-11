# Get contrast names from SEStats

Extract contrast names from the SEStats object's sedesign.

## Usage

``` r
## S7 method for class <SEStats>
contrast_names(object)
```

## Arguments

- x:

  An `SEStats` object

## Value

Character vector of contrast names

## Details

This is a read-only accessor that retrieves the contrast names from the
SEDesign object stored in the SEStats object. Values cannot be modified
through this accessor, as the design is fixed when results are created.

## See also

Other SEStats objects: [`.class_list`](dot-class_list.md),
[`factors,SEStats-method`](factors-SEStats-method.md),
[`groups,SEStats-method`](groups-SEStats-method.md),
[`print,SEStats-method`](print-SEStats-method.md),
[`samples,SEStats-method`](samples-SEStats-method.md)

## Examples

``` r
if (FALSE) {
  contrast_names(sestats_object)
}
```
