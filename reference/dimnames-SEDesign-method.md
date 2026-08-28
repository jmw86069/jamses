# Get dimnames from SEDesign

Extract the dimension names (samples, groups, contrast_names) from an
SEDesign object, similar to how
[`dimnames()`](https://rdrr.io/r/base/dimnames.html) works on arrays or
matrices.

## Usage

``` r
## S7 method for class <SEDesign>
dimnames(x)
```

## Arguments

- x:

  An `SEDesign` object

## Value

`list` with named elements: `$samples`, `$groups`, `$contrasts`

## Details

Returns a `list` with three named elements:

- `samples`: Sample names (from `samples(x)`)

- `groups`: Group names (from `groups(x)`)

- `contrasts`: Contrast names (from `contrast_names(x)`)

This provides a unified way to access all dimension names from a
SEDesign object, consistent with base R's
[`dimnames()`](https://rdrr.io/r/base/dimnames.html) function for arrays
and matrices.

## Examples

``` r
if (FALSE) {
  dimnames(sedesign)
}
```
