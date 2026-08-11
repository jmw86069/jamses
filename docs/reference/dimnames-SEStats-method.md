# Get dimnames from SEStats

Extract the dimension names from an `SEStats` object, how
[`dimnames()`](https://rdrr.io/r/base/dimnames.html) works on
'hit_array' data.

## Usage

``` r
## S7 method for class <SEStats>
dimnames(x)
```

## Arguments

- x:

  An `SEStats` object

## Value

`list` with named elements: `'Cutoffs'`, `'Contrasts'`, `'Signal'`

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
  dimnames(sestats)
}
```
