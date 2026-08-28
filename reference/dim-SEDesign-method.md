# Get dimensions for SEDesign

Summarize the dimensions using samples, groups, contrast_names from an
`SEDesign` object.

## Usage

``` r
## S7 method for class <SEDesign>
dim(x)
```

## Arguments

- x:

  An `SEDesign` object

## Value

`integer` vector with lengths of `samples`, `groups`, `contrasts`

## Details

Returns a `list` with three named elements:

- `samples`: Sample names (from `samples(x)`)

- `groups`: Group names (from `groups(x)`)

- `contrasts`: Contrast names (from `contrast_names(x)`)

This method provides a unified way to access all dimensions for an
`SEDesign` object.

## Examples

``` r
if (FALSE) {
  dim(sedesign)
}
```
