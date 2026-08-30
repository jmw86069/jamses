# Find the shortest abbrevation to retain unique values

Find the shortest abbrevation to retain unique values

## Usage

``` r
shortest_unique_abbreviation(
  x,
  retain_contig_numbers = TRUE,
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  `character` vector

- retain_contig_numbers:

  `logical`, default `TRUE`, whether numbers at the end of an
  abbreviated string should remain contiguous.

  - When `TRUE`, the goal is not to split a numeric value in the middle
    of the number.

  - When `FALSE` the string will be abbreviated at the first position of
    uniqueness.

- ...:

  additional arguments are ignored.

## Value

`character` vector named using unique values in `x`, and whose values
are the shortest abbreviated substrings which maintain consistent
uniqueness.

## Details

This function is intended to abbreviate factor levels used in
statistical contrasts to the smallest substring that uniquely represents
the unique entries provided in `x`.

For example, `c("one", "two", "three", "four")` would be converted to
`c("on", "tw", "th", "fo")`.

The default `retain_contig_numbers=TRUE` will attempt to retain numeric
values at the end of a string, to avoid splitting the number at an
intermediate position. This option only applies when the character
substring is not already unique before encountering the numeric
substring.

    * For this input:

`c("a", "p6", "p12", "p21")` the output keeps the contiguous numbers
together: `c("a", "p6", "p12", "p21")`

- For this input: `c("a", "b6", "c12", "d21")` only the first character
  is retained, because it is already unique: `c("a", "b", "c", "d")`

## Todo

- Consider some method to retain contiguous numbers at the end of a long
  string, while abbreviating the long string.

  - For this input:
    `c("adult", "prenatal6", "prenatal12", "prenatal21")` the ideal
    output would be: `c("a", "p6", "p12", "p21")`

  - To be fair, I do not know how to describe this logic. It may
    required breaking into words by character/non-character breakpoints,
    then applying substring to each?

## See also

Other jamses utilities:
[`choose_annotation_colnames()`](https://jmw86069.github.io/jamses/reference/choose_annotation_colnames.md),
[`combine_sestats()`](https://jmw86069.github.io/jamses/reference/combine_sestats.md),
[`contrast2comp_dev()`](https://jmw86069.github.io/jamses/reference/contrast2comp_dev.md),
[`fold_to_log2fold()`](https://jmw86069.github.io/jamses/reference/fold_to_log2fold.md),
[`intercalate()`](https://jmw86069.github.io/jamses/reference/intercalate.md),
[`list2im_opt()`](https://jmw86069.github.io/jamses/reference/list2im_opt.md),
[`list2im_value_internal()`](https://jmw86069.github.io/jamses/reference/list2im_value_internal.md),
[`list_to_sestats()`](https://jmw86069.github.io/jamses/reference/list_to_sestats.md),
[`log2fold_to_fold()`](https://jmw86069.github.io/jamses/reference/log2fold_to_fold.md),
[`make_block_arrow_polygon()`](https://jmw86069.github.io/jamses/reference/make_block_arrow_polygon.md),
[`mark_stat_hits()`](https://jmw86069.github.io/jamses/reference/mark_stat_hits.md),
[`matrix_normalize()`](https://jmw86069.github.io/jamses/reference/matrix_normalize.md),
[`merge_statdf_all_test()`](https://jmw86069.github.io/jamses/reference/merge_statdf_all_test.md),
[`point_handedness()`](https://jmw86069.github.io/jamses/reference/point_handedness.md),
[`point_slope_intercept()`](https://jmw86069.github.io/jamses/reference/point_slope_intercept.md),
[`shrinkDataFrame()`](https://jmw86069.github.io/jamses/reference/shrinkDataFrame.md),
[`shrink_df()`](https://jmw86069.github.io/jamses/reference/shrink_df.md),
[`shrink_matrix()`](https://jmw86069.github.io/jamses/reference/shrink_matrix.md),
[`sort_samples()`](https://jmw86069.github.io/jamses/reference/sort_samples.md),
[`strsplitOrdered()`](https://jmw86069.github.io/jamses/reference/strsplitOrdered.md),
[`sub_split_vector()`](https://jmw86069.github.io/jamses/reference/sub_split_vector.md),
[`update_function_params()`](https://jmw86069.github.io/jamses/reference/update_function_params.md),
[`update_list_elements()`](https://jmw86069.github.io/jamses/reference/update_list_elements.md)

## Examples

``` r
x <- c("a", "p6", "p12", "p21");
shortest_unique_abbreviation(x)
#>     a    p6   p12   p21 
#>   "a"  "p6" "p12" "p21" 

shortest_unique_abbreviation(x, retain_contig_numbers=TRUE)
#>     a    p6   p12   p21 
#>   "a"  "p6" "p12" "p21" 

x1 <- c("male", "female");
shortest_unique_abbreviation(x1)
#>   male female 
#>    "m"    "f" 

x2 <- c("Control", "Nicotine");
shortest_unique_abbreviation(x2)
#>  Control Nicotine 
#>      "C"      "N" 

x3 <- c("Control", "Nicotine10", "Nicotine12", "Nicotine20");
shortest_unique_abbreviation(x3)
#>      Control   Nicotine10   Nicotine12   Nicotine20 
#>    "Control" "Nicotine10" "Nicotine12" "Nicotine20" 

x4 <- c("one", "two", "three", "four");
shortest_unique_abbreviation(x4)
#>   one   two three  four 
#>  "on"  "tw"  "th"  "fo" 
```
