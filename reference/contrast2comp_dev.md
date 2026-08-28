# Convert contrast to short-form comp, convert comp to contrast (DEV)

Convert contrast to short-form comp, convert comp to contrast (DEV)

## Usage

``` r
contrast2comp_dev(contrast_names, verbose = FALSE, ...)
```

## Details

This method is developmental and is being tested as potentially faster
version of
[`contrast2comp()`](https://jmw86069.github.io/jamses/reference/contrast2comp.md).

## See also

Other jamses utilities:
[`choose_annotation_colnames()`](https://jmw86069.github.io/jamses/reference/choose_annotation_colnames.md),
[`combine_sestats()`](https://jmw86069.github.io/jamses/reference/combine_sestats.md),
[`fold_to_log2fold()`](https://jmw86069.github.io/jamses/reference/fold_to_log2fold.md),
[`intercalate()`](https://jmw86069.github.io/jamses/reference/intercalate.md),
[`list2im_opt()`](https://jmw86069.github.io/jamses/reference/list2im_opt.md),
[`list_to_sestats()`](https://jmw86069.github.io/jamses/reference/list_to_sestats.md),
[`log2fold_to_fold()`](https://jmw86069.github.io/jamses/reference/log2fold_to_fold.md),
[`make_block_arrow_polygon()`](https://jmw86069.github.io/jamses/reference/make_block_arrow_polygon.md),
[`mark_stat_hits()`](https://jmw86069.github.io/jamses/reference/mark_stat_hits.md),
[`matrix_normalize()`](https://jmw86069.github.io/jamses/reference/matrix_normalize.md),
[`merge_statdf_all_test()`](https://jmw86069.github.io/jamses/reference/merge_statdf_all_test.md),
[`point_handedness()`](https://jmw86069.github.io/jamses/reference/point_handedness.md),
[`point_slope_intercept()`](https://jmw86069.github.io/jamses/reference/point_slope_intercept.md),
[`shortest_unique_abbreviation()`](https://jmw86069.github.io/jamses/reference/shortest_unique_abbreviation.md),
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
contrast_names1 <- c("(A_c-B_c)-(A_d-B_d)",
   "A_c-B_c",
   "GroupA-GroupB",
   "((A_c_J-B_c_J)-(A_d_J-B_d_J))-((A_c_W-B_c_W)-(A_d_W-B_d_W))")

contrast2comp_dev(contrast_names1)
#> [1] "A-B:c-d"       "A-B:c"         "GroupA-GroupB" "A-B:c-d:J-W"  

# Note: this function fails when contrasts are not balanced
contrast2comp_dev(c("(A_c-B_c)-(A_d-C_d)"), verbose=TRUE)
#> (17:55:19) 28Aug2026: contrast2comp_dev(): ilist:
#> $`(A_c-B_c)-(A_d-C_d)`
#> [1] ""    "A_c" "B_c" "A_d" "C_d"
#> 
#> (17:55:19) 28Aug2026: contrast2comp_dev(): idf:
#>                                        C                   B A.1 A.2 split1
#> (A_c-B_c)-(A_d-C_d)2 (A_c-B_c)-(A_d-C_d) (A_c-B_c)-(A_d-C_d)   A   c      1
#> (A_c-B_c)-(A_d-C_d)3 (A_c-B_c)-(A_d-C_d) (A_c-B_c)-(A_d-C_d)   B   c      1
#> (A_c-B_c)-(A_d-C_d)4 (A_c-B_c)-(A_d-C_d) (A_c-B_c)-(A_d-C_d)   A   d      2
#> (A_c-B_c)-(A_d-C_d)5 (A_c-B_c)-(A_d-C_d) (A_c-B_c)-(A_d-C_d)   C   d      2
#>                      idf_B_split
#> (A_c-B_c)-(A_d-C_d)2           A
#> (A_c-B_c)-(A_d-C_d)3           A
#> (A_c-B_c)-(A_d-C_d)4           B
#> (A_c-B_c)-(A_d-C_d)5           B
#> (17:55:19) 28Aug2026: contrast2comp_dev(): idf_new:
#>                     C                   B A.1 A.2 split2
#> A (A_c-B_c)-(A_d-C_d) (A_c-B_c)-(A_d-C_d) A-B   c      1
#> B (A_c-B_c)-(A_d-C_d) (A_c-B_c)-(A_d-C_d) A-C   d      1
#> (17:55:19) 28Aug2026: contrast2comp_dev(): idf_newer:
#>                                       C                   B     A.1 A.2
#> (A_c-B_c)-(A_d-C_d) (A_c-B_c)-(A_d-C_d) (A_c-B_c)-(A_d-C_d) A-B-A-C c-d
#> (17:55:19) 28Aug2026: contrast2comp_dev(): two-way imbalance:
#>   A.1   A.2 
#>  TRUE FALSE 
#> (17:55:19) 28Aug2026: contrast2comp_dev(): idf_newer:
#>                         C                      B A.1 A.2                newB
#> (A_c-B_c)-(A_d-C_d)_v1 v1 (A_c-B_c)-(A_d-C_d)_v1 A-B   c (A_c-B_c)-(A_d-C_d)
#> (A_c-B_c)-(A_d-C_d)_v2 v2 (A_c-B_c)-(A_d-C_d)_v2 A-C   d (A_c-B_c)-(A_d-C_d)
#> (17:55:19) 28Aug2026: contrast2comp_dev(): idf_newer:
#>                         C                      B A.1 A.2                newB
#> (A_c-B_c)-(A_d-C_d)_v1 v1 (A_c-B_c)-(A_d-C_d)_v1 A-B   c (A_c-B_c)-(A_d-C_d)
#> (A_c-B_c)-(A_d-C_d)_v2 v2 (A_c-B_c)-(A_d-C_d)_v2 A-C   d (A_c-B_c)-(A_d-C_d)
#>                         comp
#> (A_c-B_c)-(A_d-C_d)_v1 A-B:c
#> (A_c-B_c)-(A_d-C_d)_v2 A-C:d
#> [1] "(A-B:c)-(A-C:d)"
contrast2comp(c("(A_c-B_c)-(A_d-C_d)"), verbose=TRUE)
#> (17:55:19) 28Aug2026: im:
#>      [,1] [,2]
#> [1,] "A"  "c" 
#> [2,] "B"  "c" 
#> [3,] "A"  "d" 
#> [4,] "C"  "d" 
#> (17:55:19) 28Aug2026: imnew:
#>      [,1]  [,2]
#> [1,] "A-B" "c" 
#> [2,] "A-C" "d" 
#> [1] "(A-B:c)-(A-C:d)"

contrast2comp_dev(c(contrast_names1[1:2], "(A_c-B_c)-(A_d-C_d)"), verbose=TRUE)
#> (17:55:19) 28Aug2026: contrast2comp_dev(): ilist:
#> $`(A_c-B_c)-(A_d-B_d)`
#> [1] ""    "A_c" "B_c" "A_d" "B_d"
#> 
#> $`A_c-B_c`
#> [1] "A_c" "B_c"
#> 
#> $`(A_c-B_c)-(A_d-C_d)`
#> [1] ""    "A_c" "B_c" "A_d" "C_d"
#> 
#> (17:55:19) 28Aug2026: contrast2comp_dev(): idf:
#>                                        C                   B A.1 A.2 split1
#> (A_c-B_c)-(A_d-B_d)2 (A_c-B_c)-(A_d-B_d) (A_c-B_c)-(A_d-B_d)   A   c      1
#> (A_c-B_c)-(A_d-B_d)3 (A_c-B_c)-(A_d-B_d) (A_c-B_c)-(A_d-B_d)   B   c      1
#> (A_c-B_c)-(A_d-B_d)4 (A_c-B_c)-(A_d-B_d) (A_c-B_c)-(A_d-B_d)   A   d      2
#> (A_c-B_c)-(A_d-B_d)5 (A_c-B_c)-(A_d-B_d) (A_c-B_c)-(A_d-B_d)   B   d      2
#> A_c-B_c1                         A_c-B_c             A_c-B_c   A   c      3
#> A_c-B_c2                         A_c-B_c             A_c-B_c   B   c      3
#> (A_c-B_c)-(A_d-C_d)2 (A_c-B_c)-(A_d-C_d) (A_c-B_c)-(A_d-C_d)   A   c      4
#> (A_c-B_c)-(A_d-C_d)3 (A_c-B_c)-(A_d-C_d) (A_c-B_c)-(A_d-C_d)   B   c      4
#> (A_c-B_c)-(A_d-C_d)4 (A_c-B_c)-(A_d-C_d) (A_c-B_c)-(A_d-C_d)   A   d      5
#> (A_c-B_c)-(A_d-C_d)5 (A_c-B_c)-(A_d-C_d) (A_c-B_c)-(A_d-C_d)   C   d      5
#>                      idf_B_split
#> (A_c-B_c)-(A_d-B_d)2           A
#> (A_c-B_c)-(A_d-B_d)3           A
#> (A_c-B_c)-(A_d-B_d)4           B
#> (A_c-B_c)-(A_d-B_d)5           B
#> A_c-B_c1                       C
#> A_c-B_c2                       C
#> (A_c-B_c)-(A_d-C_d)2           D
#> (A_c-B_c)-(A_d-C_d)3           D
#> (A_c-B_c)-(A_d-C_d)4           E
#> (A_c-B_c)-(A_d-C_d)5           E
#> (17:55:19) 28Aug2026: contrast2comp_dev(): idf_new:
#>                     C                   B A.1 A.2 split2
#> A (A_c-B_c)-(A_d-B_d) (A_c-B_c)-(A_d-B_d) A-B   c      1
#> B (A_c-B_c)-(A_d-B_d) (A_c-B_c)-(A_d-B_d) A-B   d      1
#> C             A_c-B_c             A_c-B_c A-B   c      2
#> D (A_c-B_c)-(A_d-C_d) (A_c-B_c)-(A_d-C_d) A-B   c      3
#> E (A_c-B_c)-(A_d-C_d) (A_c-B_c)-(A_d-C_d) A-C   d      3
#> (17:55:19) 28Aug2026: contrast2comp_dev(): idf_newer:
#>                                       C                   B     A.1 A.2
#> (A_c-B_c)-(A_d-B_d) (A_c-B_c)-(A_d-B_d) (A_c-B_c)-(A_d-B_d)     A-B c-d
#> A_c-B_c                         A_c-B_c             A_c-B_c     A-B   c
#> (A_c-B_c)-(A_d-C_d) (A_c-B_c)-(A_d-C_d) (A_c-B_c)-(A_d-C_d) A-B-A-C c-d
#> (17:55:19) 28Aug2026: contrast2comp_dev(): two-way imbalance:
#>        A.1   A.2
#> [1,] FALSE FALSE
#> [2,] FALSE FALSE
#> [3,]  TRUE FALSE
#> (17:55:19) 28Aug2026: contrast2comp_dev(): idf_newer:
#>                         C                      B A.1 A.2                newB
#> (A_c-B_c)-(A_d-B_d)          (A_c-B_c)-(A_d-B_d) A-B c-d (A_c-B_c)-(A_d-B_d)
#> A_c-B_c                                  A_c-B_c A-B   c             A_c-B_c
#> (A_c-B_c)-(A_d-C_d)_v1 v1 (A_c-B_c)-(A_d-C_d)_v1 A-B   c (A_c-B_c)-(A_d-C_d)
#> (A_c-B_c)-(A_d-C_d)_v2 v2 (A_c-B_c)-(A_d-C_d)_v2 A-C   d (A_c-B_c)-(A_d-C_d)
#> (17:55:19) 28Aug2026: contrast2comp_dev(): idf_newer:
#>                         C                      B A.1 A.2                newB
#> (A_c-B_c)-(A_d-B_d)          (A_c-B_c)-(A_d-B_d) A-B c-d (A_c-B_c)-(A_d-B_d)
#> A_c-B_c                                  A_c-B_c A-B   c             A_c-B_c
#> (A_c-B_c)-(A_d-C_d)_v1 v1 (A_c-B_c)-(A_d-C_d)_v1 A-B   c (A_c-B_c)-(A_d-C_d)
#> (A_c-B_c)-(A_d-C_d)_v2 v2 (A_c-B_c)-(A_d-C_d)_v2 A-C   d (A_c-B_c)-(A_d-C_d)
#>                           comp
#> (A_c-B_c)-(A_d-B_d)    A-B:c-d
#> A_c-B_c                  A-B:c
#> (A_c-B_c)-(A_d-C_d)_v1   A-B:c
#> (A_c-B_c)-(A_d-C_d)_v2   A-C:d
#> [1] "A-B:c-d"         "A-B:c"           "(A-B:c)-(A-C:d)"
```
