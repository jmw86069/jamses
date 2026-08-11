# Convert contrast to short-form comp, convert comp to contrast (DEV)

Convert contrast to short-form comp, convert comp to contrast (DEV)

## Usage

``` r
contrast2comp_dev(contrast_names, verbose = FALSE, ...)
```

## Details

This method is developmental and is being tested as potentially faster
version of [`contrast2comp()`](contrast2comp.md).

## See also

Other jamses utilities:
[`choose_annotation_colnames()`](choose_annotation_colnames.md),
[`combine_sestats()`](combine_sestats.md),
[`fold_to_log2fold()`](fold_to_log2fold.md),
[`intercalate()`](intercalate.md), [`list2im_opt()`](list2im_opt.md),
[`list_to_sestats()`](list_to_sestats.md),
[`log2fold_to_fold()`](log2fold_to_fold.md),
[`make_block_arrow_polygon()`](make_block_arrow_polygon.md),
[`mark_stat_hits()`](mark_stat_hits.md),
[`matrix_normalize()`](matrix_normalize.md),
[`merge_statdf_all_test()`](merge_statdf_all_test.md),
[`point_handedness()`](point_handedness.md),
[`point_slope_intercept()`](point_slope_intercept.md),
[`shortest_unique_abbreviation()`](shortest_unique_abbreviation.md),
[`shrinkDataFrame()`](shrinkDataFrame.md),
[`shrink_df()`](shrink_df.md), [`shrink_matrix()`](shrink_matrix.md),
[`sort_samples()`](sort_samples.md),
[`strsplitOrdered()`](strsplitOrdered.md),
[`sub_split_vector()`](sub_split_vector.md),
[`update_function_params()`](update_function_params.md),
[`update_list_elements()`](update_list_elements.md)

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
#> (16:00:42) 10Aug2026: contrast2comp_dev(): ilist:
#> $`(A_c-B_c)-(A_d-C_d)`
#> [1] ""    "A_c" "B_c" "A_d" "C_d"
#> 
#> (16:00:42) 10Aug2026: contrast2comp_dev(): idf:
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
#> (16:00:42) 10Aug2026: contrast2comp_dev(): idf_new:
#>                     C                   B A.1 A.2 split2
#> A (A_c-B_c)-(A_d-C_d) (A_c-B_c)-(A_d-C_d) A-B   c      1
#> B (A_c-B_c)-(A_d-C_d) (A_c-B_c)-(A_d-C_d) A-C   d      1
#> (16:00:42) 10Aug2026: contrast2comp_dev(): idf_newer:
#>                                       C                   B     A.1 A.2
#> (A_c-B_c)-(A_d-C_d) (A_c-B_c)-(A_d-C_d) (A_c-B_c)-(A_d-C_d) A-B-A-C c-d
#> (16:00:42) 10Aug2026: contrast2comp_dev(): two-way imbalance:
#>   A.1   A.2 
#>  TRUE FALSE 
#> (16:00:42) 10Aug2026: contrast2comp_dev(): idf_newer:
#>                         C                      B A.1 A.2                newB
#> (A_c-B_c)-(A_d-C_d)_v1 v1 (A_c-B_c)-(A_d-C_d)_v1 A-B   c (A_c-B_c)-(A_d-C_d)
#> (A_c-B_c)-(A_d-C_d)_v2 v2 (A_c-B_c)-(A_d-C_d)_v2 A-C   d (A_c-B_c)-(A_d-C_d)
#> (16:00:42) 10Aug2026: contrast2comp_dev(): idf_newer:
#>                         C                      B A.1 A.2                newB
#> (A_c-B_c)-(A_d-C_d)_v1 v1 (A_c-B_c)-(A_d-C_d)_v1 A-B   c (A_c-B_c)-(A_d-C_d)
#> (A_c-B_c)-(A_d-C_d)_v2 v2 (A_c-B_c)-(A_d-C_d)_v2 A-C   d (A_c-B_c)-(A_d-C_d)
#>                         comp
#> (A_c-B_c)-(A_d-C_d)_v1 A-B:c
#> (A_c-B_c)-(A_d-C_d)_v2 A-C:d
#> [1] "(A-B:c)-(A-C:d)"
contrast2comp(c("(A_c-B_c)-(A_d-C_d)"), verbose=TRUE)
#> (16:00:42) 10Aug2026: im:
#>      [,1] [,2]
#> [1,] "A"  "c" 
#> [2,] "B"  "c" 
#> [3,] "A"  "d" 
#> [4,] "C"  "d" 
#> (16:00:43) 10Aug2026: imnew:
#>      [,1]  [,2]
#> [1,] "A-B" "c" 
#> [2,] "A-C" "d" 
#> [1] "(A-B:c)-(A-C:d)"

contrast2comp_dev(c(contrast_names1[1:2], "(A_c-B_c)-(A_d-C_d)"), verbose=TRUE)
#> (16:00:43) 10Aug2026: contrast2comp_dev(): ilist:
#> $`(A_c-B_c)-(A_d-B_d)`
#> [1] ""    "A_c" "B_c" "A_d" "B_d"
#> 
#> $`A_c-B_c`
#> [1] "A_c" "B_c"
#> 
#> $`(A_c-B_c)-(A_d-C_d)`
#> [1] ""    "A_c" "B_c" "A_d" "C_d"
#> 
#> (16:00:43) 10Aug2026: contrast2comp_dev(): idf:
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
#> (16:00:43) 10Aug2026: contrast2comp_dev(): idf_new:
#>                     C                   B A.1 A.2 split2
#> A (A_c-B_c)-(A_d-B_d) (A_c-B_c)-(A_d-B_d) A-B   c      1
#> B (A_c-B_c)-(A_d-B_d) (A_c-B_c)-(A_d-B_d) A-B   d      1
#> C             A_c-B_c             A_c-B_c A-B   c      2
#> D (A_c-B_c)-(A_d-C_d) (A_c-B_c)-(A_d-C_d) A-B   c      3
#> E (A_c-B_c)-(A_d-C_d) (A_c-B_c)-(A_d-C_d) A-C   d      3
#> (16:00:43) 10Aug2026: contrast2comp_dev(): idf_newer:
#>                                       C                   B     A.1 A.2
#> (A_c-B_c)-(A_d-B_d) (A_c-B_c)-(A_d-B_d) (A_c-B_c)-(A_d-B_d)     A-B c-d
#> A_c-B_c                         A_c-B_c             A_c-B_c     A-B   c
#> (A_c-B_c)-(A_d-C_d) (A_c-B_c)-(A_d-C_d) (A_c-B_c)-(A_d-C_d) A-B-A-C c-d
#> (16:00:43) 10Aug2026: contrast2comp_dev(): two-way imbalance:
#>        A.1   A.2
#> [1,] FALSE FALSE
#> [2,] FALSE FALSE
#> [3,]  TRUE FALSE
#> (16:00:43) 10Aug2026: contrast2comp_dev(): idf_newer:
#>                         C                      B A.1 A.2                newB
#> (A_c-B_c)-(A_d-B_d)          (A_c-B_c)-(A_d-B_d) A-B c-d (A_c-B_c)-(A_d-B_d)
#> A_c-B_c                                  A_c-B_c A-B   c             A_c-B_c
#> (A_c-B_c)-(A_d-C_d)_v1 v1 (A_c-B_c)-(A_d-C_d)_v1 A-B   c (A_c-B_c)-(A_d-C_d)
#> (A_c-B_c)-(A_d-C_d)_v2 v2 (A_c-B_c)-(A_d-C_d)_v2 A-C   d (A_c-B_c)-(A_d-C_d)
#> (16:00:43) 10Aug2026: contrast2comp_dev(): idf_newer:
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
