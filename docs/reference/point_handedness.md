# Determine which side one point is to another, given a slope or angle

Determine which side one point is to another, given a slope or angle

## Usage

``` r
point_handedness(
  pt1,
  pt2,
  slope = NULL,
  angle = NULL,
  do_plot = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- pt1:

  `numeric` matrix of 2 columns, with x and y coordinates.

- pt2:

  `numeric` matrix of 2 columns, with x and y coordinates.

- slope:

  `numeric` slope for each point in pt1 and pt2.

- do_plot:

  `logical` indicating whether to plot the result.

- ...:

  additional arguments are ignored.

## Value

`character` vector equal to the number of points, `nrow(pt1)`:

- `"right"` indicates `pt1` is on the right side of `pt2`

- `"left"` indicates `pt1` is on the left side of `pt2`

## Details

The result describes the position of the first line relative to the
second line, assuming both lines are parallel with identical slope. For
example "left" indicates that line 1 is on the left side of line 2. Note
that the `angle` is a more accurate measure of directionality, otherwise
`slope` is always assumed to desribe an angle moving to the right.

When both lines are exactly overlapping, the result may be unstable,
however the result tends to favor "right" by default.

## See also

Other jamses utilities:
[`choose_annotation_colnames()`](choose_annotation_colnames.md),
[`combine_sestats()`](combine_sestats.md),
[`contrast2comp_dev()`](contrast2comp_dev.md),
[`fold_to_log2fold()`](fold_to_log2fold.md),
[`intercalate()`](intercalate.md), [`list2im_opt()`](list2im_opt.md),
[`list_to_sestats()`](list_to_sestats.md),
[`log2fold_to_fold()`](log2fold_to_fold.md),
[`make_block_arrow_polygon()`](make_block_arrow_polygon.md),
[`mark_stat_hits()`](mark_stat_hits.md),
[`matrix_normalize()`](matrix_normalize.md),
[`merge_statdf_all_test()`](merge_statdf_all_test.md),
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
pt1 <- matrix(ncol=2, c(1, 1))
pt2 <- matrix(ncol=2, c(2, 2))
point_handedness(pt1, pt2, angle=0, do_plot=TRUE)

#>      [,1] [,2]
#> [1,]    1    1
#> [1] ""
#> [1] "right"
point_handedness(pt1=c(1, 1), pt2=c(2, 2), angle=0, do_plot=TRUE)
#>      [,1] [,2]
#> [1,]    1    1
#> [1] ""
#> [1] "right"

pt1 <- matrix(ncol=2, c(1, 1))
pt2 <- matrix(ncol=2, c(0, 1))

point_handedness(pt1, pt2, angle=45, do_plot=TRUE)

#>      [,1] [,2]
#> [1,]    1    1
#> [1] ""
#> [1] "right"
point_handedness(pt2, pt1, angle=45, do_plot=TRUE)

#>      [,1] [,2]
#> [1,]    0    1
#> [1] ""
#> [1] "left"

point_handedness(rbind(pt2, pt2-1), rbind(pt1, pt1-1),
   angle=c(45, 45+180), do_plot=TRUE)

#>      [,1] [,2]
#> [1,]    0    1
#> [2,]   -1    0
#> [1] "A:" "B:"
#> [1] "left"  "right"

point_handedness(pt1, pt2, angle=45, do_plot=TRUE)
#>      [,1] [,2]
#> [1,]    1    1
#> [1] ""
#> [1] "right"
title(main="angle = 45,\n(slope = 1)")

point_handedness(pt1, pt2, angle=45 + 180, do_plot=TRUE)
#>      [,1] [,2]
#> [1,]    1    1
#> [1] ""
#> [1] "left"
title(main="angle = 225,\n(slope = 1)")


point_handedness(pt1, pt2, slope=Inf, do_plot=TRUE)
#>      [,1] [,2]
#> [1,]    1    1
#> [1] ""
#> [1] "right"
title(main="slope = Inf")

point_handedness(pt1, pt2, angle=90, do_plot=TRUE)
#>      [,1] [,2]
#> [1,]    1    1
#> [1] ""
#> [1] "right"
title(main="angle = 90")

point_handedness(pt2, pt1, slope=-Inf, do_plot=TRUE)
#>      [,1] [,2]
#> [1,]    0    1
#> [1] ""
#> [1] "right"
title(main="slope = -Inf")

point_handedness(pt2, pt1, angle = 270, do_plot=TRUE)
#>      [,1] [,2]
#> [1,]    0    1
#> [1] ""
#> [1] "right"
title(main="angle = 270")


point_handedness(pt1, pt2, slope=0, do_plot=TRUE)

#>      [,1] [,2]
#> [1,]    1    1
#> [1] ""
#> [1] "right"
point_handedness(pt2, pt1, slope=0, do_plot=TRUE)

#>      [,1] [,2]
#> [1,]    0    1
#> [1] ""
#> [1] "right"

pt1 <- matrix(ncol=2, c(-2, 5))
pt2 <- matrix(ncol=2, c(2, 3))
point_handedness(pt1, pt2, slope=1, do_plot=TRUE)

#>      [,1] [,2]
#> [1,]   -2    5
#> [1] ""
#> [1] "left"
point_handedness(pt1, pt2, slope=-1, do_plot=TRUE)

#>      [,1] [,2]
#> [1,]   -2    5
#> [1] ""
#> [1] "right"
point_handedness(pt1, pt2, slope=-1/3, do_plot=TRUE)

#>      [,1] [,2]
#> [1,]   -2    5
#> [1] ""
#> [1] "left"
point_handedness(pt1, pt2, slope=-1/3, do_plot=TRUE)
#>      [,1] [,2]
#> [1,]   -2    5
#> [1] ""
#> [1] "left"

point_handedness(pt1, pt2, slope=Inf, do_plot=TRUE)

#>      [,1] [,2]
#> [1,]   -2    5
#> [1] ""
#> [1] "left"
point_handedness(pt1, pt2, slope=-Inf, do_plot=TRUE)

#>      [,1] [,2]
#> [1,]   -2    5
#> [1] ""
#> [1] "right"
```
