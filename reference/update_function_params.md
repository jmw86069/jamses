# Update function default parameters

Update function default parameters

## Usage

``` r
update_function_params(
  function_name = NULL,
  param_name = NULL,
  new_values = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- function_name:

  `character` function name, passed to
  [`formals()`](https://rdrr.io/r/base/formals.html)

- param_name:

  `character` string of one parameter / argument name.

- new_values:

  any relevant R object type to use for the function parameter /
  argument.

- verbose:

  `logical` whether to print verbose output

- ...:

  additional arguments are ignored.

## Value

a `list` with class matching the R function argument.

## Details

This function is a minor extension to
[`update_list_elements()`](https://jmw86069.github.io/jamses/reference/update_list_elements.md)
intended to help update function parameters which are defined as a
nested `list`. It will update internal components of the nested `list`
without replacing the whole list.

It is very similar to
[`modifyList()`](https://rdrr.io/r/utils/modifyList.html) except for a
nested list.

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
[`shortest_unique_abbreviation()`](https://jmw86069.github.io/jamses/reference/shortest_unique_abbreviation.md),
[`shrinkDataFrame()`](https://jmw86069.github.io/jamses/reference/shrinkDataFrame.md),
[`shrink_df()`](https://jmw86069.github.io/jamses/reference/shrink_df.md),
[`shrink_matrix()`](https://jmw86069.github.io/jamses/reference/shrink_matrix.md),
[`sort_samples()`](https://jmw86069.github.io/jamses/reference/sort_samples.md),
[`strsplitOrdered()`](https://jmw86069.github.io/jamses/reference/strsplitOrdered.md),
[`sub_split_vector()`](https://jmw86069.github.io/jamses/reference/sub_split_vector.md),
[`update_list_elements()`](https://jmw86069.github.io/jamses/reference/update_list_elements.md)

## Examples

``` r
params1 <- list(
   jammanorm=list(
      controlGenes=c("HK1", "HK2")
   ),
   cyclicLoess=list(method="affy"),
   quantile=list(ties=FALSE)
)

# 1. Not correct: update arguments with modifyList()
fn_defaults <- eval(formals("matrix_normalize"))
new_fn_args <- modifyList(fn_defaults, list(params=params1))
str(new_fn_args$params)
#> List of 3
#>  $ jammanorm  :List of 1
#>   ..$ controlGenes: chr [1:2] "HK1" "HK2"
#>  $ cyclicLoess:List of 1
#>   ..$ method: chr "affy"
#>  $ quantile   :List of 1
#>   ..$ ties: logi FALSE
str(new_fn_args$params$jammanorm)
#> List of 1
#>  $ controlGenes: chr [1:2] "HK1" "HK2"

# 2. update only params argument with modifyList()
new_fn_params <- modifyList(eval(fn_defaults$params), params1)
str(new_fn_params)
#> List of 6
#>  $ cyclicLoess       :List of 5
#>   ..$ weights      : NULL
#>   ..$ span         : num 0.7
#>   ..$ adaptive.span: logi FALSE
#>   ..$ iterations   : num 3
#>   ..$ method       : chr "affy"
#>  $ quantile          :List of 1
#>   ..$ ties: logi FALSE
#>  $ jammanorm         :List of 8
#>   ..$ controlGenes     : chr [1:2] "HK1" "HK2"
#>   ..$ minimum_mean     : num 0
#>   ..$ maximum_mean     : num Inf
#>   ..$ controlSamples   : NULL
#>   ..$ centerGroups     : NULL
#>   ..$ useMedian        : logi FALSE
#>   ..$ noise_floor      : NULL
#>   ..$ noise_floor_value: NULL
#>  $ limma_batch_adjust:List of 2
#>   ..$ batch: NULL
#>   ..$ group: NULL
#>  $ TMM               :List of 5
#>   ..$ refColumn   : NULL
#>   ..$ logratioTrim: num 0.3
#>   ..$ sumTrim     : num 0.05
#>   ..$ doWeighting : logi TRUE
#>   ..$ Acutoff     : NULL
#>  $ TMMwsp            :List of 5
#>   ..$ refColumn   : NULL
#>   ..$ logratioTrim: num 0.3
#>   ..$ sumTrim     : num 0.05
#>   ..$ doWeighting : logi TRUE
#>   ..$ Acutoff     : NULL
str(new_fn_params$jammanorm)
#> List of 8
#>  $ controlGenes     : chr [1:2] "HK1" "HK2"
#>  $ minimum_mean     : num 0
#>  $ maximum_mean     : num Inf
#>  $ controlSamples   : NULL
#>  $ centerGroups     : NULL
#>  $ useMedian        : logi FALSE
#>  $ noise_floor      : NULL
#>  $ noise_floor_value: NULL

# 3. use update_function_params()
params2 <- update_function_params(
   function_name="matrix_normalize",
   param_name="params",
   new_values=params1)
str(params2)
#> List of 6
#>  $ cyclicLoess       :List of 5
#>   ..$ weights      : NULL
#>   ..$ span         : num 0.7
#>   ..$ adaptive.span: logi FALSE
#>   ..$ iterations   : num 3
#>   ..$ method       : chr "affy"
#>  $ quantile          :List of 1
#>   ..$ ties: logi FALSE
#>  $ jammanorm         :List of 8
#>   ..$ controlGenes     : chr [1:2] "HK1" "HK2"
#>   ..$ minimum_mean     : num 0
#>   ..$ maximum_mean     : num Inf
#>   ..$ controlSamples   : NULL
#>   ..$ centerGroups     : NULL
#>   ..$ useMedian        : logi FALSE
#>   ..$ noise_floor      : NULL
#>   ..$ noise_floor_value: NULL
#>  $ limma_batch_adjust:List of 2
#>   ..$ batch: NULL
#>   ..$ group: NULL
#>  $ TMM               :List of 5
#>   ..$ refColumn   : NULL
#>   ..$ logratioTrim: num 0.3
#>   ..$ sumTrim     : num 0.05
#>   ..$ doWeighting : logi TRUE
#>   ..$ Acutoff     : NULL
#>  $ TMMwsp            :List of 5
#>   ..$ refColumn   : NULL
#>   ..$ logratioTrim: num 0.3
#>   ..$ sumTrim     : num 0.05
#>   ..$ doWeighting : logi TRUE
#>   ..$ Acutoff     : NULL
str(params2$jammanorm)
#> List of 8
#>  $ controlGenes     : chr [1:2] "HK1" "HK2"
#>  $ minimum_mean     : num 0
#>  $ maximum_mean     : num Inf
#>  $ controlSamples   : NULL
#>  $ centerGroups     : NULL
#>  $ useMedian        : logi FALSE
#>  $ noise_floor      : NULL
#>  $ noise_floor_value: NULL
```
