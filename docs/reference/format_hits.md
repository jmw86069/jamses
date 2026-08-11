# Format list of hit vectors into summary counts

Format list of hit vectors into summary counts

## Usage

``` r
format_hits(hits, style = c("text", "integer", "vector"), ...)
```

## Arguments

- hits:

  `list`, `array`, or `vector` with values `c(-1, 1)` indicating hits
  down, or hits up, respectively. When any vector contains only `NA`
  values, then `NA` is also returned. This distinction is as follows:

  - `NA` values indicate the statistical contrast was not performed,
    which happens when [`se_contrast_stats()`](se_contrast_stats.md)
    arguments for interaction contrasts differ from pairwise contrasts,
    for example `int_adjp_cutoff`, `int_p_cutoff`, or `int_fold_cutoff`.

  - `length==0` indicates the statistical contrast was performed, and
    there were no statistical hits for the given cutoffs.

- style:

  `character` string indicating the output format:

  - `"text"`: `character` string with `"hits(hits up, hits down)"`.

  - `"integer"`: `integer` number of hits, or `NA` when the test was not
    performed.

  - `"vector"`: `integer` vector with names `("hit", "up", "down")`.

- ...:

  additional arguments are ignored.

## Value

the same data type as input, where the hit vector is replaced with a
single value summarizing the hits. The data types have three expected
options:

- `array`: when used with `hit_array` data from
  [`se_contrast_stats()`](se_contrast_stats.md)

- `list`: when used on a particular subset of `hit_array` data

- `numeric`: when used with a single contrast

## Details

This function is used by [`sestats_to_df()`](sestats_to_df.md), to
summarize hits for each contrast into one of these formats:

- string summary of statistical hits with `"hits,up,down"` when
  `style="text"`, for example: `"623 hits (267 up, 379 down)"`

- integer count of statistical hits when `style="integer"`, for example:
  `623`.

- `vector` of integer counts for `c("hits", "up", "down")`, for example:
  `c(hits=623, up=267, down=379)`.

The function may be useful outside of
[`sestats_to_df()`](sestats_to_df.md) so it is exported as a convenience
function.

## See also

Other jamses stats: [`ebayes2dfs()`](ebayes2dfs.md),
[`handle_na_values()`](handle_na_values.md),
[`hit_array_to_list()`](hit_array_to_list.md),
[`process_sestats_to_hitim()`](process_sestats_to_hitim.md),
[`run_limma_replicate()`](run_limma_replicate.md),
[`save_sestats()`](save_sestats.md),
[`se_contrast_stats()`](se_contrast_stats.md),
[`sestats_to_df()`](sestats_to_df.md),
[`sestats_to_dfs()`](sestats_to_dfs.md), [`voom_jam()`](voom_jam.md)

## Examples

``` r
set.seed(123)
hitlist <- list(
   `groupA-groupB`=sample(c(-1, 1), size=25, replace=TRUE),
   `groupA-groupC`=sample(c(-1, 1), size=50, replace=TRUE))
format_hits(hitlist, style="text")
#> $`groupA-groupB`
#> [1] "25 hits (10 up, 15 down)"
#> 
#> $`groupA-groupC`
#> [1] "50 hits (18 up, 32 down)"
#> 
format_hits(hitlist, style="integer")
#> $`groupA-groupB`
#> [1] 25
#> 
#> $`groupA-groupC`
#> [1] 50
#> 
```
