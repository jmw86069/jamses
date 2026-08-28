# Stack a vector of contrasts into a full string

Stack a vector of contrasts into a full string

## Usage

``` r
stack_contrasts(imgroups, contrast_delim = "-", ...)
```

## Arguments

- imgroups:

  `character` vector of group names.

- contrast_delim:

  `character` string used as a delimiter between group names to indicate
  a statistical contrast.

- ...:

  additional arguments are ignored.

## Value

`character` vector with fully written contrasts.

## Details

This is an internal function called by
[`contrast2comp()`](https://jmw86069.github.io/jamses/reference/contrast2comp.md)
and
[`comp2contrast()`](https://jmw86069.github.io/jamses/reference/contrast2comp.md).
It takes a single vector of groups, with 2^(integer) number of entries,
and combines them into a series of pairwise contrasts.

## Examples

``` r
stack_contrasts(head(LETTERS, 2))
#> [1] "A-B"
stack_contrasts(head(LETTERS, 4))
#>             A 
#> "(A-B)-(C-D)" 
stack_contrasts(head(LETTERS, 8))
#>                             A 
#> "((A-B)-(C-D))-((E-F)-(G-H))" 
contrast2comp(stack_contrasts(head(LETTERS, 16)))
#>                                                             A 
#> "(((A-B)-(C-D))-((E-F)-(G-H)))-(((I-J)-(K-L))-((M-N)-(O-P)))" 

# more typical case
stacked <- stack_contrasts(c("KO_Treated", "KO_Control",
   "WT_Treated", "WT_Control"))
stacked
#>                                                 A 
#> "(KO_Treated-KO_Control)-(WT_Treated-WT_Control)" 

# this contrast can be rewritten as a "comp"
stacked_comp <- contrast2comp(stacked)
stacked_comp
#>                       A 
#> "KO-WT:Treated-Control" 

# of course the "comp" can be converted back to full contrast name
comp2contrast(stacked_comp)
#>                                                 A 
#> "(KO_Treated-WT_Treated)-(KO_Control-WT_Control)" 

# note when the contrast cannot be stacked, an error is thrown
tryCatch({
   stack_contrasts(head(LETTERS, 5))
}, error=function(e){
   print(e)
})
#> [1] "A" "B" "C" "D" "E"
#> <simpleError in stack_contrasts(head(LETTERS, 5)): length(imgroups) must be a value in 2^(int)>
```
