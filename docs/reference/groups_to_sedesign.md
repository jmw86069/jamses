# Create SEDesign from experimental groups

Create SEDesign from experimental groups

## Usage

``` r
groups_to_sedesign(
  ifactors,
  group_colnames = NULL,
  isamples = NULL,
  idesign = NULL,
  factor_order = NULL,
  omit_grep = "[-,]",
  max_depth = 2,
  factor_sep = "_",
  contrast_sep = "-",
  remove_pairs = NULL,
  pre_control_terms = NULL,
  add_contrastdf = NULL,
  contrast_names = NULL,
  current_depth = 1,
  rename_first_depth = TRUE,
  return_sedesign = TRUE,
  default_order = c("asis", "sort_samples", "mixedSort"),
  verbose = FALSE,
  ...
)
```

## Arguments

- ifactors:

  `data.frame` or `character` vector.

  - When `data.frame` is supplied, each column is used as a design
    factor, and rownames are recognized as sample identifiers.

  - When `character` vector is supplied, it is converted to `data.frame`
    by splitting values with a delimiter `factor_sep`, and names are
    recognized as sample identifiers.

- group_colnames:

  `character` vector or `NULL`, defines a subset of columns to use when
  `ifactors` is supplied as a `data.frame`.

  - When `ifactors` is supplied as a `character` vector, this argument
    is used to define the `colnames`.

  - `group_colnames` will become the [`factors()`](factors.md) of the
    `SEDesign`, and can be edited later using `factors()<-`.

- isamples:

  `character` vector or `NULL`, optionally used to subset the sample
  identifiers used in subsequent steps. Note that only groups and
  contrasts that contain samples will be defined.

- idesign:

  `numeric` matrix or `NULL`, intended as an optional method to use an
  existing design matrix.

- factor_order:

  `integer` or `character` vector, used to define a specific order of
  factors when generating contrasts, useful when there are multiple
  experimental factors. It can be helpful to force a secondary factor to
  be compared before a primary factor especially in two-way contrasts.
  Note that `factor_order` refers to the columns (factors) and not the
  factor levels (not column values).

- omit_grep:

  `character` regular expression pattern used to exclude secondary
  factors from contrasts.

- max_depth:

  `integer` value indicating the maximum depth of statistical contrasts
  to create. For example `max_depth=2` will allow two-way contrasts, and
  `max_depth=1` will only create one-way contrasts.

- factor_sep:

  `character` string used as a delimiter to separate experimental
  factors, when recognizing or creating experimental group names.

- contrast_sep:

  `character` string used as a delimiter to separate groups within each
  contrast name.

- remove_pairs:

  `list` of `character` vectors of factors that should not be compared.
  Each `character` vector should contain two factor levels for any given
  experimental factor, where those two factor levels should not be
  compared in the same pairwise contrast. For example, consider an
  experimental factor defined
  `treatment <- c("control", "dex", "compoundx")`. To prevent a direct
  comparison of `"dex"` to `"compoundx"`, use argument
  `remove_pairs=list(c("dex", "compoundx"))`.

- pre_control_terms:

  `character` vector used to place factor levels first in the order of
  levels, so these terms will be the denominator for contrasts. This
  approach is useful when the input `ifactors` does not already contain
  a `factor` with a specific order of factor levels.

- add_contrastdf:

  `data.frame` or `character` or `NULL`, intended to include a specific
  contrast in the output. This argument is typically used during
  iterative processing, and is not usually user-defined. It must contain

- contrast_names:

  `character` optional vector of specific contrasts to use when creating
  the contrast matrix. When `contrast_names=NULL` as default, the
  function defines contrasts using its internal logic. When
  `contrast_names` is supplied, only these `contrast_names` are used,
  with no other contrasts.

- current_depth:

  `integer` value used during iterative operations of this function.

- rename_first_depth:

  `logical` value used during iterative operations of this function.

- return_sedesign:

  `logical` used during iterative operations of this function. When
  `return_sedesign=FALSE` this function returns a `list`:

  - `"contrast_df"`: a `data.frame` as used in argument
    `add_contrastdf`, which describes each unique contrast.

  - `"contrast_names"`: a `character` vector of contrast names, which
    become [`colnames()`](https://rdrr.io/r/base/colnames.html) of the
    contrast matrix.

  - `"idesign"`: a `numeric` design matrix as defined by the input data,
    suitable for debugging purposes for example.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are ignored.

- make_unique:

  `logical` indicating whether to make output contrasts unique.

## Value

`SEDesign` object with the following slots:

- `design`: `numeric` matrix with sample-to-group association

- `contrasts`: `numeric` matrix with group-to-contrast association

- `samples`: `character` vector that represents individual sample
  replicates, equivalent to
  [`rownames()`](https://rdrr.io/r/base/colnames.html) of the `design`
  matrix.

## Details

This function creates `SEDesign` with appropriate design and contrasts,
based upon experimental groups. This approach will use multiple
experimental factors to create appropriate one-way and n-way contrasts,
where each contrast represents a symmetric comparison of each
independent factor.

Input can be provided in one of two ways:

1.  `SummarizedExperiment` where experiment design is derived from
    [`SummarizedExperiment::colData()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
    of the `se` object, and uses columns defined by `group_colnames`.
    This input should be equivalent to providing a `data.frame` whose
    [`rownames()`](https://rdrr.io/r/base/colnames.html) are equal to
    `colnames(se)`.

2.  `data.frame` where each column represents a design factor.

    - An example of `data.frame` input:

        ifactors <- data.frame(
           treatment=c("Control", "Control", "Treated", "Treated"),
           genotype=c("Wildtype", "Knockout", "Wildtype", "Knockout"))

3.  `character` vector, where design factor levels are separated by a
    delimiter such as underscore `"_"`. This input will be converted to
    `data.frame` before processing.

    - An example of `character` input:

        ifactors <- c(
           "Control_Wildtype",
           "Control_Knockout",
           "Treated_Wildtype",
           "Treated_Knockout")

When rownames are provided in the `data.frame`, or names are provided
with a `character` vector, they are retained and used as sample
identifiers.

Note: This function will change any `"-"` in a factor name to `"."`
prior to detecting valid contrasts, in order to prevent confusion and
potential problems using the contrast names in downstream analyses. This
step does not call
[`base::make.names()`](https://rdrr.io/r/base/make.names.html), so that
step should be run beforehand if required.

### Troubleshooting

- When this function returns no contrasts, or returns an unexpected
  error during processing, it is most likely due to the limitation of
  comparing one factor at a time. For example, the logic will not define
  contrast `time1_treatment1-time2_treatment2`, because this contrast
  changes two factors, it will only permit either
  `time1_treatment1-time1_treatment2` or
  `time1_treatment1-time2_treatment1`.

- `max_depth` and `factor_order` are used to define the order in which
  factors are compared, but do not affect the order of factors used for
  things like group names.

## See also

Other jam experiment design: [`SEDesign()`](SEDesign.md),
[`[,SEDesign-method`](sub-SEDesign-method.md),
[`check_sedesign()`](check_sedesign.md),
[`contrast2comp()`](contrast2comp.md),
[`contrast_colors_by_group()`](contrast_colors_by_group.md),
[`contrast_names_to_sedesign()`](contrast_names_to_sedesign.md),
[`contrastnames()`](contrastnames.md), [`contrasts()`](contrasts.md),
[`contrasts<-()`](contrasts-set.md),
[`contrasts_to_factors()`](contrasts_to_factors.md),
[`contrasts_to_venn_setlists()`](contrasts_to_venn_setlists.md),
[`design,SEDesign-method`](design-SEDesign-method.md),
`design<-,SEDesign,matrix-method`,
[`draw_oneway_contrast()`](draw_oneway_contrast.md),
[`draw_twoway_contrast()`](draw_twoway_contrast.md),
[`factors()`](factors.md),
[`filter_contrast_names()`](filter_contrast_names.md),
[`groups()`](groups.md), [`plot.SEDesign()`](plot.SEDesign.md),
[`plot_sedesign()`](plot_sedesign.md),
[`print.SEDesign()`](print.SEDesign.md), [`samples()`](samples.md),
[`sedesign_to_factors()`](sedesign_to_factors.md),
[`sort_contrasts()`](sort_contrasts.md),
[`validate_sedesign()`](validate_sedesign.md)

## Examples

``` r
# first define a vector of sample groups
igroups <- jamba::nameVector(paste(rep(c("WT", "KO"), each=6),
   rep(c("Control", "Treated"), each=3),
   sep="_"),
   suffix="_rep");
igroups <- factor(igroups, levels=unique(igroups));
igroups;
#> WT_Control_rep1 WT_Control_rep2 WT_Control_rep3 WT_Treated_rep1 WT_Treated_rep2 
#>      WT_Control      WT_Control      WT_Control      WT_Treated      WT_Treated 
#> WT_Treated_rep3 KO_Control_rep1 KO_Control_rep2 KO_Control_rep3 KO_Treated_rep1 
#>      WT_Treated      KO_Control      KO_Control      KO_Control      KO_Treated 
#> KO_Treated_rep2 KO_Treated_rep3 
#>      KO_Treated      KO_Treated 
#> Levels: WT_Control WT_Treated KO_Control KO_Treated

sedesign <- groups_to_sedesign(igroups, group_colnames=c("Genotype", "Treatment"));
design(sedesign);
#>                 WT_Control WT_Treated KO_Control KO_Treated
#> WT_Control_rep1          1          0          0          0
#> WT_Control_rep2          1          0          0          0
#> WT_Control_rep3          1          0          0          0
#> WT_Treated_rep1          0          1          0          0
#> WT_Treated_rep2          0          1          0          0
#> WT_Treated_rep3          0          1          0          0
#> KO_Control_rep1          0          0          1          0
#> KO_Control_rep2          0          0          1          0
#> KO_Control_rep3          0          0          1          0
#> KO_Treated_rep1          0          0          0          1
#> KO_Treated_rep2          0          0          0          1
#> KO_Treated_rep3          0          0          0          1
contrasts(sedesign);
#>             Contrasts
#> Levels       KO_Control-WT_Control KO_Treated-WT_Treated WT_Treated-WT_Control
#>   WT_Control                    -1                     0                    -1
#>   WT_Treated                     0                    -1                     1
#>   KO_Control                     1                     0                     0
#>   KO_Treated                     0                     1                     0
#>             Contrasts
#> Levels       KO_Treated-KO_Control
#>   WT_Control                     0
#>   WT_Treated                     0
#>   KO_Control                    -1
#>   KO_Treated                     1
#>             Contrasts
#> Levels       (KO_Treated-WT_Treated)-(KO_Control-WT_Control)
#>   WT_Control                                               1
#>   WT_Treated                                              -1
#>   KO_Control                                              -1
#>   KO_Treated                                               1

# plot the design and contrasts
plot_sedesign(sedesign)


# the two-way contrasts can be visibly flipped, since they are equivalent
plot_sedesign(sedesign, flip_twoway=TRUE)


# the design can be subset by sample
all_samples <- samples(sedesign)
subset_samples1 <- all_samples[-1:-3];
plot_sedesign(sedesign[subset_samples1, ])


# the group n=# replicates are updated
subset_samples2 <- all_samples[c(-1, -6, -11)];
plot_sedesign(sedesign[subset_samples2, ])


# The design * contrast matrix can be displayed in full
design(sedesign) %*%  contrasts(sedesign);
#>                  Contrasts
#>                   KO_Control-WT_Control KO_Treated-WT_Treated
#>   WT_Control_rep1                    -1                     0
#>   WT_Control_rep2                    -1                     0
#>   WT_Control_rep3                    -1                     0
#>   WT_Treated_rep1                     0                    -1
#>   WT_Treated_rep2                     0                    -1
#>   WT_Treated_rep3                     0                    -1
#>   KO_Control_rep1                     1                     0
#>   KO_Control_rep2                     1                     0
#>   KO_Control_rep3                     1                     0
#>   KO_Treated_rep1                     0                     1
#>   KO_Treated_rep2                     0                     1
#>   KO_Treated_rep3                     0                     1
#>                  Contrasts
#>                   WT_Treated-WT_Control KO_Treated-KO_Control
#>   WT_Control_rep1                    -1                     0
#>   WT_Control_rep2                    -1                     0
#>   WT_Control_rep3                    -1                     0
#>   WT_Treated_rep1                     1                     0
#>   WT_Treated_rep2                     1                     0
#>   WT_Treated_rep3                     1                     0
#>   KO_Control_rep1                     0                    -1
#>   KO_Control_rep2                     0                    -1
#>   KO_Control_rep3                     0                    -1
#>   KO_Treated_rep1                     0                     1
#>   KO_Treated_rep2                     0                     1
#>   KO_Treated_rep3                     0                     1
#>                  Contrasts
#>                   (KO_Treated-WT_Treated)-(KO_Control-WT_Control)
#>   WT_Control_rep1                                               1
#>   WT_Control_rep2                                               1
#>   WT_Control_rep3                                               1
#>   WT_Treated_rep1                                              -1
#>   WT_Treated_rep2                                              -1
#>   WT_Treated_rep3                                              -1
#>   KO_Control_rep1                                              -1
#>   KO_Control_rep2                                              -1
#>   KO_Control_rep3                                              -1
#>   KO_Treated_rep1                                               1
#>   KO_Treated_rep2                                               1
#>   KO_Treated_rep3                                               1

# make "KO" the control term instead of "WT"
contrast_names(groups_to_sedesign(igroups, pre_control_terms=c("KO")))
#> [1] "WT_Control-KO_Control"                          
#> [2] "WT_Treated-KO_Treated"                          
#> [3] "KO_Treated-KO_Control"                          
#> [4] "WT_Treated-WT_Control"                          
#> [5] "(WT_Treated-KO_Treated)-(WT_Control-KO_Control)"

# change the order of factors compared
contrast_names(groups_to_sedesign(igroups, factor_order=2:1))
#> [1] "WT_Treated-WT_Control"                          
#> [2] "KO_Treated-KO_Control"                          
#> [3] "KO_Control-WT_Control"                          
#> [4] "KO_Treated-WT_Treated"                          
#> [5] "(KO_Treated-KO_Control)-(WT_Treated-WT_Control)"

# prevent comparisons of WT to WT, or KO to KO
sedesign_2 <- groups_to_sedesign(as.character(igroups),
   remove_pairs=list(c("WT"), c("KO")))
contrast_names(sedesign_2)
#> [1] "KO_Control-WT_Control"                          
#> [2] "KO_Treated-WT_Treated"                          
#> [3] "(KO_Treated-WT_Treated)-(KO_Control-WT_Control)"
plot_sedesign(sedesign_2)


# prevent comparisons of Treated to Treated, or Control to Control
sedesign_3 <- groups_to_sedesign(as.character(igroups),
   remove_pairs=list(c("Treated"), c("Control")))
contrast_names(sedesign_3)
#> [1] "WT_Treated-WT_Control"                          
#> [2] "KO_Treated-KO_Control"                          
#> [3] "(KO_Treated-KO_Control)-(WT_Treated-WT_Control)"
plot_sedesign(sedesign_3)


# input as a data.frame with ordered factor levels
ifactors <- data.frame(Genotype=factor(c("WT","WT","KO","KO"),
   levels=c("WT","KO")),
   Treatment=factor(c("Treated","Control"),
      levels=c("Control","Treated")))
# not necessary, but define rownames
rownames(ifactors) <- jamba::pasteByRow(ifactors);
ifactors;
#>            Genotype Treatment
#> WT_Treated       WT   Treated
#> WT_Control       WT   Control
#> KO_Treated       KO   Treated
#> KO_Control       KO   Control
contrast_names(groups_to_sedesign(ifactors))
#> [1] "KO_Control-WT_Control"                          
#> [2] "KO_Treated-WT_Treated"                          
#> [3] "WT_Treated-WT_Control"                          
#> [4] "KO_Treated-KO_Control"                          
#> [5] "(KO_Treated-WT_Treated)-(KO_Control-WT_Control)"
plot_sedesign(groups_to_sedesign(ifactors))


# you can still override factor levels with pre_control_terms
plot_sedesign(groups_to_sedesign(ifactors, pre_control_terms=c("KO")))


# input as design matrix
design_matrix <- design(groups_to_sedesign(ifactors))
design_matrix
#>            WT_Control WT_Treated KO_Control KO_Treated
#> WT_Treated          0          1          0          0
#> WT_Control          1          0          0          0
#> KO_Treated          0          0          0          1
#> KO_Control          0          0          1          0
contrast_names(groups_to_sedesign(design_matrix))
#> [1] "KO_Control-WT_Control"                          
#> [2] "KO_Treated-WT_Treated"                          
#> [3] "WT_Treated-WT_Control"                          
#> [4] "KO_Treated-KO_Control"                          
#> [5] "(KO_Treated-WT_Treated)-(KO_Control-WT_Control)"

# again the "KO" group can be the control by using pre_control_terms
contrast_names(groups_to_sedesign(design_matrix, pre_control_terms="KO"))
#> [1] "WT_Control-KO_Control"                          
#> [2] "WT_Treated-KO_Treated"                          
#> [3] "KO_Treated-KO_Control"                          
#> [4] "WT_Treated-WT_Control"                          
#> [5] "(WT_Treated-KO_Treated)-(WT_Control-KO_Control)"

# default_order="asis"
# contrasts show A-B, because B appears fist
# contrasts show Untreated-Treated because Treated appears first
df_test <- data.frame(
   set=c("B", "B", "A", "A"),
   treat=c("Treated", "Untreated"))
plot_sedesign(groups_to_sedesign(df_test))

plot_sedesign(groups_to_sedesign(jamba::pasteByRow(df_test)))

# default_order="sort_samples"
# contrasts show B-A, because A is sorted first
# contrasts show Treated-Untreated because sort_samples()
#    determines "Untreated" is a preferred control term
plot_sedesign(groups_to_sedesign(df_test,
   default_order="sort_samples"))


# default_order="mixedSort"
# contrasts show B-A, because A is sorted first
# contrasts show Untreated-Treated because Treated is sorted first
plot_sedesign(groups_to_sedesign(df_test,
   default_order="mixedSort"))

plot_sedesign(groups_to_sedesign(df_test,
   default_order="mixedSort",
   pre_control_terms=c("Untreated")))

```
