# Define contrast colors by group colors

Define contrast colors by group colors

## Usage

``` r
contrast_colors_by_group(
  sedesign,
  sample_color_list = NULL,
  group_name = "Group",
  contrast_sep = "-",
  C_min = 80,
  C_max = 100,
  L_min = 35,
  L_max = 70,
  default_color_style = c("categorical", "navy"),
  verbose = FALSE,
  ...
)
```

## Arguments

- sedesign:

  `SEDesign` object

- sample_color_list:

  `list` of colors, with one list element with name defined by
  `group_name`, default `"Group"`, so the colors use
  `sample_color_list[[group_name]]`. This element should be a
  `character` vector of colors, whose names are group names used in the
  `contrast_names(sedesign)`. The `sample_color_list` is produced by
  `platjam::design2colors()`. Alternatively, a `character` vector of
  colors, named by contrast.

- group_name:

  `character` string for the element in `sample_color_list` that
  contains group colors.

- contrast_sep:

  `character` string used as delimiter between group names in each
  contrast name, default is `"-"`.

- C_min, C_max:

  `numeric` value to impose a minimum and maximum chroma (C) value for
  the resulting contrast color.

- L_min, L_max:

  `numeric` value to impose a minimum and maximum luminance (C) value
  for the resulting contrast color.

- default_color_style:

  `character` string indicating how to assign colors to groups without a
  color assignment.

  - `"categorical"`: calls
    [`colorjam::group2colors()`](https://jmw86069.github.io/colorjam/reference/group2colors.html)
    which calls
    [`colorjam::rainbowJam()`](https://jmw86069.github.io/colorjam/reference/rainbowJam.html)
    for rainbow categorical colors.

  - any single R color: this one R color is applied to all other groups.

  - any non-R color: assumed to be the name of a color gradient or set,
    and is resolved by calling
    [`jamba::getColorRamp()`](https://jmw86069.github.io/jamba/reference/getColorRamp.html).
    It recognizes RColorBrewer and viridis color sets for example.

  - vector of multiple colors: also resolved by calling
    [`jamba::getColorRamp()`](https://jmw86069.github.io/jamba/reference/getColorRamp.html)
    which creates a gradient across the colors.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are passed to
  [`colorjam::group2colors()`](https://jmw86069.github.io/colorjam/reference/group2colors.html),
  or
  [`jamba::getColorRamp()`](https://jmw86069.github.io/jamba/reference/getColorRamp.html)
  as needed.

## Value

`character` vector of R colors, with names defined by
`contrast_names(sedesign)`.

## Details

This function takes colors assigned to groups, and applies them to
contrasts based upon the groups represented in each contrast. Group
colors are blended into one output color per contrast, then the color is
adjusted for minimum chroma (C) and maximum luminance (L) to ensure the
colors have acceptable vibrance.

When `sample_color_list` is supplied, which defines colors per group,
these colors are used as appropriate.

When `sample_color_list` is not supplied, categorical colors are defined
for all observed group names using
[`colorjam::group2colors()`](https://jmw86069.github.io/colorjam/reference/group2colors.html).

## See also

Other jam experiment design:
[`SEDesign()`](https://jmw86069.github.io/jamses/reference/SEDesign.md),
[`[,SEDesign-method`](https://jmw86069.github.io/jamses/reference/sub-SEDesign-method.md),
[`check_sedesign()`](https://jmw86069.github.io/jamses/reference/check_sedesign.md),
[`contrast2comp()`](https://jmw86069.github.io/jamses/reference/contrast2comp.md),
[`contrast_names_to_sedesign()`](https://jmw86069.github.io/jamses/reference/contrast_names_to_sedesign.md),
[`contrastnames()`](https://jmw86069.github.io/jamses/reference/contrastnames.md),
[`contrasts()`](https://jmw86069.github.io/jamses/reference/contrasts.md),
[`contrasts<-()`](https://jmw86069.github.io/jamses/reference/contrasts-set.md),
[`contrasts_to_factors()`](https://jmw86069.github.io/jamses/reference/contrasts_to_factors.md),
[`contrasts_to_venn_setlists()`](https://jmw86069.github.io/jamses/reference/contrasts_to_venn_setlists.md),
[`design,SEDesign-method`](https://jmw86069.github.io/jamses/reference/design.md),
[`draw_oneway_contrast()`](https://jmw86069.github.io/jamses/reference/draw_oneway_contrast.md),
[`draw_twoway_contrast()`](https://jmw86069.github.io/jamses/reference/draw_twoway_contrast.md),
[`factors()`](https://jmw86069.github.io/jamses/reference/factors.md),
[`filter_contrast_names()`](https://jmw86069.github.io/jamses/reference/filter_contrast_names.md),
[`groups()`](https://jmw86069.github.io/jamses/reference/groups.md),
[`groups_to_sedesign()`](https://jmw86069.github.io/jamses/reference/groups_to_sedesign.md),
[`plot.SEDesign()`](https://jmw86069.github.io/jamses/reference/plot.SEDesign.md),
[`plot_sedesign()`](https://jmw86069.github.io/jamses/reference/plot_sedesign.md),
[`print,SEDesign-method`](https://jmw86069.github.io/jamses/reference/print-SEDesign-method.md),
[`samples()`](https://jmw86069.github.io/jamses/reference/samples.md),
[`sedesign_to_factors()`](https://jmw86069.github.io/jamses/reference/sedesign_to_factors.md),
[`sort_contrasts()`](https://jmw86069.github.io/jamses/reference/sort_contrasts.md),
[`validate_sedesign()`](https://jmw86069.github.io/jamses/reference/validate_sedesign.md)
