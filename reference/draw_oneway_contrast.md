# Draw one-way contrast using block arrows

Draw one-way contrast using block arrows

## Usage

``` r
draw_oneway_contrast(
  x0,
  x1,
  y0,
  y1,
  color = "peachpuff",
  border = "black",
  plot_type = c("base", "grid"),
  label = NULL,
  label_sep = "\n",
  na.rm = TRUE,
  label_color = "white",
  label_cex = 1,
  label_font = 2,
  oneway_position = 0.5,
  do_plot = TRUE,
  verbose = FALSE,
  debug = FALSE,
  ...
)
```

## Arguments

- x0, x1, y0, y1:

  `numeric` values with the start and end coordinates, for the x and y
  axes, respectively.

- color, border:

  `character` R colors used to define color fill and border,
  respectively, for each block arrow,

- plot_type:

  `character` string indicating the type of plot output:

  - `"base"`: base R graphics

  - `"grid"`: grid graphics (not yet implemented)

- label:

  `character` vector or `list` with optional label to display atop each
  block arrow. For base R graphics, the label is drawn using
  [`jamba::shadowText()`](https://jmw86069.github.io/jamba/reference/shadowText.html)
  to render an outline around the text.

  - When `label` is a `character` vector, it is converted to a `list` in
    two ways depending upon the number of block arrows (`length(x0)`):

    - `length == 1`: `label` is converted to `list` with length == 1.

    - `length > 1`: `label` is converted to `list` using `as.list`, then
      expanded to `length(x0)`.

  - When `label` is passed as a `list`, or after `label` is converted to
    a `list`:

    - Each block arrow label uses one concatenated string after calling
      `jamba::cPaste(..., sep=label_sep)` which joins values by default
      using newline `"\n"` between each value.

- label_sep:

  `character` string used as separator, passed to
  [`jamba::cPaste()`](https://jmw86069.github.io/jamba/reference/cPaste.html),
  so that each block arrow may contain a vector which is concatenated
  using `label_sep` between each value. By default `label_sep="\n"`
  which prints each value on a new line.

- na.rm:

  `logical` passed to
  [`jamba::cPaste()`](https://jmw86069.github.io/jamba/reference/cPaste.html)
  to define how to display NA labels:

  - `na.rm=FALSE`: `"NA"`

  - `na.rm=TRUE`: `""`.

- label_color:

  `character` color used for the `label`.

- label_cex:

  `numeric` label font expansion factor, used to adjust the font size of
  the text label.

- label_font:

  `numeric` indicating the font face, defined as:

  - 1 = normal font

  - 2 = bold font

  - 3 = italic font

  - 4 = bold, italic font

- do_plot:

  `logical` indicating whether to draw the block arrow.

- verbose:

  `logical` indicating whether to print verbose output.

- debug:

  `logical` indicating whether to print additional debug info.

- ...:

  additional arguments are passed to
  [`make_block_arrow_polygon()`](https://jmw86069.github.io/jamses/reference/make_block_arrow_polygon.md),
  including `arrow_ex` the arrow size expansion factor, and `head_ex`
  the arrow head expansion factor, which is adjusted relative to the
  arrow stem width.

## See also

Other jam experiment design:
[`SEDesign()`](https://jmw86069.github.io/jamses/reference/SEDesign.md),
[`[,SEDesign-method`](https://jmw86069.github.io/jamses/reference/sub-SEDesign-method.md),
[`check_sedesign()`](https://jmw86069.github.io/jamses/reference/check_sedesign.md),
[`contrast2comp()`](https://jmw86069.github.io/jamses/reference/contrast2comp.md),
[`contrast_colors_by_group()`](https://jmw86069.github.io/jamses/reference/contrast_colors_by_group.md),
[`contrast_names_to_sedesign()`](https://jmw86069.github.io/jamses/reference/contrast_names_to_sedesign.md),
[`contrastnames()`](https://jmw86069.github.io/jamses/reference/contrastnames.md),
[`contrasts()`](https://jmw86069.github.io/jamses/reference/contrasts.md),
[`contrasts<-()`](https://jmw86069.github.io/jamses/reference/contrasts-set.md),
[`contrasts_to_factors()`](https://jmw86069.github.io/jamses/reference/contrasts_to_factors.md),
[`contrasts_to_venn_setlists()`](https://jmw86069.github.io/jamses/reference/contrasts_to_venn_setlists.md),
[`design,SEDesign-method`](https://jmw86069.github.io/jamses/reference/design-SEDesign-method.md),
`design<-,SEDesign,matrix-method`,
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

## Examples

``` r
plot(NULL, xlim=c(0, 5), ylim=c(0, 4), asp=1, xlab="", ylab="")
draw_oneway_contrast(1, 4, 1, 1, label="contrast label")
draw_oneway_contrast(1, 4, 2, 2, head_ex=2, label_cex=1, label="contrast label")
draw_oneway_contrast(1, 4, 3, 3, arrow_ex=2, label_cex=2, label="contrast label")
draw_oneway_contrast(3, 0, 4, 1, arrow_ex=2, label_cex=2, label="contrast label")

```
