# Add Heatmap column group labels

Add Heatmap column group labels, specifically for ComplexHeatmap output.

## Usage

``` r
heatmap_column_group_labels(
  hm_group_list,
  hm_drawn = NULL,
  se = NULL,
  add_group_label = TRUE,
  add_group_line = TRUE,
  add_group_box = FALSE,
  group_line_buffer = grid::unit(1, "mm"),
  group_line_lwd = 2,
  group_line_requires_label = TRUE,
  group_box_lwd = 2,
  group_box_outer = TRUE,
  font_cex = 1,
  hm_title_base = NULL,
  hm_body_base = NULL,
  hm_title_base_default = "centered\nexpression_column_title_",
  hm_body_base_default = "centered\nexpression_heatmap_body_1_",
  y_offset_lines = 0,
  endlines = 1,
  use_gridtext = TRUE,
  verbose = FALSE,
  ...
)
```

## Arguments

- hm_group_list:

  `character` or `list` with one of the following types of content:

  - `character` colnames in the column data of `se` argument,
    specifically `colnames(colData(se))` to define groupings.

  - `list` of `integer` vectors, where list names become labels for each
    group.

- hm_drawn:

  [`ComplexHeatmap::HeatmapList`](https://rdrr.io/pkg/ComplexHeatmap/man/HeatmapList.html)
  object which is returned by the function
  [`ComplexHeatmap::draw()`](https://rdrr.io/pkg/ComplexHeatmap/man/draw-dispatch.html).
  This object is a pointer to the grid graphical elements required.

- se:

  `SummarizedExperiment` object, required only when `hm_group_list` is
  supplied in the form of colnames of the
  `SummarizedExperiment::colData(se)`.

- add_group_label:

  `logical` indicating whether to draw the group `character` label above
  each group defined in `hm_group_list`.

- add_group_line:

  `logical` indicating whether to draw the group line for each group
  defined in `hm_group_list`. The line is intended to appear below the
  text label when `add_group_label=TRUE`.

- add_group_box:

  `logical` indicating whether to draw a box around the heatmap region
  defined by each group in `hm_group_list`.

- group_line_buffer:

  [`grid::unit`](https://rdrr.io/r/grid/unit.html) object indicating the
  whitespace buffer region betwen adjacent lines when
  `add_group_line=TRUE`. The default enforces 1-mm white space between
  lines so they do not touch.

- group_line_lwd:

  `numeric` line width when `add_group_line=TRUE`.

- group_line_requires_label:

  `logical` whether to require the associated group label to contain
  non-whitespace visible text.

  - `group_line_requires_label=TRUE` (default) requires group label to
    have visible characters, which means the presence of an empty labels
    will cause the group line not to be drawn.

  - `group_line_requires_label=FALSE` does not require a group label,
    therefore all group lines are drawn.

- group_box_lwd:

  `numeric` line width when `add_group_box=TRUE`.

- group_box_outer:

  `logical` indicating whether to draw the group box as an outer border,
  which means the border will be drawn only on the very outer edge of
  the heatmap body, and will not overlap the heatmap contents. This
  option is only active when `add_group_box=TRUE`.

- font_cex:

  `numeric` adjustment for font sizes overall.

- hm_title_base:

  `character` string used to search for matching heatmap column title
  grid layout regions. It is derived from the `Heatmap` argument `name`,
  usually followed by
  `"_column_title_". When `NULL`the default values are defined by`detect_heatmap_components()`element`column_title_components\`.

- hm_body_base:

  `character` string used to search for matching heatmap body grid
  layout regions. It is derived from the `Heatmap` argument `name`,
  usually followed by
  `"_body_1_". When `NULL`the default values are defined by`detect_heatmap_components()`element`heatmap_body_components\`.
  Note: Currently the group box only includes the first row_split,
  although it is intended to include the entire set of heatmap rows in
  future iterations.

- hm_title_base_default, hm_body_base_default:

  `character` string used as reference for hm_title_base,hm_body_base.

- y_offset_lines:

  `numeric` adjustment used to shift the text label and underline by
  this many lines of character height. It is mainly used internally for
  iterative calls, when `hm_group_list` includes multiple layers of
  groups. The text and underline are intended to fill 2 character height
  lines, although this argument can be used to make manual adjustments.

- endlines:

  `integer` number of blank lines to append to the end of each group
  label, which has the effect of shifting the group label upwards
  slightly.

- use_gridtext:

  `logical` (default TRUE) whether to render text using
  [`gridtext::richtext_grob()`](https://wilkelab.org/gridtext/reference/richtext_grob.html)
  in order to enable markdown and limited HTML-based formatting.

- verbose:

  `logical` indicating whether to print verbose output.

- ...:

  additional arguments are ignored.

## Details

This function is currently experimental, and is intended only for a
specific scenario, to augment a
[`ComplexHeatmap::Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
object, as produced by [`heatmap_se()`](heatmap_se.md), that used the
argument `column_split` to sub-divide the heatmap columns into
subgroups. This function draws labels above the heatmap to describe
group factor values associated with column splits.

When there are multiple layers of grouping, this function will also draw
multiple layers. For example, when columns are split by
`column_split=c("Treatment", "Time")`, it will produce a heatmap where
each column slice has one combination of Treatment and Time. This
function can be used to add a layer of labels by `"Time"`, and a layer
of labels by `"Treatment"`. The labels are shown by default with a broad
underline to indicate contiguous column slices that contain the same
`"Treatment"` or `"Time"` values.

The output can provide a cleaner visualization than the alternative of
displaying colorized annotation boxes at the top of the heatmap.

### Input Formats

There are two strategies for defining the column group data to display:

#### Input using colnames from colData(se)

This option is recommended as the "easiest" method. It requires:

- `hm_group_list` provides a `character` vector of
  `colnames(colData(se))`. Note that the columns are applied
  bottom-to-top, so it is sometimes helpful to supply columns in reverse
  order for `column_split`.

- `se` must be provided, since it supplies `colData(se)`

Each value in `hm_group_list` is applied from bottom-to-top, to define a
row of labels. By default, `add_group_line=TRUE`, so labels are also
underlined to indicate samples included in each group.

#### Input using a named list

This option is recommended when labeling a subset of column slices, or
when the column slice groups need to be customized in some way.

### Additional requirements

- The Heatmap `column_title` should (usually) be empty, so that no text
  labels are drawn which may overlap the labels drawn by this function.
  Use `column_title=" "` (with a whitespace character) to prevent
  [`ComplexHeatmap::Heatmap()`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
  from using its internal default text label.

- There usually needs to be whitespace above the heatmap, which can be
  accomplished when drawing a Heatmap object `hm` like this:
  `ComplexHeatmap::draw(hm, column_title="\n\n\n\n")`

- When also displaying a heatmap title as defined by
  [`heatmap_se()`](heatmap_se.md), the adjustment can be done like this:

    ComplexHeatmap::draw(hm,
       column_title=paste0(attr(hm, "hm_title"),
          "\n\n\n\n"))

- The arguments `hm_title_base` and `hm_body_base` should match the
  heatmap name, which is defined in [`heatmap_se()`](heatmap_se.md) with
  `data_type`, usually prefixed `"centered\n"` when the data is centered
  by that function. For example, when `data_type="abundance"`, the
  corresponding argument value should be
  `hm_title_base="centered\nabundance_column_title_"`.

- use [`detect_heatmap_components()`](detect_heatmap_components.md) to
  help identify the appropriate `grid` elements available to be used by
  this function.

When the heatmap contains `"column_title"` elements defined in
[`ComplexHeatmap::list_components()`](https://rdrr.io/pkg/ComplexHeatmap/man/list_components.html),
and there is no element `"global_column_title"`, the `y_offset_lines` is
adjusted down so that the position is inside the column_title region,
typically below the column_title when using trailing whitespace (see
argument `hm_title_buffer` in [`heatmap_se()`](heatmap_se.md)). However,
when element `"global_column_title"` is present, the `y_offset_lines`
are not adjusted, so that the position is above the column_title region.
In this case, to position the labels below the column_title, you can
adjust `y_offset_lines` down manually like this: `y_offset_lines=-9`.

### Todo

- Automate determining the column_title grid layout name.

  - Specifically when only one `heatmapname_column_title_1` is present,
    but there are multiple column groups, it should take the x-axis
    coordinate values (left and right boundaries) from the heatmap body
    instead of the column_title region.

- When `add_group_box=TRUE`, and `row_split` indicates multiple row
  groups, it should calculate the y-axis coordinate values (top and
  bottom boundaries) using the full set of row groups.

- Enable blank annotations, either by passing a subset `se`, or by
  annotations with no associated label.

## See also

Other jamses heatmaps:
[`detect_heatmap_components()`](detect_heatmap_components.md),
[`heatmap_profile_plot()`](heatmap_profile_plot.md),
[`heatmap_se()`](heatmap_se.md)
