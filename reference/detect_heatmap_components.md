# Detect ComplexHeatmap Heatmap grid layout components

Detect ComplexHeatmap Heatmap grid layout components

## Usage

``` r
detect_heatmap_components(...)
```

## Arguments

- ...:

  additional arguments are ignored.

## Value

`list` of `character` vectors, where `list` names represent different
features of the heatmap, and each `character` vector includes the grid
layout component name stem suitable for use in
[`heatmap_column_group_labels()`](https://jmw86069.github.io/jamses/reference/heatmap_column_group_labels.md).

## Details

This function is a wrapper around
[`ComplexHeatmap::list_components()`](https://rdrr.io/pkg/ComplexHeatmap/man/list_components.html)
which also creates a `list` of component names based upon common naming
conventions used in the `ComplexHeatmap` package.

## See also

Other jamses heatmaps:
[`heatmap_column_group_labels()`](https://jmw86069.github.io/jamses/reference/heatmap_column_group_labels.md),
[`heatmap_profile_plot()`](https://jmw86069.github.io/jamses/reference/heatmap_profile_plot.md),
[`heatmap_se()`](https://jmw86069.github.io/jamses/reference/heatmap_se.md)
