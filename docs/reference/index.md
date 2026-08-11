# Package index

## All functions

- [`check_sedesign()`](check_sedesign.md) : Check SEDesign object
- [`choose_annotation_colnames()`](choose_annotation_colnames.md) :
  Choose interesting annotation colnames from a data.frame
- [`combine_sestats()`](combine_sestats.md) : Combine SEStats objects by
  contrast,signal
- [`contrast2comp()`](contrast2comp.md)
  [`comp2contrast()`](contrast2comp.md)
  [`names_contrast2comp()`](contrast2comp.md)
  [`names_comp2contrast()`](contrast2comp.md) : Convert contrast to
  short-form comp, convert comp to contrast
- [`contrast2comp_dev()`](contrast2comp_dev.md) : Convert contrast to
  short-form comp, convert comp to contrast (DEV)
- [`contrastnames()`](contrastnames.md)
  [`` `contrastnames<-`() ``](contrastnames.md)
  [`contrast_names()`](contrastnames.md)
  [`` `contrast_names<-`() ``](contrastnames.md) : Contrast names for
  SEDesign objects
- [`` `contrasts<-`() ``](contrasts-set.md) : Contrast matrix setter for
  SEDesign objects
- [`contrasts()`](contrasts.md) : Contrast matrix accessor for SEDesign
  objects
- [`contrasts_to_factors()`](contrasts_to_factors.md) : Convert
  contrasts to data.frame of design factors
- [`contrasts_to_venn_setlists()`](contrasts_to_venn_setlists.md) :
  Convert contrast names to Venn setlists for visual comparison
- [`contrast_colors_by_group()`](contrast_colors_by_group.md) : Define
  contrast colors by group colors
- [`contrast_names()`](contrast_names-SEStats-method.md) : Get contrast
  names from SEStats
- [`contrast_names_to_sedesign()`](contrast_names_to_sedesign.md) :
  Convert contrast names to sedesign object
- [`design(`*`<SEDesign>`*`)`](design-SEDesign-method.md) : Design
  matrix accessor for SEDesign objects
- [`` `design<-`( ``*`<SEDesign>`*`,`*`<matrix>`*`)`](design-set-SEDesign-matrix-method.md)
  : Design matrix setter for SEDesign objects
- [`detect_heatmap_components()`](detect_heatmap_components.md) : Detect
  ComplexHeatmap Heatmap grid layout components
- [`dimnames()`](dimnames-SEDesign-method.md) : Get dimnames from
  SEDesign
- [`dimnames()`](dimnames-SEStats-method.md) : Get dimnames from SEStats
- [`.class_list`](dot-class_list.md) : SEStats S7 Object
- [`draw_oneway_contrast()`](draw_oneway_contrast.md) : Draw one-way
  contrast using block arrows
- [`draw_twoway_contrast()`](draw_twoway_contrast.md) : Draw two one-way
  contrasts using block arrows, showing a two-way connector
- [`ebayes2dfs()`](ebayes2dfs.md) : Convert limma eBayes fit to
  data.frame with annotated hits
- [`factors()`](factors-SEStats-method.md) : Get factor names from
  SEStats
- [`factors()`](factors.md) [`` `factors<-`() ``](factors.md) :
  Experimental factor labels for SEDesign objects
- [`filter_contrast_names()`](filter_contrast_names.md) : Filter
  contrast names
- [`fold_to_log2fold()`](fold_to_log2fold.md) : Convert normal signed
  fold change to log2 fold change
- [`format_hits()`](format_hits.md) : Format list of hit vectors into
  summary counts
- [`geomx_to_se()`](geomx_to_se.md) : Convert NanoStringGeoMxSet to
  SummarizedExperiment
- [`groups()`](groups-SEStats-method.md) : Get groups from SEStats
- [`groups()`](groups.md) [`` `groups<-`() ``](groups.md) : Design group
  names for SEDesign objects
- [`groups_to_sedesign()`](groups_to_sedesign.md) : Create SEDesign from
  experimental groups
- [`handle_na_values()`](handle_na_values.md) : Handle NA values in a
  numeric matrix
- [`heatmap_column_group_labels()`](heatmap_column_group_labels.md) :
  Add Heatmap column group labels
- [`heatmap_profile_plot()`](heatmap_profile_plot.md) : Profile Line
  Plot for Heatmap of SummarizedExperiment data
- [`heatmap_se()`](heatmap_se.md) : Heatmap for SummarizedExperiment
  data
- [`hit_array_to_list()`](hit_array_to_list.md) : Quick conversion of
  hit_array to hit_list
- [`intercalate()`](intercalate.md) : Intercalate two or more vectors
- [`list2im_opt()`](list2im_opt.md) : Convert list to incidence matrix
- [`list_to_sestats()`](list_to_sestats.md) : Prepare SEStats from a
  list of stat data.frame
- [`log2fold_to_fold()`](log2fold_to_fold.md) : Convert log2 fold change
  to signed fold change
- [`make_block_arrow_polygon()`](make_block_arrow_polygon.md) : Make
  block arrow polygon coordinates for line segments
- [`make_se_test()`](make_se_test.md) : Make SummarizedExperiment test
  data
- [`mark_stat_hits()`](mark_stat_hits.md) : Mark statistical hits by
  threshold cutoffs
- [`matrix_normalize()`](matrix_normalize.md) : Normalize a numeric data
  matrix
- [`merge_statdf_all_test()`](merge_statdf_all_test.md) : Merge stats
  data.frame from two sestats results
- [`plot(`*`<SEDesign>`*`)`](plot.SEDesign.md) : Plot method for
  SEDesign objects
- [`plot_sedesign()`](plot_sedesign.md) : Plot sedesign object contrasts
- [`point_handedness()`](point_handedness.md) : Determine which side one
  point is to another, given a slope or angle
- [`point_slope_intercept()`](point_slope_intercept.md) : convert
  point-slope to axis intercept
- [`print()`](print-SEStats-method.md) : Print SEStats Object
- [`print(`*`<SEDesign>`*`)`](print.SEDesign.md) : Print / show method
  for SEDesign objects
- [`process_sestats_to_hitim()`](process_sestats_to_hitim.md) : Process
  sestats input into a hit incidence matrix
- [`run_limma_replicate()`](run_limma_replicate.md) : Run limma
  contrasts with optional probe replicates
- [`samples()`](samples-SEStats-method.md) : Get samples from SEStats
- [`samples()`](samples.md) [`` `samples<-`() ``](samples.md) : Sample
  identifiers for SEDesign objects
- [`save_sestats()`](save_sestats.md) : Save SE contrast stats output
- [`SEDesign()`](SEDesign.md) : SEDesign: experiment design and
  contrasts object (S7 class)
- [`sedesign_to_factors()`](sedesign_to_factors.md) : Convert SEDesign
  to data.frame of design factors
- [`sestats_to_df()`](sestats_to_df.md) : Convert sestats to table
  summary
- [`sestats_to_dfs()`](sestats_to_dfs.md) : Extract stats as data.frame
  from SEStats results
- [`se_collapse_by_column()`](se_collapse_by_column.md) : Collapse
  SummarizedExperiment data by column
- [`se_collapse_by_row()`](se_collapse_by_row.md) : Collapse
  SummarizedExperiment data by row
- [`se_contrast_stats()`](se_contrast_stats.md) : Compute contrast
  statistics on SummarizedExperiment data
- [`se_detected_rows()`](se_detected_rows.md) : SummarizedExperiment
  heuristics to define detected rows
- [`se_normalize()`](se_normalize.md) : Normalize SummarizedExperiment
  data
- [`se_rbind()`](se_rbind.md) : Combine SummarizedExperiment objects by
  row
- [`se_to_assay_data()`](se_to_assay_data.md) : Get SE assay data
- [`se_to_assay_names()`](se_to_assay_names.md) : Get SE assay names
- [`se_to_rowcoldata()`](se_to_rowcoldata.md) : Get SE colData and
  rowData
- [`shortest_unique_abbreviation()`](shortest_unique_abbreviation.md) :
  Find the shortest abbrevation to retain unique values
- [`shrinkDataFrame()`](shrinkDataFrame.md) : Shrink data.frame by row
  groups
- [`shrink_df()`](shrink_df.md) : Shrink data.frame by row groups
- [`shrink_matrix()`](shrink_matrix.md) : Shrink a numeric matrix by
  groups of rows
- [`sort_contrasts()`](sort_contrasts.md) : Sort contrasts by factor and
  level
- [`sort_samples()`](sort_samples.md) : Sort biological sample labels
  for experimental design
- [`stack_contrasts()`](stack_contrasts.md) : Stack a vector of
  contrasts into a full string
- [`strsplitOrdered()`](strsplitOrdered.md) : Split the elements of an
  ordered factor vector
- [`` `[` ``](sub-SEDesign-method.md) : Subset a SEDesign object by
  samples and/or design groups
- [`sub_split_vector()`](sub_split_vector.md) : Sub-split a split vector
  (Internal)
- [`update_function_params()`](update_function_params.md) : Update
  function default parameters
- [`update_list_elements()`](update_list_elements.md) : Update a subset
  of list elements
- [`validate_sedesign()`](validate_sedesign.md) : Validate SEDesign
  object contents
- [`voom_jam()`](voom_jam.md) : Limma-voom customized for Jam
