
context("contrast2comp")

testthat::test_that("contrast2comp one factor", {
   contrast_names_1fac <- c(
      "Treated-Control",
      "Knockout-Control"
   );
   testthat::expect_equal(
      contrast2comp(contrast_names_1fac),
      contrast_names_1fac)
})

testthat::test_that("contrast2comp two factors", {
   contrast_names_2fac <- c(
      "CellA_Treated-CellA_Control",
      "CellB_Treated-CellB_Control",
      "CellB_Treated-CellA_Control",
      "(CellA_Treated-CellA_Control)-(CellB_Treated-CellB_Control)",
      "(CellB_Treated-CellB_Control)-(CellA_Treated-CellA_Control)",
      "(CellA_Treated-CellB_Treated)-(CellA_Control-CellB_Control)"
   );
   expected_comps_2fac <- c(
      "CellA:Treated-Control",
      "CellB:Treated-Control",
      "CellB_Treated-CellA_Control",
      "CellA-CellB:Treated-Control",
      "CellB-CellA:Treated-Control",
      "CellA-CellB:Treated-Control")
   testthat::expect_equal(
      contrast2comp(contrast_names_2fac),
      expected_comps_2fac)
})

testthat::test_that("comp2contrast two factors", {
   contrast_names_2fac <- c(
      "CellA_Treated-CellA_Control",
      "CellB_Treated-CellB_Control",
      "CellB_Treated-CellA_Control",
      "(CellA_Treated-CellA_Control)-(CellB_Treated-CellB_Control)",
      "(CellB_Treated-CellB_Control)-(CellA_Treated-CellA_Control)",
      "(CellA_Treated-CellB_Treated)-(CellA_Control-CellB_Control)"
   );
   contrast_names_2fac_adj <- c(
      "CellA_Treated-CellA_Control",
      "CellB_Treated-CellB_Control",
      "CellB_Treated-CellA_Control",
      "(CellA_Treated-CellB_Treated)-(CellA_Control-CellB_Control)",
      "(CellB_Treated-CellA_Treated)-(CellB_Control-CellA_Control)",
      "(CellA_Treated-CellB_Treated)-(CellA_Control-CellB_Control)"
   );
   expected_comps_2fac <- c(
      "CellA:Treated-Control",
      "CellB:Treated-Control",
      "CellB_Treated-CellA_Control",
      "CellA-CellB:Treated-Control",
      "CellB-CellA:Treated-Control",
      "CellA-CellB:Treated-Control")
   # test output
   testthat::expect_equal(
      comp2contrast(expected_comps_2fac),
      contrast_names_2fac_adj)
   testthat::expect_equal(
      comp2contrast(expected_comps_2fac,
         factor_order=list(1:2, 1:2, 1:2, 2:1, 2:1, 1:2)),
      contrast_names_2fac)
})


testthat::test_that("comp2contrast two factors mixed order", {
   contrast_names_2fac_mixed <- c(
      "(CellA_Treated-CellA_Control)-(CellB_Control-CellB_Treated)"
   )
   expected_comps_2fac_mixed <- c(
      "(CellA:Treated-Control)-(CellB:Control-Treated)"
   )
   # test output
   testthat::expect_equal(
      contrast2comp(contrast_names_2fac_mixed),
      expected_comps_2fac_mixed)
   testthat::expect_equal(
      comp2contrast(expected_comps_2fac_mixed),
      contrast_names_2fac_mixed)
})


testthat::test_that("comp2contrast two factors mixed order", {
   contrast_names_3fac <- c(
      paste0("((CellA_Treated_Mut-CellB_Treated_Mut)-",
         "(CellA_Control_Mut-CellB_Control_Mut))-",
         "((CellA_Treated_WT-CellB_Treated_WT)-",
         "(CellA_Control_WT-CellB_Control_WT))"),
      paste0("((CellA_Treated_Mut-CellA_Control_Mut)-",
         "(CellB_Treated_Mut-CellB_Control_Mut))-",
         "((CellA_Treated_WT-CellA_Control_WT)-",
         "(CellB_Treated_WT-CellB_Control_WT))"))
   expected_comps_3fac <- c(
      "CellA-CellB:Treated-Control:Mut-WT",
      "CellA-CellB:Treated-Control:Mut-WT"
   )
   # test output
   testthat::expect_equal(
      contrast2comp(contrast_names_3fac),
      expected_comps_3fac)
   testthat::expect_equal(
      comp2contrast(expected_comps_3fac,
         factor_order=list(c(1, 2, 3), c(2, 1, 3))),
      contrast_names_3fac)
   testthat::expect_equal(
      comp2contrast(expected_comps_3fac),
      contrast_names_3fac[c(1, 1)])
})
