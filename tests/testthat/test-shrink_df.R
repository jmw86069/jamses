
context("shrink_df")

testthat::test_that("shrink_df basic arguments", {
   testdf <- data.frame(check.names=FALSE,
      SYMBOL=rep(c("ACTB", "GAPDH", "PPIA"), c(2, 3, 1)),
      `logFC B-A`=c(1.4, 1.4, 2.3, NA, 2.5, 5.1),
      probe=paste0("probe", c(2, 1, 3:6)))
   sdf <- shrink_df(testdf, by="SYMBOL")
   # colnames should not be mangled by make.names
   testthat::expect_equal(
      colnames(sdf),
      colnames(testdf))
   # rows should have exact order
   testthat::expect_equal(
      sdf$SYMBOL,
      c("ACTB", "GAPDH", "PPIA"))
   # "probe" column should have exact content
   testthat::expect_equal(
      sdf$probe,
      c("probe1,probe2", "probe3,probe4,probe5", "probe6"))
   # "logFC B-A" column should have exact content
   testthat::expect_equal(
      sdf[["logFC B-A"]],
      c(1.4, 2.4, 5.1))

   # confirm collapse_method=1 is identical
   sdf1 <- shrink_df(testdf, by="SYMBOL", collapse_method=1)
   testthat::expect_identical(
      sdf, sdf1)

   # include_num_reps adds "num_reps" column
   sdf1nr <- shrink_df(testdf, by="SYMBOL", include_num_reps=TRUE)
   testthat::expect_equal(
      colnames(sdf1nr),
      c(colnames(testdf), "num_reps"))
   sdf1nr2 <- shrink_df(testdf, by="SYMBOL", includeNumReps=TRUE)
   testthat::expect_equal(
      colnames(sdf1nr),
      colnames(sdf1nr2))

   sdf2 <- shrink_df(testdf[c(3, 4, 5, 6, 1, 2), ], by="SYMBOL")
   # rows should have specific order
   testthat::expect_equal(
      sdf2$SYMBOL,
      c("GAPDH", "PPIA", "ACTB"))
   # "probe" column should have specific content
   testthat::expect_equal(
      sdf2$probe,
      c("probe3,probe4,probe5", "probe6", "probe1,probe2"))

   testdf3 <- testdf;
   testdf3[1:2, "SYMBOL"] <- NA;
   sdf3 <- shrink_df(testdf3, by="SYMBOL")
   # SYMBOL should include NA
   testthat::expect_equal(
      sdf3$SYMBOL,
      c(NA, "GAPDH", "PPIA"))
   # "probe" column should have specific content
   testthat::expect_equal(
      sdf3$probe,
      c("probe1,probe2", "probe3,probe4,probe5", "probe6"))

   sdf4 <- shrink_df(testdf3, by="SYMBOL", keep_na_groups=FALSE)
   # SYMBOL should not include NA
   testthat::expect_equal(
      sdf4$SYMBOL,
      c("GAPDH", "PPIA"))
   # "probe" column should have specific content
   testthat::expect_equal(
      sdf4$probe,
      c("probe3,probe4,probe5", "probe6"))

   # custom string_func
   sdf5 <- shrink_df(testdf, by="SYMBOL", string_func=jamba::cPasteU)
   # "probe" column should have exact content
   testthat::expect_equal(
      sdf5$probe,
      c("probe2,probe1", "probe3,probe4,probe5", "probe6"))

   # custom string_func
   sdf5b <- shrink_df(testdf, by="SYMBOL", stringShrinkFunc=jamba::cPasteU)
   # "probe" column should have exact content
   testthat::expect_equal(
      sdf5b$probe,
      c("probe2,probe1", "probe3,probe4,probe5", "probe6"))

   # custom num_func
   sdf6 <- shrink_df(testdf, by="SYMBOL",
      num_func=function(x)max(x, na.rm=TRUE))
   # "probe" column should have exact content
   testthat::expect_equal(
      sdf6[[2]],
      c(1.4, 2.5, 5.1))

   # custom num_func
   sdf6b <- shrink_df(testdf, by="SYMBOL",
      numShrinkFunc=function(x)max(x, na.rm=TRUE))
   # "probe" column should have exact content
   testthat::expect_equal(
      sdf6b[[2]],
      c(1.4, 2.5, 5.1))

})
