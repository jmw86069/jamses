
context("point_handedness")

testthat::test_that("point_handedness vector input", {
   testthat::expect_equal(
      point_handedness(pt1=c(1, 1), pt2=c(2, 2), angle=0),
      "right")

   testthat::expect_equal(
      point_handedness(pt1=c(0, 4, 6, 9), pt2=c(2, 9, 4, 4), angle=0),
      c("right", "left"))
})

testthat::test_that("point_handedness angle test", {
   # test points
   pt1 <- matrix(ncol=2, rep(1, 16))
   pt2 <- pt1;
   pt2[,1] <- c(0, 1, 2, 0, 2, 0, 1, 2)
   pt2[,2] <- c(2, 2, 2, 1, 1, 0, 0, 0);

   # pointing right
   testthat::expect_equal(
      point_handedness(pt1=pt1, pt2=pt2, angle=0),
      rep(c("right", "left"), c(5, 3)))
   testthat::expect_equal(
      point_handedness(pt1=pt2, pt2=pt1, angle=0),
      rep(c("left", "right"), c(3, 5)))

   # pointing up
   testthat::expect_equal(
      point_handedness(pt1=pt1, pt2=pt2, angle=90),
      c("left", "right")[c(2, 2, 1, 2, 1, 2, 2, 1)])
   testthat::expect_equal(
      point_handedness(pt1=pt2, pt2=pt1, angle=90),
      c("left", "right")[c(1, 2, 2, 1, 2, 1, 2, 2)])

   # pointing left
   testthat::expect_equal(
      point_handedness(pt1=pt1, pt2=pt2, angle=180),
      c("left", "right")[c(1, 1, 1, 2, 2, 2, 2, 2)])
   testthat::expect_equal(
      point_handedness(pt1=pt2, pt2=pt1, angle=180),
      c("left", "right")[c(2, 2, 2, 1, 2, 1, 1, 1)])

   # pointing down
   testthat::expect_equal(
      point_handedness(pt1=pt1, pt2=pt2, angle=270),
      point_handedness(pt1=pt1, pt2=pt2, angle=-90))
   testthat::expect_equal(
      point_handedness(pt1=pt1, pt2=pt2, angle=270),
      c("left", "right")[c(1, 1, 2, 1, 2, 1, 1, 2)])

})

point_handedness(pt1=c(1, 1), pt2=c(2, 2), angle=0)
