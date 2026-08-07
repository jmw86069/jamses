group_levels <- c("one", "two", "three", "four")
mm2 <- model.matrix(~0 + factor(rep(group_levels, each=3), levels=group_levels))
colnames(mm2) <- group_levels
rownames(mm2) <- paste0("sample", 1:12)
icontrasts2 <- limma::makeContrasts(
   contrasts=c("two-one", "four-three", "(four-three)-(two-one)"),
   levels=mm2)

# multi-factor design, group names encode two underlying factors: A/B, X/Y
mf_levels <- c("A_X", "A_Y", "B_X", "B_Y")
mm3 <- model.matrix(~0 + factor(rep(mf_levels, each=3), levels=mf_levels))
colnames(mm3) <- mf_levels
rownames(mm3) <- paste0("s", 1:12)

test_that("SEDesign construction works", {
   sd_design_only <- SEDesign(design=mm2)
   expect_s3_class(sd_design_only, "SEDesign")
   expect_equal(samples(sd_design_only), rownames(mm2))
   expect_equal(groups(sd_design_only), colnames(mm2))
   # group names have no "_", so there is one underlying factor,
   # auto-labeled "factor1"
   expect_equal(factors(sd_design_only), "factor1")

   sd_full <- SEDesign(design=mm2, contrasts=icontrasts2)
   expect_equal(ncol(contrasts(sd_full)), 3);
   expect_equal(contrastnames(sd_full), colnames(icontrasts2));
   expect_equal(contrast_names(sd_full), colnames(icontrasts2));
})

test_that("SEDesign validity is enforced", {
   expect_error(SEDesign(contrasts=icontrasts2));
})

test_that("design_df and contrasts_df caches are populated", {
   sd_full <- SEDesign(design=mm2, contrasts=icontrasts2)
   expect_equal(nrow(sd_full@design_df), 4);
   expect_equal(nrow(sd_full@contrasts_df), 3);

   sd_mf <- SEDesign(design=mm3)
   expect_equal(colnames(sd_mf@design_df), c("factor1", "factor2"));
   expect_equal(sd_mf@design_df[["factor1"]], c("A", "A", "B", "B"));
   expect_equal(sd_mf@design_df[["factor2"]], c("X", "Y", "X", "Y"));
})

test_that("[ subsetting works and updates caches", {
   sd_full <- SEDesign(design=mm2, contrasts=icontrasts2)
   sd_sub <- sd_full[, c("one", "two")]
   expect_equal(groups(sd_sub), c("one", "two"));
   expect_equal(ncol(contrasts(sd_sub)), 1);
   expect_equal(nrow(sd_sub@design_df), 2);

   sd_sub_samples <- sd_full[paste0("sample", 4:12), ];
   expect_equal(samples(sd_sub_samples), paste0("sample", 4:12));
   # "one" group loses all its samples and should be dropped
   expect_false("one" %in% groups(sd_sub_samples));
})

test_that("factors<- renames design_df columns only, with no other effect", {
   sd_mf <- SEDesign(design=mm3)
   factors(sd_mf) <- c("Genotype", "Treatment")
   expect_equal(factors(sd_mf), c("Genotype", "Treatment"));
   expect_equal(colnames(sd_mf@design_df), c("Genotype", "Treatment"));
   # groups(), design(), contrasts() are unaffected by factors<-
   expect_equal(groups(sd_mf), mf_levels);
   expect_equal(colnames(design(sd_mf)), mf_levels);

   # length must match the number of underlying factors
   expect_error(factors(sd_mf) <- c("OnlyOne"));
})

test_that("groups<- renames design groups and propagates to contrasts", {
   sd_full <- SEDesign(design=mm2, contrasts=icontrasts2)
   groups(sd_full) <- c("A", "B", "C", "D")
   expect_equal(groups(sd_full), c("A", "B", "C", "D"));
   expect_equal(colnames(design(sd_full)), c("A", "B", "C", "D"));
   expect_equal(rownames(contrasts(sd_full)), c("A", "B", "C", "D"));
})

test_that("design<- and contrasts<- reject mismatched group names", {
   sd_full <- SEDesign(design=mm2, contrasts=icontrasts2)
   bad_design <- matrix(1:8, nrow=2,
      dimnames=list(c("s1", "s2"), c("X", "Y", "Z", "W")));
   expect_error(design(sd_full) <- bad_design);

   bad_contrasts <- matrix(1, nrow=1, dimnames=list("bad", "c1"));
   expect_error(contrasts(sd_full) <- bad_contrasts);
})

test_that("samples<- updates samples and design rownames", {
   sd_full <- SEDesign(design=mm2, contrasts=icontrasts2)
   new_samples <- paste0("s", 1:12);
   samples(sd_full) <- new_samples;
   expect_equal(samples(sd_full), new_samples);
   expect_equal(rownames(design(sd_full)), new_samples);
})

test_that("print.SEDesign runs without error", {
   sd_full <- SEDesign(design=mm2, contrasts=icontrasts2)
   expect_output(print(sd_full), "SEDesign");
})
