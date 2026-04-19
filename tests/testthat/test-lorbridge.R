## tests/testthat/test-ccm_helpers.R

test_that("lor_ci_2x2 returns correct structure", {
  res <- lor_ci_2x2(30, 25, 28, 24)
  expect_named(res, c("lor", "se", "lo", "hi"))
  expect_true(res$lo < res$lor)
  expect_true(res$lor < res$hi)
})

test_that("lor_ci_2x2 applies Haldane-Anscombe correction for zero cells", {
  # Should not throw; 0.5 correction applied
  expect_no_error(lor_ci_2x2(0, 5, 10, 20))
})

test_that("ccm_row returns correct column names", {
  lc  <- lor_ci_2x2(30, 25, 28, 24)
  out <- ccm_row(exp(lc$lor), exp(lc$lo), exp(lc$hi),
                 lc$lor, lc$lo, lc$hi)
  expect_s3_class(out, "data.frame")
  expect_true(all(c("OR","LOR","YuleQ","YuleY","r_meta") %in% names(out)))
})

test_that("YuleQ is bounded [-1, 1]", {
  lc  <- lor_ci_2x2(100, 5, 5, 100)
  out <- ccm_row(exp(lc$lor), exp(lc$lo), exp(lc$hi),
                 lc$lor, lc$lo, lc$hi)
  expect_true(abs(out$YuleQ) <= 1)
})

test_that("r_meta is bounded [-1, 1]", {
  lc  <- lor_ci_2x2(100, 5, 5, 100)
  out <- ccm_row(exp(lc$lor), exp(lc$lo), exp(lc$hi),
                 lc$lor, lc$lo, lc$hi)
  expect_true(abs(out$r_meta) <= 1)
})


## tests/testthat/test-blr_functions.R

test_that("blr_continuous returns required list elements", {
  set.seed(1)
  y <- rbinom(100, 1, 0.4)
  x <- rnorm(100)
  res <- blr_continuous(y, x)
  expect_named(res, c("model", "summary_table", "predictor_mean", "predictor_sd"))
  expect_s3_class(res$model, "glm")
  expect_s3_class(res$summary_table, "data.frame")
})

test_that("blr_categorical returns required list elements", {
  set.seed(1)
  y   <- rbinom(60, 1, 0.4)
  grp <- rep(c("A","B","C"), each = 20)
  res <- blr_categorical(y, grp, ref_level = "A")
  expect_named(res, c("model", "results", "tab_2xJ"))
  expect_equal(nrow(res$results), 2)   # B and C vs A
})

test_that("cosine_theta_2row returns named vector of length J-1", {
  tab <- matrix(c(10,20,30,40,15,25), nrow = 2,
                dimnames = list(c("Maj","Min"), c("A","B","C")))
  res <- cosine_theta_2row(tab, ref_col = "B")
  expect_length(res, 2)
  expect_named(res, c("A","C"))
  expect_true(all(abs(res) <= 1 | is.na(res)))
})


## tests/testthat/test-sonsca_functions.R

test_that("inertia_pct sums to 100", {
  tab <- matrix(c(30,42,32,25,9,2,
                  28,38,38,24,10,9,
                   2, 7,10,12,11,11,
                  15,51,121,172,125,76),
                nrow = 4, byrow = TRUE)
  pct <- inertia_pct(tab)
  expect_equal(round(sum(pct), 6), 100)
})

test_that("sonsca_cosines returns a matrix", {
  tab <- matrix(c(30,42,32,25,9,2,
                  28,38,38,24,10,9,
                   2, 7,10,12,11,11,
                  15,51,121,172,125,76),
                nrow = 4, byrow = TRUE,
                dimnames = list(paste0("R",1:4), paste0("C",1:6)))
  sc  <- sonsca_coords(tab)
  cos <- sonsca_cosines(sc$F, sc$G, row_anchor = "R2", col_anchor = "C4")
  expect_true(is.matrix(cos))
})
