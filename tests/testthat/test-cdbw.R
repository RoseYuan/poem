set.seed(1)
d1 <- mockData()
x <- as.matrix(d1[, seq_len(2)])
lab <- as.integer(d1$class)

test_that("CDbw does not depend on label order", {
  relab <- max(lab) + 1L - lab
  expect_equal(CDbw(x, lab), CDbw(x, relab))
})

test_that("CDbw densities are scale-free", {
  # densities count points within one stdev, so they don't change with scale;
  # compactness carries 1/stdev and separation a distance
  a <- CDbw(x, lab)
  b <- CDbw(10 * x, lab)
  expect_equal(b[["compactness"]], a[["compactness"]] / 10)
  expect_equal(b[["sep"]], a[["sep"]] * 10)
})
