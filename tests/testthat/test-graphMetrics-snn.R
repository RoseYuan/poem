# Tests for SNN graphs with potential isolated nodes

set.seed(42)
d1 <- mockData()


test_that("Class-level graph metrics with SNN handle NAs by averaging with na.rm", {
  # Use small k to increase chance of isolates in SNN
  m <- getGraphMetrics(d1[,seq_len(2)], d1$class, k=3, level="class", shared=TRUE,
                       metrics=c("SI","NP","NCE","PWC"))
  expect_false(any(is.na(as.matrix(m)) | is.infinite(as.matrix(m))))
  expect_true(all(m$NP>=0 & m$NP<=1))
  expect_true(all(m$PWC>=0 & m$PWC<=1))
})


test_that("Dataset-level graph metrics with SNN aggregate with na.rm", {
  m <- getGraphMetrics(d1[,seq_len(2)], d1$class, k=3, level="dataset", shared=TRUE,
                       metrics=c("SI","NP","PWC"))
  expect_equal(nrow(m),1)
  expect_false(any(is.na(as.matrix(m)) | is.infinite(as.matrix(m))))
  expect_true(m$NP>=0 && m$NP<=1)
  expect_true(m$PWC>=0 && m$PWC<=1)
})
