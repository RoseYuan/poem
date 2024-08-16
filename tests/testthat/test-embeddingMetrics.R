# generate mock data:
set.seed(123)
d1 <- mockData()
d2 <- mockData(classDiff=4)

test_that("Embedding metrics work", {
  m1 <- getClassEmbeddingMetrics(d1[,1:2], d1$class)
  m2 <- getGlobalEmbeddingMetrics(d1[,1:2], d1$class)
  m3 <- getGlobalEmbeddingMetrics(d2[,1:2], d2$class)
  expect_false(any(is.na(as.matrix(m1)) | is.infinite(as.matrix(m1))))
  expect_false(any(is.na(as.matrix(m2)) | is.infinite(as.matrix(m2))))
  expect_true(all(m1$pnSW>=0 & m1$pnSW<=1))
  expect_true(m2$meanSW<m3$meanSW & m2$meanClassSW<m3$meanClassSW & 
                m2$pnSW>m3$pnSW & m2$cdbw<m3$cdbw)
})
