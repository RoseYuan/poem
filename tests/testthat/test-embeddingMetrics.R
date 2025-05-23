# generate mock data:
set.seed(123)
d1 <- mockData()
d2 <- mockData(classDiff=4)

test_that("Embedding metrics work", {
  m1 <- getEmbeddingMetrics(d1[,seq_len(2)], d1$class, level="class", 
                            shared=TRUE)
  m1 <- getEmbeddingMetrics(d1[,seq_len(2)], d1$class, level="class")
  m2 <- getEmbeddingMetrics(d1[,seq_len(2)], d1$class, level="dataset", 
                            metrics=c("meanSW", "meanClassSW", "pnSW", 
                                      "minClassSW", "cdbw", "cohesion", 
                                      "compactness", "sep", "dbcv"))
  m3 <- getEmbeddingMetrics(d2[,seq_len(2)], d2$class, level="dataset",
                            metrics=c("meanSW", "meanClassSW", "pnSW", 
                                      "minClassSW", "cdbw", "cohesion", 
                                      "compactness", "sep", "dbcv"))
  expect_false(any(is.na(as.matrix(m1)) | is.infinite(as.matrix(m1))))
  expect_false(any(is.na(as.matrix(m2)) | is.infinite(as.matrix(m2))))
  expect_true(all(m1$pnSW>=0 & m1$pnSW<=1))
  expect_true(m2$meanSW<m3$meanSW & m2$meanClassSW<m3$meanClassSW & 
                m2$pnSW>m3$pnSW & m2$cdbw<m3$cdbw)
})

test_that("Dist input works", {
  m1 <- getEmbeddingMetrics(dist(d1[,seq_len(2)]), d1$class, level="class", 
                            shared=TRUE)
  expect_true(!any(is.na(as.matrix(m1[,-1]))))
})
