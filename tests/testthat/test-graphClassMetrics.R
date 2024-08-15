# generate mock data:
set.seed(123)
d1 <- mockData()
set.seed(123)
d2 <- mockData(spread=c(5,1))
set.seed(123)
d3 <- mockData(classDiff=0.5)

test_that("Graph class metrics from data.frame works", {
  m1 <- getGraphClassMetrics(d1[,1:2], d1$class, k=5)
  expect_false(any(is.na(as.matrix(m1)) | is.infinite(as.matrix(m1))))
  expect_true(all(m1$NP>=0 & m1$NP<=1))
  expect_true(all(m1$PWC>=0 & m1$PWC<=1))
  expect_true(all(m1$AMSP>0))
})

test_that("Graph class metrics from knn works", {
  k1 <- .emb2knn(as.matrix(d1[,1:2]), k=5)
  m1b <- getGraphClassMetrics(k1, d1$class)
  expect_false(any(is.na(as.matrix(m1b)) | is.infinite(as.matrix(m1b))))
})

test_that("Graph class metrics from graph works", {
  k1 <- .emb2knn(as.matrix(d1[,1:2]), k=5)
  g <- bluster::neighborsToKNNGraph(k1$index)
  m1c <- getGraphClassMetrics(g, d1$class)
  expect_false(any(is.na(as.matrix(m1c)) | is.infinite(as.matrix(m1c))))
})

test_that("Relative graph class metrics are consistent with expectations", {
  m1 <- getGraphClassMetrics(d1[,1:2], d1$class, k=5)
  m2 <- getGraphClassMetrics(d2[,1:2], d2$class, k=5)
  m3 <- getGraphClassMetrics(d3[,1:2], d3$class, k=5)
  expect_gte(m2$AMSP[1], m1$AMSP[1])
  expect_true(all(m3$SI<=m1$SI))
  expect_true(all(m3$PWC>=m1$PWC))
  expect_true(all(m3$NP<=m1$NP))
  expect_true(all(m3$NCE<=m1$NCE))
})

