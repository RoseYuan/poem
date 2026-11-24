# NCE = log2((1 + same-class neighbours) / (1 + degree * class frequency))

test_that("NCE on ragged neighbour lists matches hand computation", {
  labels <- factor(c("a","a","a","b"))   # freq: a = 3/4, b = 1/4
  nn <- list(2:3, c(1L,4L), 1L, integer(0))
  exp <- log2(c((1+2)/(1+2*3/4),   # node 1: 2 a-neighbours of 2
                (1+1)/(1+2*3/4),   # node 2: 1 of 2
                (1+1)/(1+1*3/4),   # node 3: 1 of 1
                1))                # node 4: isolated -> log2(1/1) = 0
  expect_equal(unname(.nlog2Enrichment(nn, labels)), exp)
})

test_that("ragged and fixed-k NCE agree when every node has k neighbours", {
  set.seed(1)
  x <- matrix(rnorm(200*2), ncol=2)
  labels <- factor(sample(c("a","b","c"), 200, replace=TRUE))
  knn <- emb2knn(x, k=5)
  knn$nncl <- matrix(labels[knn$index], nrow=nrow(knn$index))
  ragged <- lapply(seq_len(nrow(knn$index)), function(i) knn$index[i,])
  expect_equal(.nlog2Enrichment(ragged, labels),
               .nlog2Enrichment(knn, labels))
})

test_that("NCE handles integer labels that are not 1..C", {
  labels <- c(2L, 2L, 9L)
  nn <- list(2L, 1L, 1L)
  expect_false(anyNA(.nlog2Enrichment(nn, labels)))
})

test_that("NCE on an SNN graph is ~0 for shuffled labels", {
  set.seed(1)
  x <- matrix(rnorm(300*2), ncol=2)
  labels <- sample(c("a","b","c"), 300, replace=TRUE)
  m <- getGraphMetrics(x, labels, k=10, shared=TRUE, level="dataset",
                       metrics="NCE")
  expect_lt(abs(m$NCE), 0.2)
})
