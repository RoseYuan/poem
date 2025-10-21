# generate mock data:
set.seed(123)
em <- mockData()
kn <- emb2knn(as.matrix(em[,1:2]), k=2)
g <- .nn2graph(kn)
cl <- sample(em[,3],nrow(em))

levels <- c("element", "class", "dataset")

test_that("Metrics args are passed correctly with all inputs to partition fns", {
  for(l in levels){
    met <- ifelse(l=="element", "SPC", "WH")
    o <- getPartitionMetrics(true=em[,3], pred=cl, level=l, metrics=met)
    expect_true(met %in% colnames(o))
    expect_true(length(setdiff(colnames(o), c("class","cluster",met)))==0)
  }
})

test_that("Metrics args are passed correctly with all inputs to graph fns", {
  for(i in list(em=em[,1:2], kn=kn, g)){
    for(l in levels){
      o <- getGraphMetrics(i, labels=em[,3], level=l, metrics=c("NP"), k=2)
      expect_true("NP" %in% colnames(o))
      expect_true(length(setdiff(colnames(o), c("class","NP")))==0)
    }
  }
})

test_that("Metrics args are passed correctly with all inputs to embed fns", {
  for(l in levels){
    met <- ifelse(l=="element", "SW", "meanSW")
    o <- getEmbeddingMetrics(em[,1:2], labels=em[,3], level=l, metrics=met)
    expect_true(met %in% colnames(o))
    expect_true(length(setdiff(colnames(o), c("class",met)))==0)
  }
})

data(sp_toys)
checks <- list( c(level="element", metric="NPC"),
                c(level="class", metric="nsWH"),
                c(level="dataset", metric="nsAccuracy"))

test_that("Metrics args are passed correctly with all inputs to spatial fns", {
  for(l in checks){
    o <- getSpatialExternalMetrics(true=sp_toys$label, pred=sp_toys$p1,
                                   metrics=l["metric"], level=l["level"],
                                   location=sp_toys[,c("x", "y")], k=3)
    expect_true(l["metric"] %in% colnames(o))
    expect_true(length(setdiff(colnames(o), c("class", "cluster", l["metric"])))==0)
  }
})


set.seed(123)
d2 <- mockData(spread=c(5,1))
set.seed(123)
d3 <- mockData(classDiff=0.5)

test_that("Graph class metrics from data.frame works", {
  m1 <- getGraphMetrics(em[,seq_len(2)], em$class, k=5, level="class")
  expect_false(any(is.na(as.matrix(m1)) | is.infinite(as.matrix(m1))))
  expect_true(all(m1$NP>=0 & m1$NP<=1))
  expect_true(all(m1$PWC>=0 & m1$PWC<=1))
  expect_true(all(m1$AMSP>0))
})

test_that("Graph class metrics from knn works", {
  k1 <- emb2knn(as.matrix(em[,seq_len(2)]), k=5)
  m1b <- getGraphMetrics(k1, em$class, level="class", k=5)
  expect_false(any(is.na(as.matrix(m1b)) | is.infinite(as.matrix(m1b))))
  m1e <- getGraphMetrics(k1, em$class, level="element", k=5,
                         metrics=c("SI", "NP"))
})

test_that("Graph class metrics from graph works", {
  k1 <- emb2knn(as.matrix(em[,seq_len(2)]), k=5)
  g <- bluster::neighborsToKNNGraph(k1$index)
  m1c <- getGraphMetrics(g, em$class, level="class", k=5)
  expect_false(any(is.na(as.matrix(m1c)) | is.infinite(as.matrix(m1c))))
  m1e <- getGraphMetrics(g, em$class, level="element", k=5,
                         metrics=c("SI", "NP"))
})

test_that("Relative graph class metrics are consistent with expectations", {
  m1 <- getGraphMetrics(em[,seq_len(2)], em$class, k=5, level="class")
  m2 <- getGraphMetrics(d2[,seq_len(2)], d2$class, k=5, level="class")
  m3 <- getGraphMetrics(d3[,seq_len(2)], d3$class, k=5, level="class")
  expect_gte(m2$AMSP[1], m1$AMSP[1])
  expect_true(all(m3$SI<=m1$SI))
  expect_true(all(m3$PWC>=m1$PWC))
  expect_true(all(m3$NP<=m1$NP))
  expect_true(all(m3$NCE<=m1$NCE))
})

test_that("Element-level graph metrics work", {
  m <- getGraphMetrics(em[,seq_len(2)], em$class, k=5, level="element",
                       metrics=c("SI","NP","NCE"))
  expect_equal(nrow(m),nrow(em))
  expect_true(all(m$SI>=0 & m$SI<=1))
  expect_true(all(m$NP>=0 & m$NP<=1))
})

test_that("Dataset-level graph metrics work", {
  m <- getGraphMetrics(em[,seq_len(2)], em$class, k=5, level="dataset",
                       metrics=c("SI","AMSP","PWC","cohesion"))
  expect_equal(nrow(m),1) 
  expect_equal(m$cohesion,0)
  expect_true(m$SI>0.6 & m$SI<0.9)
  expect_true(m$PWC>0.1 & m$PWC<0.3)
})
