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

