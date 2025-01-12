# generate non-fuzzy fuzzy metrics:
Q1 <- data.frame(c(1,1,0,0,0,0),c(0,0,1,1,0,0),c(0,0,0,0,1,1))
Q2 <- cbind(Q1[,1], c(0,0,1,1,1,1))

test_that("fuzzyPartitionMetrics works and match non-fuzzy ones in a binary setting", {
  fm <- fuzzyPartitionMetrics(Q1,Q2,nperms = 1000)
  fm2 <- c(RI=fm$NDC, WC=fm$fuzzyWallace2$global, WH=fm$fuzzyWallace1$global,
           ARI=fm$ACI, AWC=fm$fuzzyAdjW2$global, AWH=fm$fuzzyAdjW1$global)
  fm <- unlist(fm)
  pm <- getPartitionMetrics(as.factor(apply(Q1,1,which.max)),
                            as.factor(apply(Q2,1,which.max)),
                            metrics=names(fm2), level="dataset")
  expect_true(all(!is.na(fm) & fm>=-1 & fm<=1))  
  expect_gt(cor(unlist(pm),fm2),expected=0.999)
})

# generate fuzzy partitions:
m1 <- matrix(c(0.95, 0.025, 0.025,
               0.98, 0.01, 0.01,
               0.96, 0.02, 0.02,
               0.95, 0.04, 0.01,
               0.95, 0.01, 0.04,
               0.99, 0.005, 0.005,
               0.025, 0.95, 0.025,
               0.97, 0.02, 0.01,
               0.025, 0.025, 0.95),
               ncol = 3, byrow=TRUE)
m2 <- matrix(c(0.95, 0.025, 0.025,
               0.98, 0.01, 0.01,
               0.96, 0.02, 0.02,
               0.025, 0.95, 0.025,
               0.02, 0.96, 0.02,
               0.01, 0.98, 0.01,
               0.05, 0.05, 0.95,
               0.02, 0.02, 0.96,
               0.01, 0.01, 0.98),
               ncol = 3, byrow=TRUE)
colnames(m1) <- colnames(m2) <- LETTERS[seq_len(3)]
# a hard truth:
hardTrue <- apply(m1,1,FUN=which.max)
# some predicted labels:
hardPred <- c(1,1,1,1,1,1,2,2,2)



test_that("getFuzzyPartitionMetrics works on a fuzzy-fuzzy setting", {
  fm <- getFuzzyPartitionMetrics(fuzzyTrue=m1,fuzzyPred=m2, level="dataset")
  fm <- getFuzzyPartitionMetrics(fuzzyTrue=m1,fuzzyPred=m2, level="class")
  expect_true(!any(is.na(fm[seq_len(3),1])))
  fm2 <- getFuzzyPartitionMetrics(fuzzyTrue=m1,fuzzyPred=m2, level="element",
                                  metrics="fuzzySPC")
  expect_equal(which(fm2$fuzzySPC<0.5), 8)
})


test_that("getFuzzyPartitionMetrics works on a fuzzy-hard setting", {
  fm <- getFuzzyPartitionMetrics(hardPred=hardPred, hardTrue=hardTrue,
                                 fuzzyTrue=m1, nperms=3, level="class")
  expect_true(all(fm$fuzzyWC[2:3]==1))
  expect_equal(fm$fuzzyWH[4],1)
})
