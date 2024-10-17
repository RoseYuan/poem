# generate non-fuzzy fuzzy metrics:
Q1 <- data.frame(c(1,1,0,0,0,0),c(0,0,1,1,0,0),c(0,0,0,0,1,1))
Q2 <- cbind(Q1[,1], c(0,0,1,1,1,1))

test_that("fuzzy metrics work and match non-fuzzy ones in a binary setting", {
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
