data(sp_toys)
d <- sp_toys

test_that("Internal spatial metrics work", {
  sm <- getSpatialInternalMetrics(labels=d$label, location=d[,seq_len(2)], 
                                  level="dataset")
  expect_true(all(!is.na(as.matrix(sm))))
  sm <- getSpatialInternalMetrics(labels=d$label, location=d[,seq_len(2)], 
                                  level="class")
  sm2 <- getSpatialInternalMetrics(labels=d$label, location=d[,seq_len(2)], 
                                   level="element",
                                   metrics=c("PAS","ELSA"))
  medElsa <- sapply(split(sm2$ELSA, sm2$PAS), median)
  expect_true(all(medElsa[1]<0.2 & medElsa[2]>0.5))
})

test_that("External spatial metrics work", {
  sm <- getSpatialExternalMetrics(true=d$label, pred=d$p1, 
                                  location=d[,seq_len(2)], level="dataset")
  expect_true(all(!is.na(as.matrix(sm))))
  sm <- getSpatialExternalMetrics(true=d$label, pred=d$p1, 
                                  location=d[,seq_len(2)], level="class")
  expect_true(sm$nsAWH[3]==1)
  expect_true(sm$nsAWC[2]==1)
  sm2 <- getSpatialExternalMetrics(true=d$label, pred=d$p1, 
                                   location=d[,seq_len(2)], level="element",
                                   metrics=c("nsSPC","NPC"))
  sm2 <- getSpatialExternalMetrics(true=d$label, pred=d$p1, 
                                   location=d[,seq_len(2)], level="element", 
                                   metrics=c("nsSPC"), useNegatives=FALSE)
  medSPC <- sapply(split(sm2$nsSPC, d$label!=d$p1), median)
  expect_true(all(medSPC[1]>0.8 & medSPC[2]<0.3))
  
  sa <- spatialARI(d$label, d$p2, d[,1:2], original=TRUE)
  expect_true(round(sa,4)==c(0.9235, 0.7437))
})
