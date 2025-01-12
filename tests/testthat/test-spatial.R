data(sp_toys)
d <- sp_toys

test_that("Internal spatial metrics work", {
  sm <- getSpatialInternalMetrics(d$label, d[,seq_len(2)], level="dataset")
  expect_true(all(!is.na(as.matrix(sm))))
  sm <- getSpatialInternalMetrics(d$label, d[,seq_len(2)], level="class")
  sm2 <- getSpatialInternalMetrics(d$label, d[,seq_len(2)], level="element",
                                   metrics=c("PAS","ELSA"))
  medElsa <- sapply(split(sm2$ELSA, sm2$PAS), median)
  expect_true(all(medElsa[1]<0.2 & medElsa[2]>0.5))
})

test_that("External spatial metrics work", {
  sm <- getSpatialExternalMetrics(d$label, d$p1, d[,seq_len(2)], level="dataset")
  expect_true(all(!is.na(as.matrix(sm))))
  sm <- getSpatialExternalMetrics(d$label, d$p1, d[,seq_len(2)], level="class")
  expect_true(sm$SpatialAWH[3]==1)
  expect_true(sm$SpatialAWC[2]==1)
  sm2 <- getSpatialExternalMetrics(d$label, d$p1, d[,seq_len(2)], level="element",
                                   metrics=c("SpatialSPC","SpatialNPC"))
  sm2 <- getSpatialExternalMetrics(d$label, d$p1, d[,seq_len(2)], level="element", 
                                   metrics=c("SpatialSPC"), useNegatives=FALSE)
  medSPC <- sapply(split(sm2$SpatialSPC, d$label!=d$p1), median)
  expect_true(all(medSPC[1]>0.8 & medSPC[2]<0.3))
})
