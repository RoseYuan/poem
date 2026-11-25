set.seed(42)
d1 <- mockData()
x <- as.matrix(d1[, seq_len(2)])
y <- as.integer(d1$class)

test_that("dbcv counts noise in the denominator", {
  yn <- y
  yn[1:20] <- -1L
  with_noise <- dbcv(x, yn)$dbcv
  without <- dbcv(x[-(1:20), ], y[-(1:20)])$dbcv
  expect_equal(with_noise, without * (length(y) - 20) / length(y))
})

test_that("class-level dbcv gives NA for singleton classes", {
  lab <- factor(as.character(d1$class), levels = c(levels(d1$class), "Z"))
  lab[1] <- "Z"
  m <- getEmbeddingMetrics(x, lab, metrics = "dbcv", level = "class")
  expect_equal(nrow(m), nlevels(lab))
  expect_true(is.na(m$dbcv[m$class == "Z"]))
  ref <- getEmbeddingMetrics(x[-1, ], droplevels(lab[-1]), metrics = "dbcv",
                             level = "class")
  expect_equal(m$dbcv[m$class != "Z"], ref$dbcv)
})

test_that("non-euclidean dbcv uses the between-cluster block", {
  # in 1-D, manhattan and euclidean distances coincide
  x1 <- x[, 1, drop = FALSE]
  expect_equal(dbcv(x1, y, distance = "manhattan"),
               dbcv(x1, y, distance = "euclidean"))
})

test_that("dbcv handles clusters too small to have internal MST nodes", {
  y2 <- y
  y2[1:2] <- 99L
  expect_no_warning(r <- dbcv(x, y2))
  expect_true(all(is.finite(r$vcs)))
  y3 <- y
  y3[1:3] <- 99L
  expect_true(all(is.finite(dbcv(x, y3)$vcs)))
})
