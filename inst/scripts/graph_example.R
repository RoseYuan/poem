## code to prepare `toyExamples` dataset goes here

set.seed(42)
sd <- 1
d4 <- data.frame(x=c(rnorm(20, -2, sd=sd), rnorm(20, -2, sd=sd), rnorm(20, 2, sd=sd), rnorm(20, 2, sd=sd)),
                 y=c(rnorm(20, -2, sd=sd), rnorm(20, 2, sd=sd), rnorm(20, 2, sd=sd), rnorm(20, -2, sd=sd)),
                 class=paste0("class",rep(1:4,each=20)))
sd <- 1.1
d5 <- data.frame(x=c(rnorm(20, -2, sd=sd), rnorm(20, -2, sd=sd), rnorm(20, 2, sd=sd), rnorm(20, 2, sd=sd)),
                 y=c(rnorm(20, -2, sd=sd), rnorm(20, -2, sd=sd), rnorm(20, 2, sd=sd), rnorm(20, 2, sd=sd)),
                 class=paste0("class",rep(1:4,each=20)))
f <- 1.1
sd <- 1.3
d6 <- data.frame(x=c(rnorm(20, -f, sd=sd), rnorm(20, -f, sd=sd), rnorm(20, f, sd=sd), rnorm(20, f, sd=sd)),
                 y=c(rnorm(20, -f, sd=sd), rnorm(20, f, sd=sd), rnorm(20, -f, sd=sd), rnorm(20, f, sd=sd)),
                 class=paste0("class",rep(1:4,each=20)))


nMult <- 1
sd <- 0.5

set.seed(10)
d0 <- data.frame(x=c(rnorm(nMult*40, -0.4, sd=1.1*sd), rnorm(nMult*20, 0, sd=1.1*sd), rnorm(nMult*40, 0.4, sd=1.1*sd),
                     rnorm(nMult*60,0,sd=sd*0.8)),
                 y=c(rnorm(nMult*20, 0.5, sd=sd), rnorm(nMult*20, 1, sd=sd), rnorm(nMult*20, 1.5, sd=sd),rnorm(nMult*20, 1, sd=sd),rnorm(nMult*20, 0, sd=sd), rnorm(nMult*60,-2.35,sd=sd)),
                 class=paste0("class",rep(1:2,nMult*c(100,60))))

set.seed(10)
d1 <- data.frame(x=c(rnorm(nMult*20, -2, sd=sd), rnorm(nMult*20, -1, sd=sd), rnorm(nMult*20, 0, sd=sd),rnorm(nMult*20, 1, sd=sd),rnorm(nMult*20, 2, sd=sd), rnorm(nMult*60,0,sd=sd*0.8)),
                 y=c(rnorm(nMult*20, 0.2, sd=sd), rnorm(nMult*20, 1, sd=sd), rnorm(nMult*20, 2, sd=sd),rnorm(nMult*20, 1, sd=sd),rnorm(nMult*20, 0, sd=sd), rnorm(nMult*60,-2,sd=sd)),
                 class=paste0("class",rep(1:2,nMult*c(100,60))))


set.seed(10)
d2 <- data.frame(x=c(rnorm(20, -2, sd=sd), rnorm(20, -1.25, sd=sd), rnorm(20, 0.25, sd=sd),rnorm(20, 1, sd=sd),rnorm(20, 2, sd=sd), rnorm(60,0,sd=sd*0.8)),
                 y=c(rnorm(20, 0.2, sd=sd), rnorm(20, 0.5, sd=sd), rnorm(20, 2, sd=sd),rnorm(20, 1, sd=sd),rnorm(20, 0, sd=sd), rnorm(60,-2,sd=sd)),
                 class=paste0("class",rep(1:2,c(100,60))))

set.seed(10)
d3 <- data.frame(x=c(rnorm(20, -2, sd=sd), rnorm(20, -1.25, sd=sd), rnorm(20, 0.25, sd=sd),rnorm(20, 1, sd=sd),rnorm(20, 2.1, sd=0.75*sd), rnorm(60,0,sd=sd*0.8)),
                 y=c(rnorm(20, 0.2, sd=sd), rnorm(20, 0.5, sd=sd), rnorm(20, 2, sd=sd),rnorm(20, 1.3, sd=0.75*sd),rnorm(20, -0.5, sd=0.75*sd), rnorm(60,-2,sd=sd)),
                 class=paste0("class",rep(1:2,c(100,60))))

toyExamples <- dplyr::bind_rows(list(graph1=d4, graph2=d5, graph3=d6, graph4=d0, graph5=d1, graph6=d2, graph7=d3), .id="graph")
usethis::use_data(toyExamples, overwrite = TRUE)