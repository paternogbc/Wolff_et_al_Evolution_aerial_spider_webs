### Script to analyze spider web evolution
### Jonas Wolff & Gustavo Paterno
### Last update: 10/09/2018
### Script: Modelling the evolutionary dynamics of continuous traits

### START--------
### Packages-----
library(phytools)
library(geiger)
library(OUwie)
library(surface)
library(bayou)
library(sensiPhy)
rm(list = ls())

### Load data-----
d     <- read.csv(file = "data/data_spider_web.csv")
tree <- read.nexus(file = "data/infile.tre") 

### Add row names matching phylogeny:
rownames(d) <- d$species

#### Drop outgroup species and match data and phy
tree <- match_dataphy(species ~ 1, data = d, phy = tree)[[2]]
name.check(tree, d)

### 1. Centrality ###-----
# with geiger-----
centrality <- setNames(d[, "centrality"], rownames(d))
phylosig(tree, centrality, method = "lambda")
seSize = setNames(d[, "SE_centrality"], rownames(d))   # or use SE estimate seSize = 0.022
centr.bm <- fitContinuous(tree, centrality, model = "BM", SE = seSize)
centr.ou <- fitContinuous(tree, centrality, model = "OU", SE = seSize)
centr.eb <- fitContinuous(tree, centrality, model = "EB", SE = seSize)
centr.trend <- fitContinuous(tree, centrality, model = "trend", SE = seSize)
centr.lambda <- fitContinuous(tree, centrality, model = "lambda", SE = seSize)
centr.white <- fitContinuous(tree, centrality, model = "white", SE = seSize)

centr.bmAICC <- centr.bm$opt$aicc
centr.ouAICC <- centr.ou$opt$aicc
centr.ebAICC <- centr.eb$opt$aicc
centr.trendAICC <- centr.trend$opt$aicc
centr.lambdaAICC <- centr.lambda$opt$aicc
centr.whiteAICC <- centr.white$opt$aicc

aicc <- c(centr.bmAICC, centr.ouAICC, centr.ebAICC, centr.trendAICC, centr.whiteAICC)
aiccD <- aicc - min(aicc)
aw <- exp(-0.5 * aiccD)
aiccW <- aw/sum(aw)
aiccW

# OUwie-----
# cribellar and ecribellar as regimes
centrouwie1 <- read.csv(file = "data/centrouwie1.csv", header = T)
x <- setNames(centrouwie1[, "reg"], centrouwie1[, "species"])

mtree <- make.simmap(tree, x, model = matrix(c(0, 0, 1, 0), nrow = 2, ncol = 2 , byrow = TRUE), nsim = 1) 
plot(mtree)
OUwie(mtree, centrouwie1, model=c("BM1"), simmap.tree=T, root.station=TRUE, mserr="known")
OUwie(mtree, centrouwie1, model=c("BMS"), simmap.tree=T, root.station=TRUE, mserr="known")
OUwie(mtree, centrouwie1, model=c("OU1"), simmap.tree=T, root.station=TRUE, mserr="known")
OUwie(mtree, centrouwie1, model=c("OUM"), simmap.tree=T, root.station=TRUE, mserr="known")
OUwie(mtree, centrouwie1, model=c("OUMV"), simmap.tree=T, root.station=TRUE, mserr="known")
OUwie(mtree, centrouwie1, model=c("OUMA"), simmap.tree=T, root.station=TRUE, mserr="known")
OUwie(mtree, centrouwie1, model=c("OUMVA"), simmap.tree=T, root.station=TRUE, mserr="known")

# Web type (aerial web) as regime
centrouwie5 <- read.csv(file = "data/centrouwie5.csv")
x <- setNames(centrouwie5[, "reg"], centrouwie5[, "species"])

mtree <- make.simmap(tree, x, model = "ARD", nsim = 1) 
plot(mtree)
OUwie(mtree, centrouwie5, model=c("BM1"), simmap.tree=T, root.station=TRUE, mserr="known")
OUwie(mtree, centrouwie5, model=c("BMS"), simmap.tree=T, root.station=TRUE, mserr="known")
OUwie(mtree, centrouwie5, model=c("OU1"), simmap.tree=T, root.station=TRUE, mserr="known")
OUwie(mtree, centrouwie5, model=c("OUM"), simmap.tree=T, root.station=TRUE, mserr="known")
OUwie(mtree, centrouwie5, model=c("OUMV"), simmap.tree=T, root.station=TRUE, mserr="known")
OUwie(mtree, centrouwie5, model=c("OUMA"), simmap.tree=T, root.station=TRUE, mserr="known")
OUwie(mtree, centrouwie5, model=c("OUMVA"), simmap.tree=T, root.station=TRUE, mserr="known")

# Surface------
treeo <- nameNodes(tree)
centralo <- d["centrality"]
centrSurf <- runSurface (treeo, centralo) # takes some time
surfaceSummary(centrSurf$fwd, centrSurf$bwd)
kk <-length(centrSurf$bwd)
surfaceTreePlot(treeo, centrSurf$bwd[[kk]], labelshifts = F, show.tip.label = T, edge.width = 2, cex=0.5)
getAIC(L=205.962, np=19, n=104, AICc = TRUE)

olist <- convertTreeData (treeo, centralo)
otree <- olist[[1]] 
odata <- olist[[2]]
bm<-startingModel(otree, odata, brownian=TRUE)
ou1<-startingModel(otree, odata)
surfaceAICPlot(centrSurf$fwd, centrSurf$bwd)
abline(h=bm[[1]]$aic,lty="longdash")
text(c(6,6),c(bm[[1]]$aic, ou1[[1]]$aic)-2,c("BM","OU1"),cex=0.5)

# Bayou-----
library(bayou)
dev.off()
priorOU10 <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dk="cdpois", 
                        dtheta="dunif"), param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.01), 
                        dk=list(lambda=10, kmax=150), dsb=list(bmax=1, prob=1), dtheta=list(min=0, max=1)))
dev.off()
startpars <- priorSim(priorOU10, tree, plot=TRUE)$pars[[1]]
priorOU10(startpars)
mcmcOU10 <- bayou.makeMCMC(tree, centrality, SE=seSize, prior=priorOU10, 
                         new.dir="output/OU10", outname="modelOU10", plot.freq=NULL)
mcmcOU10$run(500000)
chainOU10 <- mcmcOU10$load()
chainOU10 <- set.burnin(chainOU10, 0.3)

priorOU15 <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dk="cdpois", 
                        dtheta="dunif"), param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.01), 
                        dk=list(lambda=15, kmax=150), dsb=list(bmax=1, prob=1), dtheta=list(min=0, max=1)))
startpars <- priorSim(priorOU15, tree, plot=TRUE)$pars[[1]]
priorOU15(startpars)
mcmcOU15 <- bayou.makeMCMC(tree, centrality, SE=seSize, prior=priorOU15, 
                           new.dir="output/OU15", outname="modelOU15", plot.freq=NULL)
mcmcOU15$run(500000)
chainOU15 <- mcmcOU15$load()
chainOU15 <- set.burnin(chainOU15, 0.3)

priorOU20 <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dk="cdpois", 
                        dtheta="dunif"), param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.01), 
                        dk=list(lambda=20, kmax=150), dsb=list(bmax=1, prob=1), dtheta=list(min=0, max=1)))
startpars <- priorSim(priorOU20, tree, plot=TRUE)$pars[[1]]
priorOU20(startpars)
mcmcOU20 <- bayou.makeMCMC(tree, centrality, SE=seSize, prior=priorOU20, 
                           new.dir="output/OU20", outname="modelOU20", plot.freq=NULL)
mcmcOU20$run(500000)
chainOU20 <- mcmcOU20$load()
chainOU20 <- set.burnin(chainOU20, 0.3)

priorOU25 <- make.prior(tree, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dk="cdpois", 
                        dtheta="dunif"), param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.01), 
                        dk=list(lambda=15, kmax=150), dsb=list(bmax=1, prob=1), dtheta=list(min=0, max=1)))
startpars <- priorSim(priorOU25, tree, plot=TRUE)$pars[[1]]
priorOU25(startpars)
mcmcOU25 <- bayou.makeMCMC(tree, centrality, SE=seSize, prior=priorOU25, 
                           new.dir="output/OU25", outname="modelOU25", plot.freq=NULL)
mcmcOU25$run(500000)
chainOU25 <- mcmcOU25$load()
chainOU25 <- set.burnin(chainOU25, 0.3)

# compare quality of priors using BayesFactor
Bk <- qbeta(seq(0,1, length.out=5), 0.3,1)
ssOU10 <- mcmcOU10$steppingstone(10000, chainOU10, Bk, burnin=0.3, plot=FALSE)
ssOU15 <- mcmcOU15$steppingstone(10000, chainOU15, Bk, burnin=0.3, plot=FALSE)
ssOU20 <- mcmcOU20$steppingstone(10000, chainOU20, Bk, burnin=0.3, plot=FALSE)
ssOU25 <- mcmcOU25$steppingstone(10000, chainOU25, Bk, burnin=0.3, plot=FALSE)

mlnL <- c("OU10"=ssOU10$lnr, "OU15"=ssOU15$lnr, "OU20"=ssOU20$lnr, "QU25"=ssOU25$lnr)
mlnL

### 2. Spinning track proportions ###-----
### Crop phylogeny to match with missing values (track proportion)
dat   <- match_dataphy(track_proportions ~ 1, data = d, phy = tree)
treex <- dat$phy
dx    <- dat$data

# with geiger-----
trackprop <- setNames(dx[, "track_proportions"], rownames(dx))
phylosig(treex, trackprop, method = "lambda")
seSize = setNames(dx[, "SE_track_proportions"], rownames(dx))   # or use SE estimate seSize = 0.083
track.bm <- fitContinuous(treex, trackprop, model = "BM", SE = seSize)
track.ou <- fitContinuous(treex, trackprop, model = "OU", SE = seSize)
track.eb <- fitContinuous(treex, trackprop, model = "EB", SE = seSize)
track.trend <- fitContinuous(treex, trackprop, model = "trend", SE = seSize)
track.lambda <- fitContinuous(treex, trackprop, model = "lambda", SE = seSize)
track.white <- fitContinuous(treex, trackprop, model = "white", SE = seSize)

track.bmAICC <- track.bm$opt$aicc
track.ouAICC <- track.ou$opt$aicc
track.ebAICC <- track.eb$opt$aicc
track.trendAICC <- track.trend$opt$aicc
track.lambdaAICC <- track.lambda$opt$aicc
track.whiteAICC <- track.white$opt$aicc

aicc <- c(track.bmAICC, track.ouAICC, track.ebAICC, track.trendAICC, track.whiteAICC)
aiccD <- aicc - min(aicc)
aw <- exp(-0.5 * aiccD)
aiccW <- aw/sum(aw)
aiccW

# OUwie-----
# cribellar and ecribellar as regimes
smuouwie1 <- read.csv(file = "data/smuouwie1.csv", header = T)
x <- setNames(smuouwie1[, "reg"], smuouwie1[, "species"])

mtree <- make.simmap(treex, x, model = matrix(c(0, 0, 1, 0), nrow = 2, ncol = 2 , byrow = TRUE), nsim = 1) 
plot(mtree)
OUwie(mtree, smuouwie1, model=c("BM1"), simmap.tree=T, root.station=TRUE, mserr="known")
OUwie(mtree, smuouwie1, model=c("BMS"), simmap.tree=T, root.station=TRUE, mserr="known")
OUwie(mtree, smuouwie1, model=c("OU1"), simmap.tree=T, root.station=TRUE, mserr="known")
OUwie(mtree, smuouwie1, model=c("OUM"), simmap.tree=T, root.station=TRUE, mserr="known")
OUwie(mtree, smuouwie1, model=c("OUMV"), simmap.tree=T, root.station=TRUE, mserr="known")
OUwie(mtree, smuouwie1, model=c("OUMA"), simmap.tree=T, root.station=TRUE, mserr="known")
OUwie(mtree, smuouwie1, model=c("OUMVA"), simmap.tree=T, root.station=TRUE, mserr="known")

# Web type (aerial web) as regime
smuouwie5 <- read.csv(file = "data/smuouwie5.csv")
x <- setNames(smuouwie5[, "reg"], smuouwie5[, "species"])

mtree <- make.simmap(treex, x, model = "ARD", nsim = 1) 
plot(mtree)
OUwie(mtree, smuouwie5, model=c("BM1"), simmap.tree=T, root.station=TRUE, mserr="known")
OUwie(mtree, smuouwie5, model=c("BMS"), simmap.tree=T, root.station=TRUE, mserr="known")
OUwie(mtree, smuouwie5, model=c("OU1"), simmap.tree=T, root.station=TRUE, mserr="known")
OUwie(mtree, smuouwie5, model=c("OUM"), simmap.tree=T, root.station=TRUE, mserr="known")
OUwie(mtree, smuouwie5, model=c("OUMV"), simmap.tree=T, root.station=TRUE, mserr="known")
OUwie(mtree, smuouwie5, model=c("OUMA"), simmap.tree=T, root.station=TRUE, mserr="known")
OUwie(mtree, smuouwie5, model=c("OUMVA"), simmap.tree=T, root.station=TRUE, mserr="known")

# Surface----
treexo <- nameNodes(treex)
tracko <- dx["track_proportions"] 
trackSurf <- runSurface (treexo, tracko) # takes 10-15 min !
surfaceSummary(trackSurf$fwd, trackSurf$bwd)
kk<-length(trackSurf$bwd)
surfaceTreePlot(treexo, centrSurf$bwd[[kk]], labelshifts = F, show.tip.label = T, edge.width = 2, cex=0.5)
getAIC(L=37.9916, np=11, n=71, AICc = TRUE)

olist <- convertTreeData (treexo, tracko)
otree <- olist[[1]] 
odata <- olist[[2]]
bm<-startingModel(otree, odata, brownian=TRUE)
ou1<-startingModel(otree, odata)
surfaceAICPlot(trackSurf$fwd, trackSurf$bwd)
abline(h=bm[[1]]$aic,lty="longdash")
text(c(6,6),c(bm[[1]]$aic, ou1[[1]]$aic)-2,c("BM","OU1"),cex=0.5)

# Bayou-----
dev.off()
priorOU10 <- make.prior(treex, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dk="cdpois", 
                                         dtheta="dunif"), param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.01), 
                                                                     dk=list(lambda=10, kmax=75), dsb=list(bmax=1, prob=1), dtheta=list(min=0, max=1.5)))
dev.off()
startpars <- priorSim(priorOU10, treex, plot=TRUE)$pars[[1]]
priorOU10(startpars)
mcmcOU10 <- bayou.makeMCMC(treex, trackprop, SE=seSize, prior=priorOU10, 
                           new.dir="output/OU10t", outname="modelOU10t", plot.freq=NULL)
mcmcOU10$run(500000)
chainOU10 <- mcmcOU10$load()
chainOU10 <- set.burnin(chainOU10, 0.3)

priorOU15 <- make.prior(treex, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dk="cdpois", 
                                         dtheta="dunif"), param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.01), 
                                                                     dk=list(lambda=15, kmax=75), dsb=list(bmax=1, prob=1), dtheta=list(min=0, max=1.5)))
startpars <- priorSim(priorOU15, treex, plot=TRUE)$pars[[1]]
priorOU15(startpars)
mcmcOU15 <- bayou.makeMCMC(treex, trackprop, SE=seSize, prior=priorOU15, 
                           new.dir="output/OU15t", outname="modelOU15t", plot.freq=NULL)
mcmcOU15$run(500000)
chainOU15 <- mcmcOU15$load()
chainOU15 <- set.burnin(chainOU15, 0.3)

priorOU20 <- make.prior(treex, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dk="cdpois", 
                                         dtheta="dunif"), param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.01), 
                                                                     dk=list(lambda=20, kmax=75), dsb=list(bmax=1, prob=1), dtheta=list(min=0, max=1.5)))
startpars <- priorSim(priorOU20, treex, plot=TRUE)$pars[[1]]
priorOU20(startpars)
mcmcOU20 <- bayou.makeMCMC(treex, trackprop, SE=seSize, prior=priorOU20, 
                           new.dir="output/OU20t", outname="modelOU20t", plot.freq=NULL)
mcmcOU20$run(500000)
chainOU20 <- mcmcOU20$load()
chainOU20 <- set.burnin(chainOU20, 0.3)

priorOU25 <- make.prior(treex, dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy", dk="cdpois", 
                                         dtheta="dunif"), param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.01), 
                                                                     dk=list(lambda=15, kmax=75), dsb=list(bmax=1, prob=1), dtheta=list(min=0, max=1.5)))
startpars <- priorSim(priorOU25, treex, plot=TRUE)$pars[[1]]
priorOU25(startpars)
mcmcOU25 <- bayou.makeMCMC(treex, trackprop, SE=seSize, prior=priorOU25, 
                           new.dir="output/OU25t", outname="modelOU25t", plot.freq=NULL)
mcmcOU25$run(500000)
chainOU25 <- mcmcOU25$load()
chainOU25 <- set.burnin(chainOU25, 0.3)

# compare quality of priors using BayesFactor
Bk <- qbeta(seq(0,1, length.out=5), 0.3,1)
ssOU10 <- mcmcOU10$steppingstone(10000, chainOU10, Bk, burnin=0.3, plot=FALSE)
ssOU15 <- mcmcOU15$steppingstone(10000, chainOU15, Bk, burnin=0.3, plot=FALSE)
ssOU20 <- mcmcOU20$steppingstone(10000, chainOU20, Bk, burnin=0.3, plot=FALSE)
ssOU25 <- mcmcOU25$steppingstone(10000, chainOU25, Bk, burnin=0.3, plot=FALSE)

mlnL <- c("OU10"=ssOU10$lnr, "OU15"=ssOU15$lnr, "OU20"=ssOU20$lnr, "QU25"=ssOU25$lnr)
mlnL