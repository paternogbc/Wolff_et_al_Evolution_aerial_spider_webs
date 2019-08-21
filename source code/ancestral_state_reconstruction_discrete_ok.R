### Script to analyze spider web evolution
### Gustavo Paterno & Jonas Wolff
### Last update: 25/10/2018
### Script: Ancestral state reconstruction (discrete)

### START--------
### Packages-----
library(phytools)
library(geiger)
library(sensiPhy)
rm(list = ls())

### Load data-----
d     <- read.csv(file = "data/data_spider_web.csv")
trees <- read.tree(file = "data/all.trees.tre") 
tree  <- read.nexus(file = "data/infile.tre")

### Add row names matching phylogeny:
rownames(d) <- d$species

#### Drop outgroup species and match data and phy
trees <- match_dataphy(species ~ 1, data = d, phy = trees)[[2]]
tree  <- match_dataphy(species ~ 1, data = d, phy = tree)[[2]]

### Ancestral state reconstruction ---- 
### 1. Cribellum ----
x <- setNames(d[, "cribellum"], rownames(d))

## Compare models
fit_ER <- fitDiscrete(tree, x, model = "ER") 
fit_ARD <- fitDiscrete(tree, x, model = "ARD")  
Dollo <- matrix(c(0, 0, 1, 0), nrow = 2, ncol = 2 , byrow = TRUE)
fit_Dollo <- fitDiscrete(tree, x, model = Dollo)  
aicc<-setNames(c(fit_ER$opt$aicc, fit_ARD$opt$aicc, fit_Dollo$opt$aicc), c("ER","ARD","Dollo"))
aic.w(aicc)

## Stochastic Character Mapping
mtree <- make.simmap(tree, x, model = Dollo, nsim = 100)
describe.simmap(mtree)

## Density Plot:
cols<-setNames(palette()[1:length(unique(x))],sort(unique(x)))
obj   <- densityMap(mtree, lwd = 3, outline = TRUE, fsize = .4, ftype = "i")
## Pie Plot:
pd <- summary(mtree, plot = F)
plot(pd, cex = .3, fsize = .3)
add.simmap.legend(colors = cols, prompt = FALSE, y = 5 ,fsize=0.8)


### 2. Web Type ----
x   <- setNames(d[, "web_type"], rownames(d))
cols<-setNames(palette()[1:length(unique(x))],sort(unique(x)))

## Compare models
fit_ER <- fitDiscrete(tree, x, model = "ER")   
fit_SYM <- fitDiscrete(tree, x, model = "SYM")  
fit_ARD <- fitDiscrete(tree, x, model = "ARD")  
aicc<-setNames(c(fit_ER$opt$aicc, fit_SYM$opt$aicc, fit_ARD$opt$aicc), c("ER","SYM","ARD"))
aic.w(aicc)

## Stochastic Character Mapping
mtree <- make.simmap(tree, x, model="ER", nsim = 100)
describe.simmap(mtree)

## Export Graph:
pd <- summary(mtree, plot=FALSE, fsize = .25)
plot(pd, cex = .4, fsize = .3, led = .7)
add.simmap.legend(colors = cols, prompt = FALSE, y = 5 ,fsize=0.3)

### 3. Considering Phylogenetic Uncertainty-----
### 3.1 Cribellum Stochastic Mapping-----
x   <- setNames(d[, "cribellum"], rownames(d))
cols<-setNames(palette()[1:length(unique(x))],sort(unique(x)))

## Density Plot:
mtree <- make.simmap(trees, x, model = Dollo, nsim = 1)
plot(mtree)

pd <- summary(mtree, plot = F)
plot(pd, cex = .3, fsize = .3)
add.simmap.legend(colors = cols, prompt = FALSE, y = 5 ,fsize=0.8)
###--------