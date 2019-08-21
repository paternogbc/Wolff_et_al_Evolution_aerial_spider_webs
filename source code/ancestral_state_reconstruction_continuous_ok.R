### Script to analyze spider web evoltuion
### Gustavo Paterno & Jonas Wolff
### Last update: 28/09/2018
### Script: Ancestral state reconstruction (continuous)

### START--------
### Packages-----
library(phytools)
library(geiger)
rm(list = ls())

### Load data-----
d <- read.csv(file = "data/data_spider_web.csv")
tree <- read.nexus(file = "data/infile.tre") 

### Add row names matching phylogeny:
rownames(d) <- d$species

#### Drop outgroup species and match data and phy
tree <- match_dataphy(species ~ 1, data = d, phy = tree)[[2]]
name.check(tree, d)

### 1. Ancestral state reconstruction ---- 
### continuous traits (by Maximum Likelihood)
### 1.1 Centrality -----
# make a contMap
centrality <- setNames(d[, "centrality"], rownames(d))
fit.centra <-fastAnc(tree, centrality, vars=TRUE,CI=TRUE)
fit.centra$CI95
obj <- contMap(tree, centrality, plot=F)
obj <- setMap(obj,colors=c("blue","yellow","red"))
plot(obj, lwd = 3, fsize = 0.5) 

# make a traitgram with overlayed density web and tip labels
x <- setNames(d[, "cribellum"], rownames(d))
mtree <- make.simmap(tree, x, model = matrix(c(0, 0, 1, 0), nrow = 2, ncol = 2 , byrow = TRUE), nsim = 100)
obj   <- densityMap(mtree, lwd = 3, 
                    outline = TRUE, fsize = .4, ftype = "i")
# Phenogram (Centrality, Web_Type)
dev.off()
phenogram(obj$tree, centrality, colors = obj$cols, lwd=2, spread.labels = F, add = F, fsize = 0, 
          xlab = "time", ylab = "centrality")
y    <- setNames(d[, "web_type"], rownames(d))
cols <-setNames(palette()[1:length(unique(y))],sort(unique(y)))
tiplabels(pie=to.matrix(y,sort(unique(y))),piecol=cols,cex=0.7)
add.simmap.legend(colors=cols,prompt=T, y = 5,fsize = 1)
