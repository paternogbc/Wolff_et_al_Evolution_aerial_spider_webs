### Script to analyze spider web evolution
### Jonas Wolff, Arie van der Meijden & Gustavo Paterno
### Last update: 10/09/2018
### Reference for analyses: 
### Script: Geometric morphometric PCMs on spinning tracks

### START--------
### Packages-----
library(geomorph)
library(phytools)
rm(list = ls())

### Load data-----
d     <- read.csv("data/data_spider_web.csv")
smus  <- read.csv(file = "data/smu-landmarks.csv", row.names = 1) # smu-landmarks.csv
trees <- read.tree("data/all.trees.tre")
tree  <- trees[[17]]

### Add row names matching phylogeny:
rownames(d) <- d$species

### Crop phylogeny to match with missing values (track proportion & landmarks)
treex   <- match_dataphy(track_proportions ~ 1, data = d, phy = tree)[[2]]
dx   <- match_dataphy(track_proportions ~ 1, data = d, phy = tree)[[1]]

A <- arrayspecs(smus, 50, 2)
classifier <- read.csv(file = "data/smu_classifier.csv", header=T, row.names = 1) # smu-classifier.csv
classifier <- classifier[match(treex$tip.label, rownames(classifier)),]

### Procrustes Alignement ###
Y.gpa<-gpagen(A) 

### Phylogenetic signal ###
physignal(Y.gpa$coords, treex, iter=999, print.progress = TRUE) # GPA-aligned
physignal(A, treex, iter=999, print.progress = TRUE) # Spinneret-aligned

### Plot Phylo-Morphospace ###
plotGMPhyloMorphoSpace(treex, A, tip.labels = F)

### Comparing track shapes between groups ###
## Cribellar vs. ecribellar
cribx <- classifier$cribellum
names(cribx) <- row.names(classifier)
procD.pgls(A~cribx, phy=treex, iter = 999)

## Web types
webx <- classifier$web_type
names(webx) <- row.names(classifier)
procD.pgls(A~webx, phy=treex, iter = 999)

## centrality
centrx <- dx$centrality
names(centrx) <- row.names(dx)
procD.pgls(A~centrx, phy=treex, iter = 999)
plot(procD.pgls(A~centrx, phy=treex, iter = 999))

### Compare evolutionary rates ###
## Cribellar vs. ecribellar
compare.evol.rates(A, treex, cribx, iter = 999)

## Web types
compare.evol.rates(A, treex, webx, iter = 999)