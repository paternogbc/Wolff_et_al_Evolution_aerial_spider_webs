### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Phylogenetic Generalized Linear Squares analysis 
### Wolff et al.: Physical optimum in anchor points as a global driver of spider web evolution 
### Electronic Supplemental Material ESM.6
### Last update: 29/09/2018
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

### START--------------------------------------------------------------------------------------
### Packages-----------------------------------------------------------------------------------
library(phylolm)
library(sensiPhy)

### CLeaning environment and loading additional functions
rm(list = ls())
source("R/utils.R")


### Load data----------------------------------------------------------------------------------
d <- read.csv("data/data_spider_web.csv")
p <- read.tree("data/all.trees.tre")

### Enforce factors to categorical variables:
d$cribellum  <- as.factor(d$cribellum)
d$back.forth <- as.factor(d$back.forth)
d$web_type   <- as.factor(d$web_type)

### Add row names matching phylogeny:
rownames(d) <- d$species

### PGLS regressions---------------------------------------------------------------------------
### 1 centrality ~ cribellum----------------------------------
### Match data and tree
dc <- match_dataphy(centrality ~ cribellum, data = d, phy = p)

### PGLS regression accounting for phylogenetic uncertainty:
sensi.mod1 <- tree_phylm.m(formula = centrality ~ cribellum, data = d, phy = p,
                           n.tree = 100, model = "lambda")  

### Summary with average model estimates across all trees:
summary.complex(sensi.mod1)

### Diagnostic graphs for phylogentic uncertainty
### A: Estimates and B: P value for Cribellum across all trees 
sensi_plot_M(sensi.mod1, param = 2)
#-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --#

### 2 track_proportions ~ cribellum ----------------------------
### Match data and tree
dc <- match_dataphy(track_proportions ~ cribellum, data = d, phy = p)

### PGLS regression accounting for phylogenetic uncertainty:
sensi.mod2 <- tree_phylm.m(formula = track_proportions ~ cribellum, data = d,
                           phy = p, n.tree = 100, model = "lambda")  

### Summary with average model estimates across all trees:
summary.complex(sensi.mod2)

### Diagnostic graphs for phylogentic uncertainty
### A: Estimates and B: P values for cribellum across all trees 
sensi_plot_M(sensi.mod2, param = 2)
#-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --#

### 3 centrality ~ web type ------------------------------------
### Match data and tree
dc <- match_dataphy(centrality ~ web_type, data = d, phy = p)

### PGLS regression accounting for phylogenetic uncertainty:
sensi.mod3 <- tree_phylm.m(formula = centrality ~ web_type, data = d, phy = p, 
                           n.tree = 100, model = "lambda")  

### Summary with average model estimates across all trees:
summary.complex(sensi.mod3)

### Diagnostic graphs for phylogentic uncertainty
### A: Estimates and B: P values for web type1 across all trees 
sensi_plot_M(sensi.mod3, param = 2)
### A: Estimates and B: P values for web type2 across all trees 
sensi_plot_M(sensi.mod3, param = 3)
#-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --#

### 4 track_proportions ~ web_type----------------------------
### Match data and tree
dc <- match_dataphy(track_proportions ~ web_type, data = d, phy = p)

### PGLS regression accounting for phylogenetic uncertainty:
sensi.mod4 <- tree_phylm.m(formula = track_proportions ~ web_type, data = d,
                           phy = p, n.tree = 100, model = "lambda")  

### Summary with average model estimates across all trees:
summary.complex(sensi.mod4)

### Diagnostic graphs for phylogentic uncertainty
### A: Estimates and B: P values for web_type1 across all trees 
sensi_plot_M(sensi.mod4, param = 2)
### A: Estimates and B: P values for web_type2 across all trees 
sensi_plot_M(sensi.mod4, param = 3)
#-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --#

### 5 centrality ~ track_proportions------------------------------------
### Match data and tree
dc <- match_dataphy(centrality ~ track_proportions, data = d, phy = p)

### PGLS regression accounting for phylogenetic uncertainty:
sensi.mod5 <- tree_phylm(formula = centrality ~ track_proportions, data = d, phy = p, 
                           n.tree = 100, model = "lambda")  

### Summary with average model estimates across all trees:
summary(sensi.mod5)

### Diagnostic graphs for phylogentic uncertainty
### A: Estimates and B: P values for back.forth across all trees 
sensi_plot(sensi.mod5)
#-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --#

