### Script to analyse spider web evolution
### Jonas Wolff & Gustavo Paterno
### Last update: 10/09/2018
### Script: OUwie - Comparing evolution between ecological regimes

### START--------
### Packages-----
library(ape)
library(OUwie)
library(qpcR)
library(phytools)
library(dplyr)
library(ggplot2)
rm(list = ls())

### Load data----------------------------------------------------------------------------------
d     <- read.csv("data/data_spider_web.csv")
trees <- read.tree("data/all.trees.tre")
centrouwie1 <- read.csv("data/centrouwie1.csv")
centrouwie5 <- read.csv("data/centrouwie5.csv")
c <- setNames(centrouwie1[, "reg"], centrouwie1[, "species"])
w <- setNames(centrouwie5[, "reg"], centrouwie5[, "species"])

### Add row names matching phylogeny:
rownames(d) <- d$species

#### Drop outgroup species and match data and phy
trees <- match_dataphy(species ~ 1, data = d, phy = trees)[[2]]
name.check(trees[[1]], d)

# This procedure fits single and multi-regime BM and OU models to find out, if cribellum and web type had an 
# effect on the evolution of variables. 
### Create storage data.frames
aic.tab     <- data.frame("model" = character(), "deltaAIC" = numeric(), "aiccw" = numeric(), tree = numeric())
estimates   <- data.frame("model" = character(), "parameters" = numeric(), "tree" = numeric())

  
### START LOOP:
for (j in 1:100) {
paste("Tree = ", j)
### Storage values:
### Select tree  
tree <- trees[[j]]


  mtree_c <- make.simmap(tree, c, model = matrix(c(0, 0, 1, 0), nrow = 2, ncol = 2 , byrow = TRUE), nsim = 1) 
  mtree_w <- make.simmap(tree, w, model = "ARD", nsim = 1)
  BM1 <- OUwie(mtree_c, centrouwie1, model=c("BM1"), simmap.tree=T, root.station=TRUE, mserr="known")
  BMSc <- OUwie(mtree_c, centrouwie1, model=c("BMS"), simmap.tree=T, root.station=TRUE, mserr="known")
  BMSw <- OUwie(mtree_w, centrouwie5, model=c("BMS"), simmap.tree=T, root.station=TRUE, mserr="known")
  OU1 <- OUwie(mtree_c, centrouwie1, model=c("OU1"), simmap.tree=T, root.station=TRUE, mserr="known")
  OUMc <- OUwie(mtree_c, centrouwie1, model=c("OUM"), simmap.tree=T, root.station=TRUE, mserr="known")
  OUMw <- OUwie(mtree_w, centrouwie5, model=c("OUM"), simmap.tree=T, root.station=TRUE, mserr="known")
  OUMVc <- OUwie(mtree_c, centrouwie1, model=c("OUMV"), simmap.tree=T, root.station=TRUE, mserr="known")
  OUMVw <- OUwie(mtree_w, centrouwie5, model=c("OUMV"), simmap.tree=T, root.station=TRUE, mserr="known")
  OUMAc <- OUwie(mtree_c, centrouwie1, model=c("OUMA"), simmap.tree=T, root.station=TRUE, mserr="known")
  OUMAw <- OUwie(mtree_w, centrouwie5, model=c("OUMA"), simmap.tree=T, root.station=TRUE, mserr="known")
  OUMVAc <- OUwie(mtree_c, centrouwie1, model=c("OUMVA"), simmap.tree=T, root.station=TRUE, mserr="known")
  OUMVAw <- OUwie(mtree_w, centrouwie5, model=c("OUMVA"), simmap.tree=T, root.station=TRUE, mserr="known")

  # get AICc
  BM1_AICc <- BM1$AICc
  BMSc_AICc <- BMSc$AICc
  BMSw_AICc <- BMSw$AICc
  OU1_AICc <- OU1$AICc
  OUMc_AICc <- OUMc$AICc
  OUMw_AICc <- OUMw$AICc
  OUMVc_AICc <- OUMVc$AICc
  OUMVw_AICc <- OUMVw$AICc
  OUMAc_AICc <- OUMAc$AICc
  OUMAw_AICc <- OUMAw$AICc
  OUMVAc_AICc <- OUMVAc$AICc
  OUMVAw_AICc <- OUMVAw$AICc

  # calculate dAICc and AICc-weights
  aicc <- c(BM1_AICc, BMSc_AICc, BMSw_AICc, OU1_AICc, OUMc_AICc, OUMw_AICc, OUMVc_AICc, OUMVw_AICc,
            OUMAc_AICc, OUMAw_AICc, OUMVAc_AICc, OUMVAw_AICc)
  aics <- akaike.weights(c(BM1_AICc, BMSc_AICc, BMSw_AICc, OU1_AICc, OUMc_AICc, OUMw_AICc, OUMVc_AICc, OUMVw_AICc,
                         OUMAc_AICc, OUMAw_AICc, OUMVAc_AICc, OUMVAw_AICc))
  daicc <-  aics$deltaAIC
  aiccw <- aics$weights

  model <- c("BM1", "BMSc" , "BMSw", "OU1", "OUMc", "OUMw", "OUMVc", "OUMVw",
            "OUMAc", "OUMAw" , "OUMVAc" , "OUMVAw") 

aic.val <- data.frame(model, AICc = aicc, deltaAICc = daicc, aiccw = aiccw, tree = j)

  # summarize parameter estimates in table 
param.est <- rbind(
  BM1_sig2 = BM1$solution[2,1],
  BMSc_sig2_0 = BMSc$solution[2,1],
  BMSc_sig2_1 =BMSc$solution[2,2],
  BMSw_sig2_1 =BMSw$solution[2,1],
  BMSw_sig2_2 =BMSw$solution[2,2],
  OU1_sig2 =OU1$solution[2,1],
  OU1_alpha =OU1$solution[1,1],
  OU1_theta =OU1$theta[1,1],
  OUMc_sig2 =OUMc$solution[2,1],
  OUMc_alpha =OUMc$solution[1,1],
  OUMc_theta_0 =OUMc$theta[1,1],
  OUMc_theta_1 =OUMc$theta[2,1],
  OUMw_sig2 =OUMw$solution[2,1],
  OUMw_alpha =OUMw$solution[1,1],
  OUMw_theta_1 =OUMw$theta[1,1],
  OUMw_theta_2 =OUMw$theta[2,1],
  OUMVc_sig2_0 =OUMVc$solution[2,1],
  OUMVc_sig2_1 =OUMVc$solution[2,2],
  OUMVc_alpha =OUMVc$solution[1,1],
  OUMVc_theta_0 =OUMVc$theta[1,1],
  OUMVc_theta_1 =OUMVc$theta[2,1],
  OUMVw_sig2_1 =OUMVw$solution[2,1],
  OUMVw_sig2_2 =OUMVw$solution[2,2],
  OUMVw_alpha =OUMVw$solution[1,1],
  OUMVw_theta_1 =OUMVw$theta[1,1],
  OUMVw_theta_2 =OUMVw$theta[2,1],
  OUMAc_sig2 =OUMAc$solution[2,1],
  OUMAc_alpha_0 =OUMAc$solution[1,1],
  OUMAc_alpha_1 =OUMAc$solution[1,2],
  OUMAc_theta_0 =OUMAc$theta[1,1],
  OUMAc_theta_1 =OUMAc$theta[2,1],
  OUMAw_sig2 =OUMAw$solution[2,1],
  OUMAw_alpha_1 =OUMAw$solution[1,1],
  OUMAw_alpha_2 =OUMAw$solution[1,2],
  OUMAw_theta_1 =OUMAw$theta[1,1],
  OUMAw_theta_2 =OUMAw$theta[2,1],
  OUMVAc_sig2_0 =OUMVAc$solution[2,1],
  OUMVAc_sig2_1 =OUMVAc$solution[2,2],
  OUMVAc_alpha_0 =OUMVAc$solution[1,1],
  OUMVAc_alpha_1 =OUMVAc$solution[1,2],
  OUMVAc_theta_0 =OUMVAc$theta[1,1],
  OUMVAc_theta_1 =OUMVAc$theta[2,1],
  OUMVAw_sig2_1 =OUMVAw$solution[2,1],
  OUMVAw_sig2_2 =OUMVAw$solution[2,2],
  OUMVAw_alpha_1 =OUMVAw$solution[1,1],
  OUMVAw_alpha_2 =OUMVAw$solution[1,2],
  OUMVAw_theta_1 =OUMVAw$theta[1,1],
  OUMVAw_theta_2 =OUMVAw$theta[2,1])

  param.est <- data.frame(model = rownames(param.est), parameters = as.numeric(param.est), tree = j)
  rownames(param.est) <- NULL

  ### Store Estimates
  aic.tab   <- rbind(aic.tab, aic.val)
  estimates <- rbind(estimates, param.est)
  
}

### Summary STATS:----------
### Summarise values across all trees:
aic.avg <- 
  summarise(group_by(aic.tab, model),
            AICc_median         = median(AICc, na.rm = T),
            AICc_sd             = sd( AICc, na.rm = T),
            deltaAICc_median    = median(deltaAICc, na.rm = T),
            deltaAICc_sd        = sd(deltaAICc, na.rm = T),
            aiccw_median        = median(aiccw, na.rm = T),
            aiccw_sd            = sd(aiccw, na.rm = T)) %>% arrange(-aiccw_median)

est.avg <- 
  summarise(group_by(estimates, model),
            estimate_median = median(parameters, na.rm = T),
            estimate_sd   = sd(parameters, na.rm = T)) %>% arrange(estimate_median)
  
## PLOTS-----
### Summary plot for AIC:
g1 <- ggplot(aic.tab, aes(y = AICc, x = reorder(model, AICc), color = model, fill = model)) +
  geom_boxplot(show.legend = F) +
  geom_jitter(width = .02, show.legend = F) +
  theme_bw(base_size = 18) 
g2 <- ggplot(aic.tab, aes(y = deltaAICc, x = reorder(model, deltaAICc), color = model, fill = model)) +
  geom_boxplot(show.legend = F) +
  geom_jitter(width = .02, show.legend = F) +
  theme_bw(base_size = 18)
g3 <- ggplot(aic.tab, aes(y = aiccw, x = reorder(model, -aiccw), color = model, fill = model)) +
  geom_boxplot(show.legend = F) +
  geom_jitter(width = .02, show.legend = F) +
  theme_bw(base_size = 18)
g3

library(ggpubr)
ggarrange(g1, g2, g3, nrow = 3)


### Summary plot for parameters estimates:
### Choose which models to plot:
m.choosed <- c("BMSc_sig2_0", "BMSc_sig2_1", "BMSw_sig2_1", "BMSw_sig2_2")
estimates.crop <- filter(estimates, model %in% m.choosed) %>% droplevels()
g4 <- ggplot(estimates.crop, aes(y = parameters, x = model, color = model, fill = model)) +
  geom_boxplot(show.legend = F) +
  geom_jitter(width = .02, show.legend = F) +
  theme_bw(base_size = 18) 
g4

m.choosed <- c("OU1_alpha","OUMc_alpha","OUMVAw_alpha_1", "OUMVAw_alpha_2", "OUMAw_alpha_1", "OUMAw_alpha_2")
estimates.crop <- filter(estimates, model %in% m.choosed) %>% droplevels()
g5 <- ggplot(estimates.crop, aes(y = parameters, x = model, color = model, fill = model)) +
  geom_boxplot(show.legend = F) + coord_cartesian(ylim = c(0, 0.015))  +
  geom_jitter(width = .02, show.legend = F) +
  theme_bw(base_size = 18) 
g5

m.choosed <- c("OU1_theta","OUMc_theta_0","OUMc_theta_1","OUMVAw_theta_1", "OUMVAw_theta_2", "OUMAw_theta_1", "OUMAw_theta_2")
estimates.crop <- filter(estimates, model %in% m.choosed) %>% droplevels()
g6 <- ggplot(estimates.crop, aes(y = parameters, x = model, color = model, fill = model)) +
  geom_boxplot(show.legend = F) + coord_cartesian(ylim = c(0, 0.4))  +
  geom_jitter(width = .02, show.legend = F) +
  theme_bw(base_size = 18) 
g6

m.choosed <- c("BM1_sig2","OU1_sig2","OUMc_sig2","OUMVAw_sig2_1", "OUMVAw_sig2_2", "OUMAw_sig2")
estimates.crop <- filter(estimates, model %in% m.choosed) %>% droplevels()
g7 <- ggplot(estimates.crop, aes(y = parameters, x = model, color = model, fill = model)) +
  geom_boxplot(show.legend = F)   +
  geom_jitter(width = .02, show.legend = F) +
  theme_bw(base_size = 18) 
g7

### Export Results:
#write.csv(x = est.avg, "output/estimates_avg.csv")
#write.csv(x = aic.avg, "output/aic_avg.csv")
#write.csv(x = aic.tab, "output/aic_tab.csv")
