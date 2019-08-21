### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Additional functions to perform sensitivity analysis 
### Wolff et al.: Physical optimum in anchor points as a global driver of spider web evolution 
### Electronic Supplemental Material ESM.6
### Last update: 26/09/2018
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

### START------------------------------------------------------------------------------------------
#' Performs complex PGLS models evaluating uncertainty across phylogenetic trees.
#' Author: Caterina Penone & Gustavo Paterno / 2018
#' tree_phylm.m-----------
tree_phylm.m <- function(formula,data,phy,
                         n.tree=2,model="lambda",track=TRUE,...){
  #Error check
  if(class(formula)!="formula") stop("formula must be class 'formula'")
  if(class(data)!="data.frame") stop("data must be class 'data.frame'")
  if(class(phy)!="multiPhylo") stop("phy must be class 'multiPhylo'")
  if(length(phy)<n.tree) stop("'n.tree' must be smaller (or equal) than the number of trees in the 'multiPhylo' object")
  if ( (model == "trend") & (ape::is.ultrametric(phy[[1]])))
    stop("Trend is unidentifiable for ultrametric trees., see ?phylolm for details")
  else
    
    
    #Matching tree and phylogeny using utils.R
    datphy    <- sensiPhy::match_dataphy(formula,data,phy)
  full.data <-datphy[[1]]
  phy       <-datphy[[2]]
  
  # If the class of tree is multiphylo pick n=n.tree random trees
  trees <- sample(length(phy), n.tree, replace = F)
  
  #Create the results list
  estimates <- list()
  
  #Model calculation
  counter=1
  errors <- NULL
  if(track==TRUE) pb <- utils::txtProgressBar(min = 0, max = n.tree, style = 3)
  for (j in trees){
    #Match data order to tip order
    full.data <- full.data[phy[[j]]$tip.label,]
    
    #phylolm model
    mod = try(phylolm::phylolm(formula, data = full.data, model = model, phy = phy[[j]]), FALSE)
    
    
    if(isTRUE(class(mod)=="try-error")) {
      error <- j
      names(error) <- rownames(c.data$full.data)[j]
      errors <- c(errors,error)
      next }
    
    if(track==TRUE) utils::setTxtProgressBar(pb, counter)
    
    #write in a table
    estimates[[counter]]  <- phylolm::summary.phylolm(mod)[[2]] 
    counter=counter+1
    
  }
  
  if(track==TRUE) on.exit(close(pb))
  
  ### Merge lists into a single data.frame
  sensi.estimates <- do.call("rbind", estimates)
  
  ### Re-organize global estimates:
  params  <- rownames(sensi.estimates)
  sensi.estimates <- data.frame(sensi.estimates, row.names = NULL)
  n.p <- length(names(mod$coefficients))
  tree.n <- rep(as.character(trees), each = n.p)
  sensi.estimates$tree <- tree.n
  sensi.estimates$params <- params
  sensi.estimates <- sensi.estimates[c("tree", "params", "Estimate", "StdErr", "t.value", "p.value")]
  
  ### output
  res <- list(   call = match.call(),
                 formula =formula,
                 data =full.data,
                 sensi.estimates = sensi.estimates,
                 N.obs = mod$n)
  return(res)
}

#' Summarise results from sensitivity analysis.
summary.complex <- function(x) { 
  
  ### Get data:  
  es <- x$sensi.estimates
  n.tree <- length(unique(es$tree))
  par <- unique(es$params)
  out <- list()
  counter <- 1
  ### summary output:
  for (p in par) {
    apply(es[es$params == p, 3:6], 2, mean)
    statresults <- data.frame(min=apply(es[es$params == p, 3:6],2,min),
                              max=apply(es[es$params == p, 3:6],2,max),
                              mean=apply(es[es$params == p, 3:6],2,mean),
                              sd_tree=apply(es[es$params == p, 3:6],2,stats::sd))
    
    statresults$CI_low  <- statresults$mean - qt(0.975, df = n.tree-1) * statresults$sd_tree / sqrt(n.tree)
    statresults$CI_high  <- statresults$mean + qt(0.975, df = n.tree-1) * statresults$sd_tree / sqrt(n.tree)
    
    out[[counter]]  <- statresults
    counter <- counter + 1
  }
  names(out) <- par
  return(out)
}

#' Plot results from sensitivity analysis.
sensi_plot_M <- function(x, param = 1, ...){
  
  estim <- x$sensi.estimates
  pa    <- unique(estim$params)[param]
  
  ### Nulling variables
  pval <- DIFestimate <- Estimate <- Significant <- Species.removed  <- NULL
  element_line <- estimate <- geom_jitter <- NULL
  
  ### Distribution of Estimated values estimated
  dd <- subset(estim, estim$params == pa)
  e1 <- ggplot2::ggplot(dd, aes(x = Estimate))+
    geom_histogram(fill="yellow",colour="black", size=.2,
                   alpha = .3) +
    geom_vline(xintercept = mean(dd$Estimate), color="red",linetype=2,size=.7)+
    xlab(paste("Estimated", pa, "values", sep = " "))+
    ylab("Frequency")+
    theme(axis.title=element_text(size=12),
          axis.text = element_text(size=12),
          panel.background = element_rect(fill="white",
                                          colour="black"))
  
  ### Distribution of Estimated P values estimated
  e2 <- ggplot2::ggplot(dd, aes(x =  p.value))+
    geom_histogram(fill="yellow",colour="black", size=.2,
                   alpha = .3) +
    geom_vline(xintercept = mean(dd$p.value), color="red",linetype=2,size=.7)+
    geom_vline(xintercept = 0.05, color="black",linetype=1,size=.5)+
    
    xlab(paste("Estimated","P values for ", pa,   sep = " "))+
    ylab("Frequency")+
    theme(axis.title=element_text(size=12),
          axis.text = element_text(size=12),
          panel.background = element_rect(fill="white",
                                          colour="black"))
  
  suppressMessages(return(multiplot(e1,e2, cols=2)))
} 

### Function to plot multiple ggplo2 graphs together:
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
    }
  }
}


