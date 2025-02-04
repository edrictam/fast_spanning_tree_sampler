library(coda)

## return the effective sample size per iteration 
## (chain that is passed in should already have burnin removed )
ESSPerIter <- function(chain, niter, burnin){
  effectiveSize(chain)/(niter-burnin)
}

## return the mean and sd of effective sample size per iteration 
## for vector parameters
## return mean and sd of ESS per iterations, averaged across the coordinates provided
## (chain that is passed in should already have burnin removed )
ESSPerIterVector <- function(chain, niter, burnin){
  list(avg = mean(effectiveSize(chain)/(niter-burnin)), 
       stdev = sd(effectiveSize(chain)/(niter-burnin)))
}


## return number of nodes in a graph
num_nodes = function(g, upper_bound){
  length(getTreeNodes(g, upper_bound))
}


## get all the leaves in the tree
# return nodes' names...need to turn back to character before use
get_leaves <- function(g){
  results = c()
  allNodes<- V(g)
  for (node in allNodes$name){
    if(length(neighbors(g,as.character(node), mode = "out"))==0){
      results = c(results, node)
    }
  }
  results
}

## get number of leaves in a tree graph
num_leaves = function(g){
  length(get_leaves(g))
}

## obtain max degree
max_degree = function(g){
  length(degree_distribution(g, mode = "all")) - 1
}

## find the root of a tree graph (name in character returned)
find_root = function(g){
  allNodes<- V(g)
  for (node in allNodes){
    if(length(neighbors(g,as.character(node), mode = "in"))==0){
      return(as.character(node))
    }
  }
  return("0")
}

## return the maximum depth of a rooted tree 
## (length of longest path from root to any vertex) 
max_depth = function(g, root = "1"){
  dists = shortest_paths(g, root)
  depth = 0
  for (i in 1:length(V(g))){
    len = length(dists[[1]][[i]])
    if (len  > depth){
      depth = len
    }
  }
  len
}

make.trace <- function(chain, name, x_off = TRUE, title_scale = 2, axis_scale = 1.5, title_off = TRUE){
  xlab_string = ifelse(x_off, "", "Iterations")
  title_string = ifelse(title_off, "", paste("Trace plot of", name)) 
  plot(1:length(chain), chain, pch = NA,  
       xlab = xlab_string, ylab = "Counts", main = title_string, 
       cex.main = title_scale, cex.lab = axis_scale, cex.axis = axis_scale)
  lines(1:length(chain), chain)
}


make.bars <- function(chain, name, x_off = TRUE, title_scale = 2, axis_scale = 1.5,  title_off = TRUE){
  xlab_string = ifelse(x_off, "", name)
  title_string = ifelse(title_off, "", paste("Bar plot of", name)) 
  barplot(table(chain), xlab = xlab_string, ylab = "Counts", 
          main = title_string, 
          cex.main = title_scale, cex.lab = axis_scale, cex.axis = axis_scale,
          cex.names= axis_scale)
}

make.acf <- function(chain, name, x_off = TRUE, title_scale = 2, axis_scale = 1.5, title_off = TRUE){
  xlab_string = ifelse(x_off, "", name)
  title_string = ifelse(title_off, "", paste("Autocorrelation plot of", name))
  plot(acf(chain,plot=F)$acf, type = "h", xlab = xlab_string, ylab = "Autocorrelation", ylim = c(0, 1),
       main = title_string, 
       cex.main = title_scale, cex.lab = axis_scale, cex.axis = axis_scale)
  abline(h = 0)
}