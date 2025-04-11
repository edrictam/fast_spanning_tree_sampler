source("helpers/helper.R")
source("helpers/helper_process_samples.R")

## *only used for processing gibbs sampler output* 
## given an undirected tree and the root, orient all edges to point away from root
## give the vertex *Name* of root as input
bfs_orient_graph = function(g, root = "1"){
  gg = make_empty_graph(n = length(V(g)), directed = TRUE)
  V(gg)$name <- 1:length(V(g))
  need_to_visit = c(root)
  visited = c()
  while (length(need_to_visit) > 0){
    parent = need_to_visit[1]
    need_to_visit = need_to_visit[-1]
    neighborhood = neighbors(g, as.character(parent), mode = "all")
    for (neighbor in neighborhood$name){
      if (!(neighbor %in% visited)){
        gg <- add_edges(gg, c(as.character(parent), as.character(neighbor)))
        need_to_visit = c(need_to_visit, as.character(neighbor))
      }
    }
    visited = c(visited, parent)
  }
  gg
}

## *only used for processing gibbs sampler output* 
## given an adjacency matrix, construct an undirected graph with appropriate names
reconstruct_graph_from_adjacency <- function(spanning){
  K_tilde = dim(spanning)[1]
  colnames(spanning) = 1:K_tilde
  g = graph_from_adjacency_matrix(spanning, mode = "undirected")
  bfs_orient_graph(g)
}

## ** not really used anymore, but saved here for reference** 
find_root = function(g){
  allNodes<- V(g)$name
  for (node in allNodes){
    if(length(neighbors(g,as.character(node), mode = "in"))==0){
      return(as.character(node))
    }
  }
  return("0")
}

## v is *name* of vertex
contract_empty_non_leaf_nodes <-function(g, v){
  parent = neighbors(g, as.character(v), mode = "in")
  children = neighbors(g, as.character(v), mode = "out")
  gg <- delete_edges(g, paste(parent$name, v, sep = "|"))
  for (child in children$name){
    gg <- delete_edges(gg, paste(v, child, sep = "|"))
  }
  for (child in children$name){
    gg <- add_edges(gg, c(as.character(parent$name), as.character(child)))
  }
  gg <- delete_vertices(gg, as.character(v))
  gg
} 

prune_tree <- function(g, z, K, root = "1"){
  K_tilde = K + 1
  empty_nodes = get_empty_nodes(z, K)
  leaves = get_leaves(g)
  nonleaves = setdiff(1:K_tilde, leaves)
  empty_leaf_nodes = intersect(leaves, empty_nodes)
  empty_non_leaf_nodes = intersect(empty_nodes, nonleaves)
  empty_non_leaf_non_root_nodes = setdiff(empty_non_leaf_nodes, as.integer(root))
  for (v in empty_non_leaf_non_root_nodes){
    g <- contract_empty_non_leaf_nodes(g, as.character(v))
  }
  for (v in empty_leaf_nodes){
    g <- delete_vertices(g, as.character(v))
  }
  g
}


evaluate_density = function(mu, Sigma, y, pi, z, lambda_sq, edge_list, K, nu, Sigma_0, a1, b1, conc){
  n = dim(y)[1]
  p = dim(y)[2]
  K_tilde = K + 1
  
  logL = 0
  eps = 1e-30
  for (i in 1:n){
    logL = logL + log(dmvnorm(y[i,], mu[z[i], ], Sigma))
    logL = logL + log(pi[z[i]] + eps) 
  }
  logL = logL + log(dinvwishart(Sigma, nu, Sigma_0))
  logL = logL + log(dinvgamma(lambda_sq, a1, b1))
  logL = logL +  ddirichlet(pi[2:K_tilde], rep(conc, K), log = TRUE)
  
  ## edgelist is E x 2 matrix, with left being source and right being sink node
  num_edge = dim(edge_list)[1]
  for (q in 1:num_edge){
    logL = logL + log(dmvnorm(mu[edge_list[q,2],], mu[edge_list[q,1],], lambda_sq * diag(p)))
  }
  
  logL
}


## given a graph g, plot it with tree layout
plot_as_tree <- function(g, root = "1"){
  plot(as.undirected(g), layout = layout_as_tree(g, root = root), vertex.label=NA)
}


