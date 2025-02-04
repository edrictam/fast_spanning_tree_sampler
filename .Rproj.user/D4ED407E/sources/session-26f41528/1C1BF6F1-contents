## Gibbs Sampler Full Conditionals
library(gtools)
library(mvtnorm)
library(invgamma)
library(LaplacesDemon)
library(igraph)

source("helpers/helper.R")


## Birth probability calculator 
p_birth <- function(g, p0,K){
  K_tilde = K + 1
  if (length(getTreeNodes(g, K_tilde)) < K_tilde){
    return(p0)
  } else{
    return(0)
  }
}

## Death probability calculator 
p_death <- function(g, p0,K){
  K_tilde = K + 1
  nodes = getTreeNodes(g, K_tilde)
  if (length(nodes) > 1){
    return(1 - p0)
  } else{
    return(0)
  }
}

## T sampler (rjmcmc version)
sample_T_rjmcmc = function(g, p0, K, z, Sigma, mu, verbose = FALSE, delta, tuning){
  returnlist = list(g = g, mu = mu)
  K_tilde = K + 1
  u = rbern(1, p0)
  
  V_T = getTreeNodes(g, K_tilde)$name
  V_T_size = length(V_T)
  V_T_C = setdiff(1:K_tilde, V_T)
  V_T_C_size = length(V_T_C)
  
  if ((u == 1) & (V_T_C_size > 0)){
    if (V_T_size == 1){
      j = V_T[1]
    } else {
      j = as.vector(sample(V_T, 1))
    }
    l = min(V_T_C)
    g_prime = add_vertices(g,1, attr = list(name = l))
    g_prime <- add_edges(g_prime, c(as.character(j), as.character(l)))
    accept_prob = min(1, tuning*p_death(g_prime, p0, K) * 1/length(emptyleaves(g_prime, z, K)) / (p_birth(g, p0, K)/V_T_size))
    q = rbern(1, accept_prob)
    if (q){
      g = g_prime
      mu[l, ] <- rmvnorm(1, mu[j, ], (1/delta)*Sigma)
      if (verbose){
        print("birth")
      }
    } else {
      if (verbose){
        print("failed birth")
        print(accept_prob)
      }
    }
    returnlist = list(g = g, mu = mu)
    return(returnlist)
  } 
  
  V_T_empty_leaves = emptyleaves(g, z, K)
  V_T_empty_leaves_size = length(V_T_empty_leaves)
  
  if ((u == 0) & (V_T_size > 1)){
    ## death
    if (V_T_empty_leaves_size == 0){
      return(returnlist)
    }
    if (V_T_empty_leaves_size == 1){
      j = V_T_empty_leaves[1]
    } else{
      j = as.vector(sample(V_T_empty_leaves, 1))
    }
    g_prime <- delete_vertices(g, as.character(j))
    accept_prob = min(1, (p_birth(g_prime, p0, K)/(V_T_size - 1))/(tuning*p_death(g, p0, K)/V_T_empty_leaves_size))
    q = rbern(1, accept_prob)
    if (q) {
      if (verbose){
        print("death")
      }
      g = g_prime
    } else {
      if (verbose){
        print("failed death")
        print(accept_prob)
      }
    }
    returnlist = list(g = g, mu = mu)
    return(returnlist)
  }
  return(returnlist)
}

## z sampler 
sample_z_rjmcmc = function(y,pi,mu, Sigma, V_T){
  n = dim(y)[1]
  p = dim(y)[2]
  
  V_T_size = length(V_T)
  p.z.given.y = matrix(rep(0, n*V_T_size), ncol = V_T_size)
  for (i in 1:n){
    for (k in 1:V_T_size){
      node = V_T[k]
      p.z.given.y[i, k] = pi[node]*dmvnorm(y[i,], mu[node,], Sigma)
    }
    denominator = sum(p.z.given.y[i,])
    p.z.given.y[i, ] = p.z.given.y[i, ]/denominator
  }
  z = rep(0, n)
  if (V_T_size == 1){
    z = rep(V_T[1], n)
    return(z)
  }
  for(i in 1:n){
    z[i] = sample(V_T, size=1,prob=p.z.given.y[i,],replace=TRUE)
  }
  return(z)
}

## pi sampler
sample_pi = function(z,K, conc = 0.1){
  K_tilde = K + 1
  counts = colSums(outer(z,2:K_tilde,FUN="=="))
  pi = rdirichlet(1,counts + conc)  
  return(c(0, pi))
}

# mu sampler
sample_mu_rjmcmc = function(y, z, K, p, curr_mu, Sigma, V_T, g, delta){
  mu = curr_mu
  
  for(k in V_T){
    sample.size = sum(z==k)
    delta_k = delta + sample.size
    sample.sum = rep(0, p)
    
    if (sample.size > 1) {
      sample.sum = as.matrix(colSums(y[z==k,], ))
    }
    
    if (sample.size == 1) {
      sample.sum =  as.matrix(y[z==k,])
    }
    neighboring = neighbors(g, as.character(k), mode = "all")$name
    neighbor_values = mu[neighboring, ]
    if (is.vector(neighbor_values)){
      neighbors.sum = as.matrix(neighbor_values)
    }else{
      neighbors.sum = as.matrix(colSums(neighbor_values))
    }
    num_neighbors = length(neighboring)
    neighbors.avg = neighbors.sum/num_neighbors
    post.mean = (delta*neighbors.avg + sample.sum)/(delta_k)
    # post.mean = (delta*neighbors.sum + sample.sum)/(delta_k)
    mu[k,] = rmvnorm(1, post.mean, Sigma/delta_k)
  }
  return(mu)
}

# Sigma sampler
sample_Sigma_rjmcmc = function(y,K, nu, Sigma_0, delta, z, mu, V_T){
  n = dim(y)[1]
  p = dim(y)[2]
  Sigma_new = Sigma_0
  nu_new = nu
  for(k in V_T){
    sample.size = sum(z==k)
    nu_new = nu_new + sample.size
    
    if (sample.size > 1) {
      mu_hat = colMeans(y[z==k,])
      sample_diff = sweep(y[z==k,], 2, mu_hat)
      mu_diff = mu[k, ] - mu_hat
      S = t(sample_diff)%*%sample_diff
      var_term = mu_diff %*% t(mu_diff)
      Sigma_new = Sigma_new + S + var_term*(delta*sample.size/(sample.size + delta))
    }
    
    if (sample.size == 1) {
      X = y[z==k,] - mu[k, ]
      vec = as.vector(X)
      sample.var = vec %*% t(vec)
      Sigma_new = Sigma_new + sample.var*(delta*sample.size/(sample.size + delta))
    }
  }
  Sigma = rinvwishart(nu_new, Sigma_new)
  
  return(Sigma)
}


## Main rjmcmc function
rjmcmc = function(y,K,niter =1000, nu, Sigma_0, delta, p0, verbose = FALSE, tuning = 0.1){

  print("In RJMCMC function")
  ## initialize variables
  n = dim(y)[1]
  p = dim(y)[2]
  
  ## set up starting points
  pi = sample_pi(rep(0, n), K)
  mu = matrix(rep(0, K * p), K, p)  
  mu_root = rep(0, p)
  mu = rbind(mu_root, mu)
  Sigma = diag(p)
  
  K_tilde = K + 1

  g <- make_empty_graph(2, directed = TRUE)
  vertex_attr(g) <- list(name = 1:2)
  g <- add_edges(g, c(as.character(1), as.character(2)))

  V_T = getTreeNodes(g, K_tilde)$name
 
  z = sample_z_rjmcmc(y,pi,mu,Sigma, V_T)
  
  ## set up array to store samples
  res = list(mu=array(0, dim = c(niter, K + 1, p)), 
             pi = matrix(nrow=niter,ncol=K + 1), 
             z = matrix(nrow=niter, ncol=n), 
             Sigma = array(0, dim = c(niter, p, p)))
  write_graph(g, paste("saved_trajectories/rjmcmc_trees/tree", 1, ".graphml", sep = ""), format = "graphml")
  
  ## initialize chain
  res$mu[1,,]=mu
  res$pi[1,]=pi
  res$z[1,]=z 
  res$Sigma[1,,]=Sigma
  
  print("Initialization finished")
  ## loop through chain
  for(i in 2:niter){
    if (i %% 100 == 0){
      print("In RJMCMC loop iteration number: ")
      print(i)
    }
    ## sample from full conditionals
    pi = sample_pi(z,K)
    mu = sample_mu_rjmcmc(y, z, K, p, mu, Sigma, V_T, g, delta)
    z = sample_z_rjmcmc(y,pi,mu, Sigma, V_T)
    Sigma = sample_Sigma_rjmcmc(y,K, nu, Sigma_0, delta, z, mu, V_T)
    mu_tree_list = sample_T_rjmcmc(g, p0, K, z, Sigma, mu, verbose, delta, tuning)
    mu = mu_tree_list$mu
    g = mu_tree_list$g

    V_T = getTreeNodes(g, K_tilde)$name

    ## save results 
    res$mu[i,,] = mu
    res$pi[i,] = pi
    res$z[i,] = z
    res$Sigma[i,,]=Sigma
    write_graph(g, paste("saved_trajectories/rjmcmc_trees/tree", i, ".graphml", sep = ""), format = "graphml")
  }
  return(res)
}

