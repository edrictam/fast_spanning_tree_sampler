## Gibbs Sampler Full Conditionals
library(gtools)
library(mvtnorm)
library(invgamma)
library(LaplacesDemon)
library(Rcpp)
library(RcppArmadillo)

## Load fastcover cpp tree sampler
sourceCpp("helpers/fastcover/fastcover_older_version.cpp")

# source("helpers/helper.R")

## z sampler
sample_z = function(y,pi,mu, Sigma){
  n = dim(y)[1]
  p = dim(y)[2]
  
  K_tilde = length(pi)
  p.z.given.y = matrix(rep(0, n*K_tilde), ncol = K_tilde)
  for (i in 1:n){
    for (k in 2:K_tilde){
      p.z.given.y[i, k] = pi[k]*dmvnorm(y[i,], mu[k,], Sigma)
    }
    denominator = sum(p.z.given.y[i,])
    p.z.given.y[i, ] = p.z.given.y[i, ]/denominator
  }
  z = rep(0, n)
  for(i in 1:n){
    z[i] = sample(1:K_tilde, size=1,prob=p.z.given.y[i,],replace=TRUE)
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

## mu sampler
sample_mu = function(y, z, K, p, curr_mu, Sigma, delta, spanning){
  mu = curr_mu
  
  for(k in 2:(K+1)){
 
    sample.size = sum(z==k)
    delta_k = delta + sample.size
    
    sample.sum = rep(0, p)
    if (sample.size > 1) {
      sample.sum = as.matrix(colSums(y[z==k,], ))
    }
    if (sample.size == 1) {
      sample.sum =  as.matrix(y[z==k,])
    }

    neighbors = mu[spanning[k,] == 1, ]
    if (is.vector(neighbors)){
      neighbors.sum = as.matrix(neighbors)
    }else{
      neighbors.sum = as.matrix(colSums(neighbors))
    }
    num_neighbors = sum(spanning[k,] == 1)
    neighbors.avg = neighbors.sum/num_neighbors
    
    post.mean = (delta*neighbors.avg + sample.sum)/(delta_k)
    # post.mean = (delta*neighbors.sum + sample.sum)/(delta_k)
    mu[k,] = rmvnorm(1, post.mean, Sigma/delta_k)
  }
  return(mu)
}

## Sigma sampler
sample_Sigma = function(y,K, nu, Sigma_0, z, mu){
  n = dim(y)[1]
  p = dim(y)[2]
  Sigma_new = Sigma_0
  nu_new = nu
  for(k in 2:(K+ 1)){
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




## sample T
sample_T = function(y, K, p, mu, delta, Sigma){
  init = 0
  n = dim(y)[1]
  A = matrix(0, K + 1, K + 1)
  for (i in 1:(K + 1)){
    for (j in 1:(K + 1)){
      if (i < j) {
        A[i, j] = dmvnorm(mu[i,], mu[j,], (1/(delta))*Sigma)
        A[j, i] = A[i, j]
      }
      
    }
  }
  spanning.tree <- fastCoverThreshold(A, init, K + 1, verbose = FALSE)
  return(spanning.tree)
}




## Gibbs sampler function
gibbs_new = function(y, K, niter =1000, nu, Sigma_0, delta){
  
  print("In Gibbs function")
  ## initialize variables
  n = dim(y)[1]
  p = dim(y)[2]
  
  ## set up starting points
  pi = sample_pi(rep(0, n), K)
  mu = matrix(rep(0, K*p), K, p)
  mu_root = rep(0, p)
  mu = rbind(mu_root, mu)
  Sigma = 0.2*0.2*diag(p)
  z = sample_z(y,pi,mu,Sigma)
  spanning = sample_T(y, K, p, mu, delta, Sigma)
  
  ## set up array to store samples
  res = list(mu=array(0, dim = c(niter, K + 1, p)), 
             pi = matrix(nrow=niter,ncol=K + 1), 
             z = matrix(nrow=niter, ncol=n), 
             Sigma = array(0, dim = c(niter, p, p)),  
             spanning = array(0, dim = c(niter, K + 1, K + 1)))
  
  ## initialize chain
  res$mu[1,,]=mu
  res$pi[1,]=pi
  res$z[1,]=z 
  res$Sigma[1,,]=Sigma
  res$spanning[1,,]=spanning
  
  print("Initialization finished")
  
  for(i in 2:niter){
    if (i %% 100 == 0){
      print("In Gibbs loop iteration number: ")
      print(i)
    }
    ## sample from full conditionals
    pi = sample_pi(z,K)
    z = sample_z(y,pi,mu, Sigma)
    Sigma = sample_Sigma(y,K, nu, Sigma_0, z, mu)
    mu = sample_mu(y, z, K, p, mu, Sigma, delta, spanning)
    spanning = sample_T(y, K, p, mu, delta, Sigma)
    
    
    ## save results 
    res$mu[i,,] = mu
    res$pi[i,] = pi
    res$z[i,] = z
    res$Sigma[i,,]=Sigma
    res$spanning[i,,]=spanning
  }
  return(res)
}
