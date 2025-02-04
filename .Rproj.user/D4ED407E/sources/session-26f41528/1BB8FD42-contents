#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
// sample from an unnormalized categorical dist
// replaced gumbelMax since the log returns an issue when the edge is 0
unsigned int sample_Categorical(const vec& probs) {
  vec normprobs = probs/accu(probs);
  NumericVector p = as<NumericVector>(wrap(normprobs));
  int k = p.size();
  IntegerVector ans(k);
  rmultinom(1, p.begin(), k, ans.begin());
  arma::vec sample = as<arma::vec>(wrap(ans));
  arma::uvec ind = find(sample);
  return (int)ind[0];
}

// create an escape probability vector from the set of visited vertices
mat create_Escape(const mat& W, const uvec& visited, 
                  const uvec& unvisited){
  
  rowvec degree = sum(W, 0);
  rowvec denominator = degree.cols(visited);
  mat boundaryEdges = W.submat(visited, unvisited);
  vec numerator = sum(boundaryEdges, 1); 
  return diagmat(numerator.t()/denominator);
  
}

// create normalized transition probability vector for transitions within the set of visited vertices
mat create_PA(const mat& W, const uvec& visited){
  mat within_transitions = W.submat(visited, visited);
  return normalise(within_transitions, 1, 1);
  
}

// create normalized transition probability vector for transitions from visited vertices to not yet visited vertices
mat create_PA_Escape(const mat& W, const uvec& visited, const uvec& unvisited){
  mat between_transitions = W.submat(visited, unvisited);
  return normalise(between_transitions, 1, 1);
  
}

// compute the distribution on the set of visited vertices that indicate 
// the probability that the walk eventually exits via that vertex
vec compute_x_Escape_A(const mat& P_A, const mat& x_escape, const uvec& ind){
  mat I = eye(ind.n_elem, ind.n_elem);
  return x_escape* inv(I - P_A.t() * (I - x_escape)) * ind;
}

// compute the distribution on the set of visited vertices that indicate 
// the probability that the walk eventually exits via that vertex
vec compute_x_Escape_A_solve(const mat& P_A, const mat& x_escape, const uvec& ind){
  mat I = eye(ind.n_elem, ind.n_elem);
  vec temp_product_solved = arma::solve((I - P_A.t() * (I - x_escape)), arma::conv_to<arma::vec>::from(ind));
  return x_escape* temp_product_solved;
}

// this function implements aldous-broder
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat randomWalkCover(const mat&W, int start) {
  // Initialize variables
  int n = W.n_rows;
  arma::mat A_T_(n, n, arma::fill::zeros);
  arma::Col<int> InTree(n, arma::fill::zeros);
  arma::Col<int> Prev(n);
  
  // Set up Next and InTree
  int r = start;
  InTree(r) = 1;
  
  int u = r;
  int NewNodeVisited = 1;
  
  int v = 0;
  
  while (NewNodeVisited<n) {
    
    //do a random walk until a new node is reached
    vec neighbors = W.col(u); 
    v = sample_Categorical(neighbors);
    
    if(!InTree(v)){
      Prev(v)= u;
      InTree(v)=1;
      NewNodeVisited += 1;
    };
    u = v;
  }
  
  // Construct adjacency matrix
  for (int u = 1; u < n; ++u)
  {
    A_T_(u, Prev(u)) = 1;
    A_T_(Prev(u), u) = 1;
  } 
  
  return A_T_;
}

// this function implements our fast-forwarded cover algorithm
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat fastCoverThreshold(const mat&W, int start, int threshold, bool verbose) {
  // Initialize variables
  int n = W.n_rows;
  
  uvec dummies(n-1, fill::ones);
  uvec nodes(n, fill::zeros);
  nodes.subvec(1, n-1) = cumsum(dummies);
  arma::mat A_T_(n, n, arma::fill::zeros);
  arma::Col<int> InTree(n, arma::fill::zeros);
  arma::Col<int> Prev(n);
  
  // Set up Next and InTree
  int r = start;
  InTree(r) = 1;
  int numVisited = 1;
  uvec visited = nodes.elem(find(InTree == 1));
  uvec unvisited = nodes.elem(find(InTree == 0));
  
  // Start Algo
  int u = r;
  int v;
  
  int revisit_counter = 0;
  while (numVisited<n) {
    
    vec neighbors = W.col(u);
    v = sample_Categorical(neighbors);
    
    if(InTree(v) == 0){
      
      Prev(v)= u;
      InTree(v)=1;
      numVisited += 1;
      visited = nodes.elem(find(InTree == 1));
      unvisited = nodes.elem(find(InTree == 0));
      
    } else if((InTree(v) == 1) && (revisit_counter < threshold)){
      revisit_counter += 1;
     
    } else {
      if (verbose) {
        cout << "Counter hit threshold. Accelerate!" <<endl;
      }
      revisit_counter = 0;

      uvec ind(visited.n_elem, arma::fill::zeros);
      uvec loc = find(visited == v);
      ind(loc[0]) = 1;

      //create escape probability vec
      mat x_escape = create_Escape(W, visited, unvisited);

      //create P_A
      mat P_A = create_PA(W, visited);
  
      //create P_{A,A^C}
      mat P_A_escape = create_PA_Escape(W, visited, unvisited);

      //Compute x_exit, A
      vec x_exit_A = compute_x_Escape_A(P_A, x_escape, ind);

      //sample m
      int m_loc = sample_Categorical(x_exit_A);
      int m = visited(m_loc);

      //sample m*
      vec ind2(visited.n_elem, arma::fill::zeros);
      ind2(m_loc) = 1;
      int m_star_loc = sample_Categorical(P_A_escape.t() * ind2);

      int m_star = unvisited(m_star_loc);
      
      u = m;
      v = m_star;
      Prev(v)= u;
      InTree(v)=1;
      numVisited += 1;
      visited = nodes.elem(find(InTree == 1));
      unvisited = nodes.elem(find(InTree == 0));
    }
    u = v;
  }
  if (verbose){
    cout << "finished" << endl;
  }

  // Construct adjacency matrix
  for (int u = 1; u < n; ++u)
  {
    A_T_(u, Prev(u)) = 1;
    A_T_(Prev(u), u) = 1;
  }
  
  return A_T_;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
int fastCoverThresholdCoverTimed(const mat&W, int start, int threshold, bool verbose) {
  // Initialize variables
  int n = W.n_rows;
  
  uvec dummies(n-1, fill::ones);
  uvec nodes(n, fill::zeros);
  nodes.subvec(1, n-1) = cumsum(dummies);
  arma::mat A_T_(n, n, arma::fill::zeros);
  arma::Col<int> InTree(n, arma::fill::zeros);
  arma::Col<int> Prev(n);
  
  // Set up Next and InTree
  int r = start;
  InTree(r) = 1;
  int numVisited = 1;
  uvec visited = nodes.elem(find(InTree == 1));
  uvec unvisited = nodes.elem(find(InTree == 0));
  
  // Start Algo
  int u = r;
  int v;
  int cover_time = 0;
  int revisit_counter = 0;
  while (numVisited<n) {
    
    vec neighbors = W.col(u);
    v = sample_Categorical(neighbors);
    cover_time += 1;
    if(InTree(v) == 0){
      
      Prev(v)= u;
      InTree(v)=1;
      numVisited += 1;
      visited = nodes.elem(find(InTree == 1));
      unvisited = nodes.elem(find(InTree == 0));
      
    } else if((InTree(v) == 1) && (revisit_counter < threshold)){
      revisit_counter += 1;
      if (verbose) {
        cout << "counter hit!" <<endl;
      }
    } else {
      revisit_counter = 0;
      uvec ind(visited.n_elem, arma::fill::zeros);
      uvec loc = find(visited == v);
      ind(loc[0]) = 1;

      //create escape probability vec
      mat x_escape = create_Escape(W, visited, unvisited);

      //create P_A
      mat P_A = create_PA(W, visited);

      //create P_{A,A^C}
      mat P_A_escape = create_PA_Escape(W, visited, unvisited);
  
      //Compute x_exit, A
      vec x_exit_A = compute_x_Escape_A_solve(P_A, x_escape, ind);

      //sample m
      int m_loc = sample_Categorical(x_exit_A);
      int m = visited(m_loc);
 
      //sample m*
      vec ind2(visited.n_elem, arma::fill::zeros);
      ind2(m_loc) = 1;
      int m_star_loc = sample_Categorical(P_A_escape.t() * ind2);

      int m_star = unvisited(m_star_loc);
      
      u = m;
      v = m_star;
      Prev(v)= u;
      InTree(v)=1;
      numVisited += 1;
      visited = nodes.elem(find(InTree == 1));
      unvisited = nodes.elem(find(InTree == 0));
    }
    u = v;
  }
  
  // Construct adjacency matrix
  for (int u = 1; u < n; ++u)
  {
    A_T_(u, Prev(u)) = 1;
    A_T_(Prev(u), u) = 1;
  }
  
  return cover_time;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
int randomWalkCoverTimed(const mat&W, int start) {
  // Initialize variables
  int n = W.n_rows;
  arma::mat A_T_(n, n, arma::fill::zeros);
  arma::Col<int> InTree(n, arma::fill::zeros);
  arma::Col<int> Prev(n);
  
  // Set up Next and InTree
  int r = start;
  InTree(r) = 1;
  
  int u = r;
  int NewNodeVisited = 1;
  
  int v = 0;
  int cover_time = 0;
  while (NewNodeVisited<n) {
    cover_time += 1;
    //do a random walk until a new node is reached
    vec neighbors = W.col(u); 
    v = sample_Categorical(neighbors);
    
    if(!InTree(v)){
      Prev(v)= u;
      InTree(v)=1;
      NewNodeVisited += 1;
    };
    u = v;
  }
  
  // Construct adjacency matrix
  for (int u = 1; u < n; ++u)
  {
    A_T_(u, Prev(u)) = 1;
    A_T_(Prev(u), u) = 1;
  } 
  
  return cover_time;
}
