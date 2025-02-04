#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
unsigned int sample_Categorical(const vec& probs) {
  // sample from an unnormalized categorical distribution
  vec normprobs = probs/accu(probs);
  NumericVector p = as<NumericVector>(wrap(normprobs));
  int k = p.size();
  IntegerVector ans(k);
  rmultinom(1, p.begin(), k, ans.begin());
  arma::vec sample = as<arma::vec>(wrap(ans));
  arma::uvec ind = find(sample);
  return (int)ind[0];
}

// create a diagonal matrix, where the diagonal entries 
// are the escape probability vector from the set of 
// visited vertices to the set of unvisited verties
// [[Rcpp::export]]
mat create_Escape(const mat& W, const uvec& visited, 
                 const uvec& unvisited){
  
  vec degree = sum(W, 1);
  vec denominator = degree.rows(visited);
  mat boundaryEdges = W.submat(visited, unvisited);
  vec numerator = sum(boundaryEdges, 1); 
  return diagmat(numerator/denominator);

}

// create transition probability matrix 
// for transitions within the set of visited vertices
// note that this is just the principal submatrix corresponding
// to the visited vertices of the overall transition matrix 
// [[Rcpp::export]]
mat create_PA(const mat& W, const uvec& visited){
  mat within_transitions = W.submat(visited, visited);
  return within_transitions;

}

// create normalized transition probability matrix 
// for transitions from visited vertices to not yet visited vertices
mat create_PA_Escape(const mat& W, const uvec& visited, const uvec& unvisited){
  mat between_transitions = W.submat(visited, unvisited);
  return normalise(between_transitions, 1, 1);
}


// compute the distribution on the set of visited vertices 
// that indicate the probability of the walk eventually exiting
// the set of visited vertices via that vertex
// NOTE: In practice we use the compute_x_Escape_A_solve version below
// this exact version is just included here for completeness
vec compute_x_Escape_A(const mat& P_A, const mat& x_escape, const uvec& ind){
  mat I = eye(ind.n_elem, ind.n_elem);
  return x_escape* inv(I - P_A.t()) * ind;
}

// does the same as compute_x_Escape_A, except using an exact numerical 
// solve routine for the inversion that is often much more efficient
vec compute_x_Escape_A_solve(const mat& P_A, const mat& x_escape, const uvec& ind){
  mat I = eye(ind.n_elem, ind.n_elem);
  vec solved = arma::solve((I - P_A.t()), arma::conv_to<arma::vec>::from(ind));
  return x_escape* solved;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat randomWalkCover(const mat&W, int start) {
  
  // this is the Aldous-Broder algorithm
  
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
    vec neighbors = W.row(u).t(); 
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


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat fastCoverThreshold(const mat&W, int start, int threshold, bool verbose) {
  
  // this is the fast forward cover algorithm, that takes in a threshold argument
  // to determine when to accelerate. Every time we revisit a visited node, a counter
  // is incremented. When the counter hits the threshold value, the fast-forwarded
  // step is deployed (and the counter is reset to 0). 
  // returns an undirected adjacency matrix representing the structure of the spanning
  // tree. 
  
  
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
        cout << "Counter hit threshold. Accelerate!" << endl;
      }
      revisit_counter = 0;
      
      // set up indicator vector for current location
      uvec ind(visited.n_elem, arma::fill::zeros);
      uvec loc = find(visited == v);
      ind(loc[0]) = 1;
      
      //create escape probability diagonal matrix 
      mat x_escape = create_Escape(W, visited, unvisited);
      
      //create transition matrix within visited vertices
      mat P_A = create_PA(W, visited);
      
      //create transition matrix from visited vertices to unvisited vertices
      mat P_A_escape = create_PA_Escape(W, visited, unvisited);
      
      //compute distribution over vertices that we exit *FROM*
      vec x_exit_A = compute_x_Escape_A_solve(P_A, x_escape, ind);
      
      //sample the vertex that we exit from
      int m_loc = sample_Categorical(x_exit_A);
      int m = visited(m_loc);
      
      //sample the vertex that we exit to
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
  
  // Construct an undirected adjacency matrix representing the
  // spanning tree skeleton structure and return as output of function
  // (for users who want a directed tree, they can further 
  // process this undirected adjacency matrix to obtain the 
  // directed tree adjacency matrix they want under their context)
  for (int u = 1; u < n; ++u)
  {
    A_T_(u, Prev(u)) = 1;
    A_T_(Prev(u), u) = 1;
  }
  
  return A_T_;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
int fastCoverThresholdTimed(const mat&W, int start, int threshold, bool verbose) {
  
  // a timed version of the fastCoverThresholdCoverTimed function. 
  // here time is measured in the number of random walk steps taken. 
  // the cover time is returned. 
  
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
  int cover_time = 0;
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
      
    } else {
      
      if (verbose) {
        cout << "Counter hit threshold. Accelerate!" << endl;
      }
      revisit_counter = 0;
      
      // set up indicator vector for current location
      uvec ind(visited.n_elem, arma::fill::zeros);
      uvec loc = find(visited == v);
      ind(loc[0]) = 1;
      
      //create escape probability diagonal matrix 
      mat x_escape = create_Escape(W, visited, unvisited);
      
      //create transition matrix within visited vertices
      mat P_A = create_PA(W, visited);
      
      //create transition matrix from visited vertices to unvisited vertices
      mat P_A_escape = create_PA_Escape(W, visited, unvisited);
      
      //compute distribution over vertices that we exit *FROM*
      vec x_exit_A = compute_x_Escape_A_solve(P_A, x_escape, ind);
      
      //sample the vertex that we exit from
      int m_loc = sample_Categorical(x_exit_A);
      int m = visited(m_loc);
      
      //sample the vertex that we exit to
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
  return cover_time;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
int randomWalkCoverTimed(const mat&W, int start) {
  
  // a timed version of the randomWalkCover function. 
  // here time is measured in the number of random walk steps taken. 
  // the cover time is returned. 
  
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
  return cover_time;
}
