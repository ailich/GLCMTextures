#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//Make GLCM (non-normalized)
// [[Rcpp::export]]
IntegerMatrix C_make_glcm_counts(IntegerMatrix x, int n_levels, IntegerVector shift, bool na_rm){
  IntegerMatrix GLCM(n_levels, n_levels);//initialize GLCM
  int nr= x.nrow();
  int nc= x.ncol();

  if((!na_rm) && (is_true(any(is_na(x))))){
    GLCM.fill(NA_INTEGER);
    return GLCM;
  }  //Return window of NA's if any vals in window are NA

  if(na_rm && (is_true(all(is_na(x))))){
    GLCM.fill(NA_INTEGER);
    return GLCM;
  }  //Return window of NA's if all values in window is NA
  IntegerVector focalval_neighborval(2);
  focalval_neighborval.names() = CharacterVector::create("focal_val", "neighborval");
  for(int i=0; i < nr; ++i){
    for(int j=0; j < nc; ++j){
      focalval_neighborval["focal_val"] = x(i,j); //focal val
      IntegerVector neighbor_idx = {i-shift[1], j+shift[0]};
      if((neighbor_idx[0] < nr) && (neighbor_idx[1] < nc) && (neighbor_idx[0] >= 0) && (neighbor_idx[1] >= 0)){
        focalval_neighborval["neighborval"] = x(neighbor_idx[0], neighbor_idx[1]); //neighbor val
        if(is_false(any(is_na(focalval_neighborval)))){
          GLCM(focalval_neighborval["focal_val"],focalval_neighborval["neighborval"])++;
          GLCM(focalval_neighborval["neighborval"],focalval_neighborval["focal_val"])++;
        }
      }}}
  return GLCM;
}

//Make GLCM (normalized)
// [[Rcpp::export]]
arma::mat C_make_glcm(IntegerMatrix x, int n_levels, IntegerVector shift, bool na_rm){
  IntegerMatrix GLCM = C_make_glcm_counts(x, n_levels, shift, na_rm); //tabulate counts
  arma::mat GLCM_Norm(n_levels, n_levels);
  if(is_true(any(is_na(GLCM)))){
    GLCM_Norm.fill(NA_REAL);
    return GLCM_Norm;
  } //If GLCM is NA return NA

  double total= sum(GLCM);
  for(int m=0; m < GLCM.size(); ++m){
    GLCM_Norm[m] = (GLCM[m]*1.0)/total;
  } //Normalize
  return GLCM_Norm;
}

//Gray Level Sum Vector
// [[Rcpp::export]]
NumericVector C_GLSV(arma::mat Pij, int n_levels) {
  NumericVector k_prob(2*n_levels-1); //Initialize at zero
  for (int i=0; i < n_levels; ++i) {
    for (int j=0; j <= i; ++j) {
      int k = i+j;
      if(i==j){
        k_prob[k]= k_prob[k] + Pij(i,j); //add probability to corresponding spot in k_prob
      } else{
        k_prob[k]= k_prob[k] + 2*Pij(i,j); //Since only looping through lower triangle have to count off diagonal elements twice
      }
    }
  }
  return k_prob;
}

//Calculate Texture Metrics
// [[Rcpp::export]]
NumericVector C_glcm_metrics(arma::mat Pij, arma::mat i_mat, arma::mat j_mat, int n_levels, NumericVector k_vals, CharacterVector metrics){

  NumericVector textures = rep(NA_REAL,metrics.length());
  textures.names() = metrics;
  if(!is_finite(Pij)){
    return textures;}
  if(in(CharacterVector::create("glcm_contrast"), metrics)){
    textures["glcm_contrast"] = accu(Pij % pow(i_mat-j_mat,2)); //Contrast= sum(P_ij*(i-j)^2)
    }
  if(in(CharacterVector::create("glcm_dissimilarity"), metrics)){
    textures["glcm_dissimilarity"]= accu(Pij % abs(i_mat-j_mat)); //Dissimilarity= sum(P_ij*|i-j|)
    }
  if(in(CharacterVector::create("glcm_homogeneity"), metrics)){
    textures["glcm_homogeneity"]= accu(Pij/(1+(pow(i_mat-j_mat,2)))); //Homogeneity= sum(P_ij / (1+(i-j)^2))
    }
  if(in(CharacterVector::create("glcm_ASM"), metrics)){
    textures["glcm_ASM"]= accu(pow(Pij,2)); //ASM= sum(P_ij^2)
    }
  if(in(CharacterVector::create("glcm_mean"), metrics)){
    textures["glcm_mean"]= accu(Pij % i_mat); //mean= sum(i*(P_ij))
    }
  if(in(CharacterVector::create("glcm_variance"), metrics)){
    textures["glcm_variance"]= accu(Pij % pow(i_mat-textures["glcm_mean"],2)); //varaince= sum(P_ij*(i-u)^2
    }
  if(in(CharacterVector::create("glcm_correlation"), metrics)){
    textures["glcm_correlation"]= accu(Pij % (((i_mat-textures["glcm_mean"]) % (j_mat-textures["glcm_mean"]))/(textures["glcm_variance"]))); //Correlation= sum(P_ij*[((i-u)*(j-u))/(var)])
    }
  if(in(CharacterVector::create("glcm_entropy"), metrics)){
    arma::mat glcm_entropy_mat = Pij % ((-1) * log(Pij));
    glcm_entropy_mat.replace(datum::nan,0.0);
    textures["glcm_entropy"]= accu(glcm_entropy_mat); //Entropy= sum(P_ij * (-ln(P_ij))) ; 0*ln(0)=0
    }
  //if(in(CharacterVector::create("glcm_SA"), metrics)){
    //NumericVector k_prob = C_GLSV(Pij, n_levels); //Gray Level Sum Vector
    //textures["glcm_SA"]= sum(k_prob*k_vals); //GLCM_SumAverage= sum(k*k_prob)
    //}
  return(textures);
}

//GLCM across matrix using sliding window (terra)
// [[Rcpp::export]]
NumericMatrix C_glcm_textures_helper(IntegerVector x, IntegerVector w2, int n_levels, IntegerVector shift, CharacterVector metrics, bool na_rm, size_t ni, size_t nw){

  NumericMatrix out(ni, metrics.length());
  out.fill(NA_REAL);
  colnames(out)= metrics;

  arma::mat i_mat(n_levels,n_levels);
  arma::mat j_mat(n_levels,n_levels);
  NumericVector k_vals((2*n_levels)-1);
  for (int k = 1; k < k_vals.length(); ++k) {
    k_vals[k] = k;
  } // All possible values of k, where k=i+j
  for(int i=0; i<n_levels; ++i){
    for(int j=0; j<n_levels; ++j){
      i_mat(i,j)=i;
      j_mat(i,j)=j;
    }} //Set up i_mat and j_mat outside of loop so they don't need to be redefined each time

  for(size_t i=0; i<ni; i++) {
    size_t start = i*nw;
    size_t end = start+nw-1;
    IntegerVector xw = x[Rcpp::Range(start,end)]; //Current window of values
    IntegerMatrix curr_window(w2[0],w2[1]);
    for(int r=0; r < w2[0]; r++){
      for(int c=0; c < w2[1]; c++){
        curr_window(r,c) = xw[r*(w2[1])+c];
      }
    } //fill in matrix by row
    arma::mat curr_GLCM = C_make_glcm(curr_window, n_levels, shift, na_rm); //Tabulate the GLCM
    out(i, _) =  C_glcm_metrics(curr_GLCM, i_mat, j_mat, n_levels, k_vals, metrics);
  }
  return(out);
}
