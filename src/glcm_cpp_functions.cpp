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
          GLCM(focalval_neighborval["focal_val"],focalval_neighborval["neighborval"]) = GLCM(focalval_neighborval["focal_val"],focalval_neighborval["neighborval"])+1;
          GLCM(focalval_neighborval["neighborval"],focalval_neighborval["focal_val"]) = GLCM(focalval_neighborval["neighborval"],focalval_neighborval["focal_val"])+1;
        }
      }}}
  return GLCM;
}

//Make GLCM (normalized)
// [[Rcpp::export]]
NumericMatrix C_make_glcm(IntegerMatrix x, int n_levels, IntegerVector shift, bool na_rm){
  IntegerMatrix GLCM = C_make_glcm_counts(x, n_levels, shift, na_rm); //tabulate counts
  NumericMatrix GLCM_Norm(n_levels, n_levels);
  if(is_true(any(is_na(GLCM)))){
    GLCM_Norm.fill(NA_REAL);
    return GLCM_Norm;
  } //If GLCM is NA return NA

  double total=0;
  for(int k=0; k < GLCM.size(); ++k){
    total = total + GLCM[k];
  }
  for(int m=0; m < GLCM.size(); ++m){
    GLCM_Norm[m] = (GLCM[m]*1.0)/total;
  } //Normalize
  return GLCM_Norm;
}

//Calculate Texture Metrics
// [[Rcpp::export]]
NumericVector C_glcm_metrics(NumericMatrix GLCM){
  NumericVector textures = rep(NA_REAL,8);
  textures.names() = CharacterVector::create("glcm_contrast", "glcm_dissimilarity", "glcm_homogeneity","glcm_ASM","glcm_entropy","glcm_mean","glcm_variance","glcm_correlation");
  if(is_true(any(is_na(GLCM)))){
    return textures;}
  int GLCM_dims = GLCM.nrow(); //Size of GLCM (they are square so nrows=ncols)
  arma::mat Pij = as<arma::mat>(GLCM);
  arma::mat i_mat(GLCM_dims,GLCM_dims);
  arma::mat j_mat(GLCM_dims,GLCM_dims);
  for(int i=0; i<GLCM_dims; ++i){
    for(int j=0; j<GLCM_dims; ++j){
      i_mat(i,j)=i;
      j_mat(i,j)=j;
    }}
  textures["glcm_contrast"] = accu(Pij % pow(i_mat-j_mat,2)); //Contrast= sum(P_ij*(i-j)^2)
  textures["glcm_dissimilarity"]= accu(Pij % abs(i_mat-j_mat)); //Dissimilarity= sum(P_ij*|i-j|)
  textures["glcm_homogeneity"]= accu(Pij/(1+(pow(i_mat-j_mat,2)))); //Homogeneity= sum(P_ij / (1+(i-j)^2))
  textures["glcm_ASM"]= accu(pow(Pij,2)); //ASM= sum(P_ij^2)
  textures["glcm_mean"]= accu(Pij % i_mat); //mean= sum(i*(P_ij))
  textures["glcm_variance"]= accu(Pij % pow(i_mat-textures["glcm_mean"],2)); //varaince= sum(P_ij*(i-u)^2
  textures["glcm_correlation"]= accu(Pij % (((i_mat-textures["glcm_mean"]) % (j_mat-textures["glcm_mean"]))/(textures["glcm_variance"]))); //Correlation= sum(P_ij*[((i-u)*(j-u))/(var)])
  arma::mat glcm_entropy_mat = Pij % ((-1) * log(Pij));
  glcm_entropy_mat.replace(datum::nan,0.0);
  textures["glcm_entropy"]= accu(glcm_entropy_mat); //Entropy= sum(P_ij * (-ln(P_ij))) ; 0*ln(0)=0
  return(textures);
}

//GLCM across matrix using sliding window (terra)
// [[Rcpp::export]]
NumericMatrix C_glcm_textures_helper(IntegerVector x, IntegerVector w2, int n_levels, IntegerVector shift, bool na_rm, size_t ni, size_t nw){

  NumericMatrix out(ni, 8);
  out.fill(NA_REAL);
  colnames(out)= CharacterVector::create("glcm_contrast", "glcm_dissimilarity", "glcm_homogeneity", "glcm_ASM", "glcm_entropy", "glcm_mean", "glcm_variance", "glcm_correlation");

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
    NumericMatrix curr_GLCM = C_make_glcm(curr_window, n_levels, shift, na_rm); //Tabulate the GLCM
    NumericVector curr_textures = C_glcm_metrics(curr_GLCM);
    out(i, _) = curr_textures;
  }
  return(out);
}
