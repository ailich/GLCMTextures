#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//Extracts relevant window from matrix based on position of central pixel and window size
// [[Rcpp::export]]
IntegerMatrix C_extract_window_int(IntegerMatrix r, IntegerVector w, IntegerVector idx){
  int nr = w(0);
  int nc = w(1);

  int rast_row_center = idx(0);
  int rast_row_top = rast_row_center-(nr-1)/2;

  int rast_col_center= idx(1);
  int rast_col_left = rast_col_center-(nc-1)/2;

  IntegerMatrix r_idx(nr, nc);
  IntegerMatrix c_idx(nr, nc);
  for(int i=0; i < nr; ++i){
    for(int j=0; j < nc; ++j){
      r_idx(i,j) = i+rast_row_top;
      c_idx(i,j)=j+rast_col_left;
    }}
  IntegerMatrix dat(nr, nc);
  for(int k=0; k < dat.size(); ++k){
    dat[k]= r(r_idx[k], c_idx[k]);
  }
  return dat; //extracted window as a matrix
}

//Make GLCM
// [[Rcpp::export]]
NumericMatrix C_make_glcm(IntegerMatrix x, int n_levels, IntegerVector shift, String na_opt){
  StringVector na_options(3);
  na_options(0) = "any";
  na_options(1) = "center";
  na_options(2) = "all";

  if((na_opt== (na_options(0))) && (is_true(any(is_na(x))))){
    NumericMatrix GLCM_Norm(n_levels, n_levels);
    GLCM_Norm.fill(NA_REAL);
    return GLCM_Norm;
  }  //Return window of NA's if any vals in window are NA

  int nr= x.nrow();
  int nc= x.ncol();
  int cr= (nr-1)/2;//row position of center
  int cc= (nc-1)/2; //column position of center
  IntegerVector center_val(1);
  center_val(0) = x(cr,cc);
  if((na_opt== (na_options(1))) && (is_true(all(is_na(center_val))))){
    NumericMatrix GLCM_Norm(n_levels, n_levels);
    GLCM_Norm.fill(NA_REAL);
    return GLCM_Norm;
  }  //Return window of NA's if center value in window is NA

  if((na_opt== (na_options(2))) && (is_true(all(is_na(x))))){
    NumericMatrix GLCM_Norm(n_levels, n_levels);
    GLCM_Norm.fill(NA_REAL);
    return GLCM_Norm;
  }  //Return window of NA's if all values in window is NA

  IntegerMatrix GLCM(n_levels, n_levels);//initialize GLCM
  for(int i=0; i < nr; ++i){
    for(int j=0; j < nc; ++j){
      int focal_val = x(i,j);
      IntegerVector neighbor_idx(2, 0);
      neighbor_idx(0) = i-shift(1);
      neighbor_idx(1) = j+shift(0);
      if((neighbor_idx(0) < nr) && (neighbor_idx(1) < nc) && (neighbor_idx(0) >= 0) && (neighbor_idx(1) >= 0)){
        int neighbor_val = x(neighbor_idx(0), neighbor_idx(1));
        IntegerVector focalval_neighborval(2);
        focalval_neighborval(0)=focal_val;
        focalval_neighborval(1)=neighbor_val;
        if(is_false(any(is_na(focalval_neighborval)))){
          GLCM(focal_val,neighbor_val) = GLCM(focal_val,neighbor_val)+1;
          GLCM(neighbor_val,focal_val) = GLCM(neighbor_val,focal_val)+1;
        }
      }}}

  double total=0;
  for(int k=0; k < GLCM.size(); ++k){
    total = total + GLCM[k];
  }
  NumericMatrix GLCM_Norm(n_levels, n_levels);
  for(int m=0; m < GLCM.size(); ++m){
    GLCM_Norm[m] = (GLCM[m]*1.0)/total;
  }
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

//GLCM across matrix using sliding window
// [[Rcpp::export]]
List C_glcm_textures_helper(IntegerMatrix rq, IntegerVector w, int n_levels, IntegerVector shift, String na_opt){
  int nr= rq.nrow();
  int nc= rq.ncol();
  int min_row = (w(0)-1)/2;
  int max_row = nr - ((w(0)-1)/2);
  int min_col = (w(1)-1)/2;
  int max_col = nc - ((w(1)-1)/2);
  NumericMatrix glcm_contrast = NumericMatrix(nr,nc);
  glcm_contrast.fill(NA_REAL);

  NumericMatrix glcm_dissimilarity = NumericMatrix(nr,nc);
  glcm_dissimilarity.fill(NA_REAL);

  NumericMatrix glcm_homogeneity = NumericMatrix(nr,nc);
  glcm_homogeneity.fill(NA_REAL);

  NumericMatrix glcm_ASM = NumericMatrix(nr,nc);
  glcm_ASM.fill(NA_REAL);

  NumericMatrix glcm_entropy = NumericMatrix(nr,nc);
  glcm_entropy.fill(NA_REAL);

  NumericMatrix glcm_mean = NumericMatrix(nr,nc);
  glcm_mean.fill(NA_REAL);

  NumericMatrix glcm_variance = NumericMatrix(nr,nc);
  glcm_variance.fill(NA_REAL);

  NumericMatrix glcm_correlation = NumericMatrix(nr,nc);
  glcm_correlation.fill(NA_REAL);

  for(int i = min_row; i< max_row; ++i) {
    for(int j = min_col; j < max_col; ++j){
      IntegerVector idx = IntegerVector(2);
      idx(0)= i;
      idx(1)=j;
      IntegerMatrix curr_window = C_extract_window_int(rq, w, idx);
      NumericMatrix curr_GLCM = C_make_glcm(curr_window, n_levels, shift, na_opt); //Tabulate the GLCM
      NumericVector curr_textures = C_glcm_metrics(curr_GLCM);
      glcm_contrast(i,j) = curr_textures["glcm_contrast"];
      glcm_dissimilarity(i,j) = curr_textures["glcm_dissimilarity"];
      glcm_homogeneity(i,j) = curr_textures["glcm_homogeneity"];
      glcm_ASM(i,j) = curr_textures["glcm_ASM"];
      glcm_entropy(i,j) = curr_textures["glcm_entropy"];
      glcm_mean(i,j) = curr_textures["glcm_mean"];
      glcm_variance(i,j) = curr_textures["glcm_variance"];
      glcm_correlation(i,j) = curr_textures["glcm_correlation"];
    }}
  List textures= List::create(_["glcm_contrast"]=glcm_contrast, _["glcm_dissimilarity"]=glcm_dissimilarity, _["glcm_homogeneity"]=glcm_homogeneity, _["glcm_ASM"]=glcm_ASM, _["glcm_entropy"]=glcm_entropy, _["glcm_mean"]=glcm_mean, _["glcm_variance"]=glcm_variance, _["glcm_correlation"]=glcm_correlation);
  return(textures);
}

//GLCM across matrix using sliding window
// [[Rcpp::export]]
NumericMatrix C_glcm_textures_helper2(IntegerMatrix rq, IntegerVector w, int n_levels, IntegerVector shift, String na_opt){
  int nr= rq.nrow();
  int nc= rq.ncol();
  int min_row = (w(0)-1)/2;
  int max_row = nr - ((w(0)-1)/2);
  int min_col = (w(1)-1)/2;
  int max_col = nc - ((w(1)-1)/2);
  int n_elem = nr * nc;

  NumericMatrix out = NumericMatrix(n_elem, 8);
  out.fill(NA_REAL);
  colnames(out)= CharacterVector::create("glcm_contrast", "glcm_dissimilarity", "glcm_homogeneity", "glcm_ASM", "glcm_entropy", "glcm_mean", "glcm_variance", "glcm_correlation");

  for(int i = min_row; i< max_row; ++i) {
    for(int j = min_col; j < max_col; ++j){
      IntegerVector idx = IntegerVector(2);
      idx(0)= i;
      idx(1)=j;
      int curr_elem_idx = i*nc + j; //rasters are indexed moving across rows
      IntegerMatrix curr_window = C_extract_window_int(rq, w, idx);
      NumericMatrix curr_GLCM = C_make_glcm(curr_window, n_levels, shift, na_opt); //Tabulate the GLCM
      NumericVector curr_textures = C_glcm_metrics(curr_GLCM);
      out(curr_elem_idx, 0) = curr_textures["glcm_contrast"];
      out(curr_elem_idx, 1) = curr_textures["glcm_dissimilarity"];
      out(curr_elem_idx, 2) = curr_textures["glcm_homogeneity"];
      out(curr_elem_idx, 3) = curr_textures["glcm_ASM"];
      out(curr_elem_idx, 4) = curr_textures["glcm_entropy"];
      out(curr_elem_idx, 5) = curr_textures["glcm_mean"];
      out(curr_elem_idx, 6) = curr_textures["glcm_variance"];
      out(curr_elem_idx, 7) = curr_textures["glcm_correlation"];
    }}
  return(out);
}
