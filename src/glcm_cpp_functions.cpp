#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Make GLCM
// [[Rcpp::export]]
arma::mat C_make_glcm(const IntegerVector& x,
                           const int n_levels,
                           const IntegerVector& shift,
                           const bool na_rm,
                           const int nrow,
                           const int ncol,
                           const bool normalize) {

  arma::mat GLCM(n_levels, n_levels, arma::fill::zeros); // single matrix from the start

  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      int idx = i * ncol + j;
      int focal = x[idx];
      int ni = i - shift[1];
      int nj = j + shift[0];
      if (ni >= 0 && ni < nrow && nj >= 0 && nj < ncol) {
        int nidx = ni * ncol + nj;
        int neighbor = x[nidx];
        if (focal != NA_INTEGER && neighbor != NA_INTEGER) {
          GLCM(focal, neighbor) += 1.0;
        } else if (!na_rm) {
          GLCM.fill(arma::datum::nan);
          return GLCM;
        }
      }
    }
  }

  // Make symmetric
  GLCM += GLCM.t();

  // Normalize in place if requested
  if (normalize) {
    double total = arma::accu(GLCM);
    if (total > 0.0) {
      GLCM /= total;
    } else {
      GLCM.fill(arma::datum::nan);
    }
  }

  return GLCM;
}

// [[Rcpp::export]]
NumericVector C_glcm_metrics(const arma::mat& Pij,
                              const arma::mat& i_mat,
                              const arma::mat& j_mat,
                              const arma::mat& i_minus_j,
                              const arma::mat& i_minus_j2,
                              const arma::mat& i_minus_j_abs,
                              const int n_levels,
                              const IntegerVector& metric_indices,
                              const bool impute_corr) {

  int n_metrics = metric_indices.length();
  NumericVector textures(n_metrics, NA_REAL);  // Initialize output

  bool all_na = true; // Check to see if all values are NA
  for (arma::uword i = 0; i < Pij.n_elem; ++i) {
    if (std::isfinite(Pij[i])) {
      all_na = false;
      break;
    }
  }
  if (all_na) {
    return textures;  // Already filled with NA_REAL
  }

  // Flags for what we need to compute
  bool need_mean = false, need_var = false;

  for (int k = 0; k < n_metrics; ++k) {
    int metric = metric_indices[k];
    if (metric == 5) need_mean = true;
    if (metric == 6 || metric == 7) {
      need_mean = true;
      need_var = true;
    }
  }

  // Mean and variance (only if needed)
  double glcm_mean = NA_REAL;
  if (need_mean) glcm_mean = arma::accu(Pij % i_mat);

  double glcm_variance = NA_REAL;
  if (need_var) glcm_variance = arma::accu(Pij % arma::square(i_mat - glcm_mean));

  // Compute requested metrics
  for (int m = 0; m < n_metrics; ++m) {
    int metric = metric_indices[m];
    switch (metric) {
    case 0: // Contrast
      textures[m] = arma::accu(Pij % i_minus_j2);
      break;
    case 1: // Dissimilarity
      textures[m] = arma::accu(Pij % i_minus_j_abs);
      break;
    case 2: // Homogeneity
      textures[m] = arma::accu(Pij / (1.0 + i_minus_j2));
      break;
    case 3: // ASM
      textures[m] = arma::accu(arma::square(Pij));
      break;
    case 4: { // Entropy (optimized)
        double entropy = 0.0;
        for (arma::uword i = 0; i < Pij.n_elem; ++i) {
          double p = Pij[i];
          if (p > 0.0) {
            entropy -= p * std::log(p);  // avoids log(0)
          }
        }
        textures[m] = entropy;
        break;
      }
    case 5: // Mean
      textures[m] = glcm_mean;
      break;
    case 6: // Variance
      textures[m] = glcm_variance;
      break;
    case 7: { // Correlation (optimized)
        if (glcm_variance == 0.0) {
        textures[m] = impute_corr ? 0.0 : R_NaN;
      } else {
        double corr_sum = 0.0;
        for (arma::uword i = 0; i < Pij.n_elem; ++i) {
          corr_sum += Pij[i] * (i_mat[i] - glcm_mean) * (j_mat[i] - glcm_mean);
        }
        textures[m] = corr_sum / glcm_variance;
      }
      break;
      }
    }
  }

  return textures;
}

// Focal function
// [[Rcpp::export]]
NumericMatrix C_glcm_textures_helper(const IntegerVector& x,
                                     const IntegerVector& w2,
                                     const int& n_levels,
                                     const List& shift_list,
                                     const IntegerVector& metric_indices,
                                     const bool na_rm,
                                     const bool impute_corr,
                                     size_t ni,
                                     size_t nw) {

  // Initialize output matrix
  NumericMatrix out(ni, metric_indices.length());
  out.fill(NA_REAL);

  // Precompute i_mat and j_mat
  arma::mat i_mat(n_levels, n_levels);
  arma::mat j_mat(n_levels, n_levels);
  for (int i = 0; i < n_levels; ++i) {
    for (int j = 0; j < n_levels; ++j) {
      i_mat(i, j) = i;
      j_mat(i, j) = j;
    }
  }

  const arma::mat i_minus_j = i_mat-j_mat;
  const arma::mat i_minus_j2 = arma::square(i_mat-j_mat);
  const arma::mat i_minus_j_abs = abs(i_minus_j);

  const int* x_ptr = x.begin();  // Pointer to the input vector

  const int nrow = w2[0];
  const int ncol = w2[1];

  for (size_t i = 0; i < ni; ++i) {
    const int* xw_ptr = x_ptr + i * nw;
    IntegerVector xw(xw_ptr, xw_ptr + nw);

    if ((!na_rm) && is_true(any(is_na(xw)))) continue;

    int n_metrics = metric_indices.length();
    NumericVector accum_metrics(n_metrics, 0.0);
    IntegerVector valid_shifts(n_metrics, 0);

    int n_shifts = shift_list.size();
    for (int shift_idx = 0; shift_idx < n_shifts; ++shift_idx) {
      IntegerVector shift = shift_list[shift_idx];

      arma::mat glcm = C_make_glcm(xw, n_levels, shift, na_rm, nrow, ncol, true);
      NumericVector metrics = C_glcm_metrics(glcm, i_mat, j_mat, i_minus_j, i_minus_j2, i_minus_j_abs, n_levels, metric_indices, impute_corr);

      for (int m = 0; m < n_metrics; ++m) {
        double val = metrics[m];
        if (na_rm) {
          if (!Rcpp::NumericVector::is_na(val)) {
            accum_metrics[m] += val;
            valid_shifts[m]++;
          }
        } else {
          accum_metrics[m] += val;
        }
      }
    }

    if (na_rm) {
      for (int m = 0; m < n_metrics; ++m) {
        if (valid_shifts[m] > 0) {
          accum_metrics[m] /= valid_shifts[m];
        } else {
          accum_metrics[m] = NA_REAL;
        }
      }
    } else {
      accum_metrics = accum_metrics / n_shifts;
    }

    out(i, _) = accum_metrics;
  }

  return out;
}
