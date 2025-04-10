# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

C_make_glcm_counts <- function(x, n_levels, shift, na_rm) {
    .Call(`_GLCMTextures_C_make_glcm_counts`, x, n_levels, shift, na_rm)
}

C_make_glcm <- function(x, n_levels, shift, na_rm) {
    .Call(`_GLCMTextures_C_make_glcm`, x, n_levels, shift, na_rm)
}

C_GLSV <- function(Pij, n_levels) {
    .Call(`_GLCMTextures_C_GLSV`, Pij, n_levels)
}

C_glcm_metrics <- function(Pij, i_mat, j_mat, n_levels, k_vals, metrics, impute_corr) {
    .Call(`_GLCMTextures_C_glcm_metrics`, Pij, i_mat, j_mat, n_levels, k_vals, metrics, impute_corr)
}

C_glcm_textures_helper <- function(x, w2, n_levels, shift, metrics, na_rm, impute_corr, ni, nw) {
    .Call(`_GLCMTextures_C_glcm_textures_helper`, x, w2, n_levels, shift, metrics, na_rm, impute_corr, ni, nw)
}

