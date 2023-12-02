#' Calculates the GLCM Texture Metrics from a GLCM
#'
#' @param GLCM A numeric matrix representing a Normalized GLCM.
#' @param metrics A vector of texture metrics to return. Valid entries include "glcm_contrast", "glcm_dissimilarity", "glcm_homogeneity", "glcm_ASM", "glcm_entropy", "glcm_mean", "glcm_variance", "glcm_correlation".
#' @return GLCM based texture measures as a numeric vector.
#' @examples
#' test_matrix<- matrix(data=c(2,0,1,3,0,0,0,3,2), nrow = 3, ncol=3)
#' horizontal_glcm<- make_glcm(test_matrix, n_levels = 4,
#' shift = c(1,0), normalize = TRUE)
#' metrics<-glcm_metrics(horizontal_glcm, metrics= c("glcm_contrast",
#' "glcm_dissimilarity", "glcm_homogeneity", "glcm_ASM",
#' "glcm_entropy", "glcm_mean", "glcm_variance", "glcm_correlation"))
#' @references
#' Hall-Beyer, M., 2017. GLCM Texture: A Tutorial v. 3.0. University of Calgary, Alberta, Canada.
#'
#' Haralick, R.M., Shanmugam, K., Dinstein, I., 1973. Textural features for image classification. IEEE Transactions on Systems, Man, and Cybernetics 610–621. https://doi.org/10.1109/TSMC.1973.4309314
#' @export
glcm_metrics<-function(GLCM, metrics= c("glcm_contrast", "glcm_dissimilarity", "glcm_homogeneity", "glcm_ASM", "glcm_entropy", "glcm_mean", "glcm_variance", "glcm_correlation", "glcm_SA")){
  all_metrics<- c("glcm_contrast", "glcm_dissimilarity", "glcm_homogeneity", "glcm_ASM", "glcm_entropy", "glcm_mean", "glcm_variance", "glcm_correlation", "glcm_SA")
  if (any(!(metrics %in% all_metrics))){
    stop("Error: Invlaid metric. Valid metrics include 'glcm_contrast', 'glcm_dissimilarity', 'glcm_homogeneity', 'glcm_ASM', 'glcm_entropy', 'glcm_mean', 'glcm_variance', 'glcm_correlation', 'glcm_SA'")
  }
  needed_metrics<- metrics
  if(("glcm_variance" %in% needed_metrics) & (!("glcm_mean" %in% needed_metrics))){
    needed_metrics<- c(needed_metrics, "glcm_mean")
  }
  if(("glcm_correlation" %in% needed_metrics) & (!("glcm_mean" %in% needed_metrics))){
    needed_metrics<- c(needed_metrics, "glcm_mean")
  }
  if(("glcm_correlation" %in% needed_metrics) & (!("glcm_variance" %in% needed_metrics))){
    needed_metrics<- c(needed_metrics, "glcm_variance")
  } #Some metrics needed

  i_mat<- matrix(data= 0:(nrow(GLCM)-1), nrow= nrow(GLCM), ncol=ncol(GLCM), byrow=FALSE)
  j_mat<- matrix(data= 0:(ncol(GLCM)-1), nrow= nrow(GLCM), ncol=ncol(GLCM), byrow=TRUE)
  n_levels=nrow(GLCM)
  k_vals<- seq((2*n_levels)-1)-1
  out<- C_glcm_metrics(Pij = GLCM, i_mat = i_mat, j_mat = j_mat, n_levels=n_levels, k_vals=k_vals, metrics = needed_metrics)
  out<- out[names(out) %in% metrics] #Subset to requsted metrics
  return(out)
  }
