#' Calculates the GLCM Texture Metrics from a GLCM
#'
#' @param GLCM A numeric matrix representing a Normalized GLCM
#' @examples
#' test_matrix<- matrix(data=c(2,0,1,3,0,0,0,3,2), nrow = 3, ncol=3)
#' horizontal_glcm<- make_glcm(test_matrix, n_levels = 4, shift = c(1,0), normalize = TRUE)
#' metrics<-glcm_metrics(horizontal_glcm)
#' @references
#' Hall-Beyer, M., 2017. GLCM Textrure: A Tutorial v. 3.0. University of Calgary, Alberta, Canada.
#'
#' Haralick, R.M., Shanmugam, K., Dinstein, I., 1973. Textural features for image classification. IEEE Transactions on Systems, Man, and Cybernetics 610–621. https://doi.org/10.1109/TSMC.1973.4309314
#' @export
glcm_metrics<-function(GLCM){
  return(C_glcm_metrics(GLCM=GLCM))
  }
