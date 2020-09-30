#' Calculates the GLCM Texture Metrics from a GLCM
#' 
#' @param GLCM A numeric matrix representing a Normalized GLCM
#' @export
glcm_metrics<-function(GLCM){
  return(C_glcm_metrics(GLCM=GLCM))
  }