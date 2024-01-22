#' Creates a symmetrical normalized GLCM for a given matrix and shift
#'
#' @param x a matrix of integers representing quantized values. The valid range of values is from 0 to n_levels-1 (e.g. a matrix with 32 grey levels would have a valid range of 0-31).
#' @param n_levels Number of grey levels used in the quantization
#' @param shift A vector of length 2 specifying the relationship between neighboring pixel to the reference pixel. The first number represents the shift in the x direction and the second number represents the shift in the y direction, where up and right are positive. For example c(1,0) is the pixel directly to the right. The GLCM is made symmetrical by counting each pair twice, once "forwards" and once "backwards" by interchanging reference and neighbor pixels. Therefore a shift directly to the right c(1,0) is equivalent to a shift directly to the left c(-1,0)
#' @param na.rm a logical value indicating whether NA values should be stripped before the computation proceeds (default=FALSE)
#' @param normalize a logical specifying whether to normalize the counts to probabilities by dividing by the sum of the GLCM (TRUE, the default) or to express the GLCM as counts (FALSE)
#' @return A symmetric GLCM as a matrix
#' @examples
#' test_matrix<- matrix(data=c(2,0,1,3,0,0,0,3,2), nrow = 3, ncol=3)
#' # Tabulate a GLCM of counts
#' horizontal_glcm_counts<- make_glcm(test_matrix, n_levels = 4, shift = c(1,0), normalize = FALSE)
#' # Calculate a normalized GLCM of probabilities
#' horizontal_glcm_norm<- make_glcm(test_matrix, n_levels = 4, shift = c(1,0), normalize = TRUE)
#' @references
#' Hall-Beyer, M., 2017. GLCM Texture: A Tutorial v. 3.0. University of Calgary, Alberta, Canada.
#'
#' Haralick, R.M., Shanmugam, K., Dinstein, I., 1973. Textural features for image classification. IEEE Transactions on Systems, Man, and Cybernetics 610–621. https://doi.org/10.1109/TSMC.1973.4309314
#' @export

make_glcm<- function(x, n_levels, shift, na.rm = FALSE, normalize=TRUE){
	if(isTRUE(any(x > (n_levels-1))) | isTRUE(any(x < 0))){
	  stop("Error: x must have values between 0 and n_levels-1")
	  }
  if(normalize){
	  return(C_make_glcm(x=x, n_levels=n_levels, shift=shift, na_rm=na.rm))
	  } else{
	    return(C_make_glcm_counts(x=x, n_levels=n_levels, shift=shift, na_rm=na.rm))
	  }
  }
