#' Creates a symmetrical normalized GLCM for a given matrix and shift
#' 
#' @param x a matrix of integers representing quantized values. The valid range of values is from 0 to n_levels-1 (e.g. a matrix with 32 grey levels would have a valid range of 0-31).
#' @param n_levels Number of grey levels used in the quantization
#' @param shift A vector of length 2 specifying the relationship between neighboring pixel to the reference pixel. The first number represents the shift in the x direction and the second number represents the shift in the y direction, where up and right are positive. For example c(1,0) is the pixel directly to the right. The GLCM is made symmetrical by counting each pair twice, once "forwards" and once "backwards" by interchanging reference and neighbor pixels. Therefore a shift directly to the right c(1,0) is equivalent to a shift directly to the left c(-1,0)
#' @export

make_glcm<- function(x,n_levels, shift){
  return(C_make_glcm(x=x, n_levels=n_levels, shift=shift))
  }