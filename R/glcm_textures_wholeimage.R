#' Calculates GLCM texture metrics of a Raster Layer
#'
#' Calculates GLCM texture metrics of a RasterLayer over a sliding rectangular window
#' @param r A single layer SpatRaster or RasterLayer. If already quantized set quant_method to "none". The valid range of values for a quantized raster is from 0 to n_levels-1 (e.g. a raster with 32 grey levels would have a valid range of 0-31).
#' @param n_levels Number of grey levels used in the quantization (Typically set to 16 or 32).
#' @param shift A vector of length 2, or a list of vectors each of length 2 specifying the relationship between neighboring pixel to the reference pixel. The first number represents the shift in the x direction and the second number represents the shift in the y direction, where up and right are positive. For example c(1,0) is the pixel directly to the right. The GLCM is made symmetrical by counting each pair twice, once "forwards" and once "backwards" by interchanging reference and neighbor pixels. Therefore a shift directly to the right c(1,0) is equivalent to a shift directly to the left c(-1,0). To average over "all directions" you can use shift=list(c(1,0), c(1,1), c(0,1), c(-1,1)), which is the default.
#' @param metrics A vector of glcm texture metrics to return. Valid entries include "glcm_contrast", "glcm_dissimilarity", "glcm_homogeneity", "glcm_ASM" (angular second moment), "glcm_entropy", "glcm_mean", "glcm_variance", "glcm_correlation".
#' @param quant_method quantization method (either "range", "prob", or "none"). "range" quantization will create bins that cover a range of equal size. "prob" performs equal probability quantization and will use quantiles to create bins with approximately equal number of samples. "none" means the layer has already been quantized.
#' @param min_val minimum value for equal range quantization (if not supplied, the minimum value of the raster is used)
#' @param max_val maximum value for equal range quantization (if not supplied, the maximum value of the raster is used)
#' @param maxcell positive integer used to take a regular sample for quantization if "prob" is used as quant_method (default is Inf)
#' @param na.rm a logical value indicating whether NA values should be stripped before the computation proceeds (default=FALSE)
#' @param impute_corr logical indicating whether glcm correlation should be filled with zero in the case where all values are the same. Strictly glcm correlation is NA in this case but the limit approaches zero.
#' @param wopt list with named options for writing files as in writeRaster
#' @return a vector of texture metrics
#' @import terra
#' @references
#' Hall-Beyer, M., 2017. GLCM Texture: A Tutorial v. 3.0. University of Calgary, Alberta, Canada.
#'
#' Haralick, R.M., Shanmugam, K., Dinstein, I., 1973. Textural features for image classification. IEEE Transactions on Systems, Man, and Cybernetics 610–621. https://doi.org/10.1109/TSMC.1973.4309314

glcm_textures_wholeimage<- function(r, n_levels, shift, metrics, quant_method, min_val, max_val, maxcell, na.rm, impute_corr, wopt){
  if(quant_method!="none"){
    r<- quantize_raster(r = r, n_levels = n_levels, quant_method = quant_method, min_val = min_val, max_val = max_val, maxcell=maxcell, wopt=wopt)
  }
  GLCM<- make_glcm(r, n_levels = n_levels, shift = shift, na.rm = na.rm, normalize = TRUE)
  out<- glcm_metrics(GLCM = GLCM, metrics = metrics, average = TRUE, impute_corr = impute_corr)
  return(out)
}

