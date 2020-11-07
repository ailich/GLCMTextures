#' Calculates GLCM texture metrics of a Raster Layer
#'
#' Calculates GLCM texture metrics of a RasterLayer over a sliding rectangular window
#' @param r A raster layer. If already quantized set quantization to "none". The valid range of values for a quantized raster is from 0 to n_levels-1 (e.g. a raster with 32 grey levels would have a valid range of 0-31).
#' @param w A vector of length 2 specifying the dimensions of the rectangular window to use where the first number is the number of rows and the second number is the number of columns. Window size must be an odd number.
#' @param n_levels Number of grey levels used in the quantization
#' @param shift A vector of length 2, or a list of vectors each of length 2 specifying the relationship between neighboring pixel to the reference pixel. The first number represents the shift in the x direction and the second number represents the shift in the y direction, where up and right are positive. For example c(1,0) is the pixel directly to the right. The GLCM is made symmetrical by counting each pair twice, once "forwards" and once "backwards" by interchanging reference and neighbor pixels. Therefore a shift directly to the right c(1,0) is equivalent to a shift directly to the left c(-1,0). To average over "all directions" you can use shift=list(c(1,0), c(1,1), c(0,1), c(-1,1)), which is the default.
#' @param metrics A vector of glcm texture metrics to return. Valid entries include "glcm_contrast", "glcm_dissimilarity", "glcm_homogeneity", "glcm_ASM", "glcm_entropy", "glcm_mean", "glcm_variance", "glcm_correlation".
#' @param quantization quantization method (either "equal range", "equal prob", or "none"). "equal range" quantization will create bins that cover a range of equal size. "equal prob" performs equal probability quantization and will use quantiles to create bins with approximately equal number of samples. "none" means the layer has already been quantized.
#' @param min_val minimum value for equal range quantization (if not supplied, the minimum value of the raster is used)
#' @param max_val maximum value for equal range quantization (if not supplied, the maximum value of the raster is used)
#' @param na_opt A character vector indicating how to consider NA values. "any" means that NA will be returned if any values in the window are NA. "center" means that NA will be returned only if the central pixel in the window is NA. "all" means that NA will be returned only if all values in the window are NA. "any" is the default.
#' @param pad logical value specifying whether rows/columns of NA's should be padded to the edge of the raster to remove edge effects (FALSE by default)
#'
#' @return a RasterBrick of Texture Metrics (or RasterLayer if just one metric is calculated)
#' @import raster
#' @export
#'
glcm_textures<- function(r, w, n_levels, shift=list(c(1,0), c(1,1), c(0,1), c(-1,1)), metrics= c("glcm_contrast", "glcm_dissimilarity", "glcm_homogeneity", "glcm_ASM", "glcm_entropy", "glcm_mean", "glcm_variance", "glcm_correlation"), quantization, min_val=NULL, max_val=NULL, na_opt= "any", pad=FALSE){
  if(length(w==1)){
    w<- rep(w,2)}
  if(pad==TRUE){
    og_extent<- raster::extent(r)
    r<- raster::extend(r, c(w[1], w[2]), value=NA)
    }
  if(quantization!="none"){
    r<- quantize_raster(r = r, n_levels = n_levels, method = quantization, min_val = min_val, max_val = max_val)
  }
  if((raster::cellStats(r, stat = max) > (n_levels-1)) | (raster::cellStats(r, stat = min) < 0)){
    stop("Error: raster layer must have values between 0 and n_levels-1")}
  run_in_blocks<- !raster::canProcessInMemory(r, n = 9)
  if(run_in_blocks==FALSE){
    output<- glcm_textures_helper(rq=r, w=w, n_levels=n_levels, shift=shift, metrics=metrics, na_opt=na_opt)
  } else{
    block_idx<- raster::blockSize(r, n = 9, minblocks = 2, minrows = w[1])
    out_blocks<- vector(mode = "list", length = block_idx$n)
    block_overlap<- w[1]-1
    for (i in 1:block_idx$n) {
      min_row<- block_idx$row[[i]]
      max_row<- min(min_row + block_idx$nrows[[i]] - 1 + block_overlap, nrow(r))
      block_extent<- raster::extent(r, min_row, max_row, 1, ncol(r))
      curr_block<- raster::crop(r, block_extent)
      out_blocks[[i]]<- glcm_textures_helper(rq=curr_block, w=w, n_levels=n_levels, shift=shift, metrics=metrics)
    }
    output<- do.call(raster::merge, out_blocks)
     names(output)<- metrics
  }
  if(pad==TRUE){
    output<- raster::crop(output, og_extent)
  }
  return(output)
}
