#' Calculates GLCM texture metrics of a Raster Layer
#'
#' Calculates GLCM texture metrics of a RasterLayer over a sliding rectangular window
#' @param r A single layer SpatRaster or RasterLayer. If already quantized set quantization to "none". The valid range of values for a quantized raster is from 0 to n_levels-1 (e.g. a raster with 32 grey levels would have a valid range of 0-31).
#' @param w A vector of length 2 specifying the dimensions of the rectangular window to use where the first number is the number of rows and the second number is the number of columns. Window size must be an odd number.
#' @param n_levels Number of grey levels used in the quantization
#' @param shift A vector of length 2, or a list of vectors each of length 2 specifying the relationship between neighboring pixel to the reference pixel. The first number represents the shift in the x direction and the second number represents the shift in the y direction, where up and right are positive. For example c(1,0) is the pixel directly to the right. The GLCM is made symmetrical by counting each pair twice, once "forwards" and once "backwards" by interchanging reference and neighbor pixels. Therefore a shift directly to the right c(1,0) is equivalent to a shift directly to the left c(-1,0). To average over "all directions" you can use shift=list(c(1,0), c(1,1), c(0,1), c(-1,1)), which is the default.
#' @param metrics A vector of glcm texture metrics to return. Valid entries include "glcm_contrast", "glcm_dissimilarity", "glcm_homogeneity", "glcm_ASM", "glcm_entropy", "glcm_mean", "glcm_variance", "glcm_correlation".
#' @param quantization quantization method (either "equal range", "equal prob", or "none"). "equal range" quantization will create bins that cover a range of equal size. "equal prob" performs equal probability quantization and will use quantiles to create bins with approximately equal number of samples. "none" means the layer has already been quantized.
#' @param min_val minimum value for equal range quantization (if not supplied, the minimum value of the raster is used)
#' @param max_val maximum value for equal range quantization (if not supplied, the maximum value of the raster is used)
#' @param na_opt A character vector indicating how to consider NA values. "any" means that NA will be returned if any values in the window are NA. "center" means that NA will be returned only if the central pixel in the window is NA. "all" means that NA will be returned only if all values in the window are NA. "any" is the default.
#' @return a SpatRaster or Raster* Object
#' @import terra
#' @importFrom raster raster
#' @importFrom raster stack
#' @references
#' Hall-Beyer, M., 2017. GLCM Textrure: A Tutorial v. 3.0. University of Calgary, Alberta, Canada.
#'
#' Haralick, R.M., Shanmugam, K., Dinstein, I., 1973. Textural features for image classification. IEEE Transactions on Systems, Man, and Cybernetics 610–621. https://doi.org/10.1109/TSMC.1973.4309314
#' @export
#'
glcm_textures2<- function(r, w = c(3,3), n_levels, shift=list(c(1,0), c(1,1), c(0,1), c(-1,1)), metrics= c("glcm_contrast", "glcm_dissimilarity", "glcm_homogeneity", "glcm_ASM", "glcm_entropy", "glcm_mean", "glcm_variance", "glcm_correlation"), quantization, min_val=NULL, max_val=NULL, na_opt= "any"){
  og_class<- class(r)[1]
  if(og_class=="RasterLayer"){
    r<- terra::rast(r) #Convert to SpatRaster
  }
  all_metrics<- c("glcm_contrast", "glcm_dissimilarity", "glcm_homogeneity", "glcm_ASM", "glcm_entropy", "glcm_mean", "glcm_variance", "glcm_correlation")
  # Input checks
  if(!(og_class %in% c("RasterLayer", "SpatRaster"))){
    stop("Error: Input must be a 'SpatRaster' or 'RasterLayer'")
  }
  if(terra::nlyr(r)!=1){
    stop("Error: Input raster must be one layer.")
  }
  if(length(w)==1){
    w<- rep(w,2)}
  if(length(w)>2){
    stop("Specified window exceeds 2 dimensions")
    }
  if(any(0 == (w %% 2))){
    stop("Error: w must be odd")}
  if(all(w<3)){
    stop("Error: w must be greater or equal to 3 in at least one dimension")
  }
  if (!any(na_opt==c("any", "center", "all"))){
    stop("na_opt must be 'any', 'center', or 'all'")
  }
  if(class(shift)!="list"){shift=list(shift)}
  if(any(sapply(shift, length)!=2)){
    stop("Error: each shift must be a vector of length 2")
  }
  if(!all((sapply(shift, class)=="numeric") | (sapply(shift, class)=="integer"))){
    stop("Error: shifts must be a numeric or integer")
  }
  if (any(!(metrics %in% all_metrics))){
    stop("Error: Invlaid metric. Valid metrics include 'glcm_contrast', 'glcm_dissimilarity', 'glcm_homogeneity', 'glcm_ASM', 'glcm_entropy', 'glcm_mean', 'glcm_variance', 'glcm_correlation'")
  }

  out_list<- vector(mode = "list", length=length(shift))
  if(quantization!="none"){
    r<- quantize_raster2(r = r, n_levels = n_levels, method = quantization, min_val = min_val, max_val = max_val)
  } else if(!terra::is.int(r)){
    r<- terra::as.int(r) #Make it an integer raster
    }
  if((terra::global(r, fun = max) > (n_levels-1)) | (terra::global(r, fun = min) < 0)){
    stop("Error: raster must have values between 0 and n_levels-1")}

  out_list<- vector(mode = "list", length=length(shift))
  for (k in 1:length(shift)) {
    out_list[[k]]<- terra::focalCpp(r, w=w, fun = C_glcm_textures_helper2, w2=w, n_levels= n_levels, shift = shift[[k]], na_opt=na_opt)
    out_list[[k]]<- terra::subset(out_list[[k]], metrics)
    }

  n_layers<- length(metrics)
  avg_shifts<- length(shift) > 1
  if(avg_shifts){
    output<- terra::rast()
    for (j in 1:n_layers) {
      out_layer<- mean(do.call(c, lapply(out_list, terra::subset,j))) #Average across all shifts
      output<- c(output, out_layer, warn=FALSE) #Create new stack of averaged values
    }} else{
      output<- out_list[[1]]
    }
  names(output)<- metrics #preserve names in case they were lost

  if(og_class=="RasterLayer"){
    if(terra::nlyr(output) > 1){
      output<- raster::stack(output) #Convert to RasterStack
    } else{
      output<- raster::raster(output)
    }
  }
  return(output)
}
