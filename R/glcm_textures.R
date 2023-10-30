#' Calculates GLCM texture metrics of a Raster Layer
#'
#' Calculates GLCM texture metrics of a RasterLayer over a sliding rectangular window
#' @param r A single layer SpatRaster or RasterLayer. If already quantized set quantization to "none". The valid range of values for a quantized raster is from 0 to n_levels-1 (e.g. a raster with 32 grey levels would have a valid range of 0-31).
#' @param w A vector of length 2 specifying the dimensions of the rectangular window to use where the first number is the number of rows and the second number is the number of columns. Window size must be an odd number.
#' @param n_levels Number of grey levels used in the quantization (Typically set to 16 or 32).
#' @param shift A vector of length 2, or a list of vectors each of length 2 specifying the relationship between neighboring pixel to the reference pixel. The first number represents the shift in the x direction and the second number represents the shift in the y direction, where up and right are positive. For example c(1,0) is the pixel directly to the right. The GLCM is made symmetrical by counting each pair twice, once "forwards" and once "backwards" by interchanging reference and neighbor pixels. Therefore a shift directly to the right c(1,0) is equivalent to a shift directly to the left c(-1,0). To average over "all directions" you can use shift=list(c(1,0), c(1,1), c(0,1), c(-1,1)), which is the default.
#' @param metrics A vector of glcm texture metrics to return. Valid entries include "glcm_contrast", "glcm_dissimilarity", "glcm_homogeneity", "glcm_ASM", "glcm_entropy", "glcm_mean", "glcm_variance", "glcm_correlation".
#' @param quantization quantization method (either "equal range", "equal prob", or "none"). "equal range" quantization will create bins that cover a range of equal size. "equal prob" performs equal probability quantization and will use quantiles to create bins with approximately equal number of samples. "none" means the layer has already been quantized.
#' @param min_val minimum value for equal range quantization (if not supplied, the minimum value of the raster is used)
#' @param max_val maximum value for equal range quantization (if not supplied, the maximum value of the raster is used)
#' @param na.rm a logical value indicating whether NA values should be stripped before the computation proceeds (default=FALSE)
#' @param include_scale Logical indicating whether to append window size to the layer names (default = FALSE).
#' @param filename character Output filename. Can be a single filename, or as many filenames as there are layers to write a file for each layer
#' @param overwrite logical. If TRUE, filename is overwritten (default is FALSE).
#' @param wopt list with named options for writing files as in writeRaster
#' @return a SpatRaster or Raster* Object
#' @examples
#' r<- rast(volcano, extent= ext(2667400, 2667400 + ncol(volcano)*10,
#' 6478700, 6478700 + nrow(volcano)*10), crs = "EPSG:27200")
#' txt <- glcm_textures(r, w = c(3,5), n_levels = 16,
#' quantization = "equal prob", shift = list(c(1, 0), c(1, 1),
#' c(0, 1), c(-1, 1)))
#' plot(txt)
#' @import terra
#' @importFrom raster raster
#' @importFrom raster stack
#' @importFrom raster writeRaster

#' @references
#' Hall-Beyer, M., 2017. GLCM Texture: A Tutorial v. 3.0. University of Calgary, Alberta, Canada.
#'
#' Haralick, R.M., Shanmugam, K., Dinstein, I., 1973. Textural features for image classification. IEEE Transactions on Systems, Man, and Cybernetics 610–621. https://doi.org/10.1109/TSMC.1973.4309314
#' @export
#'
glcm_textures<- function(r, w = c(3,3), n_levels, shift=list(c(1,0), c(1,1), c(0,1), c(-1,1)), metrics= c("glcm_contrast", "glcm_dissimilarity", "glcm_homogeneity", "glcm_ASM", "glcm_entropy", "glcm_mean", "glcm_variance", "glcm_correlation"), quantization, min_val=NULL, max_val=NULL, na.rm=FALSE, include_scale=FALSE, filename=NULL, overwrite=FALSE, wopt=list()){
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
  if(!is.list(shift)){shift=list(shift)}
  if(any(sapply(shift, length)!=2)){
    stop("Error: each shift must be a vector of length 2")
  }
  if(!all((sapply(shift, class)=="numeric") | (sapply(shift, class)=="integer"))){
    stop("Error: shifts must be a numeric or integer")
  }
  if (any(!(metrics %in% all_metrics))){
    stop("Error: Invlaid metric. Valid metrics include 'glcm_contrast', 'glcm_dissimilarity', 'glcm_homogeneity', 'glcm_ASM', 'glcm_entropy', 'glcm_mean', 'glcm_variance', 'glcm_correlation'")
  }

  if(n_levels > 32){
    warning("n_levels is > 32. This may be very computationally expensive.")
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
  } #Some metrics needed to calculate others

  out_list<- vector(mode = "list", length=length(shift))
  if(quantization!="none"){
    r<- quantize_raster(r = r, n_levels = n_levels, method = quantization, min_val = min_val, max_val = max_val, wopt=wopt)
    } else if(!terra::is.int(r)){
    r<- terra::as.int(r, wopt=wopt) #Make it an integer raster
    }
  if((unlist(terra::global(r, fun = max, na.rm=TRUE)) > (n_levels-1)) | (unlist(terra::global(r, fun = min, na.rm=TRUE)) < 0)){
    stop("Error: raster must have values between 0 and n_levels-1")}

  out_list<- vector(mode = "list", length=length(shift))
  for (k in 1:length(shift)) {
    out_list[[k]]<- terra::focalCpp(r, w=w, fun = C_glcm_textures_helper, w2=w, n_levels= n_levels, shift = shift[[k]], metrics = needed_metrics, na_rm=na.rm, fillvalue=NA, wopt=wopt)
    out_list[[k]]<- terra::subset(out_list[[k]], metrics, wopt=wopt) #Subset from needed to requested metrics
    }

  n_layers<- length(metrics)
  if(length(shift) > 1){
    output<- terra::app(terra::sds(out_list), fun=mean, na.rm=na.rm, wopt=wopt)
    } else{
      output<- out_list[[1]]
    }
  names(output)<- metrics #preserve names in case they were lost

  if(include_scale){names(output)<- paste0(names(output), "_", w[1],"x", w[2])} #Add scale to layer names

  if(og_class=="RasterLayer"){
    if(terra::nlyr(output) > 1){
      output<- raster::stack(output) #Convert to RasterStack
      if(!is.null(filename)){
        if(length(filename)==1){
          return(raster::writeRaster(output, filename=filename, overwrite=overwrite, bylayer=FALSE))
        } else{
          return(raster::writeRaster(output, filename=filename, overwrite=overwrite, bylayer=TRUE))
        }
        }
      } else{
        output<- raster::raster(output)
        if(!is.null(filename)){
          return(raster::writeRaster(output, filename=filename, overwrite=overwrite))
        }
      }
  }
  if(!is.null(filename)){
    return(terra::writeRaster(output, filename=filename, overwrite=overwrite, wopt=wopt))
    }
  return(output)
}
