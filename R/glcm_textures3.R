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
#' @param pad logical value specifying whether rows/columns of NA's should be padded to the edge of the raster to remove edge effects (FALSE by default). If pad is TRUE, na_opt musy be set to "center" or "all".
#'
#' @return a RasterBrick of Texture Metrics (or RasterLayer if just one metric is calculated)
#' @import raster
#' @export
#'
glcm_textures3<- function(r, w, n_levels, shift=list(c(1,0), c(1,1), c(0,1), c(-1,1)), metrics= c("glcm_contrast", "glcm_dissimilarity", "glcm_homogeneity", "glcm_ASM", "glcm_entropy", "glcm_mean", "glcm_variance", "glcm_correlation"), quantization, min_val=NULL, max_val=NULL, na_opt= "any", pad=FALSE){
  if(length(w==1)){
    w<- rep(w,2)}
  if(any(w<3) | any(0 == (w %% 2))){
    stop("Error: w must be odd and greater than or equal to 3")}
  if (!any(na_opt==c("any", "center", "all"))){
    stop("na_opt must be 'any', 'center', or 'all'")
  }
  if(class(shift)!="list"){shift=list(shift)}
  #Maybe put check for shift size < windowsize-1/2
  out_list<- vector(mode = "list", length=length(shift))
  if(pad==TRUE){
    if(na_opt=="any"){
      stop("if pad=TRUE, na_opt must be 'center' or 'all'")}
    og_extent<- raster::extent(r)
    r<- raster::extend(r, c((w[1]-1)/2, (w[2]-1)/2), value=NA)
    }
  if(quantization!="none"){
    r<- quantize_raster(r = r, n_levels = n_levels, method = quantization, min_val = min_val, max_val = max_val)
  }
  if((raster::cellStats(r, stat = max) > (n_levels-1)) | (raster::cellStats(r, stat = min) < 0)){
    stop("Error: raster layer must have values between 0 and n_levels-1")}

  out_list<- vector(mode = "list", length=length(shift))

  run_in_blocks<- !raster::canProcessInMemory(r, n = (8*length(shift))+1)
  if(run_in_blocks==FALSE){
    for (k in 1:length(shift)) {
      out_list[[k]]<- raster::brick(r, nl=8, values=FALSE)
      curr_vals<- C_glcm_textures_helper2(rq= as.matrix(r), w=w, n_levels=n_levels, shift=shift[[k]], na_opt=na_opt)
      values(out_list[[k]])<- curr_vals
      names(out_list[[k]])<- colnames(curr_vals)
      out_list[[k]]<- raster::subset(out_list[[k]], metrics, drop=TRUE)
    }
    n_layers<- raster::nlayers(out_list[[1]])
    avg_shifts<- length(shift) > 1
    if(avg_shifts){
      output<- stack()
      for (j in 1:n_layers) {
        out_layer<- mean(do.call(raster::stack, lapply(out_list, raster::subset,j))) #Average across all shifts
        output<- raster::stack(output, out_layer) #Create new stack of averaged values
      }} else{
        output<- out_list[[1]]
      }
    names(output)<- metrics
    } else{
      block_overlap<- (w[1]-1)/2
      nr<- nrow(r)
      nc<- ncol(r)
      f_out <- raster::rasterTmpFile()
      output<- raster::brick(r, nl=8, values=FALSE)
      output<- raster::writeStart(output, filename = f_out)
      block_idx<- raster::blockSize(r, n = 8, minblocks = 2, minrows = w[1])

      for (i in 1:block_idx$n) {
        print(paste(i, "of", block_idx$n))
        min_row<- max(c(block_idx$row[[i]] - block_overlap), 1)
        max_row<- min(c(block_idx$row[[i]] + block_idx$nrows[[i]] - 1 + block_overlap, nr))
        curr_block <- raster::getValues(r, row = min_row, nrows = max_row-min_row+1, format="matrix")
        out_array<- array(dim = c(length(curr_block), 8, length(shift)))
        for (k in 1:length(shift)) {
          print(paste("Starting Shift", k))
          out_array[,,k]<- C_glcm_textures_helper2(rq= curr_block, w=w, n_levels=n_levels, shift=shift[[k]], na_opt=na_opt)
        }

        out_block<- apply(out_array,MARGIN = c(1,2), FUN = mean, na.rm=TRUE) #average across shifts

        #out_block is a formatted as matrix where each column corresponds to a raster layer (this is how writeRaster needs it to be formatted)
        #As you go down rows in this matrix you move across rows in the raster object
        if(i==1){
          out_block<- out_block[1:(nrow(out_block)-(block_overlap*nc)),] #Trim bottom edge of raster
        } else if (i != block_idx$n){
          out_block<- out_block[(1+block_overlap*nc):(nrow(out_block)-(block_overlap*nc)),] #Trim top and bottom edge of raster
        } else {
          out_block<- out_block[(1+block_overlap*nc):nrow(out_block),] #Trim top edge of raster
        }
        raster::writeValues(output, v= out_block, start= block_idx$row[i])
        }
        output<- raster::writeStop(output)
        names(output)<- c("glcm_contrast", "glcm_dissimilarity", "glcm_homogeneity","glcm_ASM","glcm_entropy","glcm_mean","glcm_variance","glcm_correlation")
        output<- raster::subset(output, metrics, drop=TRUE)
        }
  if(n_layers>1){
    output<- raster::brick(output)}
  if(pad==TRUE){
    output<- raster::crop(output, og_extent)
    }
  return(output)
}
