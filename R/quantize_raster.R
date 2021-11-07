#' Quantizes raster to a set number of discrete levels
#'
#' Quantizes raster to a set number of discrete levels starting at 0. There are 2 methods of quantization are available: "uniform" and "equal prob"
#' @param r A raster layer
#' @param n_levels number of levels to quantize to
#' @param method quantization method (either "equal range" or "equal prob"). "equal range" quantization will create bins that cover a range of equal size. "equal prob" performs equal probability quantization and will use quantiles to create bins with approximately equal number of samples.
#' @param min_val minimum value for equal range quantization (if not supplied, the minimum value of the raster is used)
#' @param max_val maximum value for equal range quantization (if not supplied, the maximum value of the raster is used)
#' @importFrom  raster cellStats
#' @importFrom  raster quantile
#' @importFrom  raster cut
#' @export
#'
quantize_raster<- function(r, n_levels, method, min_val=NULL, max_val=NULL){
  if(method=="equal range"){
    if(is.null(min_val)){min_val<- raster::cellStats(x = r, stat="min")}
    if(is.null(max_val)){max_val<- raster::cellStats(x = r, stat="max")}
    qrules<- seq(min_val, max_val, length.out = (n_levels + 1))
  } else if (method == "equal prob") {
    qrules<- quantile(r, probs= seq(0,1,length.out = n_levels+1), na.rm=TRUE, type=8)
  } else {
    message("Error: method must be 'equal range' or 'equal prob'")
    }
  rq <- raster::cut(r, breaks = qrules, include.lowest = TRUE, right = FALSE)-1
  return(rq)
}
