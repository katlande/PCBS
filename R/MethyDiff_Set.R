MethyDiff_Set <- function(chromDictMeth, regions){
  
  apply(regions, 1, function(x){
    MethylDiff(chromDictMeth, x[[1]], as.numeric(x[[2]]), as.numeric(x[[3]]))
  }) -> regions$mean_methylDiff
  
  return(regions)
}