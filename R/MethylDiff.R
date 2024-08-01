MethylDiff <- function(chromDictMeth, chrom, start, end){
  return(mean(na.omit(chromDictMeth[[chrom]][.(start:end)])$MethylDiff))
}