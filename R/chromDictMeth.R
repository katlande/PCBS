chromDictMeth <- function(mat, IDs, filter_thresh=50){
  
  message("Removing low depth sites...")
  mat$totalCounts <- rowSums(mat[which(grepl("nCpG", colnames(mat)))])
  mat <- mat[mat$totalCounts > filter_thresh,]
  
  message("Calculating mean methylation differences...")
  mat[which(grepl("PercMeth", colnames(mat)))]
  mat <- cbind(mat[1],mat[which(grepl("PercMeth", colnames(mat)))]) 
  mat$MethylDiff <- rowMeans(mat[which(grepl(IDs[1], colnames(mat)))])-rowMeans(mat[which(grepl(IDs[2], colnames(mat)))])
  colnames(mat)[1] <- "cpgID"
  row.names(mat) <- mat$cpgID
  mat$chr <- gsub(":.*", "", row.names(mat))
  mat$pos <- gsub(".*:", "", row.names(mat))
  methdiff <- mat[, c("chr", "pos", "MethylDiff")]
  
  
  outlist <- list()
  message("Splitting by chromosome...")
  for(i in unique(methdiff$chr)){
    d <- methdiff[methdiff$chr==i,]
    d$pos <- as.numeric(d$pos)
    d <- d[order(d$pos),]
    data.table::setDT(d)
    data.table::setkey(d,pos)
    outlist <- append(outlist, list(d))
  }
  
  names(outlist) <- unique(methdiff$chr)
  return(outlist)
}