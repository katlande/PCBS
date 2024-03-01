#' @export
getPCRanks <- function(mat, IDs, PC=1, filter_thresh=50){
  mat$totalCounts <- rowSums(mat[which(grepl("nCpG", colnames(mat)))])
  mat <- mat[mat$totalCounts > filter_thresh,]
  row.names(mat) <- c()
  tibble::column_to_rownames(mat, "cpgID")->mat
  pca <- prcomp(t(mat[which(grepl("PercMeth", colnames(mat)))]))
  rotations <- data.frame(abs(pca$rotation[,paste0("PC",PC)]))
  colnames(rotations) <- "PC_Score"

  # get directions:
  mat[which(grepl("PercMeth", colnames(mat)))]->mat
  mat$mean_1 <- rowMeans(mat[which(grepl(IDs[1], colnames(mat)))])
  mat$mean_2 <- rowMeans(mat[which(grepl(IDs[2], colnames(mat)))])
  mat$direction <- sign(mat$mean_1-mat$mean_2)

  rotations <- cbind(rotations,mat[c(ncol(mat))])
  rotations$PC_Score <- rotations$PC_Score*rotations$direction
  return(rotations[c(1, drop=F)])
}
