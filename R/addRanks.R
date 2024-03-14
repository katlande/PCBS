#' @export
addRanks <- function(ranks){
  if(!tibble::has_rownames(ranks)){
    row.names(ranks) <- paste0(ranks$chr, ":", as.character(ranks$pos))
  }
  ranks[order(ranks$PC_Score, decreasing = T), , drop=F]->ranks
  ranks$order <- 1:nrow(ranks)
  ranks[order(abs(ranks$PC_Score), decreasing = T),]->ranks
  ranks$abs.order <- 1:nrow(ranks)
  ranks$chr <- gsub("\\:.*", "", row.names(ranks))
  ranks$pos <- as.numeric(gsub(".*\\:", "", row.names(ranks)))
  return(ranks)
}
