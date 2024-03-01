#' @export
chromDict <- function(ranks){
  ranks$chr <- gsub(":.*", "", row.names(ranks))
  ranks$pos <- gsub(".*:", "", row.names(ranks))
  outlist <- list()
  for(i in unique(ranks$chr)){
    d <- ranks[ranks$chr==i,]
    d$pos <- as.numeric(d$pos)
    d <- d[order(d$pos),]
    data.table::setDT(d)
    data.table::setkey(d,pos)
    outlist <- append(outlist, list(d))
  }

  names(outlist) <- unique(ranks$chr)
  return(outlist)
}
