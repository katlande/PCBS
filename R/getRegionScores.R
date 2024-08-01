#' @export
getRegionScores <- function(ranks, regions){

  if(! "abs.order" %in% colnames(ranks)){
    message("No ranks in input file; adding ranks. Save computational time by running addRanks() on your input file if you plan on running getRegionScores() more than once.")
    addRanks(ranks)->ranks
  }

  unlist(lapply(1:1000, function(x){
    mean(ranks$PC_Score[sample(1:nrow(ranks), 100)])
  })) -> n
  upop <- mean(n)
  sdpop <- sd(n)

  apply(regions, 1, FUN=function(x){
    paste0(x[1],":",as.numeric(x[2]):as.numeric(x[3]))
  }) -> r

  ranks[row.names(ranks)%in%unique(unlist(r)),]->tmp
  names(r)<-regions[[4]]

  means <- c()
  nCpGs <- c()
  zv <- c()
  pv <- c()
  for(i in 1:length(r)){
    a<-tmp$PC_Score[row.names(tmp)%in%r[[i]]]
    means <- c(means,mean(a))
    nCpGs <- c(nCpGs, length(a))

    zv <- c(zv,((mean(a)-upop)/sdpop))
    if(length(a) > 3){
      if(length(a) < length(n)){
        b <- n[sample(1:length(n), length(a))]
      } else{
        b <- n
      }
      pv <- c(pv, t.test(a, b)$p.value)
    } else{
      pv <- c(pv, NA)
    }

  }

  return(data.frame(feature=names(r), meanPC=means, nCpG=nCpGs, Z=zv, p=pv))
}
