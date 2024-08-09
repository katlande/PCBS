#' @export
getRegionScores <- function(ranks=NULL, regions, chromDictObj=NULL){
  
  if(is.null(chromDictObj) & is.null(ranks)){
    warning("Both ranks and chromDictObj are NULL. At least one must be used.")
  } else{
    if(is.null(chromDictObj)){
      message("Creating chromDict Object. If you plan to run the getRegionScores() function multiple times, it is strongly recommended to generate a chromDict Object with the chromDict() function, and specify it with getRegionScores(..., chromDictObj=OBJECT). This is the most computationally intensive part of most PCBS functions, and only needs to be done once per dataset.")
      chromDict(ranks) -> chromDictObj
    }
    
    # background distribution for each chromosome:
    message("Getting Background Sites...")
    local_list <- list()
    for(i in names(chromDictObj)){
      d <- chromDictObj[[i]]
      v <- ifelse(nrow(d)*0.01 > 10000, 10000, ceiling(nrow(d)*0.01))
      v <- ifelse(v < 100, 100, v)
      n <- d$PC_Score[sample(1:nrow(d),v)]
      local_list <- append(local_list, list(n))
    }
    names(local_list) <- names(chromDictObj)
    
    message(paste("Checking significance at", nrow(regions), "sites..."))
    means <- c()
    nCpGs <- c()
    zv <- c()
    pv <- c()
    for(i in 1:nrow(regions)){
      df <- chromDictObj[[regions[[1]][i]]]
      df.main <- na.omit(df[ .( c((regions[[2]][i]:regions[[3]][i])) ) ])
      
      u <- mean(df.main$PC_Score)
      n <- nrow(df.main)
      means <- c(means,u)
      nCpGs <- c(nCpGs, n)
      zv <- c(zv,((u-mean(local_list[[regions[[1]][i]]]))/sd(local_list[[regions[[1]][i]]])))
      
      if(length(local_list[[regions[[1]][i]]]) > 3 & length(df.main$PC_Score) > 3){
        pv <- c(pv, t.test(df.main$PC_Score, local_list[[regions[[1]][i]]])$p.value)
      } else{
        pv <- c(pv, NA)
      }
      }
     if(ncol(regions) > 3){
      return(data.frame(chrom=regions[[1]], start=regions[[2]], end=regions[[3]], 
                        feature=regions[[4]], meanPC=means, nCpG=nCpGs, Z=zv, p=pv))
    } else {
      return(data.frame(chrom=regions[[1]], start=regions[[2]], end=regions[[3]], 
                        meanPC=means, nCpG=nCpGs, Z=zv, p=pv))
    }
  }
}
