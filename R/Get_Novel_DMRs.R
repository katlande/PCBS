#' @export
Get_Novel_DMRs <- function(ranks, nSeeds, chromDictObj=NULL, DMR_resolution=NULL, QueryLimit=5000, minCpGs=15, minZ=1,
                           perms=1000){

  if(is.null(chromDictObj)){
    message("Splitting data by chromosome...")
    chromDict(ranks) -> chromDictObj
  }

  # getting a null background for all chromosomes:
  message("Bootstrapping background distributions for each chromosome...")
  nulls <- list()
  nulls_sum <- list()
  for(i in 1:length(chromDictObj)){
    d <- chromDictObj[[i]]

    if(nrow(d) < 1000){
      message(paste0("\n",names(chromDictObj)[i], " is sparse! Will use another chromosome as the background for DMRs on this chromosomes..." ))
      nulls <- append(nulls, NA)
      nulls_sum <- append(nulls_sum, NA)
    } else {
      if(nrow(d)*0.001 < 100){
        nss <- 100
      } else{
        nss <- floor(nrow(d)*0.01)
      }
      unlist(lapply(1:perms, FUN=function(x){
        mean(d$PC_Score[sample(1:nrow(d), nss)])
      })) -> tmpnull

      nulls <- append(nulls, list(tmpnull))
      nulls_sum <- append(nulls_sum, list(c(mean(tmpnull), sd(tmpnull))))
    }


  }

  if(any( !is.na(nulls))){
    i<-sample(which(!is.na(nulls)), 1)

    if(any( is.na(nulls))){
      jvec <- which(is.na(nulls))
      for(j in jvec){
        nulls[[j]] <- nulls[[i]]
        nulls_sum[[j]] <- nulls_sum[[i]]
      }
    }

    nulls[which(is.na(nulls))] <- nulls[[i]]
    nulls_sum[which(is.na(nulls_sum))] <- nulls_sum[[i]]
    names(nulls) <- names(chromDictObj)
    names(nulls_sum) <- names(chromDictObj)
  } else{
    warning("Error: Data is too sparse to calculate DMRs.")
    break
  }

  if(is.null(DMR_resolution)){
    DMR_resolution <- floor(QueryLimit/15)
  }

  if(DMR_resolution==0){
    DMR_resolution <- 1 # infinite while loop safeguard
  }

  if(! "abs.order" %in% colnames(ranks)){
    addRanks(ranks)->ranks
  } else{
    ranks[order(ranks$abs.order, decreasing = T), , drop=F]->ranks
  }

  message("Compressing nearby seeds... ")
  compressed_seeds <- c()
  sub <- row.names(ranks)[1:nSeeds]
  sc <- unique(gsub("\\:.*", "", sub))
  local_lim <- floor(QueryLimit*0.2)
  for(c in sc){
    m <- sub[grepl(paste0(c,":"), sub)]
    m1 <- as.numeric(gsub(".*:", "", m))
    m1[order(m1)] ->m1
    m2 <- m1[2:length(m1)]-m1[1:(length(m1)-1)]

    if(length(which(m2<local_lim)) <= 1){
      compressed_seeds <- c(compressed_seeds, m)
    } else{
      split(which(m2<local_lim), cumsum(c(1, diff(which(m2<local_lim)) != 1))) -> cc
      for(j in 1:length(cc)){
        compressed_seeds <- c(compressed_seeds, paste0(c,":",as.integer(median( m1[min(cc[[j]]-1):max(cc[[j]])]))))
      }
    }
  }
  message(paste("done! Collapsed", nSeeds, "seeds to", length(compressed_seeds),"seeds!"))

  message(paste0("\nExpanding DMRs from ", length(compressed_seeds)," seeds..."))
  get_all_DMRs(chromDictObj, compressed_seeds, res=DMR_resolution, max.dmr.size=QueryLimit,
               min.dmr.cpgs=minCpGs, min.absZscore=minZ, null=nulls_sum) -> DMRs

  regions_list <- DMRs[[2]]
  DMRs <- DMRs[[1]]
  cat(paste0("\nTrimming ", nrow(DMRs)," DMRs..."))
  for(i in 1:nrow(DMRs)){
    trimDMR(DMRs[i,], regions_list[[i]], minCpGs, QueryLimit, nulls_sum, nulls) -> d
    if(is.character(d)){
      DMRs$DMR_size[i] <- 0
    } else{
      DMRs[i,] <- d
    }
  }
  DMRs <- DMRs[DMRs$DMR_size > 0,]
  DMRs <- DMRs[DMRs$nCpGs >= minCpGs,]
  colnames(DMRs) <- c("Chr", "Start", "End", "DMR_Zscore", "nCpGs", "p", "DMR_size")

  #print((Sys.time()-time_s))

  DMRs$FDR <- p.adjust(DMRs$p)
  cat(" done!\n")
  return(DMRs)
}
