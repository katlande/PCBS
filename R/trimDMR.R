trimDMR <- function(df, region, min.dmr.cpgs, max.dmr.size, null_summary, null_values){
  if(nrow(region) < min.dmr.cpgs){
    return("nulldmr")
  } else {

    null_summary[[df$chrom[1]]][1] -> n
    null_summary[[df$chrom[1]]][2] -> n_sd
    mean(region$PC_Score) -> u
    if(sign(u)==1){

      quantile(region$PC_Score, 0.85) -> cut
      if(! any(region$PC_Score >= cut)){
        new_region <- region[1,]
      } else{
        tmpve <- split(which(region$PC_Score >= cut), cumsum(c(1, diff(which(region$PC_Score >= cut)) > floor(nrow(region)*0.05))))
        r <- unlist(tmpve[ which(unlist(lapply(tmpve, length)) > max(unlist(lapply(tmpve, length)))*0.66)])
        new_region <- region[min(r):max(r),]
      }
    } else{
      quantile(region$PC_Score, 0.15) -> cut
      if(! any(region$PC_Score <= cut)){
        new_region <- region[1,]
      } else{
        tmpve <- split(which(region$PC_Score <= cut), cumsum(c(1, diff(which(region$PC_Score <= cut)) > floor(nrow(region)*0.05))))
        r <- unlist(tmpve[ which(unlist(lapply(tmpve, length)) > max(unlist(lapply(tmpve, length)))*0.66)])
        new_region <- region[min(r):max(r),]
      }
    }
    if(nrow(new_region) > 3){
      p_t <- t.test(null_values[[df$chrom[1]]], new_region$PC_Score)$p.value
      data.frame(chrom=df$chrom[1],
                 start=min(new_region$pos),
                 end=max(new_region$pos),
                 Zscore=((mean(new_region$PC_Score))-n)/n_sd,
                 nCpGs=nrow(new_region),
                 DMR_score=p_t,
                 DMR_size=(max(new_region$pos)-min(new_region$pos))) -> d
      return(d)
    } else{
      return("nulldmr")
    }
  }
}
