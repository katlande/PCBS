#' @export
mergePair <- function(chromDictObj, df, min.dmr.cpgs){

  chrom <- chromDictObj[[df$chrom[1]]]
  region <- chrom[chrom$pos >= min(df$start) & chrom$pos <= max(df$end),]

  if(length(region$PC_Score[region$pos >= df$start[1] & region$pos <= df$end[1]]) < min.dmr.cpgs |
     mean(region$PC_Score[region$pos >= df$start[2] & region$pos <= df$end[2]]) < min.dmr.cpgs){
    mode <- "merge"
  } else{
    DMR_a <- mean(region$PC_Score[region$pos >= df$start[1] & region$pos <= df$end[1]])
    DMR_b <- mean(region$PC_Score[region$pos >= df$start[2] & region$pos <= df$end[2]])
    DMR_c <- mean(region$PC_Score)

    # if alt sign, split the DMR
    if(! sign(DMR_a) == sign(DMR_b)){
      mode <- "split"
    } else{
      if(abs(DMR_a) > abs(DMR_c) | abs(DMR_b) > abs(DMR_c)){
        mode <- "split"
      } else {
        mode <- "merge"
      }
    }
  }

  if(mode == "split"){
    v <- mean(c(min(df$end),  max(df$start)))
    v1 <- ifelse(floor(v)==v, v-1, floor(v))
    df$end[which(df$start==min(df$start))] <- v1
    df$start[which(!df$start==min(df$start))] <- v
    df[c(4:6,8)] <- NA
  } else {
    df <- data.frame(chrom=df$chrom[1],
                     start=min(df$end),
                     end=max(df$end),
                     Zscore=NA, nCpGs=NA, DMR_score=NA,
                     null=paste0(df$null[1], df$null[2]),
                     DMR_size=NA)
  }

  return(df)
}
