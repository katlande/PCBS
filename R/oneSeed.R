oneSeed <- function(chroms, seed, resolution, max.size, mincpgs, null_list, Zlim=1){
  seed.chrom <- gsub(":.*", "",seed)
  seed.pos <- as.numeric(gsub(".*:", "", seed))
  null <- null_list[[seed.chrom]]
  df <- chroms[[seed.chrom]]

  df.main <- na.omit(df[ .( c((seed.pos-(max.size/2)):(seed.pos+(max.size/2))) ) ])
  if(nrow(df.main) > mincpgs){
    u.pop <- null[1]
    sd.pop <- null[2]
    min_vec <- c()
    max_vec <- c()
    Z_vec <- c()
    nCpG <- c()
    score <- c()

    minval <- seed.pos-resolution
    maxval <- seed.pos+resolution
    while((maxval-minval) <= max.size){

      if(nrow(df.main[df.main$pos >= minval & df.main$pos <= maxval,]) == 0){
        maxval=maxval+resolution
        minval=minval-resolution
        next
      }

      min_vec <- c(min_vec, minval)
      max_vec <- c(max_vec, maxval)

      Z_vec <- c(Z_vec, ((mean(df.main$PC_Score[df.main$pos >= minval & df.main$pos <= maxval])-u.pop)/sd.pop))
      nCpG <- c(nCpG, length(df.main$PC_Score[df.main$pos >= minval & df.main$pos <= maxval]))

      if(Z_vec[length(Z_vec)] > 0){
        score <- c(score, length(which(df.main$PC_Score[df.main$pos >= minval & df.main$pos <= maxval] >=
                                         (mean(df.main$PC_Score[df.main$PC_Score > u.pop+sd.pop])))))
      } else{
        score <- c(score, length(which(df.main$PC_Score[df.main$pos >= minval & df.main$pos <= maxval] <=
                                         (mean(df.main$PC_Score[df.main$PC_Score < u.pop-sd.pop])))))
      }

      maxval=maxval+resolution
      minval=minval-resolution
    }

    data.frame(chrom=seed.chrom,
               start=min_vec,
               end=max_vec,
               Zscore=Z_vec,
               nCpGs=nCpG,
               DMR_score=score) -> out

    out[out$nCpGs >= mincpgs & abs(out$Zscore) >= Zlim,]->out
    if(nrow(out) > 1){
      out[out$DMR_score >= as.numeric(quantile(out$DMR_score,0.9)),]->out
      out <- out[out$nCpGs == min(out$nCpGs),]
      out <- out[out$end-out$start == max(out$end-out$start),]
    }

    return(list(out, df.main))
  } else {
    return("low sparsity")
  }


}
