#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Rscript --vanilla script.R file_path file_tsv file_out
# Reproducible PCA-based approach for methylation

# Step 1: combine bismark.cov files into a matrix
covList <- read.delim(args[2], header=F)

for(i in 1:nrow(covList)){
  name <- as.character(covList[i,2])
  cat(paste0("processing ", name, "...\n"))
  df <- read.delim(paste0(args[1],"/",as.character(covList[i,1])), header=F)
  df$cpgID <- paste0(df$V1,":",df$V2)
  df$count <- df$V5+df$V6
  df <- df[c(7,4,8)]
  colnames(df) <- c("cpgID", paste0(name,"_PercMeth"), paste0(name, "_nCpG"))

  if(i == 1){
    output <- df
  } else {
    cat(paste("merging", name, "to output...\n"))
    output <- merge(output, df, by="cpgID", all=T)
  }
}

# new:
rm_rows <- c()
cat("Removing rows with missing data...\n")
#for(id in c("trt", "ctl")){
for(id in unique(covList[[3]])){
  output[which(grepl(id, colnames(output)) & grepl("PercMeth", colnames(output)))] -> tmp
  tmp$keep <- T
  tmp$keep[c(as.numeric(which(apply(tmp, 1, function(z) sum(is.na(z))) > floor((ncol(tmp)-1)/2))))] <- F
  rm_rows <- c(rm_rows, which(tmp$keep==F))
}

output[-c(unique(rm_rows)),]->output

cat("Filling NA values...\n")
#for(id in c("trt", "ctl")){
for(id in unique(covList[[3]])){
  cols <- which(grepl(id, colnames(output)) & grepl("PercMeth", colnames(output)))
  output[cols] -> tmp
  v<-unique(as.data.frame(which(is.na(tmp), arr.ind=TRUE))$row)
  if(length(v)>0){
    tmp2 <- tmp[v,]
    tmp2[]<-t(zoo::na.aggregate(t(tmp2)))
    tmp[v,] <- tmp2
  }
  output[cols] <- tmp
}

write.table(output, args[3], sep="\t", quote=F, col.names = T, row.names = F)











