#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Rscript --vanilla script.R file_path file_tsv file_out
# Reproducible PCA-based approach for methylation

# Combine bismark.cov files into a matrix
# Here the script checks for headers, and only proceeds if the first column of the file_tsv input contains files with a ".cov" 
covList <- read.delim(args[2], header=F)
if(grepl("\\.cov", covList[1,1]) == F){
  covList <- read.delim(args[2], header=T)
}
if(grepl("\\.cov", covList[1,1])==F){
  message(paste("Error:", args[2], "column 1 does not contain .cov files!")) # break the script if there are no .cov files supplied
} else{
  for(i in 1:nrow(covList)){
    name <- as.character(covList[i,2])
    condition <- as.character(covList[i,3])
    message(paste0("processing ", name, "..."))
    df <- read.delim(paste0(args[1],"/",as.character(covList[i,1])), header=F)
    df$cpgID <- paste0(df$V1,":",df$V2)
    df$count <- df$V5+df$V6
    df <- df[c(7,4,8)]
    colnames(df) <- c("cpgID", paste0(name, "_", condition, "_PercMeth"), paste0(name, "_", condition, "_nCpG"))

    if(i == 1){
      output <- df
    } else {
      message(paste("merging", name, "to output..."))
      output <- merge(output, df, by="cpgID", all=F) # remove NA rows
    }
  }
  message("Removing NA rows...")
  output <- na.omit(output)
  message("Writing output...")
  write.table(output, args[3], sep="\t", quote=F, col.names = T, row.names = F)
}



