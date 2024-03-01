#' @export
DefineBestPC <- function(mat, IDs, filter_thresh=50, return.plot=T){
  mat$totalCounts <- rowSums(mat[which(grepl("nCpG", colnames(mat)))])
  mat <- mat[mat$totalCounts > filter_thresh,]
  row.names(mat) <- c()
  tibble::column_to_rownames(mat, "cpgID")->mat
  pca <- prcomp(t(mat[which(grepl("PercMeth", colnames(mat)))]))
  as.data.frame(pca$x[]) -> pca_info # extract PCA values

  # Logic gate to check if assigned IDs make sense
  if(! length(c(which(grepl(IDs[1], row.names(pca_info))), which(grepl(IDs[2], row.names(pca_info))))) == (ncol(mat)-1)/2){

    if(length(c(which(grepl(IDs[1], row.names(pca_info))), which(grepl(IDs[2], row.names(pca_info))))) > (ncol(mat)-1)/2){
      warning("Error: some samples are assigned both IDs!")
    } else {
      warning("Error: some samples are not assigned IDs!")
    }


  } else{
    # estimate the best PC to used based on what we are comparing:
    dist_v <- c()
    pass_v <- c()
    for(i in 1:ncol(pca_info)){

      pcs1 <- pca_info[which(grepl(IDs[1], row.names(pca_info))), i]
      pcs2 <- pca_info[which(grepl(IDs[2], row.names(pca_info))), i]
      m1 <- mean(pcs1)
      m2 <- mean(pcs2)

      if(m1 > m2){
        pass_v <- c(pass_v, ifelse(any(unlist(lapply(pcs1, FUN=function(x){any(as.numeric(unlist(x)) < pcs2)}))), F, T))
      } else{
        pass_v <- c(pass_v, ifelse(any(unlist(lapply(pcs1, FUN=function(x){any(as.numeric(unlist(x)) > pcs2)}))), F, T))
      }

      dist_v <- c(dist_v, dist(rbind(pca_info[which(grepl(IDs[1], row.names(pca_info))), i],
                                     pca_info[which(grepl(IDs[2], row.names(pca_info))), i])))

    }

    data.frame(dist=dist_v, pass=pass_v, PC=colnames(pca_info))->d

    if(nrow(d[d$pass==T,])>0){

      b <- d[d$pass==TRUE,]
      best_dist <- max(b$dist)
      best_pc <- d$PC[d$pass==T & d$dist==best_dist]

      varPC <- paste0(substr(as.character(((pca$sdev^2/sum(pca$sdev^2))*100)[which(d$PC==best_pc)[1]]),1,5),"%")
      message(paste0("\nBest PC to use is ", best_pc, " with a sample distance of ", formatC(best_dist, digits=3),
                     ", representing ", varPC, " of the total variance."))
    } else {
      message("\nThere are no optimal principle components to use to analyze this data.\n
      You may have outliers or mislabeled samples in you data.\n
      Select a principle component to use manually, or consider an alternative method.")
    }

    # draw best PCA:
    pca_info$sample <- row.names(pca_info)
    pca_info$ID <- ""
    pca_info$ID[grepl(IDs[1], pca_info$sample)] <- IDs[1]
    pca_info$ID[grepl(IDs[2], pca_info$sample)] <- IDs[2]

    if(nrow(d[d$pass==T,])>0){
      if(best_pc %in% c("PC1", "PC2")){
        plot.v <- "B"
      } else{
        plot.v <- "A"
      }
    } else{
      plot.v <- "B"
    }


    if(plot.v=="A"){

      colnames(pca_info)[which(best_pc==colnames(pca_info))] <- "best"
      ggplot2::ggplot(data = pca_info, ggplot2::aes(x=best, y=PC1, colour=ID)) +
        ggplot2::theme_linedraw()+
        ggplot2::geom_point(size=2, alpha=0.75) +
        ggplot2::ylab(paste0("PC1 - ", paste0(substr(as.character(((pca$sdev^2/sum(pca$sdev^2))*100)[1]),1,5),"%")," of Variance\n"))+
        ggplot2::xlab(paste0("\n", best_pc, " - ", paste0(substr(as.character(((pca$sdev^2/sum(pca$sdev^2))*100)[which(d$PC==best_pc)[1]]),1,5),"%")," of Variance"))+
        ggrepel::geom_text_repel(size=3, label=gsub("_PercMeth", "", pca_info$sample))+
        ggplot2::labs(title="PCA of % Methylation")+
        ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5, face="bold", size=14),
              plot.subtitle=ggplot2::element_text(hjust=0.5, face="italic", size=11)) -> pca_plot

    } else {
      ggplot2::ggplot(data = pca_info, ggplot2::aes(x=PC1, y=PC2, colour=ID)) +
        ggplot2::theme_linedraw()+
        ggplot2::geom_point(size=2, alpha=0.75) +
        ggplot2::ylab(paste0("PC2 - ", paste0(substr(as.character(((pca$sdev^2/sum(pca$sdev^2))*100)[2]),1,5),"%")," of Variance\n"))+
        ggplot2::xlab(paste0("\nPC1 - ", paste0(substr(as.character(((pca$sdev^2/sum(pca$sdev^2))*100)[1]),1,5),"%")," of Variance"))+
        ggrepel::geom_text_repel(size=3, label=gsub("_PercMeth", "", pca_info$sample))+
        ggplot2::labs(title="PCA of % Methylation")+
        ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5, face="bold", size=14),
              plot.subtitle=ggplot2::element_text(hjust=0.5, face="italic", size=11)) -> pca_plot

    }

    if(return.plot==T){
      return(pca_plot)
    } else{
      return(as.numeric(gsub("PC", "", best_pc)))
    }

  }

}
