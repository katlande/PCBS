#' @export
score_metagene <- function(ranks, regions, bin=100, title="", xaxis="Relative Position",
                           yaxis="PC Score", chromDictObj=NULL, return.data=F, linecol="red"){

  if(is.null(chromDictObj)){
    message("Creating chromDict Object. To speed up computing times, this object can be generated with the chromDict() function, and specified in with the chromDictObj option.")
    chromDict(ranks) -> chromDictObj
  }

  apply(regions, 1, function(x){
    s <- as.numeric(x[[2]])
    e <- as.numeric(x[[3]])
    tmp <- chromDictObj[[unlist(x[[1]])]]
    tmp <- tmp[tmp$pos >= s & tmp$pos <= e,]
    if(nrow(tmp) > 2){
      tmp$relpos <- ceiling((tmp$pos-s)*(bin/(e-s)))
      tmp <- tmp[tmp$relpos > 0 & tmp$relpos <= bin,]
      out <- aggregate(tmp$PC_Score, by = list(tmp$relpos), FUN="median")
    }
  }) -> all_score_file
  dplyr::bind_rows(all_score_file, .id = "column_label") -> all_score_file
  all_score_file <- all_score_file[-c(1)]

  # mean of all scores per bin
  cbind(aggregate(all_score_file$x, by = list(all_score_file$Group.1), FUN="mean"),
        aggregate(all_score_file$x, by = list(all_score_file$Group.1), FUN="se"))[-c(3)] -> metagene
  colnames(metagene) <- c("bin", "meanScore", "SE")

  # Mean score with +/- SE
  data.frame(bin=metagene$bin, se_min=loess(meanScore ~ bin, data = metagene)$fitted-metagene$SE,
             se_max=loess(meanScore ~ bin, data = metagene)$fitted+metagene$SE,
             meanScore=metagene$meanScore) -> se_fit


  if(return.data==F){
    # Save plot
    g <- ggplot2::ggplot(metagene, ggplot2::aes(x=bin, y=meanScore))+
      ggplot2::geom_hline(yintercept = 0, linetype="dashed")+
      ggplot2::geom_ribbon(data=se_fit, mapping=ggplot2::aes(ymin=se_min, ymax=se_max, x=bin, y=meanScore),
                  fill="grey", alpha=0.4, stat = "smooth", formula = y ~ x, method = "loess", colour = NA)+
      ggplot2::geom_smooth(colour=linecol, se=F, formula = y ~ x, method = "loess")+
      Ol_Reliable()+ ggplot2::xlab(xaxis)+ ggplot2::ylab(paste0(yaxis, "\n"))+ ggplot2::ggtitle(title)+
      ggplot2::scale_x_continuous(expand=c(0,0))
    return(g)
  } else {
    return(list(metagene, se_fit))
  }
}
