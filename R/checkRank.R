#' @export
checkRank <- function(ranks, cutoff){

  if(! "abs.order" %in% colnames(ranks)){
    addRanks(ranks)->ranks
  } else{
    ranks[order(ranks$abs.order, decreasing = T), , drop=F]->ranks
  }
  # reduce plot size if many points:
  if(nrow(ranks) <= 1e5){
    df <- ranks
  } else if(nrow(ranks) > 1e5 & nrow(ranks) <= 1e7){
    df <- ranks[unique(c(1:floor(nrow(ranks)*0.01),
                         sample(ceiling(nrow(ranks)*0.01):floor(nrow(ranks)*0.1), (nrow(ranks)*0.05)),
                         sample(ceiling(nrow(ranks)*0.1):nrow(ranks), (nrow(ranks)*0.01))))   ,]

  } else if(nrow(ranks) > 1e7 & nrow(ranks) <= 1e8){
    df <- ranks[unique(c(1:floor(nrow(ranks)*0.005),
                         sample(ceiling(nrow(ranks)*0.01):floor(nrow(ranks)*0.05), (nrow(ranks)*0.01)),
                         sample(ceiling(nrow(ranks)*0.05):nrow(ranks), (nrow(ranks)*0.001))))   ,]
  } else{
    df <- ranks[unique(c(1:floor(nrow(ranks)*0.001),
                         sample(ceiling(nrow(ranks)*0.001):floor(nrow(ranks)*0.01), (nrow(ranks)*0.005)),
                         sample(ceiling(nrow(ranks)*0.01):nrow(ranks), (nrow(ranks)*0.0005))))   ,]
  }

  g <- ggplot2::ggplot(df, ggplot2::aes(x=abs.order, y=abs(PC_Score)))+
    ggplot2::geom_line(linewidth=1.25)+Ol_Reliable()+ ggplot2::ylab("PC Score\n")+ ggplot2::xlab("\nPC Score Rank")+
    ggplot2::annotate("text", y=max(abs(ranks$PC_Score))*0.8, x=cutoff, hjust=-0.5,
                      label=cutoff, colour="red")+
    ggplot2::annotate("point", x=cutoff, y=abs(ranks$PC_Score[ranks$abs.order == cutoff][1]), size=2, colour="red")-> g

  return(g)
}
