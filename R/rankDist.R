#' @export
rankDist <- function(ranks, draw_intersects=T, noise_perc=0.5, mode="intersect", return.plot=T){

  if(nrow(ranks)<10000){
    warning("Data is very sparse! Cut-off estimation is likely inaccurate")
  }


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

  g <- ggplot2::ggplot(df, ggplot2::aes(x=abs.order, y=abs(PC_Score)))


  i <- as.integer(nrow(ranks)*noise_perc)
  ivec <- c()
  rsq_vec <- c()
  while(i >100){
    suppressMessages(suppressWarnings(summary(lm(abs(ranks$PC_Score)[1:i] ~ ranks$abs.order[1:i]))$adj.r.squared)) -> l
    if(! is.nan(l)){
      ivec <- c(ivec, i)
      rsq_vec <- c(rsq_vec, l)
    }
    i <- floor(i*0.9)
  }

  est <- ranks$abs.order[ivec[which(rsq_vec==max(rsq_vec))]][1]

  lm( abs(ranks$PC_Score)[1:est] ~
        ranks$abs.order[1:est]) -> l1
  lm( abs(ranks$PC_Score)[(nrow(ranks)*noise_perc):nrow(ranks)] ~
        ranks$abs.order[(nrow(ranks)*noise_perc):nrow(ranks)]) -> l2

  int <- as.integer(lmIntx(l1,l2)[1])
  strict<-abs(ranks$PC_Score[ranks$abs.order == est][1])

  if(mode == "intersect"){
    est_fin <- int # intersection of the NOISE fit and the SIGNAL fit
  } else{
    est_fin <- as.integer(mean(c(int,est)))
  }

  if(draw_intersects){

    zero_x <- nrow(ranks)*-0.05
    top_x <- nrow(ranks)*1.05
    zero_y <- max(abs(ranks$PC_Score))*-0.05
    top_y <- max(abs(ranks$PC_Score))*1.05

    data.frame(set=c(rep("A",4), rep("B",4)),
               x=c(zero_x, zero_x,
                   ((top_y-l1$coefficients[1])/l1$coefficients[2]),
                   ((zero_y-l1$coefficients[1])/l1$coefficients[2]),
                   zero_x, zero_x, top_x, top_x),
               y=c(zero_y, top_y, top_y, zero_y,
                   zero_y, lin(zero_x, l2), lin(top_x, l2), zero_y)) -> poly

    g+ggplot2::geom_polygon(data=poly, mapping=ggplot2::aes(x=x,y=y,group=set,fill=set), alpha=0.5,show.legend = F)+
      ggplot2::scale_x_continuous(expand=c(0,0), limits=c(zero_x, top_x))+
      ggplot2::scale_y_continuous(expand=c(0,0), limits=c(zero_y, top_y)) -> g

    if(mode=="intersect"){
      g+ggplot2::annotate("point", x=est_fin, y=lin(est_fin, l1), size=2, colour="red") ->g
    } else{
      g+ggplot2::annotate("point", x=est_fin, y=abs(ranks$PC_Score[ranks$abs.order == est_fin][1]), size=2, colour="red") ->g
    }

  }

  g+ggplot2::geom_line(linewidth=1.25)+Ol_Reliable()+ ggplot2::ylab("PC Score\n")+ ggplot2::xlab("\nPC Score Rank")+
    ggplot2::geom_vline(xintercept=est_fin, colour="red", linetype="dashed")+
    ggplot2::annotate("text", y=max(abs(ranks$PC_Score))*0.8, x=est, hjust=-0.5,
             label=est_fin, colour="red") -> g

  message(paste0("Estimated rank cut-off for significant CpGs is ", est_fin, "."))

  if(return.plot==T){
    return(g)
  } else{
    return(est_fin)
  }


}
