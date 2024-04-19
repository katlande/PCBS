#' @export
CheckOvercompression <- function(ranks, CpG_cutoff=NULL,
                            values=c(0.6,0.8,0.9,1,1.15,1.3,1.5,1.75,2,2.25,2.5,2.75,3),
                            max.dmr.size=5000, return.plot=T){

  message(paste0("Checking ",length(values), " seed values for best DMR calling"))

  if(! "abs.order" %in% colnames(ranks)){
    addRanks(ranks)->ranks
  } else{
    ranks[order(ranks$abs.order, decreasing = T), , drop=F]->ranks
  }

  compressed_seeds <- c()
  for(v in values){
    tmp_cut <- ifelse(is.null(CpG_cutoff), as.integer(v), as.integer(CpG_cutoff*v))
    tmp <- c()
    sub <- row.names(ranks)[1:tmp_cut]
    sc <- unique(gsub("\\:.*", "", sub))
    local_lim <- floor(max.dmr.size*0.2)
    for(c in sc){
      m1 <- sub[grepl(c, sub)]
      m1 <- as.numeric(gsub(".*:", "", m1))
      m1[order(m1)] ->m1
      m2 <- m1[2:length(m1)]-m1[1:(length(m1)-1)]
      split(which(m2<local_lim), cumsum(c(1, diff(which(m2<local_lim)) != 1))) -> cc
      tmp <- c(tmp, length(cc))
    }
    compressed_seeds <- c(compressed_seeds, sum(tmp))
  }
  message("done!")

  if(is.null(CpG_cutoff)){
    data.frame(n.seeds=values,
               compressed_seeds=compressed_seeds) -> d
  } else{
    data.frame(n.seeds=as.integer(CpG_cutoff*values),
               compressed_seeds=compressed_seeds) -> d
  }
  t <- d$n.seeds[d$compressed_seeds==max(d$compressed_seeds)][1]


  if(t == min(d$n.seeds)){
    message("\nSmallest tested nSeed number is the optimum; a smaller number of seeds should be used to avoid overcompression..")
  } else if(t == max(d$n.seeds)){
    message("\nNo overcompression detected!.")
  } else{
    message(paste0("\nOvercompression begins around nSeed=", t, ". This is an optimal seed number for DMR calling."))
  }

  if(return.plot==T){
    ggplot2::ggplot(d, ggplot2::aes(x=n.seeds, y=compressed_seeds))+
      ggplot2::geom_point()+ggplot2::geom_line()+Ol_Reliable()+
      ggplot2::xlab("\nn Seeds Input")+ ggplot2::ylab("n Seeds After Compression\n")+
      ggplot2::geom_vline(xintercept = t, linetype="dashed")-> g
    return(g)
  } else {return(t)}

}
