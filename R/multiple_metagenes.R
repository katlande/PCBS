#' @export
multiple_metagenes <- function(data_list, set_names, title="", xaxis="Relative Position",
                               yaxis="PC Score", legend.title=F, col=NULL, se_alpha=0.25){
  for(i in 1:length(data_list)){
    data_list[[i]][[1]]$Set <- set_names[i]
    data_list[[i]][[2]]$Set <- set_names[i]

    if(i==1){
      all_sets <- data_list[[i]][[1]]
      all_SEs <- data_list[[i]][[2]]
    } else{
      all_sets <- rbind(all_sets, data_list[[i]][[1]])
      all_SEs <- rbind(all_SEs, data_list[[i]][[2]])
    }
  }

  ggplot2::ggplot(all_sets, ggplot2::aes(x=bin, y=meanScore, colour=Set, fill=Set))+
    ggplot2::geom_hline(yintercept = 0, linetype="dashed")+
    ggplot2::geom_ribbon(data=all_SEs, mapping=ggplot2::aes(ymin=se_min, ymax=se_max, x=bin, y=meanScore),
                alpha=se_alpha, stat = "smooth", formula = y ~ x, method = "loess", colour = NA)+
    ggplot2::geom_smooth(se=F, formula = y ~ x, method = "loess")+
    Ol_Reliable()+ ggplot2::xlab(xaxis)+ ggplot2::ylab(paste0(yaxis, "\n"))+ ggplot2::ggtitle(title)+
    ggplot2::scale_x_continuous(expand=c(0,0)) -> g

  # Additional Aesthetic Parameters:
  if(legend.title==F){ g <- g+ggplot2::theme(legend.title = ggplot2::element_blank()) }
  if(! is.null(col)){ g <- g+ggplot2::scale_colour_manual(values = col)+ggplot2::scale_fill_manual(values = col) }
  return(g)
}
