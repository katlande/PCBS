#' @export
plot_metagene <- function(data, title="", xaxis="Relative Position",
                          yaxis="PC Score", linecol="red"){

  ggplot2::ggplot(data[[1]], ggplot2::aes(x=bin, y=meanScore))+
    ggplot2::geom_hline(yintercept = 0, linetype="dashed")+
    ggplot2::geom_ribbon(data=data[[2]], mapping=ggplot2::aes(ymin=se_min, ymax=se_max, x=bin, y=meanScore),
                fill="grey", alpha=0.4, stat = "smooth", formula = y ~ x, method = "loess", colour = NA)+
    ggplot2::geom_smooth(colour=linecol, se=F, formula = y ~ x, method = "loess")+
    Ol_Reliable()+ ggplot2::xlab(xaxis)+ ggplot2::ylab(paste0(yaxis, "\n"))+ ggplot2::ggtitle(title)+
    ggplot2::scale_x_continuous(expand=c(0,0)) -> g

  return(g)
}
