#' @export
plot_metagene <- function(data, title="", xaxis="Relative Position",
                          yaxis="PC Score", linecol="red"){

  ggplot22::ggplot(data[[1]], ggplot2::aes(x=bin, y=meanScore))+
    ggplot22::geom_hline(yintercept = 0, linetype="dashed")+
    ggplot22::geom_ribbon(data=data[[2]], mapping=ggplot22::aes(ymin=se_min, ymax=se_max, x=bin, y=meanScore),
                fill="grey", alpha=0.4, stat = "smooth", formula = y ~ x, method = "loess", colour = NA)+
    ggplot22::geom_smooth(colour=linecol, se=F, formula = y ~ x, method = "loess")+
    Ol_Reliable()+ ggplot22::xlab(xaxis)+ ggplot22::ylab(paste0(yaxis, "\n"))+ ggplot22::ggtitle(title)+
    ggplot22::scale_x_continuous(expand=c(0,0)) -> g

  return(g)
}
