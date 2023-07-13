#' Scatter plot of P_value comparision to study the comparison between results
#'  of two different analysis type for downregulated genes
#'
#' @param res1  results from first analysis
#' @param res2  results from second analysis
#' @param p_val_xlabel string variable for x axis label for log10(p_value) of res1
#' @param p_val_ylabel string variable for x axis label for log10(p_value) of res2
#' @param p_val_title string variable for title for log10(p_value) scatter comparison plot
#' @param p_val float value P_value to be filtered
#' @param FDR_  float value False discovery rate (p adjusted value) to be choosen
#' @param LFC_ float value for threshold LogFC value for  regulated genes
#' @importFrom dplyr inner_join
#' @importFrom ggplot2 ggplot geom_point ggtitle xlab ylab xlim ylim
#' @importFrom  tibble rownames_to_column
#' @return a plot comparing p_vales of downregulated genes
#' two different analysis types
#' @export
p_value_scatter_plot_comparison_down <- function(res1,res2,
                                               p_val_xlabel,
                                               p_val_ylabel,
                                               p_val_title,
                                               p_val,
                                               FDR_,
                                               LFC_)
{
  res1<-tibble::rownames_to_column(res1,"gene")
  res2<-tibble::rownames_to_column(res2,"gene")

  plot_df_1<-inner_join(res1,res2,by="gene")

  plot_x=plot_df_1[ which(plot_df_1$P_value.x < p_val &
                            plot_df_1$FDR.x < FDR_
                          & plot_df_1$LFC.x < -LFC_),]

  plot_y=plot_df_1[ which(plot_df_1$P_value.y < p_val &
                            plot_df_1$FDR.y < FDR_
                          & plot_df_1$LFC.y < -LFC_),]


  ggplot(plot_df_1, aes(x=log10(P_value.x), y=log10(P_value.y))) +
    geom_point(color="black",size=0.3)+
    geom_point(plot_x, mapping = aes(x=log10(plot_x$P_value.x),
                                     y=log10(plot_x$P_value.y)),color="red",size=1,
               shape = 1)+
    geom_point(plot_y, mapping = aes(x=log10(plot_y$P_value.x),
                                     y=log10(plot_y$P_value.y)),color="blue",size=1.5,
               shape = 1)+
    xlab(p_val_xlabel)+
    ylab(p_val_ylabel)+
    xlim(-8,0)+
    ylim(-8,0)+
    ggtitle(p_val_title)

}
