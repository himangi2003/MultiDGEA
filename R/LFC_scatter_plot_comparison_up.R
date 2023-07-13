#' Scatter plot of LFC value comparision to study the comparison between results
#'  of two different analysis type for upregulated genes
#'
#' @param res1  results from first analysis
#' @param res2  results from second analysis
#' @param LFC_xlabel string variable for x axis label for log10(p_value) of res1
#' @param LFC_ylabel string variable for x axis label for log10(p_value) of res2
#' @param LFC_title string variable for title for log10(p_value) scatter comparison plot
#' @param p_val float value P_value to be filtered
#' @param FDR_  float value False discovery rate (p adjusted value) to be choosen
#' @param LFC_ float value for threshold LogFC value for  regulated genes
#' @importFrom dplyr inner_join
#' @importFrom ggplot2 ggplot geom_point ggtitle xlab ylab xlim ylim
#' @importFrom  tibble rownames_to_column
#' @return a plot comparing LFC vales of upregulated genes
#' two different analysis types
#' @export

LFC_scatter_plot_comparison_up <- function(res1,res2,
                                                   LFC_xlabel,
                                                   LFC_ylabel,
                                                   LFC_title,
                                                   p_val,
                                                  FDR_,
                                                  LFC_)
{
  res1<-tibble::rownames_to_column(res1,"gene")
  res2<-tibble::rownames_to_column(res2,"gene")

  plot_df_1<-dplyr::inner_join(res1,res2,by="gene")

  plot_x=plot_df_1[ which(plot_df_1$P_value.x < 0.05 &
                            plot_df_1$FDR.x < 0.2
                          & plot_df_1$LFC.x > 0.5),]
  plot_y=plot_df_1[ which(plot_df_1$P_value.y < 0.05 &
                            plot_df_1$FDR.y < 0.2
                          & plot_df_1$LFC.y > 0.5),]



  ggplot(plot_df_1, aes(x=LFC.x, y=LFC.y))+
    geom_point(color="black",size=0.3)+
    geom_point(plot_x, mapping = aes(x=plot_x$LFC.x,
                                     y=plot_x$LFC.y),color="red",size=1,
               shape = 1)+
    geom_point(plot_y, mapping = aes(x=plot_y$LFC.x,
                                     y=plot_y$LFC.y),color="blue",size=1.5,
               shape = 1)+
    xlab(LFC_xlabel)+
    ylab(LFC_ylabel)+
    ggtitle(LFC_title)
}



