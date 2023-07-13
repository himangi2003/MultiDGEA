#' Volcano plot for EdgeR or DESeq2 analysis
#'

#' @param result result obtained from fitting data
#' @param plot_title string containing title of the plot
#' @importFrom EnhancedVolcano EnhancedVolcano
#' @return volcano plot
#' @export
#'

volcano_plot <- function (result,plot_title)
{

 a= EnhancedVolcano(result,
                  lab = rownames(result),
                  x = 'LFC',
                  y = 'P_value',
                  title = plot_title,
                  pointSize = 3.0,
                  labSize = 6.0)
 return(a)
}

