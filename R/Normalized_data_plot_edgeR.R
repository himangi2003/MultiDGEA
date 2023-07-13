#' TMM and voom normalized visualization
#'
#' @param result result obtained by fitting the data
#' @param count_matrix R dataframe containing gene counts with genenames as
#'  rownames and the sample subset selected for the analysis as column names
#' @importFrom edgeR cpm
#' @importFrom graphics boxplot par
#' @return boxplot of normalized EdgeR data plot
#' @export
Normalized_data_plot_edgeR<-function(result,count_matrix)
{

  #filtering and normalization of count matrix
  data1 <- data_filtering_EdgeR(count_matrix)
  d_ss <- data_normalization_edgeR(data1)
  pseudoNormCounts <- cpm(d_ss, log = TRUE, prior.count = 1)
  par(mar = c(8,4,1,2))
  boxplot(pseudoNormCounts, col = "gray", las = 3, main="TMM normalized")

}
