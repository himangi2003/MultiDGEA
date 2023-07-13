#' Function to normalize filtered count matrix for EdgeR analysis
#'
#' @param count_matrix R dataframe containing gene counts with genenames as
#'  rownames and the sample subset selected for the analysis as column names
#' @importFrom edgeR DGEList calcNormFactors
#' @importFrom graphics par boxplot abline
#' @return Normalized count matrix boxplot
#' @export
#'
plot_normalized_data_voom <- function (count_matrix)
{
  #filtering and normalization of count matrix
  data1 <- data_filtering_EdgeR(count_matrix)
  d_ss <-  data_normalization_edgeR(data1)


  pseudoNormCounts <- cpm(d_ss, log = TRUE, prior.count = 1)
  par(mar = c(8,4,1,2))
  a=boxplot(pseudoNormCounts, col = "gray", las = 3, main="TMM normalized")

  return(a)

}
