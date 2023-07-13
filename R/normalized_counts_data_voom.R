#' Function to normalize filtered count matrix for EdgeR analysis
#'
#' @param count_matrix R dataframe containing gene counts with genenames as
#'  rownames and the sample subset selected for the analysis as column names
#' @importFrom edgeR DGEList calcNormFactors
#' @return Normalized count matrix
#' @export
#'
#'
normalized_counts_data_voom <- function (count_matrix)
{
  #filtering and normalization of count matrix
  data1 <- data_filtering_EdgeR(count_matrix)
  d_ss <-  data_normalization_edgeR(data1)

  return(d_ss)

}
