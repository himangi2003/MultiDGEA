#' Function to filter EdgeR count matrix
#'
#' @param count_matrix R dataframe containing gene counts with genenames as
#'  rownames and the sample subset selected for the analysis as column names
#' @importFrom edgeR cpm
#' @return Filtered count matrix
#' @export
#'
#' @examples
#' \dontrun{
#' data1 <- data_filtering_EdgeR(count_matrix)
#' }
data_filtering_EdgeR <- function(count_matrix)
{
  cpm <- cpm(count_matrix)
  lcpm <- cpm(count_matrix, log=TRUE)
  # fitering  genes
  keep <- rowSums(cpm(count_matrix)>1) >= 15
  d0 <- count_matrix[keep,]
  nsamples <- ncol(d0)


  return(d0)
}

