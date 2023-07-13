#' Function to normalize filtered count matrix for EdgeR analysis
#'
#' @param filtered_count_matrix R dataframe containing filtered gene counts with
#'  genenames as rownames and the sample subset selected for the analysis
#'  as column names
#' @importFrom edgeR DGEList calcNormFactors
#' @return Normalized count matrix
#' @export
#'
#' @examples
#' \dontrun{
#' d_ss <- data_normalization_edgeR(data1)
#' }
data_normalization_edgeR <- function(filtered_count_matrix)
{
  # creating a DGE object
  normalized_matrix <- DGEList(filtered_count_matrix)
  normalized_matrix <- calcNormFactors(normalized_matrix,method = "TMM")
  return(normalized_matrix)

}
