#' Function to filter DeSeq object
#'
#' @param DESeqData Deseq object obtained after using Deseg function
#'
#' @return filtered Deseq object
#' @importFrom DESeq2 counts
#' @export
#' @examples
#' \dontrun{
#' dds <- data_filtering_Deseq2(dds)
#'  }
data_filtering_Deseq2 <- function(DESeqData)
{
  #filter genes
  keep <- rowSums(counts(DESeqData) >= 10) >= 10
  dds_filtered <- DESeqData[keep,]
  return(dds_filtered)

}
