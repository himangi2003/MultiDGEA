#' Function to find normalized counts for Deseq analysis
#' @param count_matrix R dataframe containing gene counts with genenames as
#'  rownames and the sample subset selected for the analysis as column names
#' @param meta_data R dataframe containing information about each sample
#' with subset samples selected for the analysis as rownames and the information
#' about samples like the day, participant ID as column names
#' @param model_design_parameter string variable used to define
#'  the model design parameter to be used to run the analysis.
#'  for example "~day+ptid", note: it should be the same variable or combination
#'  of variable used to define the contrast variable and patient id variable
#'  for the analysis
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors
#' estimateDispersions plotDispEsts DESeq results
#' @return list of Normalized count matrix
#' @export

data_normalization_deseq <- function(count_matrix,meta_data,
                                     model_design_parameter)
{

  dds <- DESeqDataSetFromMatrix(countData = round(count_matrix),
                                colData = meta_data,
                                design= model_design_parameter)


  dds <- data_filtering_Deseq2(dds)
  #estimate dispersion in DEseq obect
  dds <- estimateSizeFactors(dds)
  dds = estimateDispersions( dds )

  return(dds)

}

