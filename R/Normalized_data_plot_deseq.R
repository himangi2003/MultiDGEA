#' Deseq variance normalized visualization
#'
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
#' @importFrom DESeq2 vst
#' @importFrom graphics boxplot
#' @return boxplot of variance normalized DESeq data plot
#' @export
Normalized_data_plot_deseq<-function(count_matrix,meta_data,
                                     model_design_parameter)
{
  dds <- DESeqDataSetFromMatrix(countData = round(count_matrix),
                                colData = meta_data,
                                design= model_design_parameter)


  dds <- data_filtering_Deseq2(dds)
  variance_normalized <- vst((dds), blind = FALSE)
  boxplot(variance_normalized, las=2, main="variance_normalized")
  #colData(variance_normalized)
  return(variance_normalized)
}

