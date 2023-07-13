#' Function to plot normalized data plot
#'
#' @param method string to choose the analysis to run : limma_voom, limma_voom_duplicate_corr,
#' deseq
#' @param result result obtained from fitting data
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
#' @return normalized data plot
#' @export

Normalized_data_plot <- function(method,result,count_matrix,meta_data,
                                 model_design_parameter)
{
  if(method == "limma_voom" || "limma_voom_duplicate_correlation")

  {
    Normalized_data_plot_edgeR(result,count_matrix)
  }

  else if(method=="deseq")
  {
    Normalized_data_plot_deseq(count_matrix,meta_data,
                               model_design_parameter)

  }
}
