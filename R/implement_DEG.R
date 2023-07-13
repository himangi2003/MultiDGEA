#' Function to implement all different kind of analysis
#'
#' @param method string to choose the analysis to run : limma_voom, limma_voom_duplicate_corr,
#' deseq
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
#' @param contrast_var string variable used to define the contrast variable
#'  which is the time point variable for the analysis, for example "day"
#' @param patient_id_var string variable used to define the unique patient id
#' variable used for the analysis, for example in the vaccine study it is the
#' variable that stores participant id for each individual like "ptid"
#' @param kimma_model string variable  used to define the model design
#'  parameter used to perform the kimma analysis e.g "day+ptid"
#' @param form  string variable used to define the variable to be tested
#'  for a fixed effect e.g. "~ day + (1/ptid) "
#' @return list of Normalized count matrix and fit data result for EdgeR or
#' DESeq2
#' @export
#'
implement_DEG <- function(method,
                          count_matrix,
                          meta_data,
                          model_design_parameter,
                          form,
                          kimma_model,
                          contrast_var,
                          patient_id_var)
  {



  if(missing(form))
    {
    form <- NULL
  }
  else {
    form <- form
  }



  if(missing(kimma_model)) {
    kimma_model <- NULL
  }
  else {
    kimma_model <- kimma_model
  }


  if(missing(contrast_var)) {
    contrast_var <- NULL
  }
  else {
    contrast_var <- contrast_var
  }


  if(missing(patient_id_var)) {
    patient_id_var <- NULL
  }
  else {
    patient_id_var <- patient_id_var
  }



  if(method == "limma_voom")
  {
    fitting_data_limma_voom(count_matrix,
                            meta_data,
                            model_design_parameter,
                            contrast_var,
                            patient_id_var)
  }

  else if(method == "limma_voom_duplicate_correlation")
  {
    limma_voom_duplicate_correlation(count_matrix,
                                     meta_data,
                                     model_design_parameter,
                                     contrast_var,
                                     patient_id_var)
  }

  else if(method == "deseq")
  {
    fitting_data_Deseq(count_matrix,
                       meta_data,
                       model_design_parameter,
                       contrast_var,
                       patient_id_var)

  }

  else if(method == "dream")
  {
    fitting_data_DREAM(count_matrix,
                       meta_data,
                       model_design_parameter,
                       form,
                       contrast_var,
                       patient_id_var)

  }

  else if(method == "kimma")
  {
    fitting_data_kimma(count_matrix,
                       meta_data,
                       model_design_parameter,
                       kimma_model,
                       contrast_var,
                       patient_id_var)

  }


  else
  {
    print("wrong analysis method chosen")

  }

}
