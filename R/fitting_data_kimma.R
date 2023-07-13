#' fitting data using kimma method
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
#' @param kimma_model string variable  used to define the model design
#'  parameter used to perform the kimma analysis e.g "day+ptid"
#' @param contrast_var string variable used to define the contrast variable
#'  which is the time point variable for the analysis, for example "day"
#' @param patient_id_var string variable used to define the unique patient id
#' variable used for the analysis, for example in the vaccine study it is the
#' variable that stores participant id for each individual like "ptid"
#' @importFrom limma lmFit voom decideTests topTreat plotMD
#' @importFrom stats as.formula
#' @importFrom edgeR cpm
#' @importFrom tibble add_column
#' @importFrom kimma kmFit
#' @importFrom dplyr rename
#' @return result of kimma fitted data
#' @export
#' @examples
#'\dontrun{
#' result <- fitting_data_kimma(count_matrix = counts,
#' meta_data = sample_meta_data,
#' model_design_parameter = "~day",
#' kimma_model=  "~ day + (1|pubid)",
#' contrast_var="day",
#' patient_id_var="pubid")
#' }
#'

fitting_data_kimma <- function(count_matrix,
                               meta_data,
                               model_design_parameter,
                               kimma_model,
                               contrast_var,
                               patient_id_var)
{


  data1 <- data_filtering_EdgeR(count_matrix)
  d_ss <-  data_normalization_edgeR(data1)


  assign(paste0(contrast_var),meta_data[[contrast_var]])
  assign(paste0(patient_id_var),meta_data[[patient_id_var]])

  model_design_parameter <- as.formula(model_design_parameter)

  model_matrix <-  model.matrix(model_design_parameter, d_ss$samples)
  design <- model_matrix

  v <- voom(d_ss, design)




  v$targets = tibble::add_column(v$targets,
                                 libID=colnames(count_matrix),
                                 meta_data[[patient_id_var]],
                                 meta_data[[contrast_var]])


  colnames(v$targets)[5:6]<-c(patient_id_var,contrast_var)

  klm <- kmFit(dat = v,
               model = kimma_model,
               run.lme = TRUE,
               use.weights = TRUE,
               metrics = TRUE,
               patientID = patient_id_var,
               run.contrast = TRUE,
               contrast.var = contrast_var,
               processors=3)


  res <- data.frame(klm$lme.contrast)

  res <- res %>%
    rename( 'LFC'='estimate', 'P_value' = 'pval')


  return(res)
}
