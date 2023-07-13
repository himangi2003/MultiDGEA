#' Function to fit the analysis using limma voom duplicate correlation
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
#' @param contrast_var string variable used to define the contrast variable
#'  which is the time point variable for the analysis, for example "day"
#' @param patient_id_var string variable used to define the unique patient id
#' variable used for the analysis, for example in the vaccine study it is the
#' variable that stores participant id for each individual like "ptid"
#' @importFrom limma duplicateCorrelation lmFit voom eBayes
#' decideTests topTreat plotMD
#' @importFrom stats as.formula
#' @importFrom edgeR cpm
#' @importFrom dplyr rename
#' @return list of Normalized count matrix and EdgeR fit data result
#' @export
#' @examples
#' \dontrun{
#' result <- limma_voom_duplicate_correlation(count_matrix = counts,
#' meta_data = sample_meta_data,
#' model_design_parameter = "~day",
#' contrast_var="day",
#' patient_id_var="pubid")
#' }
#'
limma_voom_duplicate_correlation <- function(count_matrix,
                                             meta_data,
                                             model_design_parameter,
                                             contrast_var,
                                             patient_id_var)
{
  #filtering and normalization of count matrix
    data1 <- data_filtering_EdgeR(count_matrix)
    d_ss <-  data_normalization_edgeR(data1)


    assign(paste0(contrast_var),meta_data[[contrast_var]])
    assign(paste0(patient_id_var),meta_data[[patient_id_var]])

    model_design_parameter <- as.formula(model_design_parameter)

    #if limma_voom_duplicate_correlation = TRUE fit the model

    model_matrix <-  model.matrix(model_design_parameter, d_ss$samples)

    design <- model_matrix


    v <- voom(d_ss, design)
    patient_id <- meta_data[[patient_id_var]]

    corfit <- duplicateCorrelation(v, design, block = patient_id)

    v <- voom(d_ss, design, block = patient_id, correlation = corfit$consensus)
    fit <- lmFit(v, design, block = patient_id, correlation = corfit$consensus)

    fit_test <- eBayes(fit)
    dt <- decideTests(fit_test)
    plotMD(fit,column=2, status=dt[,2], main=colnames(fit)[2], xlim=c(-8,13))

    res <- topTreat(fit_test,coef=colnames(fit_test[["coefficients"]])[2],
                    number=Inf)


    res <- res %>%
      rename( 'LFC' = 'logFC' , 'P_value' ='P.Value', 'FDR' = 'adj.P.Val')


    return(res)


}
