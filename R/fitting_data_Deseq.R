#' Function to fit the Deseq analysis using two different kind of
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
#' @importFrom stats as.formula
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors
#' estimateDispersions plotDispEsts DESeq results
#' @importFrom dplyr rename
#' @return list of Normalized count matrix and EdgeR fit data result
#' @export
#' @examples
#' \dontrun{
#' result <- fitting_data_Deseq(count_matrix = counts,
#' meta_data = sample_meta_data,
#' model_design_parameter = "~day",
#' contrast_var="day",
#' patient_id_var="pubid")
#' }
#'
fitting_data_Deseq<- function(count_matrix,
                              meta_data,
                              model_design_parameter,
                              contrast_var,
                              patient_id_var)
{





    assign(paste0(contrast_var),meta_data[[contrast_var]])
    assign(paste0(patient_id_var),meta_data[[patient_id_var]])

    model_design_parameter <- as.formula(model_design_parameter)

    dds <- DESeqDataSetFromMatrix(countData = round(count_matrix),
                                  colData = meta_data,
                                  design= model_design_parameter)


    dds <- data_filtering_Deseq2(dds)
    #estimate dispersion in DEseq obect
    dds <- estimateSizeFactors(dds)
    dds = estimateDispersions( dds )
    #plot dispersion estimate
    plotDispEsts(dds)
    cds <- DESeq(dds)
    res <- results(cds)

    res<- data.frame(res)
    res <- res %>%
      rename( 'LFC' = 'log2FoldChange', 'P_value'= 'pvalue' , 'FDR' = 'padj')


    return(res)


}
