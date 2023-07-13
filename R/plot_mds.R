#' plot MDS for EdgeR data fitting
#'
#' @param method anaysis method to be chosen
#' @param result result obtained from fitting data
#' @param count_matrix R dataframe containing gene counts with genenames as
#'  rownames and the sample subset selected for the analysis as column names
#' @param meta_data R dataframe containing information about each sample
#' with subset samples selected for the analysis as rownames and the information
#' about samples like the day, participant ID as column names
#' @param contrast_var string variable used to define the contrast variable
#'  which is the time point variable for the analysis, for example "day"
#' @param mds_title title for the plot
#' @importFrom limma  plotMDS
#' @importFrom graphics title legend
#'
#' @return plot MDS
#' @export
#'
plot_mds <- function(method,result,count_matrix,meta_data,contrast_var,
                     mds_title)
{
  data1 <- data_filtering_EdgeR(count_matrix)
  d_ss <- data_normalization_edgeR(data1)

  if (method == "limma_voom" || "limma_voom_duplicate_correlation")
  {
    obj_factor<-factor(assign(paste0(contrast_var),
                                          meta_data[[contrast_var]]))

    col.cell <- c("purple","orange")[obj_factor]
    data.frame(meta_data[[contrast_var]],col.cell)

    a=plotMDS(d_ss,col=col.cell)
    a=legend("bottomleft",fill=c("purple","orange"),
             legend=levels(obj_factor))
    a= title(mds_title)
    a
  }

  else
  {
    print("EdgeR not selected...")

  }

}
