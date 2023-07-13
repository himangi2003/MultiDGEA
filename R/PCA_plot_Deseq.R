#' plot PCA for DESeq2
#'
#' @param method Deseq analysis
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
#' @param pca_group groups to be used to check the PCA plot
#' @importFrom DESeq2 vst plotPCA
#' @importFrom ggplot2 ggplot aes geom_point xlab ylab coord_fixed
#' @importFrom graphics title
#' @return pca plot
#' @export
#'
#' @examples
#' \dontrun{
#' pca_group <- c("day", "Treatment_Group")
#' PCA_plot_Deseq(Deseq,res,pca_group)
#' }
PCA_plot_Deseq<-function(method,count_matrix,meta_data,
                         model_design_parameter,pca_group)
{
  if(method=="deseq")
  {
    dds <- DESeqDataSetFromMatrix(countData = round(count_matrix),
                                  colData = meta_data,
                                  design= model_design_parameter)

    vsd <- vst(dds)
    pca_plot_Data <- plotPCA(vsd, intgroup = pca_group, returnData = TRUE)
    percent_Variance <- round(100 * attr(pca_plot_Data, "percentVar"))

    color_pca <- pca_group[1]
    shape_pca <- pca_group[2]
    a<-ggplot(pca_plot_Data, aes(x = pca_plot_Data$PC1, y =  pca_plot_Data$PC2,
                                 color = pca_plot_Data[[color_pca]],
                                 shape = pca_plot_Data[[shape_pca]])) +
      geom_point(size =3) +
      xlab(paste0("PC1: ", percent_Variance[1], "% variance")) +
      ylab(paste0("PC2: ", percent_Variance[2], "% variance")) +
      coord_fixed()
    a

  }
  else
  {
    print("DeSeq not choosen......")
  }


}
