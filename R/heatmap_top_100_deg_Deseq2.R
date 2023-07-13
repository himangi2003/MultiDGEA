#' heatmap plot for DESeq2
#'
#' @param result result obtained from fitting the Deseq data
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
#' @importFrom stats heatmap
#' @importFrom utils head
#'
#' @return heatmap plot for DESeq2
#' @export
heatmap_top_100_deg_Deseq2 <- function(result,count_matrix,meta_data,
                                       model_design_parameter)
{

  normalized_counts_deseq <- data_normalization_deseq(count_matrix,
                                                      meta_data,
                                                      model_design_parameter)


  result_s <- result[order(result$FDR), ]

  ## Merge with normalized count dat
  result_sdata <- merge(as.data.frame(result_s),
                      as.data.frame(counts(normalized_counts_deseq,normalized=T)),
                      by = "row.names",
                      sort = FALSE)

  names(result_sdata)[1] <- "Gene"
  sample_s<- normalized_counts_deseq$samplename

  #select only the columns containing the gene names and count data
  x<-subset(result_sdata, select = c("Gene",sample_s))

  #make the table a data frame with gene names then remove duplicate gene name column
  y<-(as.data.frame(x, row.names = x$Gene))
  x<-subset(y,select=-c(Gene))

  #scale rows
  xt<-t(x)
  xts<-scale(xt)
  xtst<-t(xts)

  #only grab top 100 by p-value
  h<-head(xtst, n = 100L)
  a= heatmap(h,main="Top 100 most variable genes across samples", scale="row")
  #aheatmap(h,main="Top 100 most variable genes across samples",annCol=meta_data$day, scale="row")
  a

}
