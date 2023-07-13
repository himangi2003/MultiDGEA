#' heatmap plot for EdgeR
#'
#' @param result result obtained from fitting the EdgeR data
#' @param count_matrix R dataframe containing gene counts with genenames as
#'  rownames and the sample subset selected for the analysis as column names
#' @importFrom stats heatmap
#' @importFrom utils head
#' @return heatmap plot for EdgeR
#' @export
heatmap_top_100_deg_edgeR <- function(result,count_matrix)
{


  #filtering and normalization of count matrix
  data1 <- data_filtering_EdgeR(count_matrix)
  normalized_counts_dataset  <-  data_normalization_edgeR(data1)


  # order the resultult from deseq based on pvalues
  result <- result[order(result$FDR), ]

  ## Merge with normalized count dat
  resultdata <- merge(as.data.frame(result),
                      as.data.frame(normalized_counts_dataset$counts),
                      by = "row.names",
                      sort = FALSE)

  names(resultdata)[1] <- "Gene"
  sample_s<- colnames(normalized_counts_dataset$counts)
  #select only the columns containing the gene names and count data
  x<-subset(resultdata, select = c("Gene",sample_s))

  #make the table a data frame with gene names then remove duplicate gene name column
  y<-(as.data.frame(x, row.names = x$Gene))
  x<-subset(y, select=-c(Gene))

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
