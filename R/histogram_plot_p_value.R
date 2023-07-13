#' histogram plot for p_value and adjusted p_value
#'  for result of fitting data
#'
#' @param res result obtained from fitting EdgeR or DESeq2
#' @importFrom graphics hist
#' @return histogram plot
#' @export
histogram_plot_p_value <- function (res)
{
  Gene=rownames(res)
  par(mfrow = c(1,2))
  a=hist(res$P_value, xlab = "p-value", main = "raw p-values")
  a=hist(res$FDR, xlab = "p-value", main = "adjusted p-values")
  return(a)
}
