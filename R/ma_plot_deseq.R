#' plot MA for Deseq
#'
#' @param method  deseq analysis is to be run
#' @param result  result obtained from fitting the data
#' @param thresh  log fold change threshold value
#' @param labelsig boolean for significant threshold display
#' @param textcx size of the display
#'
#' @return MA plot
#' @export
#'
#' @examples
#' \dontrun{
#' maplot_deseq("deseq",res)
#' }
maplot_deseq <- function (method,result, thresh=0.2, labelsig=TRUE, textcx=1)
{
  if(method=="deseq")
  {
    res <- result
    Gene=rownames(res)
    with(res, plot(baseMean, LFC, pch=20, cex=.5))
    with(subset(res, FDR<thresh), points(baseMean, LFC, col="red", pch=20, cex=1.5))
    if (labelsig) {
      with(subset(res, FDR<thresh), points(baseMean, LFC, labs=Gene, cex=textcx, col=2))
    }
  }
  else {
    print("DeSeq analysis not chosen....")
  }

}
