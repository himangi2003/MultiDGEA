#' Finding the significant from the results
#' @param result result obtained from fitting EdgeR or DESeq2
#' @param p_val float variable for P_value to be filtered
#' @param FDR_  float variable for False discovery rate (p adjusted value) to be choosen
#' @param LFC_ float variable for threshold LogFC value for  regulated genes
#'
#' @return list of signicant genes
#' @export

find_DEG_genes <- function(result,p_val,FDR_,LFC_)
{

    #find out significant genes
    resultSig= result[ which(result$P_value < p_val & result$FDR < FDR_), ]

    resup =  resultSig[ which( resultSig$LFC > LFC_), ]
    resdown <-   resultSig[ which( resultSig$LFC < -LFC_), ]

    # order the resultult from deseq based on pvalues
    result_ad_p_val <- resultSig[order(resultSig$FDR), ]

    result_data<-list(resup,resdown,result_ad_p_val)

    names(result_data) <- c("up_regulated","down_regulated",
                            "DEG")

  return(result_data)

}
