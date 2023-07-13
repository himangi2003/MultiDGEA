#' Print Venn diagram for Up regulated and down downregulated gene
#'
#' @param sig_res significant gene results for first comaparision
#' @param sig_res2 significant gene results for second comaparision
#' @param sig_res_category string variable for the title for the signicant gene results
#'  for both comaprision
#' @importFrom VennDiagram draw.pairwise.venn
#' @importFrom grid grid.draw grid.newpage
#' @return venn plot
#' @export
#'
#' @examples
#' \dontrun{
#'  sig_res<-sig_result_0_3
#'  sig_res2<-sig_result_0_3_trtment_grp
#'  sig_res_category <- c("0 vs 3", "0_3_trtment_grp")
#'  plot_venn(sig_res,sig_res2,sig_res_category)
#' }
#'
plot_venn <- function(sig_res,sig_res2,sig_res_category)
{
  venn.plot_1 <- draw.pairwise.venn(length(rownames(sig_res$up_regulated)),
                                    length(rownames(sig_res2$up_regulated)),
                                    # Calculate the intersection of the two sets
                                    length( intersect(rownames(sig_res$up_regulated),
                                                      rownames(sig_res2$up_regulated)) ),
                                    category = sig_res_category, scaled = F,
                                    fill = c("light blue", "pink"), alpha = rep(0.5, 2),
                                    cat.pos = c(0, 0))

  # Actually plot the plot
  grid.draw(venn.plot_1)

  grid.newpage()

  venn.plot_2 <- draw.pairwise.venn(length(rownames(sig_res$down_regulated)),
                                    length(rownames(sig_res2$down_regulated)),
                                    # Calculate the intersection of the two sets
                                    length( intersect(rownames(sig_res$down_regulated),
                                                      rownames(sig_res2$down_regulated)) ),
                                    category = sig_res_category, scaled = F,
                                    fill = c("yellow", "pink"), alpha = rep(0.5, 2),
                                    cat.pos = c(0, 0))

  # Actually plot the plot
  grid.draw(venn.plot_2)
}
