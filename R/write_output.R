#' function to write the fit data output to a folder
#'
#' @param pipeline_input piprline inputs to be used
#' @param res result obtained after fitting the dat set
#' @param method the analysis method choosen
#' @importFrom stringr str_split
#' @importFrom utils write.csv
#' @return csv file containing the output
#' @export
write_output <- function(pipeline_input,res,method)
{
  dd=str_split(pipeline_input$day_subset, pattern = ",")
  filename_ <- paste0('res_',dd[[1]][1], '_', dd[[1]][2],'_',method,".csv")

  path <- pipeline_input$output_directory
  new_path<-paste0(path,method)

  if (!dir.exists(new_path))
  {
    dir.create(new_path, showWarnings = TRUE, recursive = FALSE)

  }
  else
  {
    print("directory already exists")
  }

  return(write.csv(res, paste(new_path,"/", filename_,sep = '')))
}
