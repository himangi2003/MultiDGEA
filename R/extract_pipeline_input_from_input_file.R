#' Function to extracts input from the pipeline
#'
#' @param arguments_file file containing all arguments
#' @importFrom readr read_csv locale cols
#' @importFrom dplyr select filter
#' @importFrom rlang set_names
#' @return set of inputs for the pipeline
#' @export
#'
#' @examples
#' \dontrun{
#' pipeline_input<-extract_pipeline_input_from_input_file("arguments_file")
#' }
extract_pipeline_input_from_input_file <- function(arguments_file)
{

  # Remove non-existant configurations.
  input_file <- arguments_file %>%
    read_csv(col_types = cols(),locale = readr::locale(encoding = "latin1")) %>%
    filter(!(argument %>% is.na & value1 %>% is.na))

  # Extract argument types.
  arg_type <- input_file %>%
    dplyr::select(argument) %>%
    unlist(use.names = F)

  # Link argument to argument type in a named list data structure.
  pipeline_input <- input_file %>%
    dplyr::select(value1) %>%
    unlist(use.names = F) %>%
    as.list() %>%
    rlang::set_names(arg_type)


  return(pipeline_input)
}
