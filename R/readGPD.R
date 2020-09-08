#' @title readGPD
#' @description read GPD
#' Data from the Global Peatland Database provided by Alexandra Barthelmes.
#' The original xls file has been clean-up manually (country names). Turkey had two identical entries in the original xls file.
#' Sources:
#' "Inventory Reports and National Communications UNFCC 2014",
#' "soil and peatland science",
#' "European Mires Book" ,
#' "own estimates (incl. GIS data)",
#' @return List of magpie objects with results on cellular level, weight, unit and description.
#' @author Florian Humpenoeder
#' @examples
#'
#' \dontrun{
#'   readSource("GPD", convert="onlycorrect")
#' }
#' @importFrom readxl read_xls


readGPD <- function(){

  #read-in xls file Global Peatland Database
  x <- read_xls("Global_GPD_15_AB_1811020.xls")
  #convert to magclass object
  x <- as.magpie(x[,1:3],spatial=1)
  
  return(x)
}
