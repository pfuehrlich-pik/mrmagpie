#' @title correctLeifeld2018
#' @description correct potential peatland area from Leifeld2018
#' @param x magpie object provided by the read function
#' @return List of magpie objects with results on cellular level, weight, unit and description.
#' @author Florian Humpenoeder
#' @examples
#'
#' \dontrun{
#'   readSource("Leifeld2018", convert="onlycorrect")
#' }
#' @importFrom mrcommons toolCell2isoCell
#' @importFrom magclass getYears getNames


correctLeifeld2018 <- function(x){
  
  x <- toolCell2isoCell(x)
  getYears(x) <- NULL
  getNames(x) <- NULL
  x[is.na(x)] <- 0
  
  return(x)
}
