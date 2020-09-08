#' @title correctGPD
#' @description correct GPD 
#' @param x magpie object provided by the read function
#' @return List of magpie objects with results on cellular level, weight, unit and description.
#' @author Florian Humpenoeder
#' @examples
#'
#' \dontrun{
#'   readSource("GPD", convert="onlycorrect")
#' }
#' @importFrom readxl read_xls


correctGPD <- function(x){

  #convert to Mha
  x <- x/1000
  #replace country names with iso3 codes
  getCells(x) <- toolCountry2isocode(getCells(x))

  x <- toolCountryFill(x)
  CountryToCell   <- toolGetMapping("CountryToCellMapping.csv", type = "cell")
  x <- x[unique(CountryToCell$iso),,]
  x[is.na(x)] <- 0
  
  return(x)
}
