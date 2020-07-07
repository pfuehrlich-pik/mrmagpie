#' @title correctCO2Atmosphere
#' @description Correct atmosphere co2 concentration
#' @param x magpie object provided by the read function
#' @return List of magpie objects with results on cellular level, weight, unit and description.
#' @author Marcos Alves, Kristine Karstens
#' @seealso
#' \code{\link{readCO2Atmosphere}}
#' @examples
#'
#' \dontrun{
#' readSource("CO2Atmosphere", subtype="rcp85", convert="onlycorrect")
#' }
#'
#' @import magclass
#' @importFrom madrat toolConditionalReplace
#' @importFrom mrcommons toolCell2isoCell

correctCO2Atmosphere <- function(x){

  x <- toolConditionalReplace(x, conditions = c("is.na()"), replaceby = 0)

  return(x)
}
