#' @title correctSoilClassification
#' @description Correct soil classification
#' @param x Magpie object provided by the read function
#' @return List of magpie objects with results on cellular level, weight, unit and description.
#' @author Marcos Alves, Kristine Karstens
#' @seealso
#' \code{\link{readSoilClassification}},
#' @examples
#'
#' \dontrun{
#'   readSource("SoilClassification", subtype="HWSD.soil", convert="onlycorrect")
#' }
#'
#' @import magclass
#' @importFrom madrat toolConditionalReplace
#' @importFrom mrcommons toolCell2isoCell

correctSoilClassification <- function(x){

  x <- toolConditionalReplace(x, conditions = c("<1",">13", "!is.integer()"), replaceby = 0)
  x <- toolCell2isoCell(x)

  return(x)
}
