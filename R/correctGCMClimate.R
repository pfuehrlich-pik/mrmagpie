#' @title correctGCMClimate
#' @description Correct GCMs climate variables
#' @param x magpie object provided by the read function
#' @return Magpie objects with results on cellular level, weight, unit and description.
#' @author Marcos Alves
#' @seealso
#' \code{\link{readGCMClimate}}
#' @examples
#'
#' \dontrun{
#' readSource("GCMClimate", subtype="rcp85:HadGEM2.temperature", convert="onlycorrect")
#' }
#'
#' @import magclass
#' @importFrom madrat toolConditionalReplace
#' @importFrom mrcommons toolCell2isoCell

correctGCMClimate <- function(x){

  x <- toolConditionalReplace(x, conditions = c("is.na()"), replaceby = 0)
  x <- toolCell2isoCell(x)

  return(x)
}
