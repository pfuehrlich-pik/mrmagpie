#' @title readBendingTheCurve
#' @description Read bending the curve data
#' @param x magpie object provided by the read function
#' @return List of magpie objects with results on cellular level, weight, unit and description.
#' @author Patrick v. Jeetze, Michael Windisch
#' @examples
#'
#' \dontrun{
#'   readSource("BendingTheCurve", subtype="rr_layer", convert="onlycorrect")
#' }
#'
#' @importFrom madrat toolConditionalReplace
#' @importFrom mrcommons toolCell2isoCell


correctBendingTheCurve <- function(x){

  x <- toolConditionalReplace(x, conditions = c("is.na()","<0"), replaceby = 0)
  x <- toolCell2isoCell(x)

  return(x)
}
