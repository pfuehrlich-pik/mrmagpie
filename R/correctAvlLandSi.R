#' @title correctAvlLandSi
#' @description Read Available Land Si
#' @return List of magpie objects with results on cellular level, weight, unit and description.
#' @param x magpie object provided by the read function
#' @author David Chen
#' @seealso
#'   \code{\link{readAvlLandSi}}
#' @examples
#'
#' \dontrun{
#'   readSource("AvlLandSi", convert="onlycorrect")
#' }
#'
#' @importFrom madrat toolConditionalReplace
#' @importFrom mrcommons toolCell2isoCell

correctAvlLandSi <- function(x){

  x <- toolConditionalReplace(x, conditions = c("is.na()","<0"), replaceby = 0)
  x <- toolCell2isoCell(x)

  return(x)
}
