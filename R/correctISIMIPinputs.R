#' @title correctISIMIPinputs
#' @description Correct ISIMIP data
#' @return List of magpie objects with results on cellular level, weight, unit and description.
#' @param x MAgPIE object provided by the read function
#' @author Felicitas Beier
#' @seealso
#'   \code{\link{readISIMIPinputs}}
#' @examples
#' \dontrun{
#'   readSource("ISIMIPinputs", convert="onlycorrect")
#' }
#'
#' @import magclass
#' @importFrom madrat toolConditionalReplace
#' @importFrom mrcommons toolCell2isoCell

correctISIMIPinputs <- function(x){

  x <- toolConditionalReplace(x, conditions = c("is.na()","<0"), replaceby = 0)
  x <- toolCell2isoCell(x)

  return(x)
}
