#' @title correctIrrigationSystem
#' @description Correct Irrigation System data
#' @return List of magpie objects with results on cellular level, weight, unit and description.
#' @param x MAgPIE object provided by the read function
#' @author Felicitas Beier
#' @seealso
#'   \code{\link{readIrrigationSystem}}
#' @examples
#' \dontrun{
#'   readSource("IrrigationSystem", convert="onlycorrect")
#' }
#'
#' @importFrom mrcommons toolCell2isoCell

correctIrrigationSystem <- function(x){
  return(x)
}
