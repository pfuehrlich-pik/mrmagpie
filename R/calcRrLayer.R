#' @title calcRRLayer
#' @description Function extracts range-rarity as used for biodiversity loss
#'
#' @return magpie object in cellular resolution
#' @author Michael Windisch, Patrick v. Jeetze
#'
#' @examples
#' \dontrun{ calcOutput("RRLayer", aggregate = FALSE) }
#'
#' @importFrom magpiesets findset
#'

calcRRLayer <-function(){

x <- readSource("BendingTheCurve", subtype = "rr_layer", convert="onlycorrect")

return(list(
  x=x,
  weight=NULL,
  unit="Range-Rarity (-)",
  description="range-rarity layer provided by David Leclere from IIASA, Bending the curve on biodiversity loss",
  isocountries=FALSE))
}
