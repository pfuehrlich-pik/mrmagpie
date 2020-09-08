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

  x      <- readSource("BendingTheCurve", subtype = "rr_layer", convert="onlycorrect")
  weight <- dimSums(calcOutput("LanduseInitialisation", aggregate=FALSE, cellular=TRUE, land="fao", input_magpie=TRUE, years="y1995", round=6), dim=3)

return(list(
  x=x,
  weight=weight,
  unit="Range-Rarity (-)",
  description="range-rarity layer provided by David Leclere from IIASA, Bending the curve on biodiversity loss",
  isocountries=FALSE))
}
