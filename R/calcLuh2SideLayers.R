#' @title calcLuh2SideLayers
#' @description Function extracts biodiversity data for LUH2 land cover types
#'
#' @return magpie object in cellular resolution
#' @author Michael Windisch
#'
#' @examples
#' \dontrun{ calcOutput("Luh2SideLayers", aggregate = FALSE) }
#'
#' @importFrom magpiesets findset
#'

calcLuh2SideLayers <-function(){

  x      <- readSource("BendingTheCurve", subtype = "luh2_side_layers", convert="onlycorrect")
  weight <- dimSums(calcOutput("LanduseInitialisation", aggregate=FALSE, cellular=TRUE, land="fao", input_magpie=TRUE, years="y1995", round=6), dim=3)

return(list(
  x=x,
  weight=weight,
  unit="boolean",
  description="Data from LUH2 provided by David Leclere from IIASA, Bending the curve on biodiversity loss",
  isocountries=FALSE))
}
