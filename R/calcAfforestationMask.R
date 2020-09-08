#' @title calcAfforestationMask
#' @description Afforestation mask for where afforestation possible
#' @param subtype afforestation mask sub type
#' @return magpie object in cellular resolution
#' @author David Chen
#'
#' @examples
#' \dontrun{ calcOutput("AfforestationMask", subtype="noboreal", aggregate = FALSE) }
#'
#' @importFrom magpiesets findset
#' @importFrom madrat readSource

calcAfforestationMask <-function(subtype){

  x      <- readSource("AfforestationMask",subtype=subtype, convert="onlycorrect")
  weight <- dimSums(calcOutput("LanduseInitialisation", aggregate=FALSE, cellular=TRUE, land="fao", input_magpie=TRUE, years="y1995", round=6), dim=3)

  return(list(
    x=x,
    weight=weight,
    unit="binary",
    description="",
    isocountries=FALSE))
}
