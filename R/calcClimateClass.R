#' @title calcClimateClass
#' @description fraction of a cell belonging to a given climate classification based on Koeppen Geiger Classification. http://koeppen-geiger.vu-wien.ac.at/.
#'
#' @return Clustered MAgPIE object on requested resolution
#' @author Abhijeet Mishra
#'
#' @examples
#' \dontrun{ calcOutput("ClimateClass", aggregate = FALSE) }
#'
#' @export

calcClimateClass <-function(){

  x <- readSource("Koeppen", subtype="cellular",convert = FALSE)

  weight <- calcOutput("LanduseInitialisation", aggregate=FALSE, cellular=TRUE, land="fao", input_magpie=TRUE, years="y1995", round=6)

  return(list(
    x=x,
    weight=weight,
    unit="share",
    description="share of koeppen geiger area",
    isocountries=FALSE))
}
