#' @title calcTransportDistance
#' @description Function extracts travel time to major cities in minutes
#'
#' @return magpie object in cellular resolution
#' @author David Chen
#'
#' @examples
#' \dontrun{ calcOutput("TransportDistance", aggregate = FALSE) }
#'
#' @importFrom magpiesets findset
#'

calcTransportDistance <-function(){

  x      <- readSource("TransportDistance", convert="onlycorrect")
  weight <- dimSums(calcOutput("LanduseInitialisation", aggregate=FALSE, cellular=TRUE, land="fao", input_magpie=TRUE, years="y1995", round=6), dim=3)

  return(list(
    x=x,
    weight=weight,
    unit="Travel Time (minutes)",
    description="Travel time to major cities Nelson 2008 EC JRC, see model documentation",
    isocountries=FALSE))
}
