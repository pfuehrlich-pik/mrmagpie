#' @title calcProtectArea
#' @description Function extracts conservation protected area
#'
#' @return magpie object in cellular resolution with different protection scenarios
#' @author Felicitas Beier, David Chen
#'
#' @examples
#' \dontrun{ calcOutput("ProtectArea", aggregate = FALSE) }
#'
#' @importFrom magpiesets findset
#'

calcProtectArea <- function(){

  # Protection Area mz file (conservation priority area in Mha)
  x <- readSource("ProtectArea", convert="onlycorrect")

  # Half Earth Protection Share
  protect_share           <- readSource("HalfEarth", convert="onlycorrect")
  getNames(protect_share) <- "HalfEarth"
  # Land area (in Mha):
  magpie_land_area <- calcOutput("LanduseInitialisation", aggregate=FALSE, cellular=TRUE, land="fao", input_magpie=TRUE, years="y1995")
  magpie_land_area <- dimSums(magpie_land_area, dim=3)
  # Land area to be protected by 2050 (in Mha)
  protect_area     <- protect_share * magpie_land_area

  # Add HalfEarth scenario to Protection Area file
  x <- mbind(x, protect_area)

  return(list(
    x=x,
    weight=NULL,
    unit="Mha",
    description="conservation priority areas",
    isocountries=FALSE))
}
