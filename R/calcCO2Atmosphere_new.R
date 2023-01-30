#' @title calcCO2Atmosphere_new
#' @description Disaggregate CO2 global atmospheric concentration to cellular level
#' @param subtype specify the version and scenario eg. "ISIMIP3bv2:ssp126"
#' @param co2_evolution Define `rising` for rising CO2 according to the climate scenario selected or
#' `static` for stable CO2 at the last past time step level.
#' @return magpie object in cellular resolution
#' @author Marcos Alves, Kristine Karstens
#'
#' @examples
#' \dontrun{
#' calcOutput("CO2Atmosphere_new", aggregate = FALSE, subtype, co2_evolution)
#' }
#'
#' @import madrat
#' @import magclass
#' @importFrom utils tail
#' @importFrom magpiesets findset
calcCO2Atmosphere_new <- function(subtype = "ISIMIP3bv2:ssp126", # nolint: object_name_linter.
                                  co2_evolution = "rising") { # nolint: object_name_linter.

  x <- readSource("CO2Atmosphere_new", subtype = subtype, convert = FALSE)

  if (co2_evolution == "rising") {

    cells <- toolGetMapping("CountryToCellMapping.csv", type = "cell")
    cells$glo <- "GLO"
    x <- toolAggregate(x, rel = cells, from = "glo", to = "celliso")

  } else if (co2_evolution == "static") {
    past <- tail(findset("past"), 1)
    fYear <- match(past, getYears(x))
    for (i in fYear:length(getYears(x))) {
      x[, i, ] <- x[, past, ]
    }
    cells <- toolGetMapping("CountryToCellMapping.csv", type = "cell")
    cells$glo <- "GLO"
    x <- toolAggregate(x, rel = cells, from = "glo", to = "celliso")
  }

  return(list(
    x = x,
    weight = NULL,
    unit = "ppm",
    description = "Atmosphere CO2 concentration",
    isocountries = FALSE
  ))
}
