#' @title calcPastureSuit
#' @description Calculate glassland suitable for pasture management based on population and aridity criteria.
#' @param subtype Select version, climate model and period.
#' @param smoothPrecipitation Smooth precipitation climate data over time
#' @param smoothOut smooth the Pasture suitability areas variations over time
#' @return List of magpie object with results on cluster level
#' @author Marcos Alves
#' @examples
#' \dontrun{
#' calcOutput("PastureSuit")
#' }
#' @importFrom raster area rasterFromXYZ

calcPastureSuit <- function(subtype = "ISIMIP3bv2:IPSL-CM6A-LR:1850-2100", smoothPrecipitation = 10, smoothOut = 10) {
  x <- toolSplitSubtype(subtype, list(version = NULL, climatemodel = NULL, period = NULL))

  # pasture drivers
  population <- calcOutput("GridPop", subtype = "all", cellular = TRUE, FiveYear = TRUE,
                           harmonize_until = 2015, aggregate = FALSE)

  precipitation <- list()
  scenarios <- c("ssp126", "ssp245", "ssp370", "ssp460", "ssp585") # Current ISIMIP3bv2 scenarios
  for (scenario in scenarios) {
    subtype <- paste(x$version, x$climatemodel, scenario, x$period, "pr", "annual_mean", sep = ":")
    precipitation[[scenario]] <- setNames(calcOutput("GCMClimate_new", subtype = subtype,
                                                     smooth = smoothPrecipitation, aggregate = FALSE), scenario)
  }
  precipitation <- collapseNames(mbind(precipitation))

  evapotranspiration <- calcOutput("Evapotranspiration", subtype = "H08:mri-esm2-0", aggregate = FALSE)

  # temporary mapping of evapotranspiration RCP scenarios unavailable in ISIMIP3bv2
  evapotranspiration <- add_columns(evapotranspiration, addnm = "ssp245", dim = 3.1, fill = NA)
  evapotranspiration[, , "ssp245"] <- evapotranspiration[, , "ssp370"]
  evapotranspiration <- add_columns(evapotranspiration, addnm = "ssp460", dim = 3.1, fill = NA)
  evapotranspiration[, , "ssp460"] <- evapotranspiration[, , "ssp370"]

  evapotranspiration <- evapotranspiration[, , getItems(precipitation, dim = 3)]

  # matching available ssps scenarios
  regex <- paste0("[", paste0("+", strtrim(getItems(evapotranspiration, dim = 3), 4), collapse = "|"), "]", "{4}$")
  avlSSPs <- grep(regex, getNames(population), ignore.case = TRUE)
  population <- population[, , avlSSPs]

  # Cell area calculation
  landcoords <- as.data.frame(toolGetMapping("magpie_coord.rda", type = "cell"))
  landcoords <- cbind(landcoords, rep(1, nrow(landcoords)))
  landcoords <- raster::rasterFromXYZ(landcoords)
  crs(landcoords) <- "+proj=longlat" # outputs cell are in km2
  cellSize <- raster::area(landcoords)
  cellSize <- cellSize * landcoords
  cellSize <- as.magpie(cellSize)
  cellSize <- toolOrderCells(collapseDim(addLocation(cellSize), dim = c("x", "y")))

  # population density
  popDensity <- (population * 1e6) / cellSize # population density in number of people per km2
  popDensity[is.infinite(popDensity)] <- 0
  popDensity[is.nan(popDensity)] <- 0


  yearsCom <- intersect(getYears(popDensity), getYears(precipitation))

  # Aridity (the real aridity is measured as the ratio between evapotranspiration
  # and precipitarion (I have complete this calculation))
  aridity <- precipitation[, yearsCom, ] / (evapotranspiration[, yearsCom, ])
  aridity[is.infinite(aridity) | is.nan(aridity)] <- 0
  # 0.5 aridity threshold for managed pastures. Same from HYDE 3.2.
  aridity[aridity < 0.5] <- 0
  aridity[aridity >= 0.5] <- 1

  # pasture suitability check
  pastureSuit <- aridity
  popDensity <- popDensity[, getYears(pastureSuit), ]
  pastureSuit[popDensity < 5] <- 0 # 5 hab km2 population threshold for managed pastures. Same from HYDE 3.2.

  pastureSuitArea <- (pastureSuit * cellSize * 100) / 1e6 # (from km2 (x100) to mha (/1e6))
  pastureSuitArea <- collapseDim(pastureSuitArea)
  pastureSuitArea <- toolHoldConstantBeyondEnd(pastureSuitArea)

  # calibration to historical values

  histPastr <- setNames(calcOutput("LanduseInitialisation", cellular = TRUE,
                                   nclasses = "nine", aggregate = FALSE)[, , c("past")], "pastr")
  pastAll <- intersect(getYears(histPastr), getYears(pastureSuitArea))

  pastLy <- findset("past")
  pastLy <- pastLy[length(pastLy)] # past last year
  future <- setdiff(getYears(pastureSuitArea), pastAll)

  map <- toolGetMapping("CountryToCellMapping.csv", type = "cell")
  pastureSuitAreaReg <- toolAggregate(pastureSuitArea, rel = map, from = "celliso", to = "iso")
  histPastrReg <- toolAggregate(histPastr, rel = map, from = "celliso", to = "iso")
  corrReg <- histPastrReg[, pastLy, ] / pastureSuitAreaReg[, pastLy, ]
  pastureSuitArea[, future, ] <- toolAggregate(corrReg, rel = map, from = "iso", to = "celliso") *
                                 pastureSuitArea[, future, ]

  pastureSuitArea[is.infinite(pastureSuitArea) | is.nan(pastureSuitArea) | is.na(pastureSuitArea)] <- 0
  pastureSuitArea[pastureSuitArea < 0] <- 0
  pastureSuitArea[, pastAll, ] <- histPastr[, pastAll, ]

  if (smoothOut > 1) {
    pastureSuitArea <- toolTimeAverage(pastureSuitArea, averaging_range = smoothOut)
  }
  y <- intersect(pastAll, getYears(pastureSuitArea))
  pastureSuitArea[, y, ] <- histPastr[, intersect(pastAll, getYears(pastureSuitArea)), ]
  pastureSuitArea <- toolHoldConstant(pastureSuitArea, findset("time"))
  pastureSuitArea <- collapseNames(pastureSuitArea)
  pastureSuitArea[, pastAll, ] <- histPastr[, pastAll, ]

  return(list(
    x = pastureSuitArea,
    weight = NULL,
    unit = "Mha",
    description = "Area suitable for pasture management in mha",
    isocountries = FALSE
  ))
}
