#' @title calcGCMClimate_new
#' @description Disaggregate CO2 global atmospheric concentration to cellular level
#' @param subtype type of climate data to collect,
#'                consisting of data source, GDM, RCP, time period, variable and time resolution
#'                separated by ":"
#' @param smooth  set averaging value for smoothing trajectories
#' @return magpie object in cellular resolution
#' @author Marcos Alves, Kristine Karstens, Felicitas Beier
#'
#' @examples
#' \dontrun{
#' calcOutput("GCMClimate_new", subtype = "ISIMIP3bv2:IPSL-CM6A-LR:ssp126:1850-2100:tas:annual_mean")
#' }
#'
#' @importFrom madrat toolSplitSubtype toolTimeAverage
#' @importFrom magclass getNames
#' @importFrom magpiesets findset
#' @importFrom mstools toolHoldConstant
#'

calcGCMClimate_new <- function(subtype = "ISIMIP3bv2:IPSL-CM6A-LR:ssp126:1850-2100:tas:annual_mean", smooth = 0) {
  ###### CONFIG ######
  splitting_year <- 2014
  hist_name      <- "historical"
  time <- findset("time")
  time <- time[1:match("y2100", time)]
  ###### CONFIG ######

  x <- toolSplitSubtype(subtype,
                        list(version  = NULL, climatemodel = NULL,
                             scenario = NULL, period       = NULL,
                             variable = NULL, timeres      = NULL))

  period <- unlist(strsplit(x$period, "-|_"))

  .subtype_hist <- paste(x$version, x$climatemodel, hist_name, paste(period[1], splitting_year, sep = "-"), x$variable, sep = ":")
  .subtype_scen <- paste(x$version, x$climatemodel, x$scenario, paste((splitting_year + 1), period[2], sep = "-"), x$variable, sep = ":")

  y <- toolCoord2Isocell(
         mbind(readSource("GCMClimate_new", subtype = .subtype_hist,
                          subset = x$timeres, convert = "onlycorrect"),
               readSource("GCMClimate_new", subtype = .subtype_scen,
                          subset = x$timeres, convert = "onlycorrect")
               )
         )

  getNames(y) <- gsub("-", "_", paste(x$timeres, x$variable, x$version, x$climatemodel, x$scenario, sep = "_"))

  if (smooth > 1) {
    y <- toolTimeAverage(y, averaging_range = smooth)
    y <- toolHoldConstant(y, time)
    y <- y[, time, ]
  } else {
    y <- y[, time, ]
  }

  unit <- switch(x$variable,
    "tas"       = "Degree Celcius",
    "pr"        = "mm/day",
    "lwnet"     = "watt per m2",
    "rsds"      = "watt per m2",
    "wet"       = "number of rainy days"
  )

  description <- switch(x$variable,
    "tas"       = paste0("Average ", x$timeres, " air temperature"),
    "pr"        = "Average precipitation",
    "lwnet"     = "Long wave radiation",
    "rsds"      = "Short wave radiation",
    "wet"       = "number of rainy days"
  )

  return(list(x            = y,
              weight       = NULL,
              unit         = unit,
              description  = description,
              isocountries = FALSE))
}
