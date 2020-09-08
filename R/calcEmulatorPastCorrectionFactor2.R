#' @title calcEmuPastCorrectFactor2
#' @description Loads the ANN trained grass emulator with historical data and calculates the pasture correction factor
#' @return MAgPIE objects with pasture correction factor on a cellular level.
#' @param model grass harvest machine learning model ID
#' @param GCMModel Global climate model name
#' @param rcp RCP
#' @author Marcos Alves
#' @examples
#'
#' \dontrun{
#' calcOutput("EmuPastCorrectFactor2", model = "f41f19be67", GCMModel = "HadGEM2", rcp = "rcp85")
#' }
#'
#' @import madrat
#' @importFrom tidyr pivot_wider
#' @import magpiesets
#' @import magclass
#' @import dplyr
#' @export
#'

# #development
# library(mrcommons)
# library(stringi)
# library(tidyr)
# library(magpiesets)
# library(mrmagpie)
# setConfig(forcecache=T)
# setConfig(globalenv = T)
# setwd("C:/magpie_inputdata/sources")

calcEmuPastCorrectFactor2 <-
  function(model = "f41f19be67", GCMModel = "HadGEM2", rcp = "rcp85") {

    past <- findset("past")
    past <- past[7:length(past)]


    #############################
    ###    FAO pasture area   ###
    #############################

    MAGPasturearea <- calcOutput("LanduseInitialisation", cellular = T, aggregate = FALSE)[, past, "past"]

    #################################
    ### Loading emulator data      ##
    #################################

    max_grass              <- toolCell2isoCell(readSource("GrassYldEmu", subtype = paste(model,"max_harvest", sep = "."), convert="onlycorrect"))
    getYears(max_grass)    <- past
    max_grass[max_grass<0] <- 0

    #################################
    ### converting max harvest     ##
    #################################

    max_grass[max_grass < 0]        <- 0
    # max_grass                       <- round(max_grass, 2)
    gCm2_To_tDM                     <- (10000 * 2.21 / 1e6)
    max_grass_production            <- max_grass * gCm2_To_tDM
    max_grass_production            <- max_grass_production * (MAGPasturearea * 1e6)
    # max_grass_productivity          <- round(pred_pasture_production, 2)

    #############################
    ### Disaggregation weights###
    #############################

    potential_grass_cell         <- MAGPasturearea * max_grass
    mapping                      <- toolGetMapping(name = "CountryToCellMapping.csv", type = "cell")
    potential_grass_country      <- toolAggregate(potential_grass_cell, rel = mapping, from = "celliso", to = "iso", partrel = F)
    potential_grass_country_cell <- toolAggregate(potential_grass_country, rel = mapping, from = "iso", to = "celliso", partrel = F)
    grass_distribution           <- potential_grass_cell / potential_grass_country_cell
    grass_distribution[is.na(grass_distribution)] <- 0

    #################################
    ### FAO past. demand cellular ###
    #################################
    mapping                      <- toolGetMapping(name = "CountryToCellMapping.csv", type = "cell")
    FAO_pasture_demand                           <- calcOutput("FAOmassbalance", aggregate = F)[, , "dm"][, , "feed"][, past, "pasture"]
    FAO_pasture_demand                           <- FAO_pasture_demand[unique(mapping$iso)] * 1e6
    FAO_pasture_demand_cell                      <- toolAggregate(FAO_pasture_demand, rel = mapping, from = "iso", to = "celliso", partrel = F, weight = grass_distribution)
    FAO_pasture_demand_cell[MAGPasturearea == 0] <- 0 # Setting to zero pasture demand in countries with no pasture area

    #################################
    ### Correction factor         ###
    #################################

    correction_factor                                 <- FAO_pasture_demand_cell / max_grass_production
    correction_factor[is.nan(correction_factor)]      <- 0
    correction_factor[is.infinite(correction_factor)] <- 0
    # correction_factor[correction_factor>10] <- 10

    #Analysis
    # sum(max_grass_production) / sum(FAO_pasture_demand_cell)
    # summary(correction_factor[, 1, ])
    # q <- quantile(correction_factor, probs = 0.95);q
    # luplot::plotmap2(correction_factor[, 1, ], legend_range = c(0,q), lowcol = "#FFFFFF", midpoint = q/2)

    return(list(
      x = correction_factor,
      weight = NULL,
      unit = "",
      description = "Pasture correction factor for the historical dates",
      isocountries = FALSE
    ))
  }
