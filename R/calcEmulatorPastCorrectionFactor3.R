#' @title calcEmuPastCorrectFactor
#' @description Loads the ANN trained grass emulator with historical data and calculates the pasture correction factor
#' @return MAgPIE objects with pasture correction factor on a cellular level.
#' @param model grass harvest machine learning model ID
#' @param GCMModel Global climate model name
#' @param rcp RCP
#' @author Marcos Alves
#' @examples
#'
#' \dontrun{
#' calcOutput("EmuPastCorrectFactor", model = "f41f19be67", GCMModel = "HadGEM2", rcp = "rcp85")
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

calcEmuPastCorrectFactor3 <-
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

    inputs_vec             <- as.vector(readSource("GrassYldEmu", subtype = paste(model,"inputs", sep = "."), convert="onlycorrect"))
    max_grass              <- toolCell2isoCell(readSource("GrassYldEmu", subtype = paste(model,"max_harvest", sep = "."), convert="onlycorrect"))
    getYears(max_grass)    <- past
    max_grass[max_grass<0] <- 0
    magpie_center          <- readSource("GrassYldEmu", subtype =paste(model,"mean", sep = "."), convert="onlycorrect")
    center                 <- as.numeric(magpie_center)
    names(center)          <- getCells(magpie_center)
    magpie_scale           <- readSource("GrassYldEmu", subtype = paste(model,"stddevs", sep = "."), convert="onlycorrect")
    scale                  <- as.numeric(magpie_scale)
    names(scale)           <- getCells(magpie_scale)
    weights                <- readSource("GrassYldEmu", subtype = paste(model,"weights", sep = "."), convert="onlycorrect")

    #################################
    ### Loading environmental data ##
    #################################

    environmental_data_1 <- calcOutput("CellEnvironmentData", model = model, GCMModel = GCMModel, rcp = rcp, years = past, aggregate = F)
    environmental_data_2 <- as.data.frame(environmental_data_1, rev = 2)
    environmental_data_3 <- pivot_wider(environmental_data_2, values_from = ".value", names_from = "data")
    drops                <- grep("region+|year+", colnames(environmental_data_3))
    environmental_data   <- as.matrix(environmental_data_3[,-drops])

    #############################
    ### Disaggregation weights###
    #############################
#
#     potential_grass_cell         <- MAGPasturearea * max_grass
     mapping                      <- toolGetMapping(name = "CountryToCellMapping.csv", type = "cell")
#     potential_grass_country      <- toolAggregate(potential_grass_cell, rel = mapping, from = "celliso", to = "iso", partrel = F)
#     potential_grass_country_cell <- toolAggregate(potential_grass_country, rel = mapping, from = "iso", to = "celliso", partrel = F)
#     grass_distribution           <- potential_grass_cell / potential_grass_country_cell
#     grass_distribution[is.na(grass_distribution)] <- 0

    GLW3 <- readSource("GLW3", subtype = "Da", convert="onlycorrect")

    #############################
    ### FAO livestock Numbers ###
    #############################

    livestock_FAO <- readSource("FAO", subtype = "LiveHead")[, past, c("1746|Cattle and Buffaloes + (Total).stock", "1749|Sheep and Goats + (Total).stock")]

    #################################
    ###        Feedbaskets        ###
    #################################

    fbask                  <- calcOutput("FeedBasketsPast", aggregate = FALSE)
    fbask_rum_pasture      <- dimSums(fbask[, past, c("alias_livst_rum", "alias_livst_milk")][, , "pasture"])
    fbask_rum_total        <- dimSums(fbask[, past, c("alias_livst_rum", "alias_livst_milk")])
    fbask_pasture_fraction <- fbask_rum_pasture / fbask_rum_total

    #################################
    ### FAO past. demand cellular ###
    #################################

    FAO_pasture_demand                           <- calcOutput("FAOmassbalance", aggregate = F)[, , "dm"][, , "feed"][, past, "pasture"]
    FAO_pasture_demand                           <- FAO_pasture_demand[unique(mapping$iso)] * 1e6
    FAO_pasture_demand_cell                      <- toolAggregate(FAO_pasture_demand, rel = mapping, from = "iso", to = "celliso", partrel = F, weight = GLW3)
    FAO_pasture_demand_cell[MAGPasturearea == 0] <- 0 # Setting to zero pasture demand in countries with no pasture area

    FAO_LSU_demand                               <- FAO_pasture_demand/(8.9/1000*365)
    FAO_LSU_demand_cell                          <- toolAggregate(FAO_LSU_demand, rel = mapping, from = "iso", to = "celliso", partrel = F, weight = GLW3)
    FAO_LSU_demand_cell[MAGPasturearea == 0] <- 0 # Setting to zero pasture demand in countries with no pasture area

    #################################
    ### Calculating LSU/ha        ###
    #################################

    # conversion_rate_LSU            <- c(0.894, 0.1) # source calulation from EU statistics https://docs.google.com/spreadsheets/d/1SZAAVl1SLwrrK6j329tq5zo1VZhfFtxUHTmsGQKYCCk/edit#gid=0
    # conversion_rate_LSU            <- as.magpie(conversion_rate_LSU)
    # dimnames(conversion_rate_LSU)  <- list("region", "year", c("large", "small"))
    #
    # livestock_FAO_scaled           <- livestock_FAO * conversion_rate_LSU * fbask_pasture_fraction
    # livestock_FAO_scaled           <- dimSums(livestock_FAO_scaled[, , c(1, 4)])
    # getYears(livestock_FAO_scaled) <- past
    # livestock_FAO_scaled           <- livestock_FAO_scaled[unique(mapping$iso)]
    livestock_cell                   <- toolAggregate(FAO_LSU_demand, rel = mapping, from = "iso", to = "celliso", partrel = F, weight = GLW3)

    livestock_cell[MAGPasturearea == 0] <- 0
    lsu_ha                              <- livestock_cell / (MAGPasturearea * 1e6)
    lsu_ha[is.nan(lsu_ha)]              <- 0
    lsu_ha[lsu_ha > 2.5]                <- 2.5
    lsu_ha_scaled                       <- scale(lsu_ha, center = center["lsu"], scale = scale["lsu"])


    #################################
    ### Predicted grass yield     ###
    #################################

    prediction                            <- toolNeuralNet(cbind(environmental_data, lsu_ha_scaled),weights)
    prediction[prediction < 0]            <- 0
    prediction                            <- round(prediction, 2)
    gCm2_To_tDM                           <- (10000 * 2.21 / 1e6)
    pred_pasture_productivity             <- prediction * gCm2_To_tDM
    pred_pasture_productivity             <- cbind(environmental_data_3[, c(1, 2, 3)], pred_pasture_productivity)
    pred_pasture_productivity_mag         <- pred_pasture_productivity
    pred_pasture_productivity_mag$celliso <- paste(pred_pasture_productivity_mag$region, pred_pasture_productivity_mag$region1, sep = ".")
    pred_pasture_productivity_mag$year    <- paste0("y", pred_pasture_productivity_mag$year)
    pred_pasture_productivity_mag         <- pred_pasture_productivity_mag[,-c(1,2)]
    pred_pasture_productivity_mag         <- as.magpie(pred_pasture_productivity_mag[, c(3, 1, 2)])
    pred_pasture_production               <- pred_pasture_productivity_mag * (MAGPasturearea * 1e6)
    pred_pasture_production               <- round(pred_pasture_production, 2)

    #################################
    ### Correction factor         ###
    #################################

    correction_factor                                 <- FAO_pasture_demand_cell / pred_pasture_production
    correction_factor[is.nan(correction_factor)]      <- 0
    correction_factor[is.infinite(correction_factor)] <- 0
    # correction_factor[correction_factor>10] <- 10


    #Analysis
    # sum(pred_pasture_production) / sum(FAO_pasture_demand_cell)
    # summary(correction_factor[, 1, ])
    # q <- quantile(correction_factor, probs = 0.95);q
    # luplot::plotmap2(correction_factor[, 1, ], legend_range = c(0,q), lowcol = "#FFFFFF", midpoint = 4)

    return(list(
      x = correction_factor,
      weight = NULL,
      unit = "",
      description = "Pasture correction factor for the historical dates",
      isocountries = FALSE
    ))
  }
