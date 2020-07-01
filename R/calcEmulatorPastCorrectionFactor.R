#' @title EmulatorPastCorrectionFactor
#' @description Loads the ANN trained grass emulator with historical data and calculates the pasture correction factor
#' @param model ANN model ID
#' @return MAgPIE objects with results on cellular level.
#' @author Marcos Alves
#' @seealso
#' \code{\link{readGCMClimate}}
#' @examples
#'
#' \dontrun{
#' readSource("GCMClimate", subtype = "HadGEM2:rcp85.temperature", convert = "onlycorrect")
#' }
#'
#' @importFrom lpjclass read.LPJ_input
#' @importFrom madrat toolSubtypeSelect
#' @importFrom lucode2 path
#' @importFrom tidyr select pivot_wider
#' @import keras
#' @import magpiesets
#' @import magclass
#' @export
#'

#development
library(mrcommons)
library(magclass)
library(keras)
library(lucode2)
library(tidyr)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(dplyr)
library(magpiesets)
library(mrmagpie)
setConfig(forcecache=TRUE)

calcEmulatorPastCorrectionFactor <-
  function(model_id = "20f33a2280", folder = "GrassYieldsEmulator") {
    if (nchar(model_id) == "10" & !grepl("[^0-9, a-z]", model_id)) {
      folder_list <- list.dirs(path(folder))
      folder <- path(folder_list[grep(model_id, folder_list)])
    } else {
      stop("'model_id' has an incorrect format. The 10 first characters that identify a model must have only letter and numbers")
    }

    past <- findset("past")
    past <- past[7:length(past)]
    past <- "y1995"
    folder <- "C:/magpie_inputdata/sources/GrassYieldsEmulator/20f33a2280"

    #############################
    ###        FAO data       ###
    #############################

    FAO_pasture_demand <- calcOutput("FAOmassbalance",aggregate = F)[,,"dm"][,,"feed"][,past,"pasture"]
    MAGPasturearea <- calcOutput("LanduseInitialisation", cellular = T, aggregate = FALSE)[,past,"past"]

    #################################
    ### Loading environmental data ##
    #################################

    files <- list.files(path(folder))
    inputs_vec <- readRDS(path(folder,files[grep("inputs", files)]))
    center <- readRDS(path(folder,files[grep("means", files)]))
    scale <- readRDS(path(folder,files[grep("stddevs", files)]))

    soil_char <- c("Ks", "Sf", "w_pwp", "w_fc", "w_sat", "hsg", "tdiff_0", "tdiff_15", "tdiff_100", "cond_pwp", "cond_100", "cond_100_ice")
    ml_inputs <- inputs_vec[is.na(match(inputs_vec, "lsu"))]
    if (any(soil_char %in% ml_inputs)) {
      inputs <- ml_inputs[is.na(match(ml_inputs, soil_char))]
      if (!("soil" %in% inputs)) {
        inputs <- append(inputs, "soil")
      }
    }

    list_mag_obj <- list()

    for (input in inputs) {
      if ("temperature" %in% input) {
        list_mag_obj[[input]] <- calcOutput("GCMClimate", aggregate = FALSE, GCMModel = "HadGEM2", ClimateVariable = "temperature", rcp = "rcp85")[, past, ]
      }

      if ("precipitation" %in% input) {
        list_mag_obj[[input]] <- calcOutput("GCMClimate", aggregate = FALSE, GCMModel = "HadGEM2", ClimateVariable = "precipitation", rcp = "rcp85")[, past, ]
      }

      if ("lwnet" %in% input) {
        list_mag_obj[[input]] <- calcOutput("GCMClimate", aggregate = FALSE, GCMModel = "HadGEM2", ClimateVariable = "longwave_radiation", rcp = "rcp85")[, past, ]
      }

      if ("rsds" %in% input) {
        list_mag_obj[[input]] <- calcOutput("GCMClimate", aggregate = FALSE, GCMModel = "HadGEM2", ClimateVariable = "shortwave_radiation", rcp = "rcp85")[, past, ]
      }

      if ("wetdays" %in% input) {
        list_mag_obj[[input]] <- calcOutput("GCMClimate", aggregate = FALSE, GCMModel = "HadGEM2", ClimateVariable = "wetdays", rcp = "rcp85")[, past, ]
      }

      if ("co2" %in% input) {
        list_mag_obj[[input]] <- calcOutput("CO2Atmosphere", aggregate = FALSE, rcp = "rcp85", level = "cellular")[, past, ]
      }

      if ("soil" %in% input) {
        list_mag_obj[[input]] <- calcOutput("SoilCharacteristics", aggregate = FALSE)[, past, ]
      }
    }

    environmental_data_1 <- list_mag_obj[[1]]
    for (i in 2:length(inputs)) {
      environmental_data_1 <- mbind(environmental_data_1, list_mag_obj[[i]])
    }

    environmental_data_2 <- as.data.frame(environmental_data_1, rev = 2)
    environmental_data_3 <- pivot_wider(environmental_data_2, values_from = ".value", names_from = "data")
    environmental_data_4 <- environmental_data_3[,-c(1,2,3)]
    environmental_data <- scale(environmental_data_4, center = center[ml_inputs] , scale = scale[ml_inputs]) #final version

    #############################
    ### Disaggregation weights###
    #############################

    max_grass <- read.magpie("C:/Users/pedrosa/github/Models/LPJmL-emulator/data/mag_max_grass_lpjml_cell.mz")[,past,]
    getCells(max_grass) <- getCells(environmental_data_1)
    potential_grass_cell <- MAGPasturearea * max_grass

    mapping <- toolGetMapping(name = "CountryToCellMapping.csv", type = "cell")
    potential_grass_country <- toolAggregate(potential_grass_cell, rel=mapping, from="celliso", to="iso", partrel=F)
    potential_grass_country_cell <- toolAggregate(potential_grass_country, rel=mapping, from="iso", to="celliso", partrel=F)
    grass_distribution <- potential_grass_cell/potential_grass_country_cell
    grass_distribution[is.na(grass_distribution)] <- 0

    # getCells(grass_distribution)[is.na(grass_distribution)]
    # grass_distribution[grass_distribution>0.01] <-0.01
    # luplot::plotmap2(grass_distribution[,1,])

    #############################
    ### Disaggregation        ###
    #############################

    # cb <- as.data.frame(mrmagpie:::magpie_coord)
    # cell_area  <- (111e3*0.5) * (111e3*0.5)*cos(cb$lat/180*pi)
    # cellarea <- as.magpie(cell_area, spatial = 1)
    # weight <- new.magpie(1:59199,past, fill = cellarea)
    # getCells(weight) <- getCells(environmental_data_1)

    # mapping <- toolGetMapping(name = "CountryToCellMapping.csv", type = "cell") #Find a country maping with 251 countries
    # environmental_data_2 = toolAggregate(environmental_data_1, rel=mapping, from="celliso", to="iso", partrel=FALSE, weight = weight)

    # environmental_data_3 <- as.data.frame(environmental_data_2, rev = 2)
    # environmental_data_4 <- pivot_wider(environmental_data_3, values_from = ".value", names_from = "data")[,-c(1,2)]
    # environmental_data <- scale(environmental_data_4, center = center[ml_inputs] , scale = scale[ml_inputs]) #final version

    # a = as.magpie(environmental_data, spatial =1)
    # getCells(a) <- getCells(w)
    # country2map(a[,,1], c("iso_a3" = "Region"))


    #############################
    ### FAO livestock Numbers ###
    #############################

    livestock_FAO <- readSource("FAO", subtype = "LiveHead")[, past,c("1746|Cattle and Buffaloes + (Total).stock","1749|Sheep and Goats + (Total).stock")]

    #################################
    ###        Feedbaskets        ###
    #################################

    fbask <- calcOutput("FeedBasketsPast", aggregate = FALSE)
    fbask_rum_pasture <- dimSums(fbask[,past,c("alias_livst_rum","alias_livst_milk")][,,"pasture"])
    fbask_rum_total <- dimSums(fbask[,past,c("alias_livst_rum","alias_livst_milk")])
    fbask_pasture_fraction <- fbask_rum_pasture/fbask_rum_total

    #calculate share of all pasture as a ratio of the sum of pasture and forrage? I think we dont need that.

    #################################
    ### Calculating LSU/ha        ###
    #################################
    #conversion_rate_LSU <- 0.894 #source calulation from EU statistics https://docs.google.com/spreadsheets/d/1SZAAVl1SLwrrK6j329tq5zo1VZhfFtxUHTmsGQKYCCk/edit#gid=0
    conversion_rate_LSU <- c(0.894,0.1) #source calulation from EU statistics https://docs.google.com/spreadsheets/d/1SZAAVl1SLwrrK6j329tq5zo1VZhfFtxUHTmsGQKYCCk/edit#gid=0
    conversion_rate_LSU <- as.magpie(conversion_rate_LSU)
    dimnames(conversion_rate_LSU) <- list("region", "year",c("large", "small"))


    # mapping2 <- toolGetMapping(name = "regionmappingH12.csv", type = "cell")
    # livestock_FAO_scaled <- livestock_FAO * conversion_rate_LSU * toolAggregate(fbask_pasture_fraction, rel = mapping2, from="CountryCode", to="RegionCode")
    livestock_FAO_scaled <- livestock_FAO * conversion_rate_LSU * fbask_pasture_fraction
    livestock_FAO_scaled <- dimSums(livestock_FAO_scaled[,,c(1,4)])
    getYears(livestock_FAO_scaled) <- past
    livestock_FAO_scaled <- livestock_FAO_scaled[unique(mapping$iso)]
    livestock_cell <- toolAggregate(livestock_FAO_scaled, rel=mapping, from="iso", to="celliso", partrel=F, weight = grass_distribution)

    # livestock_cell[livestock_cell>5e5] <- NA
    # hist(livestock_cell, breaks = 1000)
    # summary(livestock_cell)
    # quantile(livestock_cell, probs = 0.95)
    # getCells(livestock_cell)[is.na(livestock_cell[,1,])]
    # luplot::plotmap2(livestock_cell[,1,])


    luplot::plotmap2(livestock_cell[,1,])

    livestock_cell[MAGPasturearea==0] <- 0 #setting to 0 the LSUs fed on grass in regions with no pasture areas.
    lsu_ha <- livestock_cell/(MAGPasturearea*1e6)
    lsu_ha[is.nan(lsu_ha)] <- 0
    # lsu_ha[lsu_ha>3] <- NA
    lsu_ha[lsu_ha>2.27] <- 2.27
    lsu_ha_scaled <- scale(lsu_ha, center = center["lsu"] , scale = scale["lsu"])

    # luplot::plotmap2(lsu_ha)

    # luplot::plotmap2(lsu_ha_scaled)
    # sum(MAGPasturearea[problem_countries]) / sum(MAGPasturearea) * 100
    # country2map(livestock_FAO_scaled, c("iso_a3" = "Region"))

    #################################
    ### Predicted grass yield     ###
    #################################

    model <- load_model_hdf5(path(folder,files[grep("model", files)]))
    prediction <- predict(model, cbind(environmental_data, lsu_ha_scaled))
    prediction[prediction<0] <- 0
    prediction <- round(prediction, 3)
    gCm2_To_tDM <- (10000 * 2.21 / 1e6)
    pred_pasture_productivity <- prediction * gCm2_To_tDM
    pred_pasture_productivity <- cbind(environmental_data_3[,c(1,2,3)],pred_pasture_productivity) # fix
    pred_pasture_productivity_mag <- mutate(pred_pasture_productivity, celliso = paste(region,region1,sep="."))[-c(1,2)]
    pred_pasture_productivity_mag <- mutate(pred_pasture_productivity_mag, years = paste0("y",year))
    pred_pasture_productivity_mag <- as.magpie(pred_pasture_productivity_mag[,c(3,4,2)])

    # luplot::plotmap2(pred_pasture_productivity_mag)

    pred_pasture_production <- pred_pasture_productivity_mag * (MAGPasturearea*1e6)

    # luplot::plotmap2(pred_pasture_production)

    # country2map(pred_pasture_production, c("iso_a3" = "Region"))

    #################################
    ### FAO past. productivity    ###
    #################################
    FAO_pasture_demand <- calcOutput("FAOmassbalance",aggregate = F)[,,"dm"][,,"feed"][,past,"pasture"]
    FAO_pasture_demand <- FAO_pasture_demand[unique(mapping$iso)]*1e6
    FAO_pasture_demand_cell <- toolAggregate(FAO_pasture_demand, rel=mapping, from="iso", to="celliso", partrel=F, weight = grass_distribution)
    FAO_pasture_demand_cell[MAGPasturearea==0] <- 0 #Setting to zero pasture demand in countries with no pasture area

    # luplot::plotmap2(FAO_pasture_demand_cell)

    # FAO_pasture_productivity <- FAO_pasture_demand_cell/(MAGPasturearea)
    # FAO_pasture_productivity[is.nan(FAO_pasture_productivity)] <- 0
    # FAO_pasture_productivity <- FAO_pasture_productivity[getCells(environmental_data_2)]
    # par(mfrow=c(1,2))
    # plot(FAO_pasture_productivity, pred_pasture_productivity)
    # plot(pred_pasture_productivity)
    # par(mfrow=c(1,1))
    # country2map(FAO_pasture_productivity, c("iso_a3" = "Region"))

    # FAO_pasture_productivity[FAO_pasture_productivity>15] <-15


    #################################
    ### Correction factor         ###
    #################################


    correction_factor <- FAO_pasture_demand_cell/pred_pasture_production
    correction_factor[is.nan(correction_factor)] <- 0
    correction_factor[is.infinite(correction_factor)] <- 0
    quantile(correction_factor, probs = 0.95)

    correction_factor[correction_factor> 17.54123] <-  17.54123

    sum(correction_factor) / sum(FAO_pasture_demand_cell)
    summary(correction_factor)
    luplot::plotmap2(correction_factor, lowcol="grey95",midcol="orange",
                     highcol="blue",midpoint=0.5)
    correction_factor[max(correction_factor)]
    View(as.data.frame(correction_factor))
    quantile(correction_factor, probs = 0.95)

    luplot::plotmap2(correction_factor)

    correction_factor[correction_factor>1] <- 1

    luplot::plotmap2(FAO_pasture_demand_cell)

    luplot::plotmap2(correction_factor)

    country2map(correction_factor, c("iso_a3" = "Region"))


    #Analisys
    plot(pred_pasture_production,FAO_pasture_demand)
    cor(pred_pasture_production,FAO_pasture_demand)^2
    plot(FAO_pasture_demand, col="red", type = "h")
    points(pred_pasture_production)
    points(correction_factor, col="blue", type = "h")

    sum(pred_pasture_production)/sum(FAO_pasture_demand)*100
    sd(correction_factor)
    mean(correction_factor)
    median(correction_factor)
    getCells(correction_factor)[correction_factor>9]

    return(list(
      x=correction_factor,
      weight=NULL,
      unit="",
      description="pasture correction factor for the historical dates",
      isocountries=FALSE))
  }


hist(MAGPasturearea,breaks = 200)
plot(MAGPasturearea)


luplot::plotmap2(cellarea)
plot(cellarea)
