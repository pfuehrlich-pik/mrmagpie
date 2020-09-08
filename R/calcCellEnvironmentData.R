#' @title calcCellEnvironmentData
#' @description Loads the ANN trained grass emulator with historical data and calculates the pasture correction factor
#' @param model grass harvest machine learning model ID
#' @param GCMModel Global climate model name
#' @param rcp RCP
#' @return MAgPIE objects with pasture correction factor on a cellular level.
#' @param model grass harvest machine learning model ID
#' @author Marcos Alves
#' @examples
#'
#' \dontrun{
#' calcOutput("CellEnvironmentData", model = "f41f19be67", GCMModel = "HadGEM2", rcp = "rcp85")
#' }
#'
#' @import madrat
#' @import magpiesets
#' @import magclass
#' @export
#'

calcCellEnvironmentData <-
  function(model, GCMModel, rcp) {

    time_steps <- findset("t_all")

    #################################
    ### Loading emulator data      ##
    #################################

    inputs_vec <- as.vector(readSource("GrassYldEmu", subtype = paste(model,"inputs", sep = "."), convert="onlycorrect"))
    magpie_center          <- readSource("GrassYldEmu", subtype =paste(model,"mean", sep = "."), convert="onlycorrect")
    center                 <- as.numeric(magpie_center)
    names(center)          <- getCells(magpie_center)
    magpie_scale           <- readSource("GrassYldEmu", subtype = paste(model,"stddevs", sep = "."), convert="onlycorrect")
    scale                  <- as.numeric(magpie_scale)
    names(scale)           <- getCells(magpie_scale)

    #################################
    ### Loading environmental data ##
    #################################

    soil_char  <- c("Ks", "Sf", "w_pwp", "w_fc", "w_sat", "hsg", "tdiff_0", "tdiff_15", "tdiff_100", "cond_pwp", "cond_100", "cond_100_ice")
    ml_inputs  <- inputs_vec[is.na(match(inputs_vec, "lsu"))]
    if (any(soil_char %in% ml_inputs)) {
      inputs   <- ml_inputs[is.na(match(ml_inputs, soil_char))]
      if (!("soil" %in% inputs)) {
        inputs <- append(inputs, "soil")
      }
    }

    list_mag_obj <- list()

    for (input in inputs) {
      if ("temperature" %in% input) {
        list_mag_obj[[input]] <- toolHoldConstantBeyondEnd(calcOutput("GCMClimate", aggregate = FALSE, GCMModel = GCMModel, ClimateVariable = "temperature", rcp = rcp))[, time_steps, ]
      }

      if ("precipitation" %in% input) {
        list_mag_obj[[input]] <- toolHoldConstantBeyondEnd(calcOutput("GCMClimate", aggregate = FALSE, GCMModel = GCMModel, ClimateVariable = "precipitation", rcp = rcp))[, time_steps, ]
      }

      if ("lwnet" %in% input) {
        list_mag_obj[[input]] <- toolHoldConstantBeyondEnd(calcOutput("GCMClimate", aggregate = FALSE, GCMModel = GCMModel, ClimateVariable = "longwave_radiation", rcp = rcp))[, time_steps, ]
      }

      if ("rsds" %in% input) {
        list_mag_obj[[input]] <- toolHoldConstantBeyondEnd(calcOutput("GCMClimate", aggregate = FALSE, GCMModel = GCMModel, ClimateVariable = "shortwave_radiation", rcp = rcp))[, time_steps, ]
      }

      if ("wetdays" %in% input) {
        list_mag_obj[[input]] <- toolHoldConstantBeyondEnd(calcOutput("GCMClimate", aggregate = FALSE, GCMModel = GCMModel, ClimateVariable = "wetdays", rcp = rcp))[, time_steps, ]
      }

      if ("co2" %in% input) {
        list_mag_obj[[input]] <- toolHoldConstantBeyondEnd(calcOutput("CO2Atmosphere", aggregate = FALSE, rcp = rcp, level = "cellular"))[, time_steps, ]
      }

      if ("soil" %in% input) {
        list_mag_obj[[input]] <- toolHoldConstantBeyondEnd(calcOutput("SoilCharacteristics", aggregate = FALSE))[, time_steps, ]
      }
    }

    environmental_data <- list_mag_obj[[1]]
    for (i in 2:length(inputs)) {
      environmental_data <- mbind(environmental_data, list_mag_obj[[i]])
    }

    translation <- gsub("rsds","shortwave", gsub("lwnet","longwave",  gsub("co2","CO2ATM", ml_inputs)))
    keep <- paste0(translation, collapse = "+|")
    keep <- grep(keep, getNames(environmental_data))
    environmental_data <- environmental_data[,,keep]
    getNames(environmental_data) <- ml_inputs
    #manually scaling data with mean and average
    environmental_data <- (environmental_data - as.magpie(center[ml_inputs])) / as.magpie(scale[ml_inputs])

    return(list(
      x = environmental_data,
      weight = NULL,
      unit = paste(GCMModel, paste(rcp, "inputs:", paste(ml_inputs, collapse = " ")), collapse = " "),
      description = paste("Environmental contidions used for the model", model),
      isocountries = FALSE
    ))
  }
