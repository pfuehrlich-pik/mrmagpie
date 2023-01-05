#' @title calcRangeSoilCarbonHist
#' @description calculates soil carbon for rangelands
#'
#' @param subtype subtypes
#' @param model trained model ID
#' @param lpjml lpjml version
#'
#' @return List of magpie objects with results on country level, weight on country level, unit and description.
#' @author Marcos Alves
#'
#' @examples
#' \dontrun{
#' calcOutput("GrassSoilCarbonHist ", subtype, model)
#' }
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr left_join

calcRangeSoilCarbonHist <- function(subtype = "ISIMIP3bv2:IPSL-CM6A-LR:ssp126:1965-2100", lpjml ,model = "9eaf9b") {

  x <- toolSplitSubtype(subtype, list(ismip = NULL, climatemodel = NULL, scenario = NULL, years = NULL))
  Data1 <- NULL
  Value <- NULL

  environment_data <- calcOutput("CollectEnvironmentData_new", subtype = subtype, sar = 1, aggregate = F, sel_feat = c("tas", "pr", "lwnet", "rsds", "CO2", "Ks", "Sf", "w_pwp", "w_fc", "w_sat", "hsg", "wet"))
  weights <- toolRefoldWeights(readSource("GrassSoilEmu", subtype = paste(subtype, model, "weights", sep = ":"), convert = F))
  mean_col <- readSource("GrassSoilEmu", subtype = paste(subtype, model, "mean_col", sep = ":"), convert = F)
  stddevs_col <- readSource("GrassSoilEmu", subtype = paste(subtype, model, "stddevs_col", sep = ":"), convert = F)
  mean_lab <- readSource("GrassSoilEmu", subtype = paste(subtype, model, "mean_lab", sep = ":"), convert = F)
  stddevs_lab <- readSource("GrassSoilEmu", subtype = paste(subtype, model, "stddevs_lab", sep = ":"), convert = F)
  inputs <- getItems(readSource("GrassSoilEmu", subtype = paste(subtype, model, "inputs", sep = ":"), convert = F),dim  = 1)
  hist_lsu_ha <- calcOutput("LsuDensityHist", disagg_type = "grassland", aggregate = F)
  land_ini_LUH2v2 <- calcOutput("LUH2v2", aggregate = F, landuse_types = "LUH2v2", cellular = TRUE)

  past <- intersect(getYears(environment_data), getYears(hist_lsu_ha))

  if (any(grepl("memory", inputs))) {
    hist_lsu_ha[hist_lsu_ha > 2.5] <- 2.5
    bin_hist_lsu_ha <- as.data.frame(hist_lsu_ha[, past[1], "range"])
    breaks <- c(seq(0, 2, 0.1), 2.25, 2.5)
    labels <- c(0.0, 0.2, 0.2, 0.4, 0.4, 0.6, 0.6, 0.8, 0.8, 1.0, 1.0, 1.2, 1.2, 1.4, 1.4, 1.6, 1.6, 1.8, 1.8, 2.0, 2.0, 2.5)
    bins <- cut(bin_hist_lsu_ha$Value, breaks = breaks, labels = labels, include.lowest = TRUE, right = F)
    bin_hist_lsu_ha$Value <- as.numeric(levels(bins)[bins])

    sc <- calcOutput("CollectSoilCarbonLSU", lsu_levels = c(seq(0, 2, 0.2), 2.5), lpjml = lpjml, climatemodel = x$climatemodel, scenario = paste0(x$scenario, "/co2/Nreturn0p5/limN"), sar = 1, aggregate = F, years = seq(1965, 2100, by = 5))
    sc_df <- as.data.frame(sc[, past[1], ])
    sc_df$Data1 <- gsub("p", ".", as.character.factor(sc_df$Data1))
    sc_df$Data1 <- as.numeric(sc_df$Data1)

    sc_df <- left_join(bin_hist_lsu_ha, sc_df, by = c("Cell" = "Cell", "Region" = "Region", "Year" = "Year", "Value" = "Data1"))
    sc_df_mg <- as.magpie(sc_df[, "Value.y"], spatial = 1)

    sc_start <- sc_df_mg

    past_num <- gsub("[^0-9.-]", "", past[-1])
    environment_data_past <- environment_data[, past[-1], ]
    hist_lsu_ha <- hist_lsu_ha[, past[-1], ]
    input_past <- mbind(setNames(hist_lsu_ha[, , "range"], grep("lsu", inputs, value = T)), environment_data_past)
    input_past <- as.data.frame(input_past)
    input_past_df <- pivot_wider(input_past, names_from = Data1, values_from = Value)
    memory_col <- grep("memory", inputs, value = T)

    input_past_df[[memory_col]] <- 0
    input_past_df[input_past_df$Year == past_num[1], memory_col] <- as.data.frame(sc_start)$Value
    input_df_scaled_past <- scale(input_past_df[, inputs], center = as.matrix(mean_col[,,inputs]), scale = as.matrix(stddevs_col[,,inputs]))

    soilc_range_past <- NULL
    for (i in 1:length(past_num)) {
      message(paste("Current timestep:", past_num[i]))
      current_tmp <- input_past_df$Year == past_num[i]
      soilc_range_past_tmp <- toolNeuralNet(input_df_scaled_past[current_tmp, ], weights, "softplus")
      soilc_range_past <- append(soilc_range_past, soilc_range_past_tmp)
      if (length(past_num) == i) {
        break
      }
      next_tmp <- input_past_df$Year == past_num[i + 1]
      input_df_scaled_past[next_tmp, memory_col] <- soilc_range_past_tmp
    }

    soilc_range_past <- soilc_range_past * as.numeric(stddevs_lab) + as.numeric(mean_lab)
    soilc_range_past <- cbind(input_past_df[, c("Cell", "Year")], soilc_range_past)
    soilc_range_past <- as.magpie(soilc_range_past, spatial = 1)
    soilc_range_past <- toolCell2isoCell(soilc_range_past)
    soilc_range_past <- add_columns(soilc_range_past, addnm = past[1], dim = 2.1)
    soilc_range_past[,past[1],] <- as.data.frame(sc_start)$Value
    soilc_range_past <- soilc_range_past[,order(getYears(soilc_range_past, as.integer = T)),]

  } else {
    environment_data_past <- environment_data[, past, ]
    hist_lsu_ha <- hist_lsu_ha[, past, ]
    input_past <- mbind(setNames(hist_lsu_ha[, , "range"], grep("lsu", inputs, value = T)), environment_data_past)
    input_past <- as.data.frame(input_past)
    input_past_df <- pivot_wider(input_past, names_from = Data1, values_from = Value)
    input_df_scaled_past <- scale(input_past_df[, inputs], center = mean_col[inputs], scale = stddevs_col[inputs])

    soilc_range_past <- toolNeuralNet(input_df_scaled_past, weights, "softplus")
    soilc_range_past <- soilc_range_past * as.numeric(stddevs_lab) + as.numeric(mean_lab)
    soilc_range_past <- cbind(input_past_df[, c("Cell", "Year")], soilc_range_past)
    soilc_range_past <- as.magpie(soilc_range_past, spatial = 1)
    soilc_range_past <- toolCell2isoCell(soilc_range_past)
  }

  return(list(
    x = soilc_range_past,
    weight = NULL,
    unit = "tC",
    description = "Historical rangelands carbon stocks",
    isocountries = FALSE
  ))
}
