#' @title calcLivestockDistribution
#' @description Disaggregate Livestock estimates based on the GLW3 dataset.
#' @return MAgPIE objects with livestock numbers on a cellular level.
#' @author Marcos Alves
#' @examples
#'
#' \dontrun{
#' calcOutput("LivestockDistribution")
#' }
#'
#' @import madrat
#' @importFrom tidyr pivot_wider
#' @import magpiesets
#' @import magclass
#' @import dplyr
#' @export
#'

#development
# library(mrcommons)
# library(stringi)
# library(tidyr)
# library(magpiesets)
# library(mrmagpie)
# setConfig(forcecache=T)
# setConfig(globalenv = T)
# setwd("C:/magpie_inputdata/sources")

calcLivestockDistribution <-
  function() {

    past <- findset("past")
    past <- past[7:length(past)]

    #############################
    ### Disaggregation weights###
    #############################

    mapping <- toolGetMapping(name = "CountryToCellMapping.csv", type = "cell")
    GLW3    <- readSource("GLW3", subtype = "Da", convert="onlycorrect")

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
    ### Calculating LSU/ha        ###
    #################################

    conversion_rate_LSU            <- c(0.7, 0.1) # source calulation from EU statistics https://docs.google.com/spreadsheets/d/1SZAAVl1SLwrrK6j329tq5zo1VZhfFtxUHTmsGQKYCCk/edit#gid=0
    conversion_rate_LSU            <- as.magpie(conversion_rate_LSU)
    dimnames(conversion_rate_LSU)  <- list("region", "year", c("large", "small"))

    livestock_FAO_scaled           <- livestock_FAO * conversion_rate_LSU * fbask_pasture_fraction
    livestock_FAO_scaled           <- dimSums(livestock_FAO_scaled[, , c(1, 4)])
    getYears(livestock_FAO_scaled) <- past
    livestock_FAO_scaled           <- livestock_FAO_scaled[unique(mapping$iso)]
    livestock_cell                 <- toolAggregate(livestock_FAO_scaled, rel = mapping, from = "iso", to = "celliso", weight = GLW3)
    livestock_cell                 <- livestock_cell/1e6

    #Analysis
    # sum(pred_pasture_production) / sum(FAO_pasture_demand_cell)
    # summary(correction_factor[, 1, ])
    # q <- quantile(correction_factor, probs = 0.95);q


    return(list(
      x = livestock_cell,
      weight = NULL,
      unit = "Total Livestock numbers per cell (mio)",
      description = "Pasture correction factor for the historical dates",
      isocountries = FALSE
    ))
  }
