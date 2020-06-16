#' @title calcCommittedAgWaterUse
#' @description This function calculates committed agricultural water uses that are used in the river routing algorithm for distributing available water across the basin
#'
#' @param version Switch between LPJmL4 and LPJmL5
#' @param climatetype Switch between different climate scenarios (default: "CRU_4")
#' @param time Time smoothing: average or spline (default)
#' @param averaging_range only specify if time=="average": number of time steps to average
#' @param dof only specify if time=="spline": degrees of freedom needed for spline
#'
#' @import magclass
#' @import madrat
#' @importFrom stats quantile
#'
#' @return magpie object in cellular resolution
#' @author Felicitas Beier, Jens Heinke
#'
#' @examples
#' \dontrun{ calcOutput("CommittedAgWaterUse", aggregate = FALSE) }
#'

### Note: Rewrite calcEnvmtlFlow so that it is only based on LPJmL discharge, not calcAvlWater (because calcAvlWater will be rewritten and being based on EFRs)

calcCommittedAgWaterUse <- function(version="LPJmL5", climatetype="HadGEM2_ES:rcp2p6:co2",
                           time="spline", dof=4, averaging_range=NULL){

  # Read in Irrigation Water Withdrawals (in m^3 per hectar per year) [smoothed]
  # NOTE: Currently: airrig as placeholder until replaced by efficiency calculation
  airrig <- calcOutput("Irrigation", version="LPJmL5", climatetype=climatetype, harmonize_baseline=FALSE, time=time, dof=dof, aggregate=FALSE)
  ## read in only 1995 (initialization year)

  # Only initialization year needed
  airrig <- airrig[,"y1995",]
  # Only irrigated needed
  airrig <- collapseNames(airrig[,,"irrigated"])

  # Read in area irrigated from land use initialization (in mio. ha) [NEEDED?]
  #area_irrigated <- calcOutput("AreaEquippedForIrrigation", aggregate=FALSE, cellular=TRUE, source="LUH2v2")

  # Read in cropland area (by crop) from crop area initialization (in mio. ha)
  ## read in only 1995 (initialization year)
  crops_grown <- calcOutput("Croparea", sectoral="kcr", physical=TRUE, cellular=TRUE, irrigation=TRUE, aggregate = FALSE)
  # Only initialization year needed
  crops_grown <- crops_grown[,"y1995",]
  # Only irrigated needed
  crops_grown <- collapseNames(crops_grown[,,"irrigated"])

  ## For now: delete pasture from airrig (later: add pasture area)
  #airrig <- airrig[,,"pasture",invert=T]
  land_area <- calcOutput("LanduseInitialisation", aggregate=FALSE, cellular=TRUE, land="fao", input_magpie=TRUE)
  land_area <- land_area[,"y1995",]



  # Committed agricultural uses (in mio. m^3 per year) [in initialization year]
  CAU <- airrig * crops_grown

  return(list(
    x=CAU,
    weight=NULL,
    unit="mio. m^3 per year",
    description="water withdrawn for irrigation in each cell", ### Note: currently not water withdrawn, but water applied to field
    isocountries=FALSE))
}
