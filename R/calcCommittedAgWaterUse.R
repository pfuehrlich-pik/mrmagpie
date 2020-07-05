#' @title calcCommittedAgWaterUse
#' @description This function calculates committed agricultural water uses that are used in the river routing algorithm for distributing available water across the basin
#'
#' @param version     Switch between LPJmL4 and LPJmL5
#' @param climatetype Switch between different climate scenarios (default: "CRU_4")
#' @param time            Time smoothing: average or spline (default)
#' @param averaging_range Only specify if time=="average": number of time steps to average
#' @param dof             Only specify if time=="spline": degrees of freedom needed for spline
#' @param iniyear Year of initialization
#'
#' @import magclass
#'
#' @return magpie object in cellular resolution
#' @author Felicitas Beier, Jens Heinke
#'
#' @examples
#' \dontrun{ calcOutput("CommittedAgWaterUse", aggregate = FALSE) }
#'

calcCommittedAgWaterUse <- function(version="LPJmL5", iniyear=1995, climatetype="HadGEM2_ES:rcp2p6:co2",
                           time="spline", dof=4, averaging_range=NULL){

  ## Read in Irrigation Water Withdrawals (in m^3 per hectar per year) [smoothed]
##### NOTE: Currently: airrig as placeholder until replaced by efficiency calculation
  airrig <- calcOutput("Irrigation2", version="LPJmL5", cells="lpjcell", selectyears=iniyear, climatetype=climatetype, harmonize_baseline=FALSE, time=time, dof=dof, aggregate=FALSE)
  # Pasture is not irrigated in MAgPIE
  airrig <- airrig[,,"pasture",invert=T]

  ## Read in cropland area (by crop) from crop area initialization (in mio. ha)
  crops_grown <- calcOutput("Croparea", sectoral="kcr", cells="lpjcell", physical=TRUE, cellular=TRUE, irrigation=TRUE, aggregate = FALSE)
  # Only initialization year needed
  crops_grown <- crops_grown[,paste0("y",iniyear),]
  # Only irrigated needed
  crops_grown <- collapseNames(crops_grown[,,"irrigated"])

  ## Committed agricultural uses (in mio. m^3 per year) [in initialization year]
  CAU <- airrig * crops_grown
  CAU <- dimSums(CAU,dim=3)

  return(list(
    x=CAU,
    weight=NULL,
    unit="mio. m^3 per year",
    description="water withdrawn for irrigation in each cell",
    isocountries=FALSE))
}
