#' @title calcCommittedAgWaterUse
#' @description This function calculates committed agricultural water uses that are used in the river routing algorithm for distributing available water across the basin
#'
#' @param version     Switch between LPJmL4 and LPJmL5
#' @param climatetype Switch between different climate scenarios (default: "CRU_4")
#' @param time            Time smoothing: average or spline (default)
#' @param averaging_range Only specify if time=="average": number of time steps to average
#' @param dof             Only specify if time=="spline": degrees of freedom needed for spline
#' @param iniyear  Year of initialization for cropland area
#' @param irrigini Initialization data set for irrigation system initialization ("Jaegermeyr_lpjcell", "LPJmL_lpjcell")
#'
#' @import magclass
#'
#' @return magpie object in cellular resolution
#' @author Felicitas Beier, Jens Heinke
#'
#' @examples
#' \dontrun{ calcOutput("CommittedAgWaterUse", aggregate = FALSE) }
#'

calcCommittedAgWaterUse <- function(version="LPJmL5", climatetype="HadGEM2_ES:rcp2p6:co2",
                           time="spline", dof=4, averaging_range=NULL, iniyear=1995, irrigini="Jaegermeyr_lpjcell"){

  ##############################
  ######## Read in data ########
  ##############################
  ## Irrigation system area initialization
  irrigation_system <- calcOutput("IrrigationSystem", source=irrigini, aggregate=FALSE)

  ## Read in Irrigation Water Withdrawals (in m^3 per hectar per year) [smoothed]
  irrig_withdrawal  <- calcOutput("Irrigation2", version="LPJmL5", cells="lpjcell", selectyears=iniyear, climatetype=climatetype, harmonize_baseline=FALSE, time=time, dof=dof, irrig_requirement="withdrawal", aggregate=FALSE)
  # Pasture is not irrigated in MAgPIE
  irrig_withdrawal  <- irrig_withdrawal[,,"pasture",invert=T]

  ## Read in Irrigation Water Consumption (in m^3 per hectar per year) [smoothed]
  irrig_consumption <- calcOutput("Irrigation2", version="LPJmL5", cells="lpjcell", selectyears=iniyear, climatetype=climatetype, harmonize_baseline=FALSE, time=time, dof=dof, irrig_requirement="consumption", aggregate=FALSE)
  # Pasture is not irrigated in MAgPIE
  irrig_consumption <- irrig_consumption[,,"pasture",invert=T]

  ## Read in cropland area (by crop) from crop area initialization (in mio. ha)
  crops_grown <- calcOutput("Croparea", sectoral="kcr", cells="lpjcell", physical=TRUE, cellular=TRUE, irrigation=TRUE, aggregate=FALSE)
  # Only initialization year needed
  crops_grown <- crops_grown[,paste0("y",iniyear),]
  # Only irrigated needed
  crops_grown <- collapseNames(crops_grown[,,"irrigated"])

  ##############################
  ######## Calculations ########
  ##############################
  ## Committed agricultural uses (in mio. m^3 per year) [in initialization year]
  # withdrawal
  CAW <- (irrigation_system[,,]*irrig_withdrawal[,,]) * crops_grown
  CAW <- dimSums(CAW,dim=3.1)
  getNames(CAW) <- paste(rep("withdrawal",3),getNames(CAW),sep=".")
  # consumption
  CAU <- (irrigation_system[,,]*irrig_consumption[,,]) * crops_grown
  CAU <- dimSums(CAU,dim=3.1)
  getNames(CAU) <- paste(rep("consumption",3),getNames(CAU),sep=".")
  # combine
  CAD <- mbind(CAW, CAU)
  names(dimnames(CAD))[1] <- "iso.cell"
  names(dimnames(CAD))[3] <- "wateruse.crop"

  return(list(
    x=CAD,
    weight=NULL,
    unit="mio. m^3 per year",
    description="committed agricultural water demands",
    isocountries=FALSE))
}
