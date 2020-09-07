#' @title calcFullIrrigationRequirement
#' @description This function calculates the requirements for full irrigation per cell per crop
#'
#' @param selectyears years to be returned
#' @param version     switch between LPJmL4 and LPJmL5
#' @param climatetype switch between different climate scenarios (default: "CRU_4")
#' @param cells       switch between "lpjcell" (67420) and "magpiecell" (59199)
#' @param time            time smoothing: average, spline or raw (default)
#' @param averaging_range only specify if time=="average": number of time steps to average
#' @param dof             only specify if time=="spline": degrees of freedom needed for spline
#' @param harmonize_baseline FALSE (default): no harmonization, TRUE: if a baseline is specified here data is harmonized to that baseline (from ref_year on)
#' @param ref_year           reference year for harmonization baseline (just specify when harmonize_baseline=TRUE)
#' @param irrig_requirement consumptive (consumption) or non-consumptive (withdrawals) irrigation water requirements
#'
#' @return magpie object in cellular resolution
#' @author Felicitas Beier
#'
#' @examples
#' \dontrun{ calcOutput("FullIrrigationRequirement", aggregate=FALSE) }
#'
#' @import madrat
#' @import magclass

calcFullIrrigationRequirement <- function(version="LPJmL5", selectyears=seq(1995, 2095,by=5), climatetype="HadGEM2_ES:rcp2p6:co2", harmonize_baseline=FALSE, time="spline", dof=4, irrig_requirement="withdrawal"){

  # read irrigation water requirements (in m^3 per hectar per year)
  irrig_wat <- calcOutput("Irrigation2", version=version, cells="magpiecell", selectyears=selectyears, climatetype=climatetype, harmonize_baseline=harmonize_baseline, time=time, dof=dof, irrig_requirement=irrig_requirement, aggregate=FALSE)
  # pasture is not irrigated in MAgPIE
  irrig_wat <- irrig_wat[,,"pasture",invert=T]
  irrig_wat <- toolCell2isoCell(irrig_wat)

  # read in land available for agricultural use (in mio. ha) and transform to ha
  land      <- collapseNames(calcOutput("AvlLandSi", aggregate=FALSE)[,,"si0"])*1000000

  # requirements for full irrigation in cell per crop (in m^3)
  x <- irrig_wat*land

  # Checks
  if(any(is.na(x))){
    stop("produced NA full irrigation requirements")
  }

  return(list(
    x=x,
    weight=NULL,
    unit="m^3",
    description="full irrigation requirements per cell per crop per irrigation system",
    isocountries=FALSE))
}
