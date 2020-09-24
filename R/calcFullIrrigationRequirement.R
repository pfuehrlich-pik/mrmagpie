#' @title calcFullIrrigationRequirement
#' @description This function calculates the water requirements for full irrigation per cell per crop given potentially available land
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
#' @param iniarea if TRUE (default): already irrigated area is subracted, if FALSE: total potential land area is used
#' @param iniyear  Year of initialization for cropland area
#' @param irrigini Initialization data set for irrigation system initialization ("Jaegermeyr_lpjcell", "LPJmL_lpjcell")
#'
#' @return magpie object in cellular resolution
#' @author Felicitas Beier
#'
#' @examples
#' \dontrun{ calcOutput("FullIrrigationRequirement", aggregate=FALSE) }
#'
#' @import madrat
#' @import magclass

calcFullIrrigationRequirement <- function(version="LPJmL5", climatetype="HadGEM2_ES:rcp2p6:co2", harmonize_baseline=FALSE, time="spline", dof=4, cells="lpjcell", selectyears=seq(1995, 2095,by=5), iniyear=1995, iniarea=TRUE, irrigini="Jaegermeyr_lpjcell", irrig_requirement="withdrawal"){

  # read in irrigation water requirements [in m^3 per hectar per year]
  irrig_wat <- calcOutput("Irrigation2", version=version, cells="magpiecell", selectyears=selectyears, climatetype=climatetype, harmonize_baseline=harmonize_baseline, time=time, dof=dof, irrig_requirement=irrig_requirement, aggregate=FALSE)
  # pasture is not irrigated in MAgPIE
  irrig_wat <- irrig_wat[,,"pasture",invert=T]
  irrig_wat <- toolCell2isoCell(irrig_wat)

  # read in land available for agricultural use (in mio. ha)
  land <- collapseNames(calcOutput("AvlLandSi", aggregate=FALSE)[,,"si0"])
  if (iniarea) {
    # subtract area already irrigated in initialization (in mio. ha)
    crops_grown <- calcOutput("Croparea", sectoral="kcr", cells="magpiecell", physical=TRUE, cellular=TRUE, irrigation=TRUE, aggregate=FALSE)[,,"irrigated"]
    crops_grown <- collapseNames(dimSums(crops_grown[,paste0("y",iniyear),],dim=3))
    land        <- land - crops_grown
  }
  # transform units (from mio. ha to ha)
  land <- land*1000000

  # water requirements for full irrigation in cell per crop (in m^3)
  tmp <- irrig_wat*land

  # cellular dimension
  if (cells=="magpiecell") {
    out <- tmp
  } else if (cells=="lpjcell") {
    lpj_cells_map  <- toolGetMapping("LPJ_CellBelongingsToCountries.csv", type="cell")
    getCells(tmp)  <- paste("GLO",magclassdata$cellbelongings$LPJ_input.Index,sep=".")
    out     <- new.magpie(1:67420,getYears(tmp),getNames(tmp))
    out[,,] <- 0
    out[paste("GLO",magclassdata$cellbelongings$LPJ_input.Index,sep="."),,] <- tmp[,,]
    getCells(out) <- paste(lpj_cells_map$ISO,1:67420,sep=".")
  } else {
    stop("Cells argument not supported. Please select lpjcell for 67420 cells or magpiecell for 59199 cells")
  }

  # dimension names
  getSets(out)[c(1,2)] <- c("iso","cell")
  getSets(out)[c(4,5)] <- c("crop","system")

  # Checks
  if(any(is.na(out))){
    stop("produced NA full irrigation requirements")
  }

  return(list(
    x=out,
    weight=NULL,
    unit="m^3",
    description="full irrigation requirements per cell per crop per irrigation system",
    isocountries=FALSE))
}
