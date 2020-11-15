#' @title calcYieldImprovementPotential
#' @description This function calculates the yield improvement potential of irrigation for different crops
#'
#' @param version     switch between LPJmL4 and LPJmL5 of calcYields function
#' @param climatetype switch between different climate scenarios (default: "CRU_4") of calcYields function
#' @param time            time smoothing of calcYields function: average, spline (default) or raw
#' @param averaging_range only specify if time=="average": number of time steps to average
#' @param dof             only specify if time=="spline": degrees of freedom needed for spline
#' @param harmonize_baseline harmonization in calcYields function: FALSE (default): no harmonization, TRUE: if a baseline is specified here data is harmonized to that baseline (from ref_year onwards)
#' @param ref_year           reference year for harmonization baseline (just specify when harmonize_baseline=TRUE)
#' @param selectyears years selected in calcYields function
#' @param cells       switch between "lpjcell" (67420) and "magpiecell" (59199)
#' @param crops       switch between "magpie" and "lpjml" (default) crops
#'
#' @return magpie object in cellular resolution
#' @author Felicitas Beier
#'
#' @examples
#' \dontrun{ calcOutput("YieldImprovementPotential", aggregate=FALSE) }
#'
#' @import madrat
#' @import magclass

calcYieldImprovementPotential <- function(version="LPJmL5", climatetype="HadGEM2_ES:rcp2p6:co2", time="spline", averaging_range=NULL, dof=4,
                       harmonize_baseline=FALSE, selectyears=seq(1995, 2095,by=5), cells="magpiecell", crops="lpjml"){

  # read in yields [in tons/ha]
  yields <- calcOutput("Yields", version=version, climatetype=climatetype, time=time, dof=dof,
                          harmonize_baseline=harmonize_baseline, aggregate=FALSE, selectyears=selectyears, crops=crops)

  # yield gap (irrigated vs. rainfed) [in tons/ha]
  tmp <- collapseNames(yields[,,"irrigated"])-collapseNames(yields[,,"rainfed"])
  # (Note: under N-stress, irrigation may lead to lower yields; also: irrigation may lead to shift in growing period -> tmp can have negative values)

  # cellular dimension
  if (cells=="magpiecell") {
    yield_gain <- tmp
  } else if (cells=="lpjcell") {
    lpj_cells_map  <- toolGetMapping("LPJ_CellBelongingsToCountries.csv", type="cell")
    getCells(tmp)  <- paste("GLO",magclassdata$cellbelongings$LPJ_input.Index,sep=".")
    yield_gain     <- new.magpie(1:67420,getYears(tmp),getNames(tmp))
    yield_gain[,,] <- 0
    yield_gain[paste("GLO",magclassdata$cellbelongings$LPJ_input.Index,sep="."),,] <- tmp[,,]
    getCells(yield_gain) <- paste(lpj_cells_map$ISO,1:67420,sep=".")
  } else {
    stop("Cells argument not supported. Please select lpjcell for 67420 cells or magpiecell for 59199 cells")
  }

  # Check for NAs
  if(any(is.na(yield_gain))){
    stop("Function YieldImprovementPotential produced NAs")
  }

  return(list(
    x=yield_gain,
    weight=NULL,
    unit="tons per ha",
    description="Yield improvement potential by irrigation for different crop types.",
    isocountries=FALSE))
}
