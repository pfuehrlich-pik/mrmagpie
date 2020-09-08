#' @title calcYieldImprovementPotential
#' @description This calculates the yield improvement potential of irrigation for different crops
#'
#' @param version     switch between LPJmL4 and LPJmL5
#' @param climatetype switch between different climate scenarios (default: "CRU_4")
#' @param time            time smoothing: average, spline or raw (default)
#' @param averaging_range only specify if time=="average": number of time steps to average
#' @param dof             only specify if time=="spline": degrees of freedom needed for spline
#' @param harmonize_baseline FALSE (default): no harmonization, TRUE: if a baseline is specified here data is harmonized to that baseline (from ref_year on)
#' @param ref_year           reference year for harmonization baseline (just specify when harmonize_baseline=TRUE)
#' @param selectyears defaults to all years available
#' @param cells       switch between "lpjcell" (67420) and "magpiecell" (59199)
#'
#' @return magpie object in cellular resolution
#' @author Felicitas Beier
#'
#' @examples
#' \dontrun{ calcOutput("YieldImprovementPotential", aggregate=FALSE) }
#'
#' @import madrat
#' @import magclass

calcYieldImprovementPotential <- function(version="LPJmL5", climatetype="CRU_4", time="spline", averaging_range=NULL, dof=4,
                       harmonize_baseline=FALSE, selectyears=seq(1995, 2095,by=5), cells="magpiecell"){

  # read in land available for agricultural use (in mio. ha) and transform to ha
  land <- collapseNames(calcOutput("AvlLandSi", aggregate=FALSE)[,,"si0"])*1000000

  # read in yields (in tons/ha)
  yields    <- calcOutput("Yields", version=version, climatetype=climatetype, time=time, dof=dof,
                          harmonize_baseline=harmonize_baseline, aggregate=FALSE, years=selectyears)

  # yield gap (irrigated vs. rainfed) [tons per ha]
  tmp <- collapseNames(yields[,,"irrigated"])-collapseNames(yields[,,"rainfed"])
  # calculate yield gap in tons (considering available land)
  tmp <- tmp*land
  ### cap to 0
  tmp <- pmax(tmp,0)

  # cellular dimension
  if (cells=="magpiecell") {
    yield_gain <- tmp
  } else if (cells=="lpjcell") {
    getCells(tmp)  <- paste("GLO",magclassdata$cellbelongings$LPJ_input.Index,sep=".")
    yield_gain     <- new.magpie(1:NCELLS,getYears(tmp),getNames(tmp))
    yield_gain[,,] <- 0
    yield_gain[paste("GLO",magclassdata$cellbelongings$LPJ_input.Index,sep="."),,] <- tmp[,,]
    getCells(yield_gain) <- paste(lpj_cells_map$ISO,1:67420,sep=".")
    yield_gain           <- as.array(collapseNames(yield_gain))
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
    unit="tons",
    description="Yield improvement potential by irrigation for different crop types.",
    isocountries=FALSE))
}
