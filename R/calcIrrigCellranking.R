#' @title calcIrrigCellranking
#' @description This function calculates a cellranking for the river basin discharge allocation based on yield improvement potential
#'
#' @param version     switch between LPJmL version for yields
#' @param climatetype switch between different climate scenarios for yields
#' @param time            average, spline or raw (default)
#' @param averaging_range just specify for time=="average": number of time steps to average
#' @param dof             just specify for time=="spline": degrees of freedom
#' @param harmonize_baseline FALSE (default) nothing happens, if a baseline is specified here data is harmonized to that baseline (from ref_year on)
#' @param ref_year           just specify for harmonize_baseline != FALSE : Reference year
#' @param selectyears years selected for yield gain
#' @param cells       switch between "lpjcell" (67420) and "magpiecell" (59199)
#' @param crops       switch between "magpie" and "lpjml" (default) crops
#' @param method      method of calculating the rank: "meancellrank" (default): mean over cellrank of proxy crops, "meancroprank": rank over mean of proxy crops, "cropcellrank": rank of one selected proxy crop
#' @param proxycrop   proxycrop(s) selected for rank calculation
#'
#' @return magpie object in cellular resolution
#' @author Felicitas Beier
#'
#' @examples
#' \dontrun{ calcOutput("calcIrrigCellranking", aggregate=FALSE) }

calcIrrigCellranking <- function(version="LPJmL5", climatetype="HadGEM2_ES:rcp2p6:co2", time="spline", averaging_range=NULL, dof=4, harmonize_baseline=FALSE, ref_year="y2015",
                                 selectyears="y1995", cells="lpjcell", crops="magpie", method="meancellrank", proxycrop=c("maiz", "rapeseed", "puls_pro")){

  # Read in potential yield gain per cell in initialization year
  yield_gain <- calcOutput("YieldImprovementPotential", version=version, climatetype=climatetype, harmonize_baseline=harmonize_baseline,
                           time=time, averaging_range=averaging_range, dof=dof, selectyears=selectyears,
                           cells=cells, crops=crops, aggregate=FALSE)
  # select proxy crops
  yield_gain <- yield_gain[,,proxycrop]
  # transform to array
  yield_gain <- as.array(yield_gain)

  # Calculate global cell rank
  if (method=="meancellrank"){

    # calculate rank of proxy crops
    cellrank <-  NULL
    for (crop in getNames(yield_gain)) {
      # cell ranking for crop (from highest yield gain (rank=1) to lowest yield gain (rank=1+x))
      cropcellrank     <- floor(rank(-yield_gain[,,crop]))
      cellrank[[crop]] <- cropcellrank
      rm(cropcellrank)
    }

    # calculate mean over cropcellranks
    cellrank    <- as.list(cellrank)
    glocellrank <- floor(rank(rowMeans(do.call(cbind,cellrank))))

  } else if (method=="meancroprank"){

    # normalize yield gains of proxy crops
    ##### STILL MISSING

    # calculate average yield gain over normalized proxy crops
    yield_gain <- dimSums(yield_gain,dim=3)/length(getNames(yield_gain))

    # calculate rank
    glocellrank <- floor(rank(-yield_gain))

  } else if (method=="cropcellrank"){

    # calculate yield gain of proxy crop
    glocellrank <- floor(rank(-yield_gain[,,proxycrop]))

  } else {
    stop("Please select a method for rank calculation")
  }

  glocellrank <- as.magpie(glocellrank,spatial=1)

  # Check for NAs
  if(any(is.na(glocellrank))){
    stop("Function YieldImprovementPotential produced NAs")
  }

  return(list(
    x=glocellrank,
    weight=NULL,
    unit="1",
    description="Rank of cell according to yield gain potential by irrigation",
    isocountries=FALSE))
}
