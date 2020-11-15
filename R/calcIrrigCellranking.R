#' @title calcIrrigCellranking
#' @description This function calculates a cellranking for the river basin discharge allocation based on yield improvement potential through irrigation
#'
#' @param version     switch between LPJmL version for yields
#' @param climatetype switch between different climate scenarios for yields
#' @param time            time smoothing: average, spline (default) or raw
#' @param averaging_range just specify for time=="average": number of time steps to average
#' @param dof             just specify for time=="spline": degrees of freedom
#' @param harmonize_baseline FALSE (default) no harmonization, harmonization: if a baseline is specified here data is harmonized to that baseline (from ref_year onwards)
#' @param ref_year           just specify for harmonize_baseline != FALSE : Reference year
#' @param cellrankyear year for which cell rank is calculated
#' @param cells       switch between "lpjcell" (67420) and "magpiecell" (59199)
#' @param crops       switch between "magpie" and "lpjml" (default) crops
#' @param method      method of calculating the rank: "meancellrank" (default): mean over cellrank of proxy crops, "meancroprank": rank over mean of proxy crops
#' @param proxycrop   proxycrop(s) selected for rank calculation
#'
#' @return magpie object in cellular resolution
#' @author Felicitas Beier
#'
#' @examples
#' \dontrun{ calcOutput("calcIrrigCellranking", aggregate=FALSE) }

calcIrrigCellranking <- function(version="LPJmL5", climatetype="HadGEM2_ES:rcp2p6:co2", time="spline", averaging_range=NULL, dof=4, harmonize_baseline=FALSE, ref_year="y2015",
                                 cellrankyear="y1995", cells="lpjcell", crops="magpie", method="meancellrank", proxycrop=c("maiz", "rapeseed", "puls_pro")){

  ### Read in potential yield gain per cell
  yield_gain <- calcOutput("YieldImprovementPotential", version=version, climatetype=climatetype, harmonize_baseline=harmonize_baseline,
                           time=time, averaging_range=averaging_range, dof=dof, selectyears=cellrankyear,
                           cells=cells, crops=crops, aggregate=FALSE)
  # select proxy crops
  yield_gain <- yield_gain[,,proxycrop]

  ### Calculate global cell rank
  # Def. "meancellrank": ranking of cells or proxy crops, then: average over ranks
  if (method=="meancellrank"){

    # cell ranking for crop (from highest yield gain (rank=1) to lowest yield gain (rank=1+x))
    cropcellrank <- apply(-yield_gain,c(2,3),rank)

    # calculate mean over cropcellranks
    glocellrank  <- dimSums(cropcellrank,dim=3)/length(getNames(cropcellrank))
    # ties are solved by first occurrence
    glocellrank <- apply(glocellrank,2,rank,ties.method="first")

  # Def. "meancroprank": average over yield gain of proxycrops, then: ranking of resulting average yield gain
  } else if (method=="meancroprank"){

    # normalize yield gains of proxy crops (unity-based normalization)
    min_yield  <- as.magpie(apply(yield_gain,c(2,3),min))
    max_yield  <- as.magpie(apply(yield_gain,c(2,3),max))
    yield_gain <- (yield_gain-min_yield)/(max_yield-min_yield)

    # calculate average yield gain over normalized proxy crops
    yield_gain <- dimSums(yield_gain,dim=3)/length(getNames(yield_gain))

    # calculate rank (ties are solved by first occurence)
    glocellrank <- apply(-yield_gain,c(2,3),rank,ties.method="first")

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
