#' @title calcYields
#' @description This function extracts yields from LPJmL to MAgPIE
#'
#' @param version Switch between LPJmL4 and LPJmL5
#' @param climatetype Switch between different climate scenarios (default: "CRU_4")
#' @param time average, spline or raw (default)
#' @param averaging_range just specify for time=="average": number of time steps to average
#' @param dof             just specify for time=="spline": degrees of freedom
#' @param harmonize_baseline FALSE (default) nothing happens, if a baseline is specified here data is harmonized to that baseline (from ref_year on)
#' @param ref_year just specify for harmonize_baseline != FALSE : Reference year
#' @param calib_proxy calibrated proxy to FAO values if set to TRUE
#' @param split_cropcalc split calculation for different crop types (e.g. for parallelization)
#' @param crops magpie crops or lpjml crops
#' @param selectyears defaults to all years available
#'
#' @return magpie object in cellular resolution
#' @author Kristine Karstens
#'
#' @examples
#' \dontrun{ calcOutput("Yields", aggregate = FALSE) }
#'
#' @importFrom magpiesets findset
#' @importFrom magclass getYears add_columns dimSums time_interpolate
#' @importFrom madrat toolFillYears

calcYields <- function(version="LPJmL5", climatetype="CRU_4", time="spline", averaging_range=NULL, dof=4,
                       harmonize_baseline=FALSE, ref_year="y2015", calib_proxy=TRUE, split_cropcalc=TRUE, crops="magpie", selectyears="all"){

  sizelimit <- getOption("magclass_sizeLimit")
  options(magclass_sizeLimit=1e+10)
  on.exit(options(magclass_sizeLimit=sizelimit))

  LPJ2MAG      <- toolGetMapping( "MAgPIE_LPJmL.csv", type = "sectoral", where = "mappingfolder")

  if(split_cropcalc){

    lpjml_crops  <- unique(LPJ2MAG$LPJmL)
    irrig_types  <- c("irrigated","rainfed")
    yields       <- NULL

    for(crop in lpjml_crops){

      subdata <- as.vector(outer(crop, irrig_types, paste, sep="."))
      tmp     <- calcOutput("LPJmL", version=version, climatetype=climatetype, subtype="harvest", subdata=subdata, time=time, averaging_range=averaging_range, dof=dof,
                            harmonize_baseline=harmonize_baseline, ref_year=ref_year, limited=TRUE, hard_cut=FALSE, selectyears=selectyears, aggregate=FALSE)

      yields  <- mbind(yields, tmp)
    }

  } else {

    yields    <- calcOutput("LPJmL", version=version, climatetype=climatetype, subtype="harvest", time=time, averaging_range=averaging_range, dof=dof,
                          harmonize_baseline=harmonize_baseline, ref_year=ref_year, limited=TRUE, hard_cut=FALSE, selectyears=selectyears, aggregate=FALSE)
  }
  # Aggregate to MAgPIE crops
  if (crops!="lpjml"){
    yields    <- toolAggregate(yields, LPJ2MAG, from = "LPJmL", to = "MAgPIE", dim=3.1, partrel=TRUE)

    # Check for NAs
    if(any(is.na(yields))){
      stop("produced NA yields")
    }

    if(calib_proxy){

      FAOproduction     <- collapseNames(calcOutput("FAOmassbalance_pre", aggregate=FALSE)[,,"production"][,,"dm"])
      MAGarea           <- calcOutput("Croparea", sectoral="kcr", physical=TRUE, aggregate=FALSE)


      MAGcroptypes  <- findset("kcr")
      missing       <- c("betr","begr")
      MAGcroptypes  <- setdiff(MAGcroptypes, missing)
      FAOproduction <- add_columns(FAOproduction[,,MAGcroptypes],addnm = missing,dim = 3.1)
      FAOproduction[,,missing] <- 0

      FAOYields         <- dimSums(FAOproduction,dim=1)/dimSums(MAGarea, dim=1)

      matchingFAOyears <- intersect(getYears(yields),getYears(FAOYields))
      FAOYields        <- FAOYields[,matchingFAOyears,]
      Calib            <- new.magpie("GLO", getYears(yields), c(getNames(FAOYields), "pasture"), fill=1)
      Calib[,matchingFAOyears,"oilpalm"]   <- FAOYields[,,"oilpalm"]/FAOYields[,,"groundnut"]      # LPJmL proxy for oil palm is groundnut
      Calib[,matchingFAOyears,"cottn_pro"] <- FAOYields[,,"cottn_pro"]/FAOYields[,,"groundnut"]    # LPJmL proxy for cotton is groundnut
      Calib[,matchingFAOyears,"foddr"]     <- FAOYields[,,"foddr"]/FAOYields[,,"maiz"]             # LPJmL proxy for fodder is maize
      Calib[,matchingFAOyears,"others"]    <- FAOYields[,,"others"]/FAOYields[,,"maiz"]            # LPJmL proxy for others is maize
      Calib[,matchingFAOyears,"potato"]    <- FAOYields[,,"potato"]/FAOYields[,,"sugr_beet"]       # LPJmL proxy for potato is sugar beet

      # interpolate between FAO years
      Calib <- toolFillYears(Calib, getYears(yields))

      # recalibrate yields for proxys
      yields <- Calib[,,getNames(yields, dim=1)] * yields
    }

    #check again, what makes sense irrigation=FALSE/TRUE?
    crop_area_weight <- dimSums(calcOutput("Croparea", sectoral="kcr", physical=TRUE, cellular=TRUE, irrigation=FALSE, aggregate = FALSE, years="y1995", round=6), dim=3)
  } else {
    crop_area_weight     <- yields
    crop_area_weight[,,] <- 1
  }

  return(list(
    x=yields,
    weight=crop_area_weight,
    unit="t per ha",
    description="Yields in tons per hectar for different crop types.",
    isocountries=FALSE))
}
