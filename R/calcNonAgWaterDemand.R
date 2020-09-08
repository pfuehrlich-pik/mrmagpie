#' @title calcNonAgWaterDemand
#' @description This function extracts non-agricultural water demand
#' @param selectyears years to be returned
#' @param source data source to be used
#' @param time Time smoothing: average, spline or raw (default)
#' @param averaging_range only specify if time=="average": number of time steps to average
#' @param dof             only specify if time=="spline": degrees of freedom needed for spline
#' @param waterusetype withdrawal (default) or consumption
#' @param seasonality total (default) or grper (growing period)
#'
#' @return magpie object in cellular resolution
#' @author Felicitas Beier
#'
#' @examples
#' \dontrun{ calcOutput("NonAgWaterDemand", aggregate = FALSE) }
#'
#' @import madrat
#' @import magclass
#' @importFrom mrcommons toolHarmonize2Baseline

calcNonAgWaterDemand <- function(selectyears="all", source="WATCH_ISIMIP_WATERGAP",
                                 time="raw", averaging_range=NULL, dof=NULL,
                                 waterusetype="withdrawal",seasonality="total"){

  # Old Non-Agricultural Waterdemand data (current default, will be deleted soon):
  if(source=="WATCH_ISIMIP_WATERGAP"){
    # Read in nonagricultural water demand:
    watdem_nonagr <- readSource("WATERGAP", convert="onlycorrect", subtype=source)
    # Add year 2100
    watdem_nonagr <- toolFillYears(watdem_nonagr, seq(getYears(watdem_nonagr, as.integer=TRUE)[1],2100,by=5))
  }

  # New Non-Agricultural Waterdemand data (will be new default)
  if(source=="WATERGAP2020"){

    if(time=="raw"){
      # Read in nonagricultural water demand:
      watdem_nonagr_WATERGAP      <- readSource("WATERGAP", convert="onlycorrect", subtype="WATERGAP2020")
      watdem_nonagr_ISIMIP_hist   <- readSource("ISIMIPinputs",subtype="ISIMIP3b:water:histsoc.waterabstraction",convert="onlycorrect")
      watdem_nonagr_ISIMIP_future <- readSource("ISIMIPinputs",subtype="ISIMIP3b:water:2015soc.waterabstraction",convert="onlycorrect")

      ### Combine datasets from different sources:
      # historical and future ISIMIP data:
      watdem_ISIMIP <- mbind(watdem_nonagr_ISIMIP_hist, watdem_nonagr_ISIMIP_future)

      # empty magpie object
      cells <- getCells(watdem_nonagr_WATERGAP)
      years <- getYears(watdem_ISIMIP)
      names <- c(getNames(watdem_nonagr_WATERGAP),paste0("ISIMIP.",getNames(watdem_ISIMIP)))
      watdem_nonagr <- new.magpie(cells,years,names)

      # historical and future ISIMIP data
      watdem_nonagr[,getYears(watdem_ISIMIP),paste0("ISIMIP.",getNames(watdem_ISIMIP))] <- watdem_ISIMIP[,getYears(watdem_ISIMIP),getNames(watdem_ISIMIP)]

      # historical data provided by ISIMIP (same for all scenarios)
      watdem_nonagr[,getYears(watdem_nonagr_ISIMIP_hist),] <- watdem_nonagr_ISIMIP_hist[,getYears(watdem_nonagr_ISIMIP_hist),]

      # future WATERGAP scenarios (adjusted to transition year of historical data)
      watdem_nonagr_WATERGAP_adjusted <- watdem_nonagr_WATERGAP
      watdem_nonagr_WATERGAP_adjusted[,,] <- NA

      # Calibration of WATERGAP data to ISIMIP baseline:
      # historical ISIMIP data in correct format for harmonization
      cells <- getCells(watdem_nonagr_WATERGAP)
      years <- getYears(watdem_nonagr_ISIMIP_hist)
      names <- getNames(watdem_nonagr_WATERGAP)
      tmp     <- new.magpie(cells,years,names)
      tmp[,,] <- 1
      tmp     <- tmp * watdem_nonagr_ISIMIP_hist

      # Harmonization
      watdem_nonagr_WATERGAP_adjusted <- toolHarmonize2Baseline(x=watdem_nonagr_WATERGAP, base=tmp, ref_year="y2005", limited=TRUE, hard_cut=FALSE)

      # WATERGAP adjusted future scenario data
      watdem_nonagr[,getYears(watdem_nonagr_WATERGAP),getNames(watdem_nonagr_WATERGAP)] <- watdem_nonagr_WATERGAP_adjusted[,getYears(watdem_nonagr_WATERGAP),getNames(watdem_nonagr_WATERGAP)]

      # Water consumption or water withdrawal:
      watdem_nonagr <- watdem_nonagr[,,waterusetype]

    } else {
      # Time smoothing:
      x     <- calcOutput("NonAgWaterDemand", selectyears=selectyears, source=source,
                          waterusetype=waterusetype, time="raw", averaging_range=NULL, dof=NULL, aggregate=FALSE)

      if(time=="average"){
        # Smoothing data through average:
        watdem_nonagr   <- toolTimeAverage(x, averaging_range=averaging_range)

      } else if(time=="spline"){
        # Smoothing data with spline method:
        watdem_nonagr   <- toolTimeSpline(x, dof=dof)
        # Replace value in 2100 with value from 2099
        if ("y2099" %in% getYears(watdem_nonagr)) {
          watdem_nonagr <- toolFillYears(watdem_nonagr, c(getYears(watdem_nonagr, as.integer=TRUE)[1]:2100))
        }

      } else if(time!="raw"){
        stop("Time argument not supported!")
      }
    }
  }

  # Select years to be returned
  if(selectyears!="all"){
    years         <- sort(findset(selectyears, noset="original"))
    watdem_nonagr <- watdem_nonagr[,years,]
  }

  ### Non-agricultural water demands in Growing Period
  if(seasonality=="grper"){
    # Get growing days per month
    grow_days <- calcOutput("GrowingPeriod", version="LPJmL5", climatetype=climatetype, time="spline", dof=4,
                            harmonize_baseline=harmonize_baseline, ref_year=ref_year, yield_ratio=0.1, aggregate=FALSE)
    # Growing days per year
    grow_days <- dimSums(grow_days,dim=3)

    # Adjust years
    years_watdem <- getYears(watdem_nonagr)
    years_grper  <- getYears(grow_days)
    if(length(years_watdem)>=length(years_grper)){
      years <- years_grper
    } else {
      years <- years_watdem
    }
    rm(years_grper, years_watdem)

    # Calculate non-agricultural water demand in growing period
    out         <- watdem_nonagr[,years,]*grow_days[,years,]/365
    description <- "Non-agricultural water demand (industry, electiricty, domestic) in growing period"
  } else if(seasonality=="total"){
    ### Total non-agricultural water demands per year
    out         <- watdem_nonagr[,,]
    description <- "Total non-agricultural water demand (industry, electiricty, domestic)"
  } else {
    stop("Specify seasonality! grper or total")
  }

  # Check for NAs
  if(any(is.na(watdem_nonagr))){
    stop("produced NA watdem_nonagr")
  }

  return(list(
    x=out,
    weight=NULL,
    unit="mio. m^3",
    description=description,
    isocountries=FALSE))
}
