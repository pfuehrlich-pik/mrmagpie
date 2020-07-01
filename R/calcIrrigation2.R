#' @title calcIrrigation2
#' @description This function calculates irrigation water withdrawal based on LPJmL plant blue water consumption and different irrigation efficiencies
#'
#' @param selectyears years to be returned
#' @param version Switch between LPJmL4 and LPJmL5
#' @param climatetype Switch between different climate scenarios (default: "CRU_4")
#' @param cells Switch between "lpjcell" (67420) and "magpiecell" (59199)
#' @param time Time smoothing: average, spline or raw (default)
#' @param averaging_range only specify if time=="average": number of time steps to average
#' @param dof             only specify if time=="spline": degrees of freedom needed for spline
#' @param harmonize_baseline FALSE (default): no harmonization, TRUE: if a baseline is specified here data is harmonized to that baseline (from ref_year on)
#' @param ref_year Reference year for harmonization baseline (just specify when harmonize_baseline=TRUE)
#'
#' @return magpie object in cellular resolution
#' @author Felicitas Beier, Jens Heinke
#'
#' @examples
#' \dontrun{ calcOutput("Irrigation", aggregate = FALSE) }
#'
#' @importFrom magpiesets findset
#' @importFrom magclass setNames

#### Rewrite this function!!!!
### Get: crop water requirement (per crop) from LPJmL (e.g. NIR, blue water consumption, blue water transpiration?)
  # blue water consumption: crop-specific and climate-specific and spatially explicit
### Get: Field and conveyance efficiencies from literature (or: LPJmL)
  # conveyance efficiency: irrigation-system-specific
  # field efficiency: crop-specific and irrigation-system-specific

### Calculations:
# Withdrawal – Conveyance_losses – Field_losses = Consumption
# Conveyance_lossess = Non_consumptive_conveyance_loss + consumptive_conveyance loss
# Field_losses = Non_consumptive_field_loss + (Consumptive_field_loss + Blue_transpiration)

### Report irrig_ww (irrigation water withdrawal) instead of airrig (water applied on field)


calcIrrigation2 <- function(selectyears="all", cells="lpjcell",
                           version="LPJmL5", climatetype="CRU_4", time="raw", averaging_range=NULL, dof=NULL,
                           harmonize_baseline=FALSE, ref_year=NULL){

  if(harmonize_baseline==FALSE){

    if(time=="raw"){

      ##############################
      ######## Read in data ########
      ##############################
      ### Read in blue water consumption for irrigated crops (in m^3 per ha per yr):
      blue_water_consumption <- collapseNames(calcOutput("LPJmL", version=version, climatetype=climatetype, subtype="cwater_b_lpjcell", aggregate=FALSE,
                                                                 harmonize_baseline=FALSE,
                                                                 time="raw")[,,"irrigated"])
      years       <- getYears(blue_water_consumption)
      cropnames   <- getNames(blue_water_consumption)
      systemnames <- c("drip","sprinkler","surface")

      ### Field efficiencies from Jägermeyr et al. (global values) [placeholder!]
      field_efficiency                <- new.magpie(1:67420,years,sort(paste(systemnames, rep(cropnames,3), sep=".")),sets=c("cell","year","system.crop"))
      dimnames(field_efficiency)[[1]] <- gsub("GLO","LPJ",getCells(field_efficiency))
      field_efficiency[,,"drip"]      <- 0.88
      field_efficiency[,,"sprinkler"] <- 0.78
      field_efficiency[,,"surface"]   <- 0.52
      field_loss_shr <- 1-field_efficiency
      ### Use field efficiency from LPJmL here (by system, by crop, on 0.5 degree) [Does it vary by year?]

      ### Conveyance efficiency proxy [placeholder]
      conveyance_efficiency                <- new.magpie(1:67420,years,sort(paste(systemnames, rep(cropnames,3), sep=".")),sets=c("cell","year","system.crop"))
      dimnames(conveyance_efficiency)[[1]] <- gsub("GLO","LPJ",getCells(conveyance_efficiency))
      conveyance_efficiency[,,"surface"]   <- 0.7
      conveyance_efficiency[,,"drip"]      <- 0.95
      conveyance_efficiency[,,"sprinkler"] <- 0.95
      conveyance_loss_shr <- 1-conveyance_efficiency
      ### Use field efficiency from LPJmL here (by system, on 0.5 degree) [Does it vary by year?]

      # Load LPJmL to MAgPIE mapping to aggregate to MAgPIE crops
      LPJ2MAG      <- toolGetMapping( "MAgPIE_LPJmL.csv", type = "sectoral", where = "mappingfolder")

      ##############################
      ######## Calculations ########
      ##############################
      # Water Withdrawal – Conveyance_losses – Field_losses = Water Consumption
      water_withdrawal <- new.magpie(1:67420,years,sort(paste(systemnames, rep(LPJ2MAG$MAgPIE,3), sep=".")),sets=c("cell.region","year","system.crop"),fill=0)
      dimnames(water_withdrawal)[[1]] <- gsub("GLO","LPJ",getCells(water_withdrawal))

      for (s in systemnames){
        # Water applied to field = Non-consumptive field loss + consumptive loss
        water_applied_to_field <- (collapseNames(field_loss_shr[,,s]) * blue_water_consumption) + blue_water_consumption
        # Conveyance loss
        conveyance_loss <- (collapseNames(conveyance_loss_shr[,,s]) * water_applied_to_field)
        # Withdrawal from source
        tmp <- water_applied_to_field + conveyance_loss
        # Aggregate to MAgPIE crops
        tmp <- toolAggregate(tmp, LPJ2MAG, from="LPJmL", to="MAgPIE", dim=3, partrel=TRUE)
        getNames(tmp) <- paste0(s,".",getNames(tmp))
        names(dimnames(tmp))[3] <- paste0("system.",names(dimnames(tmp))[3])
        water_withdrawal[,,s] <- tmp[,,s]
      }

      # Remove negative values
      water_withdrawal[water_withdrawal<0] <- 0

    } else {
      # Time smoothing:
      x     <- calcOutput("Irrigation", version=version, climatetype=climatetype, aggregate=FALSE,
                          harmonize_baseline=FALSE, time="raw")

      # Smoothing data through average:
      if(time=="average"){
        water_withdrawal <- toolTimeAverage(x, averaging_range=averaging_range)

      # Smoothing data with spline method:
      } else if(time=="spline"){
        water_withdrawal <- toolTimeSpline(x, dof=dof)
        # Replace value in 2100 with value from 2099 (LPJmL output ends in 2099)
        if ("y2099" %in% getYears(water_withdrawal)) {
          water_withdrawal <- toolFillYears(water_withdrawal, c(getYears(water_withdrawal, as.integer=TRUE)[1]:2100))
        }
      } else if(time!="raw"){
        stop("Time argument not supported!")
      }
    }

  } else {
    # Harmonization
    if(time=="raw"){
      stop("Harmonization with raw data not possible. Select time='spline' when applying harmonize_baseline=TRUE")
    } else {
      # Load smoothed data
      baseline   <- calcOutput("Irrigation", version=version, climatetype=harmonize_baseline, aggregate=FALSE,
                             harmonize_baseline=FALSE, time=time, dof=dof, averaging_range=averaging_range)
      x          <- calcOutput("Irrigation", version=version, climatetype=climatetype, aggregate=FALSE,
                             harmonize_baseline=FALSE, time=time, dof=dof, averaging_range=averaging_range)
      # Harmonize to baseline
      water_withdrawal <- toolHarmonize2Baseline(x=x, base=baseline, ref_year=ref_year, limited=TRUE, hard_cut=FALSE)
    }
  }

  if(selectyears!="all"){
    years       <- sort(findset(selectyears,noset="original"))
    water_withdrawal  <- water_withdrawal[,years,]
  }

  ### Correct number of cells and transform to magpie object
  if (cells=="lpjcell"){
    out <- water_withdrawal
  } else if (cells=="magpiecell"){
    water_withdrawal <- water_withdrawal[magclassdata$cellbelongings$LPJ_input.Index,,]
    dimnames(water_withdrawal)[[1]] <- paste(magclassdata$half_deg$region,1:59199,sep='.')
    out <- water_withdrawal
  } else {
    stop("Cell argument not supported. Select lpjcell for 67420 cells or magpiecell for 59199 cells")
  }

  # Check for NAs
  if(any(is.na(water_withdrawal))){
    stop("produced NA airrig")
  }

  return(list(
    x=out,
    weight=NULL,
    unit="m^3 per ha per yr",
    description="Irrigation water (water applied in addition to rainfall) for different crop types following LPJmL irrigation system assumptions",
    isocountries=FALSE))
}
