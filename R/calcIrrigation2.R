#' @title calcIrrigation2
#' @description This function calculates irrigation water requirements based on LPJmL plant blue water consumption and different irrigation efficiencies
#'
#' @param selectyears years to be returned
#' @param version Switch between LPJmL4 and LPJmL5
#' @param climatetype Switch between different climate scenarios (default: "CRU_4")
#' @param cells Switch between "lpjcell" (67420) and "magpiecell" (59199)
#' @param crops selects "magpie" (default) or "lpjml" crops
#' @param time Time smoothing: average, spline or raw (default)
#' @param averaging_range only specify if time=="average": number of time steps to average
#' @param dof             only specify if time=="spline": degrees of freedom needed for spline
#' @param harmonize_baseline FALSE (default): no harmonization, TRUE: if a baseline is specified here data is harmonized to that baseline (from ref_year on)
#' @param ref_year Reference year for harmonization baseline (just specify when harmonize_baseline=TRUE)
#' @param irrig_requirement consumptive (consumption) or non-consumptive (withdrawals) irrigation water requirements
#'
#' @return magpie object in cellular resolution
#' @author Felicitas Beier, Jens Heinke
#'
#' @examples
#' \dontrun{ calcOutput("Irrigation2", aggregate = FALSE) }
#'
#' @importFrom magpiesets findset
#' @importFrom magclass setNames

#### Rewrite this function!!!!
### Get: crop water requirement (per crop) from LPJmL (e.g. NIR, blue water consumption, blue water transpiration?)
  # blue water consumption: crop-specific and climate-specific and spatially explicit
### Get: Field and conveyance efficiencies from literature (or: LPJmL)
  # conveyance efficiency: irrigation-system-specific
  # field efficiency: crop-specific and irrigation-system-specific


calcIrrigation2 <- function(selectyears="all", cells="lpjcell", crops="magpie",
                           version="LPJmL5", climatetype="HadGEM2_ES:rcp2p6:co2", time="raw", averaging_range=NULL, dof=NULL,
                           harmonize_baseline=FALSE, ref_year=NULL, irrig_requirement="withdrawal"){

  sizelimit <- getOption("magclass_sizeLimit")
  options(magclass_sizeLimit=1e+12)
  on.exit(options(magclass_sizeLimit=sizelimit))

  if(harmonize_baseline==FALSE){

    if(time=="raw"){

      ##############################
      ######## Read in data ########
      ##############################
      ### Mappings
      lpj_cells_map <- toolGetMapping("LPJ_CellBelongingsToCountries.csv", type="cell")
      LPJ2MAG       <- toolGetMapping( "MAgPIE_LPJmL.csv", type = "sectoral", where = "mappingfolder")

      ### Read in blue water consumption for irrigated crops (in m^3 per ha per yr):
      blue_water_consumption <- collapseNames(calcOutput("LPJmL", version=version, climatetype=climatetype, subtype="cwater_b_lpjcell", aggregate=FALSE,
                                                                 harmonize_baseline=FALSE,
                                                                 time="raw")[,,"irrigated"])
      names(dimnames(blue_water_consumption))[1] <- "iso.cell"
      names(dimnames(blue_water_consumption))[3] <- "crop"
      years       <- getYears(blue_water_consumption)
      cropnames   <- getNames(blue_water_consumption)
      systemnames <- c("drip","sprinkler","surface")

      ### Field efficiencies from Jägermeyr et al. (global values) [placeholder!]
      field_efficiency                <- new.magpie(1:67420,years,sort(paste(systemnames, rep(cropnames,3), sep=".")),sets=c("iso.cell","year","system.crop"))
      getCells(field_efficiency)      <- paste(lpj_cells_map$ISO,1:67420,sep=".")
      field_efficiency[,,"drip"]      <- 0.88
      field_efficiency[,,"sprinkler"] <- 0.78
      field_efficiency[,,"surface"]   <- 0.52
      field_loss_shr <- 1-field_efficiency
      ### Use field efficiency from LPJmL here (by system, by crop, on 0.5 degree) [Does it vary by year?]

      ### Conveyance efficiency proxy [placeholder]
      conveyance_efficiency                <- new.magpie(1:67420,years,sort(paste(systemnames, rep(cropnames,3), sep=".")),sets=c("iso.cell","year","system.crop"))
      getCells(conveyance_efficiency)      <- paste(lpj_cells_map$ISO,1:67420,sep=".")
      conveyance_efficiency[,,"drip"]      <- 0.95
      conveyance_efficiency[,,"sprinkler"] <- 0.95
      conveyance_efficiency[,,"surface"]   <- 0.7
      conveyance_loss_shr <- 1-conveyance_efficiency
      ### Use field efficiency from LPJmL here (by system, on 0.5 degree) [Does it vary by year?]

      ##############################
      ######## Calculations ########
      ##############################
      # Withdrawal – Conveyance_losses – Field_losses = Consumption
      # Conveyance_lossess = Non_consumptive_conveyance_loss + consumptive_conveyance loss
      # Field_losses       = Non_consumptive_field_loss + (Consumptive_field_loss + Blue_transpiration)

      # Water Withdrawal – Conveyance_losses – Field_losses = Water Consumption
      water_applied_to_field <- blue_water_consumption*field_loss_shr + blue_water_consumption

      # Distinguish between evaporative and non-evaporative conveyance loss
      # (Note: According to ... et al. 50% of conveyance loss are evaporative)
      conveyance_loss        <- water_applied_to_field*conveyance_loss_shr
      water_withdrawal       <- conveyance_loss + water_applied_to_field
      # consumptive irrigation water = consumptive plant transpiration + evaporative conveyance loss
      water_consumption      <- blue_water_consumption + 0.5*conveyance_loss

      # Irrigation water requirements (consumption or withdrawals)
      if (irrig_requirement=="consumption"){
        irrig_requirements <- water_consumption
      } else if (irrig_requirement=="withdrawal"){
        irrig_requirements <- water_withdrawal
      } else {
        stop("Specify consumption or withdrawal in irrig_requirement")
      }

      # Aggregate to MAgPIE crops
      if (crops=="magpie") {
        irrig_requirements <- toolAggregate(irrig_requirements, LPJ2MAG, from="LPJmL", to="MAgPIE", dim=3.1, partrel=TRUE)
      }

    } else {
      # Time smoothing:
      x     <- calcOutput("Irrigation2", version=version, climatetype=climatetype, crops=crops, aggregate=FALSE,
                          harmonize_baseline=FALSE, time="raw", irrig_requirement=irrig_requirement)

      # Smoothing data through average:
      if(time=="average"){
        irrig_requirements <- toolTimeAverage(x, averaging_range=averaging_range)

      # Smoothing data with spline method:
      } else if(time=="spline"){
        irrig_requirements <- toolTimeSpline(x, dof=dof)
        # Replace value in 2100 with value from 2099 (LPJmL output ends in 2099)
        if ("y2099" %in% getYears(irrig_requirements)) {
          irrig_requirements <- toolFillYears(irrig_requirements, c(getYears(irrig_requirements, as.integer=TRUE)[1]:2100))
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
      baseline   <- calcOutput("Irrigation2", version=version, climatetype=harmonize_baseline, crops=crops, aggregate=FALSE,
                             harmonize_baseline=FALSE, time=time, dof=dof, averaging_range=averaging_range, irrig_requirement=irrig_requirement)
      x          <- calcOutput("Irrigation2", version=version, climatetype=climatetype, crops=crops, aggregate=FALSE,
                             harmonize_baseline=FALSE, time=time, dof=dof, averaging_range=averaging_range, irrig_requirement=irrig_requirement)
      # Harmonize to baseline
      irrig_requirements <- toolHarmonize2Baseline(x=x, base=baseline, ref_year=ref_year, limited=TRUE, hard_cut=FALSE)
    }
  }

  if(selectyears!="all"){
    years               <- sort(findset(selectyears,noset="original"))
    irrig_requirements  <- irrig_requirements[,years,]
  }

  ### Correct number of cells
  if (cells=="lpjcell"){
    out <- irrig_requirements
  } else if (cells=="magpiecell"){
    irrig_requirements                <- irrig_requirements[magclassdata$cellbelongings$LPJ_input.Index,,]
    dimnames(irrig_requirements)[[1]] <- paste(magclassdata$half_deg$region,1:59199,sep='.')
    out <- irrig_requirements
  } else {
    stop("Cell argument not supported. Select lpjcell for 67420 cells or magpiecell for 59199 cells")
  }

  # Check for NAs and negative values
  if(any(is.na(out))){
    stop("produced NA irrigation water requirements")
  }
  if(any(out<0)){
    stop("produced negative irrigation water requirements")
  }

  return(list(
    x=out,
    weight=NULL,
    unit="m^3 per ha per yr",
    description="Irrigation water consumption and withdrawal requirements for irrigation for different crop types under different irrigation systems",
    isocountries=FALSE))
}
