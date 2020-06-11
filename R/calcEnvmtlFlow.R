#' @title calcEnvmtlFlow
#' @description This function calculates environmental flow requirements (EFR) for MAgPIE retrieved from LPJmL monthly discharge and water availability
#'
#' @param version Switch between LPJmL4 and LPJmL5
#' @param climatetype Switch between different climate scenarios (default: "CRU_4")
#' @param time Time smoothing: average or spline (default)
#' @param averaging_range only specify if time=="average": number of time steps to average
#' @param dof only specify if time=="spline": degrees of freedom needed for spline
#' @param harmonize_baseline FALSE (default): no harmonization, TRUE: if a baseline is specified here data is harmonized to that baseline (from ref_year on)
#' @param ref_year Reference year for harmonization baseline (just specify when harmonize_baseline=TRUE)
#' @param selectyears Years to be returned
#' @param LFR_val Strictness of environmental flow requirements
#' @param HFR_LFR_less10 High flow requirements (share of total water for cells) with LFR<10percent of total water
#' @param HFR_LFR_10_20 High flow requirements (share of total water for cells) with 10percent < LFR < 20percent of total water
#' @param HFR_LFR_20_30 High flow requirements (share of total water for cells) with 20percent < LFR < 30percent of total water
#' @param HFR_LFR_more30 High flow requirements (share of total water for cells) with LFR>30percent of total water
#' @param seasonality grper (default): EFR in growing period per year; total: EFR throughout the year; monthly: monthly EFRs
#'
#' @import magclass
#' @import madrat
#' @importFrom stats quantile
#'
#' @return magpie object in cellular resolution
#' @author Felicitas Beier, Abhijeet Mishra
#'
#' @examples
#' \dontrun{ calcOutput("EnvmtlFlow", aggregate = FALSE) }
#'

### Note: Rewrite calcEnvmtlFlow so that it is only based on LPJmL discharge, not calcAvlWater (because calcAvlWater will be rewritten and being based on EFRs)

calcEnvmtlFlow <- function(selectyears="all",
                           version="LPJmL4", climatetype="CRU_4", time="spline", averaging_range=NULL, dof=4,
                           harmonize_baseline=FALSE, ref_year="y2015",
                           LFR_val=0.1, HFR_LFR_less10=0.2, HFR_LFR_10_20=0.15, HFR_LFR_20_30=0.07, HFR_LFR_more30=0.00,
                           seasonality="grper"){

  if(harmonize_baseline==FALSE){

    ############################################################
    # Step 1 Determine monthly discharge low flow requirements #
    #        (LFR_monthly_discharge)                           #
    ############################################################

    ### Monthly Discharge
    monthly_discharge_magpie <- calcOutput("LPJmL", version=version, climatetype=climatetype, subtype="mdischarge", aggregate=FALSE,
                                           harmonize_baseline=FALSE, time="raw")
    # Extract years for quantile calculation
    years <- getYears(monthly_discharge_magpie, as.integer = TRUE)
    # Transform to array (faster calculation)
    monthly_discharge_magpie <-  as.array(collapseNames(monthly_discharge_magpie))

    ### Calculate LFR_quant
    ## Note: LFRs correspond to the Q90-value (i.e. to the discharge that is exceeded in nine out of ten months)
    ## (Bonsch et al. 2015). This is calculated via the 10%-quantile of monthly discharge.

    # Empty array with magpie object names
    LFR_quant <- array(NA,dim=c(dim(monthly_discharge_magpie)[1],length(years)),dimnames=list(dimnames(monthly_discharge_magpie)[[1]],paste("y",years,sep="")))

    # Quantile calculation: Yearly LFR quantile value
    for(year in years){
      # get the LFR_val quantile for each year for all cells
      LFR_quant[,paste("y",year,sep="")] <- apply(monthly_discharge_magpie[,paste("y",year,sep=""),],MARGIN=c(1),quantile,probs=LFR_val)
    }

    # Time-smooth LFR_quant
    LFR_quant <- as.magpie(LFR_quant)
    if(time=="average"){
      # Smoothing data through average:
      LFR_quant <- toolTimeAverage(LFR_quant, averaging_range=averaging_range)
    } else if(time=="spline"){
      # Smoothing data with spline method:
      LFR_quant <- toolTimeSpline(LFR_quant, dof=dof)
      # Replace value in 2100 with value from 2099 (LPJmL output ends in 2099)
      if ("y2099" %in% getYears(LFR_quant)) {
        LFR_quant <- toolFillYears(LFR_quant, c(getYears(LFR_quant, as.integer=TRUE)[1]:2100))
      }
    } else {
      stop("Time argument not supported!")
    }

    # Raw monthly discharge no longer needed at this point
    rm(monthly_discharge_magpie)

    ### Read in smoothed monthly discharge
    monthly_discharge_magpie <- calcOutput("LPJmL", version=version, climatetype=climatetype, subtype="mdischarge", aggregate=FALSE,
                                           harmonize_baseline=FALSE,
                                           time=time, dof=dof, averaging_range=averaging_range)

    # Transform to array (faster calculation)
    LFR_quant <- as.array(collapseNames(LFR_quant))
    monthly_discharge_magpie <- as.array(collapseNames(monthly_discharge_magpie))

    ### Calculate LFR discharge values for each month
    # If LFR_quant < magpie_discharge: take LFR_quant
    # Else: take magpie_discharge
    LFR_monthly_discharge <- monthly_discharge_magpie
    for (month in 1:12){
      tmp1 <- as.vector(LFR_quant)
      tmp2 <- as.vector(monthly_discharge_magpie[,,month])
      LFR_monthly_discharge[,,month] <- pmin(tmp1,tmp2)
    }
    # Remove no longer needed objects
    rm(LFR_quant)


    ################################################
    # Step 2 Determine low flow requirements (LFR) #
    #        from available water per month        #
    ################################################
    ### Available water per month (smoothed)
    avl_water_month <- calcOutput("AvlWater", version=version, climatetype=climatetype, seasonality="monthly", aggregate=FALSE,
                                  harmonize_baseline=FALSE,
                                  time=time, dof=dof, averaging_range=averaging_range)

    # Transform to array for faster calculation
    avl_water_month <- as.array(collapseNames(avl_water_month))

    # Empty array
    LFR     <- avl_water_month
    LFR[,,] <- NA

    ### Calculate LFRs
    LFR <- avl_water_month * (LFR_monthly_discharge/monthly_discharge_magpie)
    # There are na's where monthly_discharge_magpie was 0, replace with 0:
    LFR[is.nan(LFR)] <- 0

    ###################################################################
    # Step 3 Determie monthly high flow requirements (HFR)            #
    #        based on the ratio between LFR_month and avl_water_month #
    ###################################################################
    ## Note: "For rivers with low Q90 values, high-flow events are important
    ## for river channel maintenance, wetland flooding, and riparian vegetation.
    ## HFRs of 20% of available water are therefore assigned to rivers with a
    ## low fraction of Q90 in total discharge. Rivers with a more stable flow
    ## regime receive a lower HFR." (Bonsch et al. 2015)
    HFR     <- LFR
    HFR[,,] <- NA

    HFR[LFR<0.1*avl_water_month]  <- HFR_LFR_less10 * avl_water_month[LFR<0.1*avl_water_month]
    HFR[LFR>=0.1*avl_water_month] <- HFR_LFR_10_20  * avl_water_month[LFR>=0.1*avl_water_month]
    HFR[LFR>=0.2*avl_water_month] <- HFR_LFR_20_30  * avl_water_month[LFR>=0.2*avl_water_month]
    HFR[LFR>=0.3*avl_water_month] <- HFR_LFR_more30 * avl_water_month[LFR>=0.3*avl_water_month]
    HFR[avl_water_month<=0]       <- 0

    EFR <- LFR+HFR
    EFR <- as.magpie(EFR)

  } else {
    # Load baseline and climate EFR:
    baseline <- calcOutput("EnvmtlFlow", version=version, climatetype=harmonize_baseline, seasonality="monthly", aggregate=FALSE,
                           harmonize_baseline=FALSE, time=time, dof=dof, averaging_range=averaging_range)
    x        <- calcOutput("EnvmtlFlow", version=version, climatetype=climatetype, seasonality="monthly", aggregate=FALSE,
                           harmonize_baseline=FALSE, time=time, dof=dof, averaging_range=averaging_range)
    # Harmonize to baseline
    EFR <- toolHarmonize2Baseline(x=x, base=baseline, ref_year=ref_year, limited=TRUE, hard_cut=FALSE)
  }

  if(selectyears!="all"){
    years   <- sort(findset(selectyears, noset="original"))
    EFR     <- EFR[,years,]
  }

  ###########################################
  ############ FUNCTION OUTPUT  #############
  ###########################################

  ### EFR per cell per month
  if(seasonality=="monthly"){

    # Check for NAs
    if(any(is.na(EFR))){
      stop("produced NA EFR")
    }
    out=EFR
    description="Environmental flow requirements per cell per month"
  }

  ### Total water available per cell per year
  if(seasonality=="total"){

    # Sum up over all month:
    EFR_total <- dimSums(EFR, dim=3)

    # Read in available water (for Smakthin calculation)
    avl_water_total <- calcOutput("AvlWater", version=version, climatetype=climatetype, seasonality="total", aggregate=FALSE,
                                  harmonize_baseline=harmonize_baseline,
                                  time=time, dof=dof, averaging_range=averaging_range)

    # Reduce EFR to 50% of available water where it exceeds this threshold (according to Smakhtin 2004)
    EFR_total[which(EFR_total/avl_water_total>0.5)] <- 0.5*avl_water_total[which(EFR_total/avl_water_total>0.5)]

    # Check for NAs
    if(any(is.na(EFR_total))){
      stop("produced NA EFR_total")
    }
    out=EFR_total
    description="Total EFR per year"
  }

  ### Water available in growing period per cell per year
  if(seasonality=="grper"){

    # magpie object with days per month with same dimension as EFR
    tmp <- c(31,28,31,30,31,30,31,31,30,31,30,31)
    month_days <- new.magpie(names=dimnames(EFR)[[3]])
    month_days[,,] <- tmp
    month_day_magpie <- as.magpie(EFR)
    month_day_magpie[,,] <- 1
    month_day_magpie <- month_day_magpie * month_days

    # Daily water availability
    EFR_day   <- EFR/month_day_magpie

    # Growing days per month
    grow_days <- calcOutput("GrowingPeriod", version="LPJmL5", climatetype=climatetype, time=time, dof=dof, averaging_range=averaging_range,
                            harmonize_baseline=harmonize_baseline, ref_year=ref_year, yield_ratio=0.1, aggregate=FALSE)

    # Available water in growing period
    EFR_grper <- EFR_day*grow_days
    # Available water in growing period per year
    EFR_grper <- dimSums(EFR_grper, dim=3)
    # Read in available water (for Smakthin calculation)
    avl_water_grper <- calcOutput("AvlWater", version=version, climatetype=climatetype, seasonality="grper", aggregate=FALSE,
                                  harmonize_baseline=harmonize_baseline,
                                  time=time, dof=dof, averaging_range=averaging_range)

    # Reduce EFR to 50% of available water where it exceeds this threshold (according to smakhtin 2004)
    EFR_grper[which(EFR_grper/avl_water_grper>0.5)] <- 0.5*avl_water_grper[which(EFR_grper/avl_water_grper>0.5)]

    # Check for NAs
    if(any(is.na(EFR_grper))){
      stop("produced NA EFR_grper")
    }
    out=EFR_grper
    description="EFR in growing period per year"
  }

  return(list(
    x=out,
    weight=NULL,
    unit="mio. m^3",
    description=description,
    isocountries=FALSE))
}
