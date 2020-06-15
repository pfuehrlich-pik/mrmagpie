#' @title calcCommittedAgWaterUse
#' @description This function calculates committed agricultural water uses that are used in the river routing algorithm for distributing available water across the basin
#'
#' @param version Switch between LPJmL4 and LPJmL5
#' @param climatetype Switch between different climate scenarios (default: "CRU_4")
#' @param time Time smoothing: average or spline (default)
#' @param averaging_range only specify if time=="average": number of time steps to average
#' @param dof only specify if time=="spline": degrees of freedom needed for spline
#' @param selectyears Years to be returned
#'
#' @import magclass
#' @import madrat
#' @importFrom stats quantile
#'
#' @return magpie object in cellular resolution
#' @author Felicitas Beier, Jens Heinke
#'
#' @examples
#' \dontrun{ calcOutput("EnvmtlFlow", aggregate = FALSE) }
#'

### Note: Rewrite calcEnvmtlFlow so that it is only based on LPJmL discharge, not calcAvlWater (because calcAvlWater will be rewritten and being based on EFRs)

calcCommittedAgWaterUse <- function(selectyears="all",
                           version="LPJmL5", climatetype="HadGEM2_ES:rcp2p6:co2",
                           time="spline", dof=4, averaging_range=NULL){

  # Read in Irrigation Water Withdrawals (in m^3 per hectar per year) [smoothed]
  # NOTE: Currently: airrig as placeholder until replaced by efficiency calculation
  airrig <- calcOutput("Irrigation", version="LPJmL5", climatetype=climatetype, harmonize_baseline=FALSE, time=time, dof=dof)

  # Read in area irrigated from land use initialization
  area_irrigated <-

  # Read in cropland area (by crop) from land use initialization
  crops_grown <-


  # Output: mio. m^3 of water used for irrigation (withdrawal) in each cell


    ### Monthly Discharge
    monthly_discharge_magpie <- calcOutput("LPJmL", selectyears=selectyears, version=version, climatetype=climatetype, subtype="mdischarge", aggregate=FALSE,
                                           harmonize_baseline=FALSE, time="raw")
    # Time frame for EFR calculation
    monthly_discharge_magpie <- monthly_discharge_magpie[,selectyears,]

    # Extract years for quantile calculation
    years <- getYears(monthly_discharge_magpie, as.integer = TRUE)
    # Transform to array (faster calculation)
    monthly_discharge_magpie <-  as.array(collapseNames(monthly_discharge_magpie))

    ### Calculate LFR_quant
    ## Note: LFRs correspond to the Q90-value (i.e. to the discharge that is exceeded in nine out of ten months)
    ## (Bonsch et al. 2015). This is calculated via the 10%-quantile of monthly discharge.

    # Empty array with magpie object names
    #LFR_quant_yearly <- array(NA,dim=c(dim(monthly_discharge_magpie)[1],dim(monthly_discharge_magpie)[2]),dimnames=list(dimnames(monthly_discharge_magpie)[[1]], dimnames(monthly_discharge_magpie)[[2]]))

    # Quantile calculation (according to Smakthin):
    # Get the monthly LFR_val quantile for all cells (across selected time frame)
    LFR_quant <- apply(monthly_discharge_magpie, MARGIN=c(1), quantile, probs=LFR_val)
    # Yearly LFRs
    LFR <- LFR_quant*12

    ### Mean annual discharge
    mean_annual_discharge <- apply(monthly_discharge_magpie, MARGIN=c(1), sum)/length(years)


    ### Calculate HFR

    ###################################################################
    # Step 3 Determie monthly high flow requirements (HFR)            #
    #        based on the ratio between LFR_month and avl_water_month #
    ###################################################################
    ## Note: "For rivers with low Q90 values, high-flow events are important
    ## for river channel maintenance, wetland flooding, and riparian vegetation.
    ## HFRs of 20% of available water are therefore assigned to rivers with a
    ## low fraction of Q90 in total discharge. Rivers with a more stable flow
    ## regime receive a lower HFR." (Bonsch et al. 2015)
    HFR <- LFR
    HFR <- NA

    HFR[LFR<0.1*mean_annual_discharge]  <- HFR_LFR_less10 * mean_annual_discharge[LFR<0.1*mean_annual_discharge]
    HFR[LFR>=0.1*mean_annual_discharge] <- HFR_LFR_10_20  * mean_annual_discharge[LFR>=0.1*mean_annual_discharge]
    HFR[LFR>=0.2*mean_annual_discharge] <- HFR_LFR_20_30  * mean_annual_discharge[LFR>=0.2*mean_annual_discharge]
    HFR[LFR>=0.3*mean_annual_discharge] <- HFR_LFR_more30 * mean_annual_discharge[LFR>=0.3*mean_annual_discharge]
    HFR[mean_annual_discharge<=0]       <- 0

    EFR <- LFR+HFR

    # Reduce EFR to 50% of available water where it exceeds this threshold (according to Smakhtin 2004)
    EFR <- pmin(EFR, 0.5*mean_annual_discharge)
    EFR <- as.magpie(EFR)

    out=EFR
    description="Total EFR per year"

  return(list(
    x=out,
    weight=NULL,
    unit="mio. m^3",
    description=description,
    isocountries=FALSE))
}
