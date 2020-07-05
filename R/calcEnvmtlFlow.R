#' @title calcEnvmtlFlow
#' @description This function calculates environmental flow requirements (EFR) for MAgPIE based on LPJmL monthly discharge following Smakthin et al. (2004)
#'
#' @param version Switch between LPJmL4 and LPJmL5
#' @param climatetype Switch between different climate scenarios (default: "CRU_4")
#' @param cells Number of cells to be reported: lpjcell (67420, default) or magpiecell (59199)
#' @param LFR_val Strictness of environmental flow requirements
#' @param HFR_LFR_less10 High flow requirements (share of total water for cells) with LFR<10percent of total water
#' @param HFR_LFR_10_20 High flow requirements (share of total water for cells) with 10percent < LFR < 20percent of total water
#' @param HFR_LFR_20_30 High flow requirements (share of total water for cells) with 20percent < LFR < 30percent of total water
#' @param HFR_LFR_more30 High flow requirements (share of total water for cells) with LFR>30percent of total water
#' @param EFRyears Long-term reference time frame for EFR calculation
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

calcEnvmtlFlow <- function(version="LPJmL4", climatetype="HadGEM2_ES:rcp2p6:co2", cells="lpjcell",
                           LFR_val=0.1, HFR_LFR_less10=0.2, HFR_LFR_10_20=0.15, HFR_LFR_20_30=0.07, HFR_LFR_more30=0.00,
                           EFRyears=c(1985:2015)){

  # Long-term reference period for EFR calculation
  EFRyears <- paste0("y",EFRyears)

  ### Monthly Discharge from LPJmL
  monthly_discharge_magpie <- calcOutput("LPJmL", selectyears=EFRyears, version=version, climatetype=climatetype, subtype="mdischarge_lpjcell", aggregate=FALSE,
                                         harmonize_baseline=FALSE, time="raw")

  # Extract years
  years <- getYears(monthly_discharge_magpie, as.integer = TRUE)
  # Transform to array (faster calculation)
  monthly_discharge_magpie <-  as.array(collapseNames(monthly_discharge_magpie))

  ### Calculate LFRs
  ## Note: LFRs correspond to the Q90-value (i.e. to the discharge that is exceeded in nine out of ten months)
  ## (Bonsch et al. 2015). This is calculated via the 10%-quantile of monthly discharge.

  # Get the monthly LFR_val quantile for all cells (across selected long-term reference time period)
  LFR_quant <- apply(monthly_discharge_magpie, MARGIN=c(1), quantile, probs=LFR_val)

  # Yearly LFRs
  LFR <- LFR_quant*12

  ### Mean annual discharge
  mean_annual_discharge <- apply(monthly_discharge_magpie, MARGIN=c(1), sum)/length(years)

  ### Calculate HFR
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

  ### EFR
  EFR <- LFR+HFR

  # Reduce EFR to 50% of available water where it exceeds this threshold (according to Smakhtin 2004)
  EFR <- pmin(EFR, 0.5*mean_annual_discharge)

  ### Correct number of cells and transform to magpie object
  if (cells=="lpjcell"){
    EFR <- as.magpie(EFR, spatial=1)
  } else if (cells=="magpiecell"){
    EFR <- EFR[magclassdata$cellbelongings$LPJ_input.Index]
    EFR <- as.magpie(EFR, spatial=1)
    dimnames(EFR)[[1]] <- paste(magclassdata$half_deg$region,1:59199,sep='.')
  } else {
    stop("Cell argument not supported. Select lpjcell for 67420 cells or magpiecell for 59199 cells")
  }


  # Check for NAs
  if(any(is.na(EFR))){
    stop("produced NA EFR")
  }

  return(list(
    x=EFR,
    weight=NULL,
    unit="mio. m^3",
    description="Total EFR per year",
    isocountries=FALSE))
}
