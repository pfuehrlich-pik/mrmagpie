library(mrcommons)
library(mrMAgPIE)


#####################
###### STEP 1 #######
#####################
###### STEP 1.1 #####
#####################

###############################################################
# River Routing with EFR + Non-Agricultural Water Consumption #
###############################################################
# Settings
climatetype <- "HadGEM2_ES:rcp2p6:co2"
time="spline"
averaging_range=NULL
dof=4
harmonize_baseline="CRU_4"
ref_year="y2015"


# Read in mdischarge
### Monthly Discharge (unit (after calcLPJmL): mio. m^3/month)
monthly_discharge_magpie <- calcOutput("LPJmL", version="LPJmL4", climatetype=climatetype, subtype="mdischarge", aggregate=FALSE,
                                       harmonize_baseline=FALSE,
                                       time="raw")
#monthly_discharge_magpie <- as.array(collapseNames(monthly_discharge_magpie))

# Read in mrunoff
### Monthly Runoff (unit (after calcLPJmL): mio. m^3/month)
monthly_runoff_magpie    <- calcOutput("LPJmL", version="LPJmL4", climatetype=climatetype, subtype="mrunoff", aggregate=FALSE,
                                       harmonize_baseline=FALSE,
                                       time="raw")
#monthly_runoff_magpie    <- as.array(collapseNames(monthly_runoff_magpie))

# Read in mevap
### Monthly lake evapotranspiration (unit (after calcLPJmL): m^3/ha)
monthly_evap_magpie    <- calcOutput("LPJmL", version="LPJmL4", climatetype=climatetype, subtype="mevaporation", aggregate=FALSE,
                                     harmonize_baseline=FALSE,
                                     time="raw")

# Read in non-agricultural water consumption
nonag_wc_magpie <- calcOutput("NonAgWaterDemand", source="WATERGAP2020", seasonality="total", waterusetype="consumption", time="raw", calibration_approach="harmonizefunction", aggregate=FALSE)

# Read in EFRs
# (NOTE: The calculation currently still includes avl_water from calcAvlWater
# --> needs to be replaced by calculation that is only based on discharge!!!)
EFR_magpie <- calcOutput("EnvmtlFlow", version="LPJmL4", climatetype=climatetype, harmonize_baseline=harmonize_baseline, ref_year=ref_year, time="spline", dof=dof, aggregate=FALSE, seasonality="monthly")


### River Routing
# distribute available water across the river basin

### Outcome: available water per cell (optionally: per month, per year, per grper)


#####################
###### STEP 1.2 #####
#####################


#########################################################
# River Routing with non-agricultural water withdrawals #
#########################################################


# Read in non-agricultural water withdrawal
nonag_ww_magpie <- calcOutput("NonAgWaterDemand", source="WATERGAP2020", seasonality="total", waterusetype="withdrawal", time="raw", calibration_approach="harmonizefunction", aggregate=FALSE)



#####################
###### STEP 2 #######
#####################

#########################################
# Including committed agricultural uses #
#########################################

##### Notes:
### Preconditions
## Agricultural water use efficiency calculation
# switch from airrig to blue water consumption
# take non-consumptive and consumptive water losses into account
# NECESSARY DATA:
# from LPJmL: consumptive blue water use
# from literature (e.g. Jonas JÃ¤germeyr?): values for system efficiencies (later: from LPJmL)
# from MAgPIE: land use in initialization year (irrigated agriculture by croptype)


