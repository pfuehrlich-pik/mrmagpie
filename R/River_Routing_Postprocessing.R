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
#time="spline"
#averaging_range=NULL
dof=4
harmonize_baseline="CRU_4"
#ref_year="y2015"


months <- c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec")

# Read in mdischarge
### Monthly Discharge (unit (after calcLPJmL): mio. m^3/month)
monthly_discharge_magpie <- calcOutput("LPJmL", version="LPJmL4", climatetype=climatetype, subtype="mdischarge", aggregate=FALSE,
                                       harmonize_baseline=FALSE,
                                       time="raw")
yearly_discharge_magpie <- dimSums(monthly_discharge_magpie, dim=3)
#monthly_discharge_magpie <- as.array(collapseNames(monthly_discharge_magpie))

# Read in mrunoff
### Monthly Runoff (unit (after calcLPJmL): mio. m^3/month)
monthly_runoff_magpie    <- calcOutput("LPJmL", version="LPJmL4", climatetype=climatetype, subtype="mrunoff", aggregate=FALSE,
                                       harmonize_baseline=FALSE,
                                       time="raw")
yearly_runoff_magpie <- dimSums(monthly_runoff_magpie, dim=3)
#monthly_runoff_magpie    <- as.array(collapseNames(monthly_runoff_magpie))

# Read in EFRs
# (NOTE: The calculation currently still includes avl_water from calcAvlWater
# --> needs to be replaced by calculation that is only based on discharge!!!)
# Also: splined (due to current form of function)
EFR_magpie <- calcOutput("EnvmtlFlow", version="LPJmL4", climatetype=climatetype, harmonize_baseline=FALSE, time="spline", dof=dof, aggregate=FALSE, seasonality="monthly")
EFR_magpie <- EFR_magpie[,"y2100",,invert=TRUE]


# Read in mevap
### Monthly lake evapotranspiration (unit (after calcLPJmL): mio. m^3/ha)
# Will be imported from LPJmL at later stage
# For now: place-holder variable (set to 0)
monthly_evap_magpie     <- monthly_runoff_magpie
monthly_evap_magpie[,,] <- 0
getNames(monthly_evap_magpie) <- months

# Read in non-agricultural water consumption (mio. m^3/yr)
nonag_wc_magpie <- calcOutput("NonAgWaterDemand", source="WATERGAP2020", seasonality="total", waterusetype="consumption", time="raw", aggregate=FALSE)


common_yrs <- intersect(getYears(monthly_runoff_magpie),getYears(nonag_wc_magpie))
nonag_wc_magpie <- nonag_wc_magpie[,common_yrs,]


### River Routing
# distribute available water across the river basin

### Outcome: available water per cell (optionally: per month, per year, per grper)




### Rivers from Jens

rm(list=ls(all=TRUE))
gc()

NCELL <- 67420
NCRUCELL <- 67420



zz <- file("/data/biosx/LPJ/input_longheader/drainage.bin","rb")
seek(zz,where=43,origin="start")
x <- readBin(zz, integer(), n=2*NCRUCELL, size=4)
nextcell <- x[c(1:NCRUCELL)*2-1]
dist <- x[c(1:NCRUCELL)*2]
close(zz)

nextcell[which(nextcell<0)] <- -1
nextcell[which(nextcell>=0)] <- nextcell[which(nextcell>=0)] + 1

# determine downstream cell list
dummy <- array(data=-9999,dim=c(NCRUCELL))
c <- 1
i <- 1
dummy[i] <- nextcell[c]
while(dummy[i]>0)
{
  i <- i + 1
  dummy[i] <- nextcell[dummy[i-1]]
}
dsclist <- list(dummy[0:(i-1)])

for(c in 2:NCRUCELL)
{
  dummy[] <- -9999
  i <- 1
  dummy[i] <- nextcell[c]
  while(dummy[i]>0)
  {
    i <- i + 1
    dummy[i] <- nextcell[dummy[i-1]]
  }
  dsclist[[length(dsclist)+1]] <- dummy[0:(i-1)]
}

# determine endcell and dictance to endcell
cellstoend <- array(data=0,dim=c(NCRUCELL))
endcell <- array(data=0,dim=c(NCRUCELL))
for(i in 1:NCRUCELL)
{
  endcell[i] <- i
  while(nextcell[endcell[i]]>0)
  {
    endcell[i] <- nextcell[endcell[i]]
    cellstoend[i] <- cellstoend[i] + 1
  }
}

# determine calcorder
basinids <- unique(endcell)
calcorder <- array(data=0,dim=c(NCRUCELL))
for(b in 1:length(basinids))
{
  basincells <- which(endcell==basinids[b])
  calcorder[basincells] <- (cellstoend[basincells] - max(cellstoend[basincells]) - 1)*(-1)
}

# determine upstream cell list
usclist <- list()
for(c in 1:NCRUCELL)
{
  if(c%%1000 == 0) print(c)
  basincells <- which(endcell==endcell[c])
  dummy <- numeric()
  for(cell in basincells)
  {
    if(is.element(c,dsclist[[cell]])) dummy <- c(dummy,cell)
  }
  usclist[[c]] <- dummy
}

#save(nextcell,dsclist,usclist,endcell,calcorder,file="/data/open/Jens/WaterCC/analysis/river_routing.RData")


# River basin information
basin_code <- toolGetMapping("rivermapping.csv",type="cell")
basin_code <- basin_code$basincode







#####################
###### STEP 1.2 #####
#####################


#########################################################
# River Routing with non-agricultural water withdrawals #
#########################################################


# Read in non-agricultural water withdrawal (mio. m^3/yr)
nonag_ww_magpie <- calcOutput("NonAgWaterDemand", source="WATERGAP2020", seasonality="total", waterusetype="withdrawal", time="raw", aggregate=FALSE)










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


