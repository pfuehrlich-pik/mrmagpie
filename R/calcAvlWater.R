#' @title calcAvlWater
#' @description This function calculates water availability for MAgPIE retrieved from LPJmL using a river routing algorithm for distribution of discharge across the river basin
#'
#' @param version Switch between LPJmL4 and LPJmL5
#' @param climatetype Switch between different climate scenarios (default: "CRU_4")
#' @param time Time smoothing: average, spline or raw (default)
#' @param averaging_range only specify if time=="average": number of time steps to average
#' @param dof             only specify if time=="spline": degrees of freedom needed for spline
#' @param harmonize_baseline FALSE (default): no harmonization, TRUE: if a baseline is specified here data is harmonized to that baseline (from ref_year on)
#' @param ref_year Reference year for harmonization baseline (just specify when harmonize_baseline=TRUE)
#' @param selectyears Years to be returned
#'
#' @import magclass
#' @import madrat
#'
#' @return magpie object in cellular resolution
#' @author Jens Heinke, Felicitas Beier
#'
#' @examples
#' \dontrun{ calcOutput("AvlWater", aggregate = FALSE) }
#'

calcAvlWater <- function(selectyears="all",
                         version="LPJmL4", climatetype="HadGEM2_ES:rcp2p6:co2", time="raw", averaging_range=NULL, dof=NULL,
                         harmonize_baseline=FALSE, ref_year="y2015"){
  load("C:/Users/beier/Documents/doktorarbeit/MAgPIE_Water/River_Routing_Postprocessing/river_routing_stn.RData")
  load("C:/Users/beier/Documents/doktorarbeit/MAgPIE_Water/River_Routing_Postprocessing/cells_magpie2lpj.Rda")




  ## Required inputs:
  # Yearly runoff (mio. m^3 / yr) [smoothed]
  yearly_runoff_magpie <- calcOutput("LPJmL", version="LPJmL4", climatetype=climatetype, subtype="runoff", aggregate=FALSE,
                                        harmonize_baseline=FALSE, time=time, dof=dof, averaging_range=averaging_range)
  yearly_runoff_magpie  <- as.array(collapseNames(yearly_runoff_magpie))
  years <- getYears(monthly_runoff_magpie)

  # Environmental Flow Requirements (in mio. m^3 / yr) [long-term average]
  EFR_magpie <- calcOutput("EnvmtlFlow", version="LPJmL4", climatetype=climatetype, aggregate=FALSE,
             LFR_val=0.1, HFR_LFR_less10=0.2, HFR_LFR_10_20=0.15, HFR_LFR_20_30=0.07, HFR_LFR_more30=0.00,
             EFRyears=c(1985:2015))
  EFR_magpie <- as.array(collapseNames(EFR_magpie))

  # Yearly lake evapotranspiration (in mio. m^3/ha) [place holder]
  lake_evap_magpie <- new.magpie(1:59199,years)
  lake_evap_magpie[,,] <- 0

  # Non-Agricultural Water Withdrawals (in mio. m^3 / yr) [smoothed]
  NAg_ww_magpie <- calcOutput("NonAgWaterDemand", source="WATERGAP2020", time=time, dof=dof, averaging_range=averaging_range, waterusetype="withdrawal", aggregate=FALSE)
  NAg_ww_magpie <- as.array(collapseNames(NAg_ww_magpie))

  # Non-Agricultural Water Consumption (in mio. m^3 / yr) [smoothed]
  NAg_wc_magpie <- calcOutput("NonAgWaterDemand", source="WATERGAP2020", time=time, dof=dof, averaging_range=averaging_range, waterusetype="consumption", aggregate=FALSE)
  NAg_wc_magpie <- as.array(collapseNames(NAg_wc_magpie))

  # Committed agricultural uses (in mio. m^3 / yr) [for initialization year]
  CAD_magpie <- calcOutput("CommittedAgWaterUse",aggregate=FALSE)
  CAD_magpie <- as.array(collapseNames(CAD_magpie))


  #### River routing

  EFR <- numeric(length(calcorder))
  EFR[magpie2lpj] <- EFR_magpie

  CAD <- numeric(length(calcorder))
  CAD[magpie2lpj] <- CAD_magpie

  CAD_w <- CAD/0.4
  CAD_c <- CAD

  NCELLS <- length(calcorder)



  for (y in "y1995"){
    runoff <- numeric(length(calcorder))
    runoff[magpie2lpj] <- yearly_runoff_magpie[,y,]

    lake_evap <- numeric(length(calcorder))
    lake_evap[magpie2lpj] <- lake_evap_magpie[,y,]

    ### River Routing ###
    # Naturalized discharge
    discharge_nat <- numeric(NCELLS)
    test <- numeric(NCELLS)

    for (c in 1:NCELLS){
      discharge_nat[c] <- sum(runoff[c(usclist[[c]],c)]) - sum(lake_evap[c(usclist[[c]],c)])
    }


    for (scenario in "ssp2"){
      NAg_wc <- numeric(length(calcorder))
      NAg_wc[magpie2lpj] <- NAg_wc_magpie[,y,scenario]

      NAg_ww <- numeric(length(calcorder))
      NAg_ww[magpie2lpj] <- NAg_ww_magpie[,y,scenario]

      ### River Routing Algorithm ###


    }
  }





  ### Rivers from Jens

 # rm(list=ls(all=TRUE))
  #gc()
#
#   NCELL <- 67420
#   NCRUCELL <- 67420
#
#
#   zz <- file("C:/Users/beier/Documents/doktorarbeit/MAgPIE_Water/River_Routing_Postprocessing/drainagestn.bin","rb")
#   seek(zz,where=43,origin="start")
#   x <- readBin(zz, integer(), n=2*NCRUCELL, size=4)
#   nextcell <- x[c(1:NCRUCELL)*2-1]
#   dist <- x[c(1:NCRUCELL)*2]
#   close(zz)
#
#   nextcell[which(nextcell<0)] <- -1
#   nextcell[which(nextcell>=0)] <- nextcell[which(nextcell>=0)] + 1
#
#   # determine downstream cell list
#   dummy <- array(data=-9999,dim=c(NCRUCELL))
#   c <- 1
#   i <- 1
#   dummy[i] <- nextcell[c]
#   while(dummy[i]>0)
#   {
#    i <- i + 1
#    dummy[i] <- nextcell[dummy[i-1]]
#   }
#   dsclist <- list(dummy[0:(i-1)])
#
#   for(c in 2:NCRUCELL)
#   {
#    dummy[] <- -9999
#    i <- 1
#    dummy[i] <- nextcell[c]
#    while(dummy[i]>0)
#    {
#      i <- i + 1
#      dummy[i] <- nextcell[dummy[i-1]]
#    }
#    dsclist[[length(dsclist)+1]] <- dummy[0:(i-1)]
#   }
#
#   # determine endcell and distance to endcell
#   cellstoend <- array(data=0,dim=c(NCRUCELL))
#   endcell <- array(data=0,dim=c(NCRUCELL))
#   for(i in 1:NCRUCELL)
#   {
#    endcell[i] <- i
#    while(nextcell[endcell[i]]>0)
#    {
#      endcell[i] <- nextcell[endcell[i]]
#      cellstoend[i] <- cellstoend[i] + 1
#    }
#   }
#
#   # determine calcorder
#   basinids <- unique(endcell)
#   calcorder <- array(data=0,dim=c(NCRUCELL))
#   for(b in 1:length(basinids))
#   {
#     basincells <- which(endcell==basinids[b])
#     calcorder[basincells] <- (cellstoend[basincells] - max(cellstoend[basincells]) - 1)*(-1)
#   }
#
#   # determine upstream cell list
#   usclist <- list()
#   for(c in 1:NCRUCELL)
#   {
#     if(c%%1000 == 0) print(c)
#     basincells <- which(endcell==endcell[c])
#     dummy <- numeric()
#     for(cell in basincells)
#     {
#       if(is.element(c,dsclist[[cell]])) dummy <- c(dummy,cell)
#     }
#     usclist[[c]] <- dummy
#   }
#
#   #save(nextcell,dsclist,usclist,endcell,calcorder,file="/data/open/Jens/WaterCC/analysis/river_routing.RData")
#
#
# ################# River Routing
#
#
#
#
#
#









#####################################################
  #### Old function

  if(harmonize_baseline==FALSE){

    if(time=="raw"){

      ### Monthly Discharge (unit (after calcLPJmL): mio. m^3/month)
      monthly_discharge_magpie <- calcOutput("LPJmL", version=version, climatetype=climatetype, subtype="mdischarge", aggregate=FALSE,
                                             harmonize_baseline=FALSE,
                                             time="raw")
      # Transform to array (faster calculation)
      monthly_discharge_magpie <- as.array(collapseNames(monthly_discharge_magpie))

      ### Monthly Runoff (unit (after calcLPJmL): mio. m^3/month)
      monthly_runoff_magpie    <- calcOutput("LPJmL", version=version, climatetype=climatetype, subtype="mrunoff", aggregate=FALSE,
                                             harmonize_baseline=FALSE,
                                             time="raw")
      # Transform to array (faster calculation)
      monthly_runoff_magpie    <- as.array(collapseNames(monthly_runoff_magpie))

      ### Calculate available water per month (avl_water_month)
      # Empty array
      avl_water_month     <- monthly_runoff_magpie
      avl_water_month[,,] <- NA

      ## River basin water allocation algorithm:
      # River basin information
      basin_code <- toolGetMapping("rivermapping.csv",type="cell")
      basin_code <- basin_code$basincode

      # Sum the runoff in all basins and allocate it to the basin cells with discharge as weight
      for(basin in unique(basin_code)){
        basin_cells     <- which(basin_code==basin)
        basin_runoff    <- colSums(monthly_runoff_magpie[basin_cells,,,drop=FALSE])
        basin_discharge <- colSums(monthly_discharge_magpie[basin_cells,,,drop=FALSE])
        for(month in dimnames(avl_water_month)[[3]]){
          avl_water_month[basin_cells,,month] <- t(basin_runoff[,month]*t(monthly_discharge_magpie[basin_cells,,month])/basin_discharge[,month])
        }
      }
      # Remove no longer needed objects
      rm(basin_discharge,basin_runoff)

      # avl_water_month contain NA's wherever basin_discharge was 0 -> Replace NA's by 0
      avl_water_month[is.nan(avl_water_month)] <- 0
      avl_water_month <- as.magpie(avl_water_month)

    } else {
      # Time smoothing:
      x     <- calcOutput("AvlWater", version=version, climatetype=climatetype, seasonality="monthly", aggregate=FALSE,
                          harmonize_baseline=FALSE, time="raw")

      if(time=="average"){

        # Smoothing data through average:
        avl_water_month <- toolTimeAverage(x, averaging_range=averaging_range)

      } else if(time=="spline"){

        # Smoothing data with spline method:
        avl_water_month <- toolTimeSpline(x, dof=dof)
        # Replace value in 2100 with value from 2099 (LPJmL output ends in 2099)
        if ("y2099" %in% getYears(avl_water_month)) {
          avl_water_month <- toolFillYears(avl_water_month, c(getYears(avl_water_month, as.integer=TRUE)[1]:2100))
        }

      } else if(time!="raw"){
        stop("Time argument not supported!")
      }
    }

  } else {

    if(time=="raw"){
      stop("Harmonization with raw data not possible. Select time='spline' when applying harmonize_baseline=TRUE")
    } else {
      # Load smoothed data
      baseline <- calcOutput("AvlWater", version=version, climatetype=harmonize_baseline, seasonality="monthly", aggregate=FALSE,
                             harmonize_baseline=FALSE, time=time, dof=dof, averaging_range=averaging_range)
      x        <- calcOutput("AvlWater", version=version, climatetype=climatetype, seasonality="monthly", aggregate=FALSE,
                             harmonize_baseline=FALSE, time=time, dof=dof, averaging_range=averaging_range)
      # Harmonize to baseline
      avl_water_month <- toolHarmonize2Baseline(x=x, base=baseline, ref_year=ref_year, limited=TRUE, hard_cut=FALSE)
    }
  }

  if(selectyears!="all"){
    years           <- sort(findset(selectyears,noset = "original"))
    avl_water_month <- avl_water_month[,years,]
  }

  ###########################################
  ######### RETURN FUNCTION OUTPUT ##########
  ###########################################

  ### Available water per cell per month
  if(seasonality=="monthly"){
    # Check for NAs
    if(any(is.na(avl_water_month))){
      stop("produced NA water availability")
    }
    out=avl_water_month
    description="Available water per cell per month (based on runoff and discharge from LPJmL)"
  }

  ### Total water available per cell per year
  if(seasonality=="total"){
    # Sum up over all month:
    avl_water_total <- dimSums(avl_water_month, dim=3)
    # Check for NAs
    if(any(is.na(avl_water_total))){
      stop("produced NA water availability")
    }
    out=avl_water_total
    description="Total available water per year"
  }

  ### Water available in growing period per cell per year
  if(seasonality=="grper"){
    # magpie object with days per month with same dimension as avl_water_month
    tmp <- c(31,28,31,30,31,30,31,31,30,31,30,31)
    month_days     <- new.magpie(names=dimnames(avl_water_month)[[3]])
    month_days[,,] <- tmp
    month_day_magpie     <- as.magpie(avl_water_month)
    month_day_magpie[,,] <- 1
    month_day_magpie     <- month_day_magpie * month_days

    # Daily water availability
    avl_water_day <- avl_water_month/month_day_magpie

    # Growing days per month
    grow_days <- calcOutput("GrowingPeriod", version="LPJmL5", climatetype=climatetype, time=time, dof=dof, averaging_range=averaging_range,
                            harmonize_baseline=harmonize_baseline, ref_year=ref_year, yield_ratio=0.1, aggregate=FALSE)

    # Adjust years
    years_wat <- getYears(avl_water_day)
    years_grper  <- getYears(grow_days)
    if(length(years_wat)>=length(years_grper)){
      years <- years_grper
    } else {
      years <- years_wat
    }
    rm(years_grper, years_wat)

    # Available water in growing period per month
    avl_water_grper <- avl_water_day[,years,]*grow_days[,years,]
    # Available water in growing period per year
    avl_water_grper <- dimSums(avl_water_grper, dim=3)

    # Check for NAs
    if(any(is.na(avl_water_grper))){
      stop("produced NA water availability")
    }
    out=avl_water_grper
    description="Available water in growing period per year"
  }

  return(list(
    x=out,
    weight=NULL,
    unit="mio. m^3",
    description=description,
    isocountries=FALSE))
}
