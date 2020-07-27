#' @title calcAvlWater
#' @description This function calculates water availability for MAgPIE retrieved from LPJmL using a river routing algorithm for distribution of discharge across the river basin
#'
#' @param version Switch between LPJmL4 and LPJmL5
#' @param climatetype Switch between different climate scenarios (default: "CRU_4")
#' @param time            Time smoothing: average, spline or raw (default)
#' @param averaging_range only specify if time=="average": number of time steps to average
#' @param dof             only specify if time=="spline": degrees of freedom needed for spline
#' @param harmonize_baseline FALSE (default): no harmonization, TRUE: if a baseline is specified here data is harmonized to that baseline (from ref_year on)
#' @param ref_year           Reference year for harmonization baseline (just specify when harmonize_baseline=TRUE)
#' @param selectyears Years to be returned
#' @param EFR Environmental flow requirements activated (TRUE) or not (FALSE)
#'
#' @import magclass
#' @import madrat
#' @importFrom mrcommons toolHarmonize2Baseline
#'
#' @return magpie object in cellular resolution
#' @author Jens Heinke, Felicitas Beier
#'
#' @examples
#' \dontrun{ calcOutput("AvlWater", aggregate = FALSE) }
#'

calcAvlWater <- function(selectyears="all",
                         version="LPJmL4", climatetype="HadGEM2_ES:rcp2p6:co2", time="raw", averaging_range=NULL, dof=NULL,
                         harmonize_baseline=FALSE, ref_year="y2015", EFR=TRUE){

  #############################
  ####### Read in Data ########
  #############################

  ### Read in river structure
  # Note: river structure derived from LPJmL input (drainage) [Later: implement readDrainage function]
  data <- toolGetMapping("River_structure.rda",type="cell")
  for (i in 1:length(data)){
    assign(paste(names(data[[i]])), data[[i]][[1]])
  }
  rm(data,i)

  # Number of cells to be used for calculation
  NCELLS <- length(calcorder)

  ### LPJ-MAgPIE cell mapping
  #load("C:/Users/beier/Documents/doktorarbeit/MAgPIE_Water/River_Routing_Postprocessing/cells_magpie2lpj.Rda")
  ### Question: Use LPJ_input.Index or LPJ.Index? What is the difference?
  magpie2lpj    <- magclassdata$cellbelongings$LPJ_input.Index
  lpj_cells_map <- toolGetMapping("LPJ_CellBelongingsToCountries.csv", type="cell")

  ### Required inputs:
  # Yearly runoff (mio. m^3 / yr) [smoothed]
  yearly_runoff <- calcOutput("LPJmL", version="LPJmL4", climatetype=climatetype, subtype="runoff_lpjcell", aggregate=FALSE,
                                        harmonize_baseline=FALSE, time="spline", dof=4, averaging_range=NULL)
  yearly_runoff <- as.array(collapseNames(yearly_runoff))
  years <- getYears(yearly_runoff)

  # Environmental Flow Requirements (in mio. m^3 / yr) [long-term average]
  if (EFR==TRUE){
    EFR_magpie <- calcOutput("EnvmtlFlow", version="LPJmL4", climatetype=climatetype, aggregate=FALSE, cells="lpjcell",
                      LFR_val=0.1, HFR_LFR_less10=0.2, HFR_LFR_10_20=0.15, HFR_LFR_20_30=0.07, HFR_LFR_more30=0.00,
                      EFRyears=c(1985:2015))
  } else if (EFR==FALSE){
    EFR_magpie <- new.magpie(1:NCELLS,fill=0)
    getCells(EFR_magpie) <- paste(lpj_cells_map$ISO,1:67420,sep=".")
  } else {
    stop("Specify whether environmental flows are activated or not via argument EFR")
  }
  EFR_magpie <- as.array(collapseNames(EFR_magpie))

  # Yearly lake evapotranspiration (in mio. m^3/ha) [place holder]
  lake_evap     <- new.magpie(1:NCELLS,years)
  lake_evap[,,] <- 0
  lake_evap     <- as.array(collapseNames(lake_evap))

  # Non-Agricultural Water Withdrawals (in mio. m^3 / yr) [smoothed]
  NAg_ww_magpie           <- calcOutput("NonAgWaterDemand", source="WATERGAP2020", time="spline", dof=4, averaging_range=NULL, waterusetype="withdrawal", seasonality="total", aggregate=FALSE)
  getCells(NAg_ww_magpie) <- paste("GLO",magclassdata$cellbelongings$LPJ_input.Index,sep=".")
  NAg_ww           <- new.magpie(1:NCELLS,getYears(NAg_ww_magpie),getNames(NAg_ww_magpie))
  NAg_ww[,,]       <- 0
  NAg_ww[paste("GLO",magclassdata$cellbelongings$LPJ_input.Index,sep="."),,] <- NAg_ww_magpie[,,]
  getCells(NAg_ww) <- paste(lpj_cells_map$ISO,1:67420,sep=".")
  NAg_ww           <- as.array(collapseNames(NAg_ww))
  rm(NAg_ww_magpie)

  # Non-Agricultural Water Consumption (in mio. m^3 / yr) [smoothed]
  NAg_wc_magpie           <- calcOutput("NonAgWaterDemand", source="WATERGAP2020", time="spline", dof=4, averaging_range=NULL, waterusetype="consumption", seasonality="total", aggregate=FALSE)
  getCells(NAg_wc_magpie) <- paste("GLO",magclassdata$cellbelongings$LPJ_input.Index,sep=".")
  NAg_wc           <- new.magpie(1:NCELLS,getYears(NAg_wc_magpie),getNames(NAg_wc_magpie))
  NAg_wc[,,]       <- 0
  NAg_wc[paste("GLO",magclassdata$cellbelongings$LPJ_input.Index,sep="."),,] <- NAg_wc_magpie[,,]
  getCells(NAg_wc) <- paste(lpj_cells_map$ISO,1:67420,sep=".")
  NAg_wc           <- as.array(collapseNames(NAg_wc))
  rm(NAg_wc_magpie)

  # Committed agricultural uses (in mio. m^3 / yr) [for initialization year]
  CAU_magpie <- calcOutput("CommittedAgWaterUse",iniyear=1995,irrigini="Jaegermeyr_lpjcell",time="raw",dof=NULL,aggregate=FALSE)
  CAW_magpie <- as.array(collapseNames(dimSums(CAU_magpie[,,"withdrawal"],dim=3)))
  CAC_magpie <- as.array(collapseNames(dimSums(CAU_magpie[,,"consumption"],dim=3)))
  rm(CAU_magpie)


  #############################
  ####### River routing #######
  #############################

  for (y in "y1995"){
    # Naturalized discharge
    discharge_nat <- numeric(NCELLS)
    # Discharge considering human uses
    discharge     <- numeric(NCELLS)
    # Actual withdrawals considering availability
    actual_withdrawal <- numeric(NCELLS)
    # Water available in cell
    avl_wat_nat   <- numeric(NCELLS) # naturally available water
    avl_ag_wat    <- numeric(NCELLS) # water available for agricultural withdrawal
    avl_ag_cons   <- numeric(NCELLS) # water available for agricultural consumption
    # Water not available for consumption
    NAC_water     <- numeric(NCELLS)
    # Water requirement from current cell for downstreamcell
    wat_dem_exceeding_runoff       <- numeric(NCELLS)
    shr_downstreamcell_requirement <- numeric(NCELLS)
    # Water reserved for particular cell
    res_water[c] <- numeric(NCELLS)

    for (c in 1:NCELLS){
      # nat.discharge  = (inflow from upstream + runoff on cell) -  lake evaporation (of cell and upstream)
      discharge_nat[c] <- sum(runoff[c(upstreamcells[[c]],c)])   -  sum(lake_evap[c(upstreamcells[[c]],c)])
    }

    ### River Routing 1: Non-agricultural uses ###
    for (o in 1:max(calcorder)) {
      # Note: the calcorder ensures that the upstreamcells are calculated first
      cells <- which(calcorder==o)

      for (scen in "ssp2"){
        for (c in cells){
          ## Water withdrawals cannot exceed availability (considering EFRs)
          # av. water in cell = nat.discharge  - environmental flow requirements
          avl_wat_nat[c]      <- discharge_nat[c] - EFR[c]
          # actual withdrawal (non-agriculture) < av. water in cell for withdrawal
          actual_withdrawal[c] <- avl_wat_nat[c] - NAg_ww[c,y,scen]
          print(paste("Non-agricultural water withdrawals exceed availability in", length(which(is.na(actual_withdrawal))) ,"cells. Non-agricultural withdrawals and consumption reduced accordingly.",sep=" "))
          if (actual_withdrawal[c]<0) {
            NAg_ww[c,y,scen] <- NAg_ww[c,y,scen] + (actual_withdrawal[c]*(-1))
            NAg_wc[c,y,scen] <- NAg_wc[c,y,scen] + (actual_withdrawal[c]*(-1))*(NAg_wc[c,y,scen]/NAg_ww[c,y,scen])
          }

          ## Outflow from one cell to the next
          # (Subtract local water consumption in current cell (committed ag. and non-ag. consumption))
          discharge[c] <- discharge_nat[c] - NAg_wc[c,y,scen]

          ## Water withdrawal accounting
          # (Water withdrawn downstream can be withdrawn upstream, but not consumed)
          # NAC_c = ((withdrawal-runoff)/sum(Inflow))_(c+1) * outflow_c
          # Note: which(nextcell==c): cells that go into current cell

          # Water demand (withdrawal) in current cell that cannot be fulfilled by runoff on that cell, i.e. that needs to come from upstream cell:
          wat_dem_exceeding_runoff[c] <- NAg_ww[c,y,scen] - yearly_runoff[c]
         ### wat_dem_exceeding_runoff[c] <- (NAg_ww[c,y,scen] + CAW_magpie[c]) - yearly_runoff[c] ###?????

          # Share of water that is needed for downstream withdrawal (cannot be consumed in current cell)
          for (i in 1:length(which(nextcell==c))){
            # water requirement from
            shr_downstreamcell_requirement[which(nextcell==c)[i]] <- wat_dem_exceeding_runoff[c]/sum(discharge[which(nextcell==c)])
          }
          # Water that is needed for downstream withdrawal (not available for consumption in current cell)
          shr_downstreamcell_requirement[c] <- max(shr_downstreamcell_requirement[c])
          NAC_water[c] <- shr_downstreamcell_requirement[c==nextcell] * discharge[c]

          # Water available for agricultural consumption in current cell
          # avl.wat.cons. = avl.water     - non-ag.consumption  - not-available for consumption - envmtl. flows
          avl_ag_cons[c]  <- avl_wat_nat[c] - NAg_wc[c,y,scen] - NAC_water[c] - EFR[c]
          # Water available for agricultural withdrawal in current cell
          avl_ag_wat[c]   <- avl_wat_nat[c] - NAg_wc[c,y,scen]

          # Water reserved for each cell (from above calculations):
          res_water[c] <- lake_evap[c] + EFR_magpie[c] + NAg_wc[c,y,scen] + NAC_water[c]
        }
      }
    }

    ### River Routing 2: Committed agricultural uses ###

    # Actual agricultural withdrawals considering availability and non-agricultural consumption
    actual_withdrawal_ag <- numeric(NCELLS)

    for (o in 1:max(calcorder)) {
      cells <- which(calcorder==o)

      for (scen in "ssp2"){
        for (c in cells){

          ## Agricultural water withdrawals cannot exceed availability (considering EFRs and non-agricultural uses)
          # actual withdrawal (agriculture) < av. water in cell for ag. withdrawal
          actual_withdrawal_ag[c] <- avl_ag_wat[c] - CAW_magpie[c]
          print(paste("Agricultural water withdrawals exceed availability in", length(which(is.na(actual_withdrawal_ag))) ,"cells. Agricultural withdrawals and consumption reduced accordingly.",sep=" "))
          if (actual_withdrawal_ag[c]<0) {
            CAW_magpie[c] <- CAW_magpie[c] + (actual_withdrawal_ag[c]*(-1))
            CAC_magpie[c] <- CAC_magpie[c] + (actual_withdrawal_ag[c]*(-1))*(CAC_magpie[c]/CAW_magpie[c])
          }

          ## Outflow from one cell to the next
          # (Subtract local water consumption in current cell (committed ag. consumption))
          discharge[c] <- discharge[c] - CAC_magpie[c]

          ## Water withdrawal accounting
          # (Water withdrawn downstream can be withdrawn upstream, but not consumed)
          # NAC_c = ((withdrawal-runoff)/sum(Inflow))_(c+1) * outflow_c
          # Note: which(nextcell==c): cells that go into current cell

          # Water demand (withdrawal) in current cell that cannot be fulfilled by runoff on that cell, i.e. that needs to come from upstream cell:
          wat_dem_exceeding_runoff[c] <- (NAg_ww[c,y,scen] + CAW_magpie[c]) - yearly_runoff[c]

          # Share of water that is needed for downstream withdrawal (cannot be consumed in current cell)
          for (i in 1:length(which(nextcell==c))) {
            # water requirement from
            shr_downstreamcell_requirement[which(nextcell==c)[i]] <- wat_dem_exceeding_runoff[c]/sum(discharge[which(nextcell==c)])
          }

          # Water that is needed for downstream withdrawal (not available for consumption in current cell)
          shr_downstreamcell_requirement[c] <- max(shr_downstreamcell_requirement[c])
          NAC_water[c] <- shr_downstreamcell_requirement[c==nextcell] * discharge[c]

          # Water reserved for each cell (from above calculations):
          res_water[c] <- lake_evap[c] + EFR_magpie[c] + NAg_wc[c,y,scen] + NAC_water[c]

        }
      }
    }

    ### Water allocation algorithm for "surplus water" across the river basin ###

    # River basin runoff to be distributed across cells of river basin by algorithm (tba)
    basin_runoff     <- dimSums(yearly_runoff,dim=1)
    res_water_basin  <- dimSums(res_water, dim=1)
    surlus_wat_basin <- basin_runoff + res_water_basin

    for (o in 1:max(calcorder)) {
      cells <- which(calcorder==o)

      for (scen in "ssp2"){
        for (c in cells){

          ## Discharge-weighted distribution
          if (algorithm=="discharge") {
            # Available water per cell
            avl_water[c] <- surplus_wat_basin * discharge[c]/sum(discharge[c]) + res_water[c]
          } else if (algorithm=="yieldimprovement") {
            ## Potential yield improvement maximization
            ### ( Not yet implemented ) ###
          } else if (algorithm=="potentialcropland") {
            ## Yield improvement threshold & potential cropland
            ### ( Not yet implemented ) ###
          }
        }
      }
    }

  }


#####################################################
  #### Old function

  # if(harmonize_baseline==FALSE){
  #
  #   if(time=="raw"){
  #
  #     ### Monthly Discharge (unit (after calcLPJmL): mio. m^3/month)
  #     monthly_discharge_magpie <- calcOutput("LPJmL", version=version, climatetype=climatetype, subtype="mdischarge", aggregate=FALSE,
  #                                            harmonize_baseline=FALSE,
  #                                            time="raw")
  #     # Transform to array (faster calculation)
  #     monthly_discharge_magpie <- as.array(collapseNames(monthly_discharge_magpie))
  #
  #     ### Monthly Runoff (unit (after calcLPJmL): mio. m^3/month)
  #     monthly_runoff_magpie    <- calcOutput("LPJmL", version=version, climatetype=climatetype, subtype="mrunoff", aggregate=FALSE,
  #                                            harmonize_baseline=FALSE,
  #                                            time="raw")
  #     # Transform to array (faster calculation)
  #     monthly_runoff_magpie    <- as.array(collapseNames(monthly_runoff_magpie))
  #
  #     ### Calculate available water per month (avl_water_month)
  #     # Empty array
  #     avl_water_month     <- monthly_runoff_magpie
  #     avl_water_month[,,] <- NA
  #
  #     ## River basin water allocation algorithm:
  #     # River basin information
  #     basin_code <- toolGetMapping("rivermapping.csv",type="cell")
  #     basin_code <- basin_code$basincode
  #
  #     # Sum the runoff in all basins and allocate it to the basin cells with discharge as weight
  #     for(basin in unique(basin_code)){
  #       basin_cells     <- which(basin_code==basin)
  #       basin_runoff    <- colSums(monthly_runoff_magpie[basin_cells,,,drop=FALSE])
  #       basin_discharge <- colSums(monthly_discharge_magpie[basin_cells,,,drop=FALSE])
  #       for(month in dimnames(avl_water_month)[[3]]){
  #         avl_water_month[basin_cells,,month] <- t(basin_runoff[,month]*t(monthly_discharge_magpie[basin_cells,,month])/basin_discharge[,month])
  #       }
  #     }
  #     # Remove no longer needed objects
  #     rm(basin_discharge,basin_runoff)
  #
  #     # avl_water_month contain NA's wherever basin_discharge was 0 -> Replace NA's by 0
  #     avl_water_month[is.nan(avl_water_month)] <- 0
  #     avl_water_month <- as.magpie(avl_water_month)
  #
  #   } else {
  #     # Time smoothing:
  #     x     <- calcOutput("AvlWater", version=version, climatetype=climatetype, seasonality="monthly", aggregate=FALSE,
  #                         harmonize_baseline=FALSE, time="raw")
  #
  #     if(time=="average"){
  #
  #       # Smoothing data through average:
  #       avl_water_month <- toolTimeAverage(x, averaging_range=averaging_range)
  #
  #     } else if(time=="spline"){
  #
  #       # Smoothing data with spline method:
  #       avl_water_month <- toolTimeSpline(x, dof=dof)
  #       # Replace value in 2100 with value from 2099 (LPJmL output ends in 2099)
  #       if ("y2099" %in% getYears(avl_water_month)) {
  #         avl_water_month <- toolFillYears(avl_water_month, c(getYears(avl_water_month, as.integer=TRUE)[1]:2100))
  #       }
  #
  #     } else if(time!="raw"){
  #       stop("Time argument not supported!")
  #     }
  #   }
  #
  # } else {
  #
  #   if(time=="raw"){
  #     stop("Harmonization with raw data not possible. Select time='spline' when applying harmonize_baseline=TRUE")
  #   } else {
  #     # Load smoothed data
  #     baseline <- calcOutput("AvlWater", version=version, climatetype=harmonize_baseline, seasonality="monthly", aggregate=FALSE,
  #                            harmonize_baseline=FALSE, time=time, dof=dof, averaging_range=averaging_range)
  #     x        <- calcOutput("AvlWater", version=version, climatetype=climatetype, seasonality="monthly", aggregate=FALSE,
  #                            harmonize_baseline=FALSE, time=time, dof=dof, averaging_range=averaging_range)
  #     # Harmonize to baseline
  #     avl_water_month <- toolHarmonize2Baseline(x=x, base=baseline, ref_year=ref_year, limited=TRUE, hard_cut=FALSE)
  #   }
  # }
  #
  # if(selectyears!="all"){
  #   years           <- sort(findset(selectyears,noset = "original"))
  #   avl_water_month <- avl_water_month[,years,]
  # }
  #
  # ###########################################
  # ######### RETURN FUNCTION OUTPUT ##########
  # ###########################################
  #
  # ### Available water per cell per month
  # if(seasonality=="monthly"){
  #   # Check for NAs
  #   if(any(is.na(avl_water_month))){
  #     stop("produced NA water availability")
  #   }
  #   out=avl_water_month
  #   description="Available water per cell per month (based on runoff and discharge from LPJmL)"
  # }
  #
  # ### Total water available per cell per year
  # if(seasonality=="total"){
  #   # Sum up over all month:
  #   avl_water_total <- dimSums(avl_water_month, dim=3)
  #   # Check for NAs
  #   if(any(is.na(avl_water_total))){
  #     stop("produced NA water availability")
  #   }
  #   out=avl_water_total
  #   description="Total available water per year"
  # }
  #
  # ### Water available in growing period per cell per year
  # if(seasonality=="grper"){
  #   # magpie object with days per month with same dimension as avl_water_month
  #   tmp <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  #   month_days     <- new.magpie(names=dimnames(avl_water_month)[[3]])
  #   month_days[,,] <- tmp
  #   month_day_magpie     <- as.magpie(avl_water_month)
  #   month_day_magpie[,,] <- 1
  #   month_day_magpie     <- month_day_magpie * month_days
  #
  #   # Daily water availability
  #   avl_water_day <- avl_water_month/month_day_magpie
  #
  #   # Growing days per month
  #   grow_days <- calcOutput("GrowingPeriod", version="LPJmL5", climatetype=climatetype, time=time, dof=dof, averaging_range=averaging_range,
  #                           harmonize_baseline=harmonize_baseline, ref_year=ref_year, yield_ratio=0.1, aggregate=FALSE)
  #
  #   # Adjust years
  #   years_wat <- getYears(avl_water_day)
  #   years_grper  <- getYears(grow_days)
  #   if(length(years_wat)>=length(years_grper)){
  #     years <- years_grper
  #   } else {
  #     years <- years_wat
  #   }
  #   rm(years_grper, years_wat)
  #
  #   # Available water in growing period per month
  #   avl_water_grper <- avl_water_day[,years,]*grow_days[,years,]
  #   # Available water in growing period per year
  #   avl_water_grper <- dimSums(avl_water_grper, dim=3)
  #
  #   # Check for NAs
  #   if(any(is.na(avl_water_grper))){
  #     stop("produced NA water availability")
  #   }
  #   out=avl_water_grper
  #   description="Available water in growing period per year"
  # }

  return(list(
    x=out,
    weight=NULL,
    unit="mio. m^3",
    description=description,
    isocountries=FALSE))
}
