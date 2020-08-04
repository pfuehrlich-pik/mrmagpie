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
  data <- toolGetMapping("River_structure_stn.rda",type="cell")
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
  yearly_runoff     <- yearly_runoff[,,1]
  years <- getYears(yearly_runoff)

  # Yearly lake evapotranspiration (in mio. m^3 per year) [place holder]
  lake_evap     <- input_lake <- calcOutput("LPJmL", version="LPJmL4", climatetype=climatetype, subtype="evap_lake_lpjcell", aggregate=FALSE,
                                            harmonize_baseline=FALSE, time="spline", dof=4, averaging_range=NULL)
  lake_evap     <- as.array(collapseNames(lake_evap))
  lake_evap     <- lake_evap[,,1]

  # Precipitation/Runoff on lakes and rivers from LPJmL (in mio. m^3 per year) [place holder]
  input_lake     <- input_lake <- calcOutput("LPJmL", version="LPJmL4", climatetype=climatetype, subtype="input_lake_lpjcell", aggregate=FALSE,
                                             harmonize_baseline=FALSE, time="spline", dof=4, averaging_range=NULL)
  input_lake     <- as.array(collapseNames(input_lake))
  input_lake     <- input_lake[,,1]

  # runoff (on land and water)
  yearly_runoff <- yearly_runoff + input_lake

  # Environmental Flow Requirements (in mio. m^3 / yr) [long-term average]
  if (EFR==TRUE){
    EFR_magpie <- calcOutput("EnvmtlFlow", version="LPJmL4", climatetype=climatetype, aggregate=FALSE, cells="lpjcell",
                      LFR_val=0.1, HFR_LFR_less10=0.2, HFR_LFR_10_20=0.15, HFR_LFR_20_30=0.07, HFR_LFR_more30=0.00,
                      EFRyears=c(1980:2010))
  } else if (EFR==FALSE){
    EFR_magpie <- new.magpie(1:NCELLS,fill=0)
    getCells(EFR_magpie) <- paste(lpj_cells_map$ISO,1:67420,sep=".")
  } else {
    stop("Specify whether environmental flows are activated or not via argument EFR")
  }
  EFR_magpie <- as.array(collapseNames(EFR_magpie))
  EFR_magpie <- EFR_magpie[,1,1]

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

  # Harmonize non-agricultural consumption and withdrawals
  NAg_ww <- pmax(NAg_ww, NAg_wc)
  NAg_wc <- pmax(NAg_wc, 0.01*NAg_ww)

  # Committed agricultural uses (in mio. m^3 / yr) [for initialization year]
  CAU_magpie <- calcOutput("CommittedAgWaterUse",iniyear=1995,irrigini="Jaegermeyr_lpjcell",time="raw",dof=NULL,aggregate=FALSE)
  CAW_magpie <- as.array(collapseNames(dimSums(CAU_magpie[,,"withdrawal"],dim=3)))
  CAC_magpie <- as.array(collapseNames(dimSums(CAU_magpie[,,"consumption"],dim=3)))
  rm(CAU_magpie)

  #############################
  ####### River routing #######
  #############################
  for (y in "y1995"){
    for (scen in "ssp2"){
      ## Global river routing variables
      # Naturalized discharge
      discharge_nat <- array(data=0,dim=NCELLS,dimnames=list(names(EFR_magpie)))
      inflow_nat    <- array(data=0,dim=NCELLS,dimnames=list(names(EFR_magpie)))
      lake_evap_new <- array(data=0,dim=NCELLS,dimnames=list(names(EFR_magpie)))
      # Discharge considering human uses
      discharge   <- array(data=0,dim=NCELLS,dimnames=list(names(EFR_magpie)))
      inflow      <- array(data=0,dim=NCELLS,dimnames=list(names(EFR_magpie)))
      avl_wat_act <- array(data=0,dim=NCELLS,dimnames=list(names(EFR_magpie)))
      # Water not available for consumption
      frac_NAg_fulfilled <- array(data=0,dim=NCELLS,dimnames=list(names(EFR_magpie)))
      frac_CAg_fulfilled <- array(data=0,dim=NCELLS,dimnames=list(names(EFR_magpie)))

      ### River Routing 1: Natural flows ###
      # Determine natural discharge
      for (o in 1:max(calcorder)){
        # Note: the calcorder ensures that the upstreamcells are calculated first
        cells <- which(calcorder==o)

        for (c in cells){
          ### Natural water balance
          # lake evap that can be fulfilled:
          # (if water available: lake evaporation considered; if not: lake evap is reduced respectively)
          lake_evap_new[c] <- min(lake_evap[c,y], inflow_nat[c]+yearly_runoff[c,y])
          # discharge
          discharge_nat[c] <- inflow_nat[c] + yearly_runoff[c,y] - lake_evap_new[c]
          # inflow into nextcell
          if (nextcell[c]>0){
            inflow_nat[nextcell[c]] <- inflow_nat[nextcell[c]] + discharge_nat[c]
          }
        }
      }

      # ### River Routing 1.2: Upstreamrouting - Environment: EFR & lake_evap [Basin Closure Check] ###
      # for (u in max(calcorder):1){
      #   # Note: the calcorder ensures that the upstreamcells are calculated first
      #   cells <- which(calcorder==u)
      #
      #   for (c in cells){
      #     # Tributary cells
      #     cells_trib <- which(nextcell==c)
      #     # EFRs that come from tributary cells
      #     inflow_guaranteed <- sum(discharge_min[cells_trib])
      #     # EFR and lake evap that cannot be fulfilled by tributary inflows (need to come from other sources, see below)
      #     EFR_add  <- discharge_min[c] - inflow_guaranteed + lake_evap_new[c]
      #
      #     # EFRs and lake evap needed fulfilled by other sources than tributary inflows:
      #     if (EFR_add>0){
      #       # "fair" contribution of runoff to EFRs (proportional to share of EFR in discharge)
      #       if (discharge_nat[c]>EFR_magpie[c]){
      #         EFR_runoff <- yearly_runoff[c,y] * EFR_magpie[c]/discharge_nat[c]
      #       } else {
      #         EFR_runoff <- yearly_runoff[c,y]
      #       }
      #       if (EFR_runoff>=EFR_add){
      #         # runoff that needs to stay untouched to fulfill EFRs
      #         runoff_allocated[c] <- EFR_add
      #       } else {
      #         # runoff that needs to stay untouched to fulfill EFRs
      #         runoff_allocated[c] <- EFR_runoff
      #         # local EFR_runoff not sufficient to fulfill additional EFR
      #         EFR_add <- EFR_add-EFR_runoff
      #         # additional water available that can be used to fulfill EFR_add
      #         water_avail <- sum(discharge_nat[cells_trib]-discharge_min[cells_trib]) + yearly_runoff[c,y] - EFR_runoff
      #
      #         # fraction of EFR_add that can be fulfilled by additionally available water
      #         if (water_avail>EFR_add){
      #           frac_add_EFR <- EFR_add/water_avail
      #         } else {
      #           frac_add_EFR <- 1
      #         }
      #
      #         # Runoff reserved for EFRs and lake evap
      #         runoff_allocated[c] <- runoff_allocated[c] + (yearly_runoff[c,y]-EFR_runoff)*frac_add_EFR
      #         # Discharge reserved for EFRs and lake evap
      #         if (length(cells_trib)>0){
      #           discharge_min[cells_trib] <- discharge_min[cells_trib] + (discharge_nat[cells_trib]-discharge_min[cells_trib])*frac_add_EFR
      #         }
      #       }
      #     }
      #     # Inflow reserved for EFR and lake evap
      #     inflow_allocated[c] <- EFR_magpie[c] + lake_evap_new[c] - runoff_allocated[c]
      #   }
      # }
      #
      # ## Outputs:
      # discharge_min_EFR <- discharge_min
      # #inflow_allocated_nat    <- inflow_allocated
      # #runoff_allocated_nat    <- runoff_allocated


      # Minimum availability of water in river to fulfill local EFRs
      discharge_min <- pmin(discharge_nat, EFR_magpie)
      # Discharge from cell to next considering human uses
      discharge <- discharge_nat
      # Environmental flows cannot exceed natural discharge
      EFR_magpie <- pmin(discharge_nat, EFR_magpie)

      ### River Routing 2: Non-agricultural uses considering local EFRs ###
      for (o in 1:max(calcorder)) {
        # Note: the calcorder ensures that the upstreamcells are calculated first
        cells <- which(calcorder==o)

        for (c in cells){

          if (discharge[c]<discharge_min[c]){
            # discharge in cell not sufficient to fulfill requirements: upstreamcells must release more water
            # -> necessary to restrict withdrawal upstream....
            if (upstreamcells[c]>0){
              for (i in upstreamcells[c]){
                discharge[i] <- discharge[i] + (discharge[c]-discharge_min[c])/length(which(nextcell==c))
              }
            } else {
              discharge[c] <- discharge_min[c]
              NAg_ww[c,y,scen] <- 0
            }
          } else {
            discharge[c] <- discharge_min[c]
          }

          # available water in cell
          avl_wat_act[c] <- max(inflow[c]+yearly_runoff[c,y]-lake_evap_new[c], 0)

          if (NAg_ww[c,y,scen]>0){
            ## Water withdrawals must not exceed local EFRs
            frac_NAg_fulfilled[c] <- min(max(avl_wat_act[c]-EFR_magpie[c], 0)/NAg_ww[c,y,scen], 1)
            if (frac_NAg_fulfilled[c]>0){
              ## Water consumption must not exceed availability
              frac_NAg_fulfilled[c] <- min(max(avl_wat_act[c]-discharge_min[c], 0)/(NAg_wc[c,y,scen]*frac_NAg_fulfilled[c]), frac_NAg_fulfilled[c])
            }
          }

          # Minimum river water requirement in cell: EFR + non-ag. withdrawal
          discharge_min[c] <- EFR_magpie[c] + NAg_ww[c,y,scen]*frac_NAg_fulfilled[c]

          ## Outflow from one cell to the next
          # (Subtract local water consumption in current cell (non-ag. consumption))
          discharge[c] <- avl_wat_act[c] - NAg_wc[c,y,scen]*frac_NAg_fulfilled[c]

          if (nextcell[c]>0){
            inflow[nextcell[c]] <- inflow[nextcell[c]] + discharge[c]
          }
        }
      }


#
#       ## Water withdrawal accounting
#       # (Water withdrawn downstream can be withdrawn upstream, but not consumed)
#       # NAC_c = ((withdrawal-runoff)/sum(Inflow))_(c+1) * outflow_c
#       # Note: which(nextcell==c): cells that go into current cell
#       # Water demand (withdrawal) in current cell that cannot be fulfilled by runoff on that cell, i.e. that needs to come from upstream cell:
#       wat_dem_exceeding_runoff[c] <- NAg_ww[c,y,scen] - yearly_runoff[c]
#
#       # Share of water that is needed for downstream withdrawal (cannot be consumed in current cell)
#       if (nextcell!=-1) {
#         for (i in 1:length(which(nextcell==c))) {
#           # water requirement from
#           shr_downstreamcell_requirement[which(nextcell==c)[i]] <- wat_dem_exceeding_runoff[c]/sum(discharge[which(nextcell==c)])
#         }
#       } else if (nextcell==-1) {
#         shr_downstreamcell_requirement[c] <- 0
#       }
#
#       # Water that is needed for downstream withdrawal (not available for consumption in current cell)
#       #shr_downstreamcell_requirement[c] <- max(shr_downstreamcell_requirement[c])
#       NAC_water[c] <- shr_downstreamcell_requirement[c==nextcell] * discharge[c]




      # ### River Routing 2.2: Upstreamrouting - Water withdrawal accounting ###
      # for (u in max(calcorder):1) {
      #   # (Water withdrawn downstream can be withdrawn upstream, but not consumed)
      #   # NAC_c = ((withdrawal-runoff)/sum(Inflow))_(c+1) * outflow_c
      #   # Note: which(nextcell==c): cells that go into current cell
      #   # Note: the calcorder ensures that the upstreamcells are calculated first
      #   cells <- which(calcorder==u)
      #
      #   for (c in cells){
      #     water_needed <- EFR_magpie[c] + NAg_ww_magpie[c,y,scen]*frac_NAg_fulfilled[c] + lake_evap_new[c]
      #     # Water demand (withdrawal) in current cell coming from runoff on that cell:
      #     withdrawal_from_runoff <- min(NAg_ww[c,y,scen]*frac_NAg_fulfilled[c], yearly_runoff[c,y]-runoff_allocated[c])
      #     #runoff_allocated[c]    <- runoff_allocated[c] + withdrawal_from_runoff
      #     # Water demand (withdrawal) in current cell that cannot be fulfilled by runoff on that cell, i.e. that needs to come from upstream cell:
      #     withdrawal_from_inflow <- NAg_ww[c,y,scen]*frac_NAg_fulfilled[c] - withdrawal_from_runoff
      #     #inflow_allocated[c]    <- inflow_allocated[c] + withdrawal_from_inflow
      #
      #     # Tributary cells
      #     cells_trib <- which(nextcell==c)
      #     # Water reserved in tributary cells
      #     inflow_guaranteed <- sum(discharge_min[cells_trib])
      #
      #     if (inflow_allocated[c]>inflow_guaranteed){
      #       # Inflow to cell that can be allocated
      #       inflow_avl        <- sum(discharge[cells_trib]-discharge_min[cells_trib])
      #       if (inflow_avl>0){
      #         frac_add_NAgww <- (inflow_allocated[c]-inflow_guaranteed)/inflow_avl
      #       } else {
      #         frac_add_NAgww <- 1
      #       }
      #       if (length(cells_trib)>0){
      #         discharge_min[cells_trib] <- discharge_min[cells_trib] + (discharge[cells_trib]-discharge_min[cells_trib])*frac_add_NAgww
      #       }
      #     }
      #   }
      # }
      #
      # ## Outputs:
      # discharge_min_NAg <- discharge_min
      # inflow_allocated_NAg    <- inflow_allocated
      # runoff_allocated_NAg    <- runoff_allocated


      ### River Routing 3: Routing of "committed agricultural uses" taking previous withdrawals and EFRs into account ###
      for (o in 1:max(calcorder)) {
        # Note: the calcorder ensures that the upstreamcells are calculated first
        cells <- which(calcorder==o)

          for (c in cells){
            # available water in cell (already reduced by non-ag withdrawals)
            avl_wat_act[c] <- max(inflow[c]+yearly_runoff[c,y]-lake_evap_new[c], 0)

            ## Water withdrawals must not exceed availability (considering EFRs)
            if (CAW_magpie[c]>0){
              frac_CAg_fulfilled[c] <- min(max(avl_wat_act[c]-discharge_min[c], 0)/CAW_magpie[c], 1)
            }

            ## Outflow from one cell to the next
            # (Subtract local water consumption in current cell (non-ag. consumption))
            discharge[c]        <- avl_wat_act[c] - CAC_magpie[c]*frac_CAg_fulfilled[c]
            if (nextcell[c]>0){
              inflow[nextcell[c]] <- inflow[nextcell[c]] + discharge[c]
            }
            # Discharge resulting from return flows of committed ag. water consumption
            discharge_min[c] <- discharge_min[c] + (CAW_magpie[c]-CAC_magpie[c])*frac_CAg_fulfilled[c]
          }
      }


      # #### Maybe;
      # ## Water withdrawal accounting
      # # (Water withdrawn downstream can be withdrawn upstream, but not consumed)
      # # NAC_c = ((withdrawal-runoff)/sum(Inflow))_(c+1) * outflow_c
      # # Note: which(nextcell==c): cells that go into current cell
      # # Water demand (withdrawal) in current cell that cannot be fulfilled by runoff on that cell, i.e. that needs to come from upstream cell:
      # wat_dem_exceeding_runoff[c] <- (NAg_ww[c,y,scen] + CAW_magpie[c]) - yearly_runoff[c]
      #
      # # Share of water that is needed for downstream withdrawal (cannot be consumed in current cell)
      # if (nextcell!=-1) {
      #   for (i in 1:length(which(nextcell==c))) {
      #     # water requirement from
      #     shr_downstreamcell_requirement[which(nextcell==c)[i]] <- wat_dem_exceeding_runoff[c]/sum(discharge[which(nextcell==c)])
      #   }
      # } else if (nextcell==-1) {
      #   shr_downstreamcell_requirement[c] <- 0
      # }
      #
      # # Water that is needed for downstream withdrawal (not available for consumption in current cell)
      # shr_downstreamcell_requirement[c] <- max(shr_downstreamcell_requirement[c])
      # NAC_water[c] <- shr_downstreamcell_requirement[c==nextcell] * discharge[c]

      # ### River Routing 2.2: Upstreamrouting - Water withdrawal accounting ###
      # for (u in max(calcorder):1) {
      #   # (Water withdrawn downstream can be withdrawn upstream, but not consumed)
      #
      #   # Water demand (withdrawal) in current cell coming from runoff on that cell:
      #   withdrawal_from_runoff <- min(CAW_magpie[c,y,scen]*frac_CAg_fulfilled[c], yearly_runoff[c,y]-runoff_allocated[c])
      #   runoff_allocated[c]    <- runoff_allocated[c] + withdrawal_from_runoff
      #   # Water demand (withdrawal) in current cell that cannot be fulfilled by runoff on that cell, i.e. that needs to come from upstream cell:
      #   withdrawal_from_inflow <- CAW_magpie[c,y,scen]*frac_CAg_fulfilled[c] - withdrawal_from_runoff
      #   inflow_allocated[c]    <- inflow_allocated[c] + withdrawal_from_inflow
      #
      #   # Tributary cells
      #   cells_trib <- which(nextcell==c)
      #   # Water reserved in tributary cells
      #   inflow_guaranteed <- sum(discharge_min[cells_trib])
      #
      #   if (inflow_allocated[c]>inflow_guaranteed){
      #     # Inflow to cell that can be allocated
      #     inflow_avl        <- sum(discharge[cells_trib]-discharge_min[cells_trib])
      #     if (inflow_avl>0){
      #       frac_add_CAgww <- (inflow_allocated[c]-inflow_guaranteed)/inflow_avl
      #     } else {
      #       frac_add_CAgww <- 1
      #     }
      #     if (length(cell_trib)>0){
      #       discharge_min[cells_trib] <- discharge_min[cells_trib] + (discharge[cells_trib]-discharge_min[cells_trib])*frac_add_CAgww
      #     }
      #   }
      # }
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
