#' @title calcAvlWater
#' @description This function calculates water availability for MAgPIE retrieved from LPJmL using a river routing algorithm for distribution of discharge across the river basin
#'
#' @param version     Switch between LPJmL4 and LPJmL5
#' @param climatetype Switch between different climate scenarios (default: "CRU_4")
#' @param time            Time smoothing: average, spline or raw (default)
#' @param averaging_range only specify if time=="average": number of time steps to average
#' @param dof             only specify if time=="spline": degrees of freedom needed for spline
#' @param harmonize_baseline FALSE (default): no harmonization, TRUE: if a baseline is specified here data is harmonized to that baseline (from ref_year on)
#' @param ref_year           Reference year for harmonization baseline (just specify when harmonize_baseline=TRUE)
#' @param selectyears Years to be returned
#' @param EFR Switch for activation of environmental flow requirements (TRUE) or not (FALSE)
#' @param allocationrule rule to be applied for river basin discharge allocation across cells of river basin ("optimization" (default), "upstreamfirst", "equality")
#' @param allocationshare share of water to be allocated to cell (only needs to be selected in case of allocationrule=="equality")
#' @param gainthreshold threshold of yield improvement potential required for water allocation in upstreamfirst algorithm (in tons per ha)
#' @param irrigationsystem irrigation system to be used for river basin discharge allocation algorithm ("surface", "sprinkler", "drip", "initialization")
#' @param irrigini when "initialization" selected for irrigation system: choose initialization data set for irrigation system initialization ("Jaegermeyr_lpjcell", "LPJmL_lpjcell")
#'
#' @import magclass
#' @import madrat
#' @importFrom mrcommons toolHarmonize2Baseline
#'
#' @return magpie object in cellular resolution
#' @author Felicitas Beier, Jens Heinke
#'
#' @examples
#' \dontrun{ calcOutput("AvlWater", aggregate = FALSE) }
#'

calcAvlWater <- function(selectyears="all",
                         version="LPJmL4", climatetype="HadGEM2_ES:rcp2p6:co2", time="raw", averaging_range=NULL, dof=NULL,
                         harmonize_baseline=FALSE, ref_year="y2015", EFR=TRUE,
                         allocationrule="optimization", allocationshare=NULL, gainthreshold=1,
                         irrigationsystem="initialization", irrigini="Jaegermeyr_lpjcell"){

  #############################
  ####### Read in Data ########
  #############################

  ### Read in river structure
  # Note: river structure derived from LPJmL input (drainage) [maybe later: implement readDrainage function]
  data <- toolGetMapping("River_structure_stn.rda",type="cell")
  for (i in 1:length(data)){
    assign(paste(names(data[[i]])), data[[i]][[1]])
  }
  rm(data,i)

  # Number of cells to be used for calculation
  NCELLS <- length(calcorder)

  ### LPJ-MAgPIE cell mapping
  #load("C:/Users/beier/Documents/doktorarbeit/MAgPIE_Water/River_Routing_Postprocessing/cells_magpie2lpj.Rda")
  magpie2lpj    <- magclassdata$cellbelongings$LPJ_input.Index
  lpj_cells_map <- toolGetMapping("LPJ_CellBelongingsToCountries.csv", type="cell")

  ### Required inputs:
  # Yearly runoff (mio. m^3 / yr) [smoothed]
  yearly_runoff <- calcOutput("LPJmL", version="LPJmL4", climatetype=climatetype, subtype="runoff_lpjcell", aggregate=FALSE,
                                        harmonize_baseline=FALSE, time="spline", dof=4, averaging_range=NULL)
  yearly_runoff <- as.array(collapseNames(yearly_runoff))
  yearly_runoff <- yearly_runoff[,,1]
  years <- getYears(yearly_runoff)

  # Yearly lake evapotranspiration (in mio. m^3 per year)
  lake_evap     <- calcOutput("LPJmL", version="LPJmL4", climatetype=climatetype, subtype="evap_lake_lpjcell", aggregate=FALSE,
                                            harmonize_baseline=FALSE, time="spline", dof=4, averaging_range=NULL)
  lake_evap     <- as.array(collapseNames(lake_evap))
  lake_evap     <- lake_evap[,,1]

  # Precipitation/Runoff on lakes and rivers from LPJmL (in mio. m^3 per year)
  input_lake     <- calcOutput("LPJmL", version="LPJmL4", climatetype=climatetype, subtype="input_lake_lpjcell", aggregate=FALSE,
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

  # Harmonize non-agricultural consumption and withdrawals (withdrawals > consumption)
  NAg_ww <- pmax(NAg_ww, NAg_wc)
  NAg_wc <- pmax(NAg_wc, 0.01*NAg_ww)

  # Committed agricultural uses (in mio. m^3 / yr) [for initialization year]
  CAU_magpie <- calcOutput("CommittedAgWaterUse",iniyear=1995,irrigini="Jaegermeyr_lpjcell",time="raw",dof=NULL,aggregate=FALSE)
  CAW_magpie <- as.array(collapseNames(dimSums(CAU_magpie[,,"withdrawal"],dim=3)))
  CAC_magpie <- as.array(collapseNames(dimSums(CAU_magpie[,,"consumption"],dim=3)))
  rm(CAU_magpie)
  CAW_magpie <- as.array(collapseNames(CAW_magpie))
  CAW_magpie <- CAW_magpie[,1,1]
  CAC_magpie <- as.array(collapseNames(CAC_magpie))
  CAC_magpie <- CAC_magpie[,1,1]

  #############################
  ####### River routing #######
  #############################
  for (y in "y1995"){
    for (scen in "ssp2"){
     # for (EFR_scen in c("EFR_on", "EFR_off"))

      ## Global river routing variables
      # Naturalized discharge
      discharge_nat <- array(data=0,dim=NCELLS,dimnames=list(names(EFR_magpie)))
      inflow_nat    <- array(data=0,dim=NCELLS,dimnames=list(names(EFR_magpie)))
      lake_evap_new <- array(data=0,dim=NCELLS,dimnames=list(names(EFR_magpie)))
      # Discharge considering human uses
      discharge   <- array(data=0,dim=NCELLS,dimnames=list(names(EFR_magpie)))
      inflow      <- array(data=0,dim=NCELLS,dimnames=list(names(EFR_magpie)))
      avl_wat_act <- array(data=0,dim=NCELLS,dimnames=list(names(EFR_magpie)))
      # Water fractions reserved for certain uses
      frac_NAg_fulfilled <- array(data=0,dim=NCELLS,dimnames=list(names(EFR_magpie)))
      frac_CAg_fulfilled <- array(data=0,dim=NCELLS,dimnames=list(names(EFR_magpie)))
      frac_fullirrig     <- array(data=0,dim=NCELLS,dimnames=list(names(EFR_magpie)))
      required_wat_min   <- array(data=0,dim=NCELLS,dimnames=list(names(EFR_magpie)))

      ### River Routing 1.1: Natural flows ###
      # Determine natural discharge
      for (o in 1:max(calcorder)){
        # Note: the calcorder ensures that upstreamcells are calculated first
        cells <- which(calcorder==o)

        for (c in cells){
          ### Natural water balance
          # lake evap that can be fulfilled (if water available: lake evaporation considered; if not: lake evap is reduced respectively):
          lake_evap_new[c] <- min(lake_evap[c,y], inflow_nat[c]+yearly_runoff[c,y])
          # natural discharge
          discharge_nat[c] <- inflow_nat[c] + yearly_runoff[c,y] - lake_evap_new[c]
          # inflow into nextcell
          if (nextcell[c]>0){
            inflow_nat[nextcell[c]] <- inflow_nat[nextcell[c]] + discharge_nat[c]
          }
        }
      }

      # Minimum availability of water in river to fulfill local EFRs
      required_wat_min <- EFR_magpie

      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      #-#-#-#-# Output reporting #-#-#-#-#
      # Check where water requirements exceed natural discharge
      summary(required_wat_min-discharge_nat)        #### Question: how do we explain this?
      check1 <- sum(required_wat_min-discharge_nat>1e-3)

      # Output reporting
      ratio_routing1 <- required_wat_min/discharge_nat
      ratio_routing1[which(required_wat_min==0 & discharge_nat==0)] <- NA
      ratio_routing1[which(ratio_routing1>1)] <- 1

      plotmap2(mrmagpie:::toolLPJarrayToMAgPIEmap(ratio_routing1))

      # Downstream constraint
      tmp                  <- pmax(discharge_nat - required_wat_min,0)
      wat_avl_consumption1 <- numeric(NCELLS)

      for (c in 1:NCELLS){
        # available for consumption in current cell considering downstream cells
        wat_avl_consumption1[c] <- min(tmp[c(downstreamcells[[c]],c)])
      }

      plotmap2(mrmagpie:::toolLPJarrayToMAgPIEmap(1-(wat_avl_consumption1/discharge_nat)))
      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

      ### River Routing 2: Non-agricultural uses considering local EFRs ###
      for (o in 1:max(calcorder)) {
        # Note: the calcorder ensures that the upstreamcells are calculated first
        cells <- which(calcorder==o)

        for (c in cells){
          # available water in cell
          avl_wat_act[c]  <- inflow[c]+yearly_runoff[c,y]-lake_evap_new[c]

          # available water in cell not sufficient to fulfill requirements
          # -> no more water can be withdrawn
          if (avl_wat_act[c]<required_wat_min[c]){
            # if cell has upstreamcells: upstreamcells must release missing water (cannot be consumed upstream)
            # -> reduce non-agricultural water consumption in upstream cells
            # -> locally: cannot withdraw
            if (length(upstreamcells[c])>0){
              # upstream non-agricultural water consumption
              upstream_cons <- sum(NAg_wc[upstreamcells[[c]],y,scen]*frac_NAg_fulfilled[upstreamcells[[c]]])
              if (upstream_cons>required_wat_min[c]-avl_wat_act[c]){
                # if missing water (difference) can be fulfilled by upstream consumption: reduce upstream consumption
                frac_NAg_fulfilled[upstreamcells[[c]]] <- (1-(required_wat_min[c]-avl_wat_act[c])/upstream_cons)*frac_NAg_fulfilled[upstreamcells[[c]]]
                discharge[c] <- required_wat_min[c]
              } else {
                # if missing water (difference) cannot be fulfilled by upstream consumption: no upstream consumption
                frac_NAg_fulfilled[upstreamcells[[c]]] <- 0
                discharge[c] <- avl_wat_act[c]+upstream_cons
              }
            }

          # available water in cell is sufficient to fulfill requirements
          # -> further withdrawals are possible
          } else {
            # Non-agricultural withdrawals
            if (NAg_ww[c,y,scen]>0){
              ## Water withdrawal constraint:
              frac_NAg_fulfilled[c] <- min((avl_wat_act[c]-required_wat_min[c])/NAg_ww[c,y,scen], 1)
            }

            ## Outflow from one cell to the next
            # (Subtract local water consumption in current cell (non-ag. consumption))
            discharge[c] <- avl_wat_act[c] - NAg_wc[c,y,scen]*frac_NAg_fulfilled[c]
          }

          if (nextcell[c]>0){
            inflow[nextcell[c]] <- inflow[nextcell[c]] + discharge[c]
          }
        }
      }

      # Update minimum water required in cell:
      required_wat_min <- required_wat_min + NAg_ww[,y,scen]*frac_NAg_fulfilled

      ### Interim routing: Update discharge and inflow considering known non-agricultural uses of river routing 2 ###
      inflow[] <- 0

      for (o in 1:max(calcorder)){
        # Note: the calcorder ensures that the upstreamcells are calculated first
        cells <- which(calcorder==o)

        for (c in cells){
          # available water
          avl_wat_act[c] <- inflow[c] + yearly_runoff[c,y] - lake_evap_new[c]

          # discharge
          discharge[c]   <- avl_wat_act[c] - NAg_wc[c,y,scen]*frac_NAg_fulfilled[c]

          # inflow into nextcell
          if (nextcell[c]>0){
            inflow[nextcell[c]] <- inflow[nextcell[c]] + discharge[c]
          }
        }
      }

      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      #-#-#-#-# Output reporting #-#-#-#-#
      # Check whether number of cells where water requirements exceed availability has increased
      if (sum(required_wat_min-avl_wat_act>1e-3)>check1) warning("River routing violation of water availability.")

      # Output reporting
      ratio_routing2 <- required_wat_min/avl_wat_act
      ratio_routing2[which(required_wat_min==0 & avl_wat_act==0)] <- NA
      ratio_routing2[which(ratio_routing2>1)] <- 1

      plotmap2(mrmagpie:::toolLPJarrayToMAgPIEmap(ratio_routing2))

      # Downstream consideration
      tmp                  <-  pmax(avl_wat_act - required_wat_min,0)
      wat_avl_consumption2 <- numeric(NCELLS)

      for (c in 1:NCELLS){
        # available for consumption in current cell considering downstream cells
        wat_avl_consumption2[c] <- min(tmp[c(downstreamcells[[c]],c)])
      }

      plotmap2(mrmagpie:::toolLPJarrayToMAgPIEmap(1-(wat_avl_consumption2/discharge_nat)))
      plotmap2(mrmagpie:::toolLPJarrayToMAgPIEmap(1-(wat_avl_consumption2/avl_wat_act)))
      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

      # inflow needs to be set to 0 prior to every river routing (is recalculated by the routing)
      inflow[] <- 0

      ### River Routing 3: Committed agricultural uses considering local EFRs and non-agricultural uses ###
      for (o in 1:max(calcorder)) {
        # Note: the calcorder ensures that the upstreamcells are calculated first
        cells <- which(calcorder==o)

        for (c in cells){
          # available water in cell
          avl_wat_act[c]  <- inflow[c]+yearly_runoff[c,y]-lake_evap_new[c]

          # available water in cell not sufficient to fulfill requirements
          # -> no more water can be withdrawn
          if (avl_wat_act[c]<required_wat_min[c]){
            # if cell has upstreamcells: upstreamcells must release missing water (cannot be consumed upstream)
            # -> reduce committed agricultural water consumption in upstream cells
            # -> locally: cannot withdraw
            if (length(upstreamcells[c])>0){
              # upstream committed agricultural water consumption:
              upstream_cons <- sum(CAC_magpie[upstreamcells[[c]]]*frac_CAg_fulfilled[upstreamcells[[c]]])
              if (upstream_cons>required_wat_min[c]-avl_wat_act[c]){
                # if upstream_cons high enough to account for difference: reduce upstream consumption respectively
                frac_CAg_fulfilled[upstreamcells[[c]]] <- (1-(required_wat_min[c]-avl_wat_act[c])/upstream_cons)*frac_CAg_fulfilled[upstreamcells[[c]]]
                discharge[c] <- required_wat_min[c]
              } else {
                # if upstream_cons not sufficient to account for difference: no water can be used upstream
                frac_CAg_fulfilled[upstreamcells[[c]]] <- 0
                discharge[c] <- avl_wat_act[c]+upstream_cons
              }
            }

          # available water in cell is sufficient to fulfill requirements
          # -> further withdrawal possible
          } else {
            # Committed agricultural withdrawals
            if (CAW_magpie[c]>0){
              ## Water withdrawal constraint:
              frac_CAg_fulfilled[c] <- min((avl_wat_act[c]-required_wat_min[c])/CAW_magpie[c], 1)
            }

            ## Outflow from one cell to the next
            # (Subtract local water consumption in current cell (committed ag. & non-agricultural consumption))
            discharge[c] <- avl_wat_act[c] - CAC_magpie[c]*frac_CAg_fulfilled[c] - NAg_wc[c,y,scen]*frac_NAg_fulfilled[c]
          }

          if (nextcell[c]>0){
            inflow[nextcell[c]] <- inflow[nextcell[c]] + discharge[c]
          }
        }
      }

      # Update minimum water required in cell:
      required_wat_min <- required_wat_min + CAW_magpie*frac_CAg_fulfilled

      ### Interim routing: Update discharge and inflow considering known non-agricultural and committed agricultural uses of river routing 3 ###
      inflow[] <- 0

      for (o in 1:max(calcorder)){
        # Note: the calcorder ensures that the upstreamcells are calculated first
        cells <- which(calcorder==o)

        for (c in cells){
          # available water
          avl_wat_act[c] <- inflow[c] + yearly_runoff[c,y] - lake_evap_new[c]
          # discharge
          discharge[c]   <- avl_wat_act[c] - NAg_wc[c,y,scen]*frac_NAg_fulfilled[c] - CAC_magpie[c]*frac_CAg_fulfilled[c]
          # inflow into nextcell
          if (nextcell[c]>0){
            inflow[nextcell[c]] <- inflow[nextcell[c]] + discharge[c]
          }
        }
      }

      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      #-#-#-#-# Output reporting #-#-#-#-#
      # Check whether number of cells where water requirements exceed availability has increased
      if (sum(required_wat_min-avl_wat_act>1e-3)>check1) warning("River routing violation of water availability")

      # Output reporting
      ratio_routing3 <- required_wat_min/avl_wat_act
      ratio_routing3[which(required_wat_min==0 & avl_wat_act==0)]<-NA
      ratio_routing3[which(ratio_routing3>1)]<-1

      plotmap2(mrmagpie:::toolLPJarrayToMAgPIEmap(ratio_routing3))

      # Downstream consideration
      tmp                  <- pmax(avl_wat_act - required_wat_min,0)
      wat_avl_consumption3 <- numeric(NCELLS)

      for (c in 1:NCELLS){
        # available for consumption in current cell considering downstream cells
        wat_avl_consumption3[c] <- min(tmp[c(downstreamcells[[c]],c)])
      }

      plotmap2(mrmagpie:::toolLPJarrayToMAgPIEmap(1-(wat_avl_consumption3/discharge_nat)))
      plotmap2(mrmagpie:::toolLPJarrayToMAgPIEmap(1-(wat_avl_consumption3/discharge)))
      plotmap2(mrmagpie:::toolLPJarrayToMAgPIEmap(1-(wat_avl_consumption3/avl_wat_act)))
      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#


      ##### MAKE SURE THAT DISCHARGE IS >0 (note: or make sure in loop...)
      #discharge <- pmax(discharge, 0)

      ################################################
      ####### River basin discharge allocation #######
      ################################################
      ### River basin discharge (that can be allocated to cells within the basin)
      # cell_basin_mapping <- array(data=0,dim=NCELLS)
      # #basin_discharge    <- array(data=0,dim=length(unique(endcell)))
      #
      # basin_code <- 1
      # for (b in unique(endcell)){
      #   # mapping of cells to basins
      #   cell_basin_mapping[which(endcell==b)] <- basin_code
      #   # river basin discharge: discharge of endcell
      #  # basin_discharge[basin_code]           <- discharge[b]
      #
      #   basin_code <- basin_code+1
      # }
      # rm(basin_code)

      ### Required water for full irrigation per cell (in m^3)
      #?? with taking into account already irrigated area in intialization PROBLEM: negative values!!!!
      required_wat_fullirrig_ww <- calcOutput("FullIrrigationRequirement", version="LPJmL5", selectyears="y1995", climatetype="HadGEM2_ES:rcp2p6:co2", harmonize_baseline=FALSE, time="spline", dof=4, iniyear=1995, iniarea=TRUE, irrig_requirement="withdrawal", cells="lpjcell", aggregate=FALSE)[,,c("maiz","rapeseed","puls_pro")]
      required_wat_fullirrig_ww <- pmax(required_wat_fullirrig_ww,0)
      required_wat_fullirrig_wc <- calcOutput("FullIrrigationRequirement", version="LPJmL5", selectyears="y1995", climatetype="HadGEM2_ES:rcp2p6:co2", harmonize_baseline=FALSE, time="spline", dof=4, iniyear=1995, iniarea=TRUE, irrig_requirement="consumption", cells="lpjcell", aggregate=FALSE)[,,c("maiz","rapeseed","puls_pro")]
      required_wat_fullirrig_wc <- pmax(required_wat_fullirrig_wc,0)
      #?? without already irrigated area in initialization
      #... note: CAW is anyway subtracted later-on. or is it? which one makes more sense?
      #required_wat_fullirrig_ww <- calcOutput("FullIrrigationRequirement", version="LPJmL5", selectyears="y1995", climatetype="HadGEM2_ES:rcp2p6:co2", harmonize_baseline=FALSE, time="spline", dof=4, iniyear=1995, iniarea=FALSE, irrig_requirement="withdrawal", cells="lpjcell", aggregate=FALSE)[,,c("maiz","rapeseed","puls_pro")]
      #required_wat_fullirrig_wc <- calcOutput("FullIrrigationRequirement", version="LPJmL5", selectyears="y1995", climatetype="HadGEM2_ES:rcp2p6:co2", harmonize_baseline=FALSE, time="spline", dof=4, iniyear=1995, iniarea=FALSE, irrig_requirement="consumption", cells="lpjcell", aggregate=FALSE)[,,c("maiz","rapeseed","puls_pro")]

      # full irrigation water requirement depending on irrigation system in use
      if (irrigationsystem=="initialization") {
        # read in irrigation system area initialization [share of AEI by system]
        irrigation_system           <- calcOutput("IrrigationSystem", source=irrigini, aggregate=FALSE)
        getYears(irrigation_system) <- getYears(required_wat_fullirrig_ww)
        required_wat_fullirrig_ww   <- dimSums(irrigation_system*required_wat_fullirrig_ww,dim=3.1)
        required_wat_fullirrig_wc   <- dimSums(irrigation_system*required_wat_fullirrig_wc,dim=3.1)
      } else {
        # whole area irrigated by one system as selected in argument "irrigationsystem"
        required_wat_fullirrig_ww <- collapseNames(required_wat_fullirrig_ww[,,irrigationsystem])
        required_wat_fullirrig_wc <- collapseNames(required_wat_fullirrig_wc[,,irrigationsystem])
      }
      # Average required water for full irrigation across selected proxy crops
      required_wat_fullirrig_ww <- dimSums(required_wat_fullirrig_ww,dim=3)/length(getNames(required_wat_fullirrig_ww))
      required_wat_fullirrig_wc <- dimSums(required_wat_fullirrig_wc,dim=3)/length(getNames(required_wat_fullirrig_wc))

      # transform to array for further calculations
      required_wat_fullirrig_ww <- as.array(collapseNames(required_wat_fullirrig_ww))[,1,1]
      required_wat_fullirrig_wc <- as.array(collapseNames(required_wat_fullirrig_wc))[,1,1]

      # full irrigation requirement that is not already fulfilled by committed agricultural use
      #required_wat_fullirrig_ww <- required_wat_fullirrig_ww - CAW_magpie*frac_CAg_fulfilled
      #required_wat_fullirrig_wc <- required_wat_fullirrig_wc - CAC_magpie*frac_CAg_fulfilled
      ####???? necessary again? or when subtract area as above not necessary? -> which approach makes more sense?


      # Global cell rank based on yield gain potential by irrigation of proxy crops: maize, rapeseed, pulses
      meancellrank <- calcOutput("IrrigCellranking", version="LPJmL5", climatetype="HadGEM2_ES:rcp2p6:co2", time="spline", averaging_range=NULL, dof=4, harmonize_baseline=FALSE, ref_year="y2015",
                                 cellrankyear="y1995", cells="lpjcell", crops="magpie", method="meancroprank", proxycrop=c("maiz", "rapeseed", "puls_pro"), aggregate=FALSE)
      meancellrank <- as.array(meancellrank)[,1,1]

        ### !!!! make rank algorithm work for meancellrank again!

        # # Solve ties in cell ranking
         # proxycrop=c("maiz", "rapeseed", "puls_pro")
         # for (crop in proxycrop) {
         #
         #   if (length(meancellrank)!=length(unique(meancellrank))){
         #     cropcellrank <- calcOutput("IrrigCellranking", version="LPJmL5", climatetype="HadGEM2_ES:rcp2p6:co2", time="spline", averaging_range=NULL, dof=4, harmonize_baseline=FALSE, ref_year="y2015",
         #                                selectyears="y1995", cells="lpjcell", crops="magpie", method="cropcellrank", proxycrop=crop, aggregate=FALSE)
         #     cropcellrank <- as.array(cropcellrank)[,1,1]
         #
         #     for (i in unique(meancellrank[duplicated(meancellrank)])) {
         #
         #       cells <- which(meancellrank==i)
         #       l     <- length(meancellrank[meancellrank==i])
         #       n     <- i - floor(l/2)
         #       meancellrank[cells] <- rank(cropcellrank[cells])+n-1
         #     }
         #   }
         # }
#
#         # Solve ties in cell ranking
#         # 1.1) in case of tie use maize cellrank
#         tmp <- meancellrank
#         meancellrank[] <- NA
#         for (i in sort(tmp)) {
#           n <- ifelse(length(meancellrank[tmp==i])>2, 2, 1)
#           meancellrank[tmp==i] <- tmp[tmp==i] + (rank(-cellrank$maiz[tmp==i])-n)
#         }
#         # 1.2) in case of still tie: use rapeseed cellrank
#         tmp <- meancellrank
#         for (i in sort(tmp)) {
#           n <- ifelse(length(meancellrank[tmp==i])>2, 2, 1)
#           meancellrank[tmp==i] <- tmp[tmp==i] + (rank(-cellrank$rapeseed[tmp==i])-n)
#         }
#         # 1.3) in case of still tie: use puls_pro cellrank
#         tmp <- meancellrank
#         for (i in sort(tmp)) {
#           n <- ifelse(length(meancellrank[tmp==i])>2, 2, 1)
#           meancellrank[tmp==i] <- tmp[tmp==i] + (rank(-cellrank$puls_pro[tmp==i])-n)
#         }
#         rm(tmp)

        # rank ties by cell-id
        if (unique(meancellrank[duplicated(meancellrank)])>0) {
          for (i in unique(meancellrank[duplicated(meancellrank)])) {
            l <- length(meancellrank[meancellrank==i])
            n <- i - floor(l/2)
            for (k in (1:l)) {
              meancellrank[meancellrank==i][[1]] <- n-1+k
            }
          }
        }

        # rank left-over ties by cell-id
        # tmp <- meancellrank
        # for (i in unique(sort(tmp))) {
        #   c <- which(meancellrank==i)
        #   if (length(meancellrank[tmp==i])>1) {
        #     if (length(meancellrank[tmp==i])>2) {
        #       n <- -1
        #     } else {
        #       n <- 0
        #     }
        #     for (k in (1:length(meancellrank[which(meancellrank==i)]))){
        #       meancellrank[c] <- meancellrank[c]+n
        #       n <- n+1
        #     }
        #   }
        # }


        ############################
        ### Allocation Algorithm ###
        ############################
        # Minimum water required
        required_wat_min_allocation <- required_wat_min

        # Allocate water for full irrigation to cell with highest yield improvement through irrigation
        if (allocationrule=="optimization") {

          for (c in (1:max(meancellrank,na.rm=T))){
            # available water for additional irrigation withdrawals
            avl_wat_ww <- max(discharge[c]-required_wat_min_allocation[c],0)

            # withdrawal constraint
            if (required_wat_fullirrig_ww[c]>0) {
              # how much withdrawals can be fulfilled by available water
              frac_fullirrig[c] <- min(avl_wat_ww/required_wat_fullirrig_ww[c],1)

              # consumption constraint
              if (required_wat_fullirrig_wc[c]>0 & length(downstreamcells[[c]])>0) {
                # available water for additional irrigation consumption (considering downstream availability)
                avl_wat_wc          <- max(min(discharge[downstreamcells[[c]]] - required_wat_min_allocation[downstreamcells[[c]]]),0)
                # how much consumption can be fulfilled by available water
                frac_fullirrig[c]   <- min(avl_wat_wc/required_wat_fullirrig_wc[c],frac_fullirrig[c])
                # adjust discharge in current cell and downstream cells (subtract irrigation water consumption)
                discharge[c(downstreamcells[[c]],c)] <- discharge[c(downstreamcells[[c]],c)] - required_wat_fullirrig_wc[c]*frac_fullirrig[c]
              }
              # update minimum water required in cell:
              required_wat_min_allocation[c] <- required_wat_min_allocation[c] + frac_fullirrig[c]*required_wat_fullirrig_ww[c]
            }
          }
        } else if (allocationrule=="upstreamfirst") {
          # Allocate full irrigation requirements to most upstream cell first (calcorder)

          # Only consider cells where irrigation potential > 0
          # (or even: above certain threshold (e.g. 1 t/ha), maybe flexible (set in argument))
          # STILL MISSING

          # loop over basin
          for (b in unique(endcell) ) {

            # allocation to upstream first (calcorder=1)
            for (o in (1:max(calcorder[which(endcell==b)],na.rm=T))){
              c <- which(endcell==b & calcorder==o)

              # several cells with same calcorder in one basin
              for (k in (1:length(which(endcell==b & calcorder==o)))){
                # available water for additional irrigation withdrawals
                avl_wat_ww <- max(discharge[c[k]]-required_wat_min_allocation[c[k]],0)

                # withdrawal constraint
                if (required_wat_fullirrig_ww[c[k]]>0) {
                  # how much withdrawals can be fulfilled by available water
                  frac_fullirrig[c[k]] <- min(avl_wat_ww/required_wat_fullirrig_ww[c[k]],1)

                  # consumption constraint
                  if (required_wat_fullirrig_wc[c[k]]>0 & length(downstreamcells[[c[k]]])>0) {
                    # available water for additional irrigation consumption (considering downstream availability)
                    avl_wat_wc          <- max(min(discharge[downstreamcells[[c[k]]]] - required_wat_min_allocation[downstreamcells[[c[k]]]]),0)
                    # how much consumption can be fulfilled by available water
                    frac_fullirrig[c[k]]   <- min(avl_wat_wc/required_wat_fullirrig_wc[c[k]],frac_fullirrig[c[k]])
                    # adjust discharge in current cell and downstream cells (subtract irrigation water consumption)
                    discharge[c(downstreamcells[[c[k]]],c[k])] <- discharge[c(downstreamcells[[c[k]]],c[k])] - required_wat_fullirrig_wc[c[k]]*frac_fullirrig[c[k]]
                  }
                  # update minimum water required in cell:
                  required_wat_min_allocation[c[k]] <- required_wat_min_allocation[c[k]] + frac_fullirrig[c[k]]*required_wat_fullirrig_ww[c[k]]
                }
              }
            }
          }
        } else if (allocationrule=="equality") {

          # Repeat optimization algorithm several times
          # Instead of full irrigation, only up to x% (e.g.20%) are allocated to most efficient cell
          # Repeat until all river basin discharge is allocated ---> STILL MISSING

          for (c in (1:max(meancellrank,na.rm=T))){
            # available water for additional irrigation withdrawals
            avl_wat_ww <- max(discharge[c]-required_wat_min_allocation[c],0)

            # withdrawal constraint
            if (required_wat_fullirrig_ww[c]>0) {
              # how much withdrawals can be fulfilled by available water
              frac_fullirrig[c] <- min(avl_wat_ww/required_wat_fullirrig_ww[c],1)

              # consumption constraint
              if (required_wat_fullirrig_wc[c]>0 & length(downstreamcells[[c]])>0) {
                # available water for additional irrigation consumption (considering downstream availability)
                avl_wat_wc          <- max(min(discharge[downstreamcells[[c]]] - required_wat_min_allocation[downstreamcells[[c]]]),0)
                # how much consumption can be fulfilled by available water
                frac_fullirrig[c]   <- min(avl_wat_wc/required_wat_fullirrig_wc[c],frac_fullirrig[c])
                # adjust discharge in current cell and downstream cells (subtract irrigation water consumption)
                discharge[c(downstreamcells[[c]],c)] <- discharge[c(downstreamcells[[c]],c)] - required_wat_fullirrig_wc[c]*frac_fullirrig[c]
              }

              frac_fullirrig[c] <- frac_fullirrig[c]*allocationshare
              # update minimum water required in cell:
              required_wat_min_allocation[c] <- required_wat_min_allocation[c] + frac_fullirrig[c]*required_wat_fullirrig_ww[c]
            }
          }
      } else {
        stop("Please choose allocation rule for river basin discharge allocation algorithm")
      }
    #}

    # update minimum water required in cell:
    #required_wat_min <- required_wat_min + frac_fullirrig*required_wat_fullirrig_ww

      ## Is the following necessary too???
      # update discharge and inflow considering known non-agricultural and committed agricultural uses and full irrigation requirements ###
      # inflow[] <- 0
      # for (o in 1:max(calcorder)){
      #   # Note: the calcorder ensures that the upstreamcells are calculated first
      #   cells <- which(calcorder==o)
      #
      #   for (c in cells){
      #     # available water
      #     avl_wat_act[c] <- inflow[c] + yearly_runoff[c,y] - lake_evap_new[c]
      #     # discharge
      #     discharge[c]   <- avl_wat_act[c] - NAg_wc[c,y,scen]*frac_NAg_fulfilled[c] - CAC_magpie[c]*frac_CAg_fulfilled[c] - required_wat_fullirrig_wc[c]*frac_fullirrig[c]
      #     # inflow into nextcell
      #     if (nextcell[c]>0){
      #       inflow[nextcell[c]] <- inflow[nextcell[c]] + discharge[c]
      #     }
      #   }
      # }

      # # Downstream consideration
      # tmp          <- pmax(avl_wat_act - required_wat_min,0)
      # wat_avl_cons <- numeric(NCELLS)
      #
      # for (c in 1:NCELLS){
      #   # available for consumption in current cell considering downstream cells
      #   wat_avl_cons[c] <- min(tmp[c(downstreamcells[[c]],c)])
      # }


      #################
      #### OUTPUTS ####
      #################
      ### MAIN OUTPUT VARIABLE: water available for irrigation (consumptive agricultural use)
      wat_avl_irrig_c[y,EFR_scen,scen] <- CAC_magpie*frac_CAg_fulfilled + frac_fullirrig*required_wat_fullirrig_wc
      wat_avl_irrig_w[y,EFR_scen,scen] <- CAW_magpie*frac_CAg_fulfilled + frac_fullirrig*required_wat_fullirrig_ww


    }
  }

  ######## REQUIRED OUTPUT FROM CALCAVLWATER:
  # water available for agricultural consumption
  # water available for agricultural withdrawal
  ## each one csv with scenarios:
  # EFP on, off
  # all SSPs of non-agricultural water (currently: 3)
  # for all years
  # option to return cellular or cluster


    # Visualization:
    # frac_NAg_fulfilled_out <- array(data=NA,dim=NCELLS,dimnames=list(names(EFR_magpie)))
    # frac_CAg_fulfilled_out <- array(data=NA,dim=NCELLS,dimnames=list(names(EFR_magpie)))
    #
    # plotmap2(mrmagpie:::toolLPJarrayToMAgPIEmap(frac_NAg_fulfilled),lowcol="red",highcol="white")
    # plotmap2(mrmagpie:::toolLPJarrayToMAgPIEmap(frac_CAg_fulfilled),lowcol="red",highcol="white")
    # plotmap2(mrmagpie:::toolLPJarrayToMAgPIEmap(ratio_routing1))
    # plotmap2(mrmagpie:::toolLPJarrayToMAgPIEmap(ratio_routing2))
    # plotmap2(mrmagpie:::toolLPJarrayToMAgPIEmap(ratio_routing3))


  # How much water available for withdrawals per cell?
  # wat_reserved_withdrawal <- NAg_ww[,y,scen]*frac_NAg_fulfilled + CAW_magpie*frac_CAg_fulfilled
  #wat_avl_withdrawal      <- discharge + CAC_magpie*frac_CAg_fulfilled - EFR_magpie


  ## Diff: wat_avl_consumption3-wat_avl_consumption2; wat_avl_consumption2-wat_avl_consumption1; wat_avl_consumption1-1)
  # (as fraction to discharge_nat)
  # always relative to discharge_nat

  return(list(
    x=out,
    weight=NULL,
    unit="mio. m^3",
    description=description,
    isocountries=FALSE))
}
