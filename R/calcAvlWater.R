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
  CAW_magpie <- as.array(collapseNames(CAW_magpie))
  CAW_magpie <- CAW_magpie[,1,1]
  CAC_magpie <- as.array(collapseNames(CAC_magpie))
  CAC_magpie <- CAC_magpie[,1,1]

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
      required_wat_min   <- array(data=0,dim=NCELLS,dimnames=list(names(EFR_magpie)))

      ### River Routing 1.1: Natural flows ###
      # Determine natural discharge
      for (o in 1:max(calcorder)){
        # Note: the calcorder ensures that the upstreamcells are calculated first
        cells <- which(calcorder==o)

        for (c in cells){
          ### Natural water balance
          # lake evap that can be fulfilled (if water available: lake evaporation considered; if not: lake evap is reduced respectively):
          lake_evap_new[c] <- min(lake_evap[c,y], inflow_nat[c]+yearly_runoff[c,y])
          # discharge
          discharge_nat[c] <- inflow_nat[c] + yearly_runoff[c,y] - lake_evap_new[c]
          # inflow into nextcell
          if (nextcell[c]>0){
            inflow_nat[nextcell[c]] <- inflow_nat[nextcell[c]] + discharge_nat[c]
          }
        }
      }

      # Minimum availability of water in river to fulfill local EFRs
      required_wat_min <- EFR_magpie

      # Output reporting
      ratio_routing1 <- required_wat_min/discharge_nat
      ratio_routing1[which(required_wat_min==0 & discharge_nat==0)]<-NA
      ratio_routing1[which(ratio_routing1>1)]<-1

      plotmap2(mrmagpie:::toolLPJarrayToMAgPIEmap(ratio_routing1))

      # Downstream consideration
      tmp <-  pmax(discharge_nat - required_wat_min,0)
      wat_avl_consumption1 <- numeric(NCELLS)

      for (c in 1:NCELLS){
        # available for consumption in current cell considering downstream cells
        wat_avl_consumption1[c] <- min(tmp[c(downstreamcells[[c]],c)])
      }

      plotmap2(mrmagpie:::toolLPJarrayToMAgPIEmap(1-(wat_avl_consumption1/discharge_nat)))

      # Check
      sum(required_wat_min-discharge_nat>1e-3)

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
            # -> reduce non-agricultural water use in upstream cells
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

      ### Interim routing: Update discharge and inflow ###
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

      # Check
      sum(required_wat_min-avl_wat_act>1e-3)

      # Output reporting
      ratio_routing2 <- required_wat_min/avl_wat_act
      ratio_routing2[which(required_wat_min==0 & avl_wat_act==0)]<-NA
      ratio_routing2[which(ratio_routing2>1)]<-1

      plotmap2(mrmagpie:::toolLPJarrayToMAgPIEmap(ratio_routing2))

      # Downstream consideration
      tmp <-  pmax(avl_wat_act - required_wat_min,0)
      wat_avl_consumption2 <- numeric(NCELLS)

      for (c in 1:NCELLS){
        # available for consumption in current cell considering downstream cells
        wat_avl_consumption2[c] <- min(tmp[c(downstreamcells[[c]],c)])
      }

      plotmap2(mrmagpie:::toolLPJarrayToMAgPIEmap(1-(wat_avl_consumption2/discharge_nat)))

      # inflow needs to be set to 0 before every river routing
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
            # -> reduce committed agricultural water use in upstream cells
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

      ### Update discharge and inflow ###
      inflow[] <- 0

      for (o in 1:max(calcorder)){
        # Note: the calcorder ensures that the upstreamcells are calculated first
        cells <- which(calcorder==o)

        for (c in cells){
          # available water
          avl_wat_act[c] <- inflow[c] + yearly_runoff[c,y] - lake_evap_new[c]
          # discharge
          discharge[c] <- avl_wat_act[c] - NAg_wc[c,y,scen]*frac_NAg_fulfilled[c] - CAC_magpie[c]*frac_CAg_fulfilled[c]
          # inflow into nextcell
          if (nextcell[c]>0){
            inflow[nextcell[c]] <- inflow[nextcell[c]] + discharge[c]
          }
        }
      }

      # Check
      sum(required_wat_min-avl_wat_act>1e-3)

      # Output reporting
      ratio_routing3 <- required_wat_min/avl_wat_act
      ratio_routing3[which(required_wat_min==0 & avl_wat_act==0)]<-NA
      ratio_routing3[which(ratio_routing3>1)]<-1

      plotmap2(mrmagpie:::toolLPJarrayToMAgPIEmap(ratio_routing3))

      # Downstream consideration
      tmp <-  pmax(avl_wat_act - required_wat_min,0)
      wat_avl_consumption3 <- numeric(NCELLS)

      for (c in 1:NCELLS){
        # available for consumption in current cell considering downstream cells
        wat_avl_consumption3[c] <- min(tmp[c(downstreamcells[[c]],c)])
      }

      plotmap2(mrmagpie:::toolLPJarrayToMAgPIEmap(1-(wat_avl_consumption3/discharge_nat)))

      # How much water available for withdrawals per cell?
     # wat_reserved_withdrawal <- NAg_ww[,y,scen]*frac_NAg_fulfilled + CAW_magpie*frac_CAg_fulfilled
      #wat_avl_withdrawal      <- discharge + CAC_magpie*frac_CAg_fulfilled - EFR_magpie


      ## Diff: wat_avl_consumption3-wat_avl_consumption2; wat_avl_consumption2-wat_avl_consumption1; wat_avl_consumption1-1)
      # (as fraction to discharge_nat)
      # always relative to discharge_nat

      plotmap2(mrmagpie:::toolLPJarrayToMAgPIEmap(wat_avl_consumption/discharge))
      plotmap2(mrmagpie:::toolLPJarrayToMAgPIEmap(wat_avl_consumption))
      plotmap2(mrmagpie:::toolLPJarrayToMAgPIEmap(log(wat_avl_consumption)))




    }

    frac_NAg_fulfilled_out <- array(data=NA,dim=NCELLS,dimnames=list(names(EFR_magpie)))
    frac_CAg_fulfilled_out <- array(data=NA,dim=NCELLS,dimnames=list(names(EFR_magpie)))

    plotmap2(mrmagpie:::toolLPJarrayToMAgPIEmap(frac_NAg_fulfilled),lowcol="red",highcol="white")
    plotmap2(mrmagpie:::toolLPJarrayToMAgPIEmap(frac_CAg_fulfilled),lowcol="red",highcol="white")
    plotmap2(mrmagpie:::toolLPJarrayToMAgPIEmap(ratio_routing1))
    plotmap2(mrmagpie:::toolLPJarrayToMAgPIEmap(ratio_routing2))
    plotmap2(mrmagpie:::toolLPJarrayToMAgPIEmap(ratio_routing3))


    ### Water allocation algorithm for "surplus water" across the river basin ###

    ## basin discharge (discharge at Mündung ist distributed... not river basin runoff)
    ## River basin runoff to be distributed across cells of river basin by algorithm
    # initialize variables
    cell_basin_mapping   <- array(data=0,dim=NCELLS)
    basin_code           <- 1
    basin_runoff         <- array(data=0,dim=length(unique(endcell)))
    reserved_water_basin <- array(data=0,dim=length(unique(endcell)))

    # river basin loop
    for (b in unique(endcell)){
      cell_basin_mapping[which(endcell==b)] <- basin_code
      basin_runoff[basin_code]                  <- sum(yearly_runoff[which(endcell==b),y])
      reserved_water_basin[basin_code]          <- sum(wat_reserved_consumption[which(endcell==b)])
      #surplus_wat_basin[basin_code] <-  basin_runoff - reserved_water_basin

      basin_code=basin_code+1

      # for (c in 1:NCELLS){
      #   avl_water[c] <- surplus_wat_basin * discharge[c]/sum(discharge[c]) + wat_reserved_consumption[c]
      #
      # }
    }

    range(basin_runoff)
    range(basin_runoff-reserved_water_basin)


    for (o in 1:max(calcorder)) {
      cells <- which(calcorder==o)

      for (scen in "ssp2"){
        for (c in cells){

          ## Discharge-weighted distribution
          if (algorithm=="discharge") {
            # Available water per cell
            avl_water[c] <- surplus_wat_basin * discharge[c]/sum(discharge[c]) + reserved_water[c]

            # for(basin in unique(basin_code)){
            #   basin_cells     <- which(basin_code==basin)
            #   basin_runoff    <- colSums(monthly_runoff_magpie[basin_cells,,,drop=FALSE])
            #   basin_discharge <- colSums(monthly_discharge_magpie[basin_cells,,,drop=FALSE])
            #   for(month in dimnames(avl_water_month)[[3]]){
            #     avl_water_month[basin_cells,,month] <- t(basin_runoff[,month]*t(monthly_discharge_magpie[basin_cells,,month])/basin_discharge[,month])
            #   }
            # }
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
