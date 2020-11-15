#' @title calcWaterAllocation
#' @description This function calculates water availability for MAgPIE retrieved from LPJmL using a river routing and allocation algorithm for distribution of discharge within the river basin
#'
#' @param version     Switch between LPJmL4 and LPJmL5
#' @param climatetype Switch between different climate scenarios (default: "CRU_4")
#' @param time            Time smoothing: average, spline or raw (default)
#' @param averaging_range only specify if time=="average": number of time steps to average
#' @param dof             only specify if time=="spline": degrees of freedom needed for spline
#' @param harmonize_baseline FALSE (default): no harmonization, TRUE: if a baseline is specified here data is harmonized to that baseline (from ref_year on)
#' @param ref_year           Reference year for harmonization baseline (just specify when harmonize_baseline=TRUE)
#' @param selectyears Years to be returned
#' @param output      Water availability output to be returned: withdrawal or consumption
#' @param allocationrule  Rule to be applied for river basin discharge allocation across cells of river basin ("optimization" (default), "upstreamfirst", "equality")
#' @param allocationshare Share of water to be allocated to cell (only needs to be selected in case of allocationrule=="equality")
#' @param gainthreshold   Threshold of yield improvement potential required for water allocation in upstreamfirst algorithm (in tons per ha)
#' @param irrigationsystem Irrigation system to be used for river basin discharge allocation algorithm ("surface", "sprinkler", "drip", "initialization")
#' @param irrigini         When "initialization" selected for irrigation system: choose initialization data set for irrigation system initialization ("Jaegermeyr_lpjcell", "LPJmL_lpjcell")
#'
#' @import magclass
#' @import madrat
#' @importFrom mrcommons toolHarmonize2Baseline
#'
#' @return magpie object in cellular resolution
#' @author Felicitas Beier, Jens Heinke
#'
#' @examples
#' \dontrun{ calcOutput("WaterAllocation", aggregate = FALSE) }
#'

calcWaterAllocation <- function(selectyears="all", output="consumption",
                         version="LPJmL4", climatetype="HadGEM2_ES:rcp2p6:co2", time="raw", averaging_range=NULL, dof=NULL,
                         harmonize_baseline=FALSE, ref_year="y2015",
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
  magpie2lpj    <- magclassdata$cellbelongings$LPJ_input.Index
  lpj_cells_map <- toolGetMapping("LPJ_CellBelongingsToCountries.csv", type="cell")

  ### Required inputs for River Routing:
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

  ### Required inputs for Allocation Algorithm:
  # Required water for full irrigation per cell (in mio. m^3)
  required_wat_fullirrig_ww <- calcOutput("FullIrrigationRequirement", version="LPJmL5", selectyears=seq(1995, 2095,by=5), climatetype="HadGEM2_ES:rcp2p6:co2", harmonize_baseline=FALSE, time="spline", dof=4, iniyear=1995, iniarea=TRUE, irrig_requirement="withdrawal", cells="lpjcell", aggregate=FALSE)[,,c("maiz","rapeseed","puls_pro")]
  required_wat_fullirrig_ww <- pmax(required_wat_fullirrig_ww,0)
  required_wat_fullirrig_wc <- calcOutput("FullIrrigationRequirement", version="LPJmL5", selectyears=seq(1995, 2095,by=5), climatetype="HadGEM2_ES:rcp2p6:co2", harmonize_baseline=FALSE, time="spline", dof=4, iniyear=1995, iniarea=TRUE, irrig_requirement="consumption", cells="lpjcell", aggregate=FALSE)[,,c("maiz","rapeseed","puls_pro")]
  required_wat_fullirrig_wc <- pmax(required_wat_fullirrig_wc,0)

  # Full irrigation water requirement depending on irrigation system in use
  if (irrigationsystem=="initialization") {
    # read in irrigation system area initialization [share of AEI by system]
    tmp               <- calcOutput("IrrigationSystem", source=irrigini, aggregate=FALSE)
    irrigation_system <- new.magpie(getCells(tmp), getYears(required_wat_fullirrig_ww), getNames(tmp))
    getYears(irrigation_system) <- getYears(required_wat_fullirrig_ww)
    for (y in getYears(required_wat_fullirrig_ww)) {
      irrigation_system[,y,] <- tmp
    }
    # full irrigation water requirements (in mio. m^3)
    required_wat_fullirrig_ww   <- dimSums(irrigation_system*required_wat_fullirrig_ww,dim=3.1)
    required_wat_fullirrig_wc   <- dimSums(irrigation_system*required_wat_fullirrig_wc,dim=3.1)
  } else {
    # whole area irrigated by one system as selected in argument "irrigationsystem"
    required_wat_fullirrig_ww <- collapseNames(required_wat_fullirrig_ww[,,irrigationsystem])
    required_wat_fullirrig_wc <- collapseNames(required_wat_fullirrig_wc[,,irrigationsystem])
  }
  # average required water for full irrigation across selected proxy crops
  required_wat_fullirrig_ww <- dimSums(required_wat_fullirrig_ww,dim=3)/length(getNames(required_wat_fullirrig_ww))
  required_wat_fullirrig_wc <- dimSums(required_wat_fullirrig_wc,dim=3)/length(getNames(required_wat_fullirrig_wc))

  # transform to array for further calculations
  required_wat_fullirrig_ww <- as.array(collapseNames(required_wat_fullirrig_ww))[,,1]
  required_wat_fullirrig_wc <- as.array(collapseNames(required_wat_fullirrig_wc))[,,1]

  # Global cell rank based on yield gain potential by irrigation of proxy crops: maize, rapeseed, pulses
  meancellrank <- calcOutput("IrrigCellranking", version="LPJmL5", climatetype="HadGEM2_ES:rcp2p6:co2", time="spline", averaging_range=NULL, dof=4, harmonize_baseline=FALSE, ref_year="y2015",
                             cellrankyear=seq(1995, 2095,by=5), cells="lpjcell", crops="magpie", method="meancroprank", proxycrop=c("maiz", "rapeseed", "puls_pro"), aggregate=FALSE)
  meancellrank <- as.array(meancellrank)[,,1]


  ############################################
  ####### Routing and Allocation Loop ########
  ############################################
  out_tmp1 <- NULL
  out_tmp2 <- NULL
  out      <- NULL

  for (EFP in c("on", "off")) {

    # Environmental Flow Requirements (in mio. m^3 / yr) [long-term average]
    if (EFP=="on"){
      EFR_magpie <- calcOutput("EnvmtlFlowRequirements", version="LPJmL4", climatetype=climatetype, aggregate=FALSE, cells="lpjcell",
                               LFR_val=0.1, HFR_LFR_less10=0.2, HFR_LFR_10_20=0.15, HFR_LFR_20_30=0.07, HFR_LFR_more30=0.00,
                               EFRyears=c(1980:2010))
    } else if (EFP=="off"){
      EFR_magpie <- new.magpie(1:NCELLS,fill=0)
      getCells(EFR_magpie) <- paste(lpj_cells_map$ISO,1:67420,sep=".")
    }
    EFR_magpie <- as.array(collapseNames(EFR_magpie))
    EFR_magpie <- EFR_magpie[,1,1]

    for (scen in c("ssp1","ssp2","ssp3")){
      for (y in selectyears){

        #############################
        ####### River routing #######
        #############################

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

        ################################################
        ####### River basin discharge allocation #######
        ################################################
        # Minimum water required
        required_wat_min_allocation <- required_wat_min

        # Allocate water for full irrigation to cell with highest yield improvement through irrigation
        if (allocationrule=="optimization") {

          for (c in (1:max(meancellrank[,y],na.rm=T))){
            # available water for additional irrigation withdrawals
            avl_wat_ww <- max(discharge[c]-required_wat_min_allocation[c],0)

            # withdrawal constraint
            if (required_wat_fullirrig_ww[c,y]>0) {
              # how much withdrawals can be fulfilled by available water
              frac_fullirrig[c] <- min(avl_wat_ww/required_wat_fullirrig_ww[c,y],1)

              # consumption constraint
              if (required_wat_fullirrig_wc[c,y]>0 & length(downstreamcells[[c]])>0) {
                # available water for additional irrigation consumption (considering downstream availability)
                avl_wat_wc          <- max(min(discharge[downstreamcells[[c]]] - required_wat_min_allocation[downstreamcells[[c]]]),0)
                # how much consumption can be fulfilled by available water
                frac_fullirrig[c]   <- min(avl_wat_wc/required_wat_fullirrig_wc[c,y],frac_fullirrig[c])
                # adjust discharge in current cell and downstream cells (subtract irrigation water consumption)
                discharge[c(downstreamcells[[c]],c)] <- discharge[c(downstreamcells[[c]],c)] - required_wat_fullirrig_wc[c,y]*frac_fullirrig[c]
              }
              # update minimum water required in cell:
              required_wat_min_allocation[c] <- required_wat_min_allocation[c] + frac_fullirrig[c]*required_wat_fullirrig_ww[c,y]
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
                if (required_wat_fullirrig_ww[c[k],y]>0) {
                  # how much withdrawals can be fulfilled by available water
                  frac_fullirrig[c[k]] <- min(avl_wat_ww/required_wat_fullirrig_ww[c[k],y],1)

                  # consumption constraint
                  if (required_wat_fullirrig_wc[c[k],y]>0 & length(downstreamcells[[c[k]]])>0) {
                    # available water for additional irrigation consumption (considering downstream availability)
                    avl_wat_wc          <- max(min(discharge[downstreamcells[[c[k]]]] - required_wat_min_allocation[downstreamcells[[c[k]]]]),0)
                    # how much consumption can be fulfilled by available water
                    frac_fullirrig[c[k]]   <- min(avl_wat_wc/required_wat_fullirrig_wc[c[k],y],frac_fullirrig[c[k]])
                    # adjust discharge in current cell and downstream cells (subtract irrigation water consumption)
                    discharge[c(downstreamcells[[c[k]]],c[k])] <- discharge[c(downstreamcells[[c[k]]],c[k])] - required_wat_fullirrig_wc[c[k],y]*frac_fullirrig[c[k]]
                  }
                  # update minimum water required in cell:
                  required_wat_min_allocation[c[k]] <- required_wat_min_allocation[c[k]] + frac_fullirrig[c[k]]*required_wat_fullirrig_ww[c[k],y]
                }
              }
            }
          }
        } else if (allocationrule=="equality") {

          # Repeat optimization algorithm several times
          # Instead of full irrigation, only up to x% (e.g.20%) are allocated to most efficient cell
          # Repeat until all river basin discharge is allocated ---> STILL MISSING

          for (c in (1:max(meancellrank[,y],na.rm=T))){
            # available water for additional irrigation withdrawals
            avl_wat_ww <- max(discharge[c]-required_wat_min_allocation[c],0)

            # withdrawal constraint
            if (required_wat_fullirrig_ww[c,y]>0) {
              # how much withdrawals can be fulfilled by available water
              frac_fullirrig[c] <- min(avl_wat_ww/required_wat_fullirrig_ww[c,y],1)

              # consumption constraint
              if (required_wat_fullirrig_wc[c,y]>0 & length(downstreamcells[[c]])>0) {
                # available water for additional irrigation consumption (considering downstream availability)
                avl_wat_wc          <- max(min(discharge[downstreamcells[[c]]] - required_wat_min_allocation[downstreamcells[[c]]]),0)
                # how much consumption can be fulfilled by available water
                frac_fullirrig[c]   <- min(avl_wat_wc/required_wat_fullirrig_wc[c,y],frac_fullirrig[c])
                # adjust discharge in current cell and downstream cells (subtract irrigation water consumption)
                discharge[c(downstreamcells[[c]],c)] <- discharge[c(downstreamcells[[c]],c)] - required_wat_fullirrig_wc[c,y]*frac_fullirrig[c]
              }

              frac_fullirrig[c] <- frac_fullirrig[c]*allocationshare
              # update minimum water required in cell:
              required_wat_min_allocation[c] <- required_wat_min_allocation[c] + frac_fullirrig[c]*required_wat_fullirrig_ww[c,y]
            }
          }
        } else {
        stop("Please choose allocation rule for river basin discharge allocation algorithm")
        }

        #################
        #### OUTPUTS ####
        #################
        ### MAIN OUTPUT VARIABLE: water available for irrigation (consumptive agricultural use)
        wat_avl_irrig_c <- CAC_magpie*frac_CAg_fulfilled + frac_fullirrig*required_wat_fullirrig_wc
        wat_avl_irrig_w <- CAW_magpie*frac_CAg_fulfilled + frac_fullirrig*required_wat_fullirrig_ww

        if (output=="consumption") {
          wat_avl_irrig <- wat_avl_irrig_c
          dataname <- "wat_avl_irrig_c"
        } else if (output=="withdrawal") {
          wat_avl_irrig <- wat_avl_irrig_w
          dataname <- "wat_avl_irrig_w"
        } else {
          stop("specify type of water availability output: withdrawal or consumption")
        }

        wat_avl_irrig <- setNames(setYears(as.magpie(wat_avl_irrig,spatial=1),y), dataname)
        wat_avl_irrig <- add_dimension(wat_avl_irrig, dim=3.1, add="nonag_scen", nm=scen)
        wat_avl_irrig <- add_dimension(wat_avl_irrig, dim=3.1, add="EFP", nm=EFP)
        out_tmp1      <- mbind(out_tmp1, wat_avl_irrig)
      }
      out_tmp2 <- mbind(out_tmp2, out_tmp1)
      out_tmp1 <- NULL
    }
    out      <- mbind(out, out_tmp2)
    out_tmp2 <- NULL
  }

  out <- out[magclassdata$cellbelongings$LPJ_input.Index,,]
  out <- toolCell2isoCell(out)
  description="Available water per year"

  return(list(
    x=out,
    weight=NULL,
    unit="mio. m^3",
    description=description,
    isocountries=FALSE))
}
