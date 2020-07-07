#' @title calcGrowingPeriod
#' @description This function determines a mean sowing date and a mean growing period for each cell
#' in order to determine when irrigation can take place.
#'
#' @param version Switch between LPJmL4 and LPJmL5
#' @param climatetype Switch between different climate scenarios (default: "CRU_4")
#' @param time average, spline or raw (default)
#' @param averaging_range just specify for time=="average": number of time steps to average
#' @param dof             just specify for time=="spline": degrees of freedom
#' @param harmonize_baseline FALSE (default) nothing happens, if a baseline is specified here data is harmonized to that baseline (from ref_year on)
#' @param ref_year just specify for harmonize_baseline != FALSE : Reference year
#' @param yield_ratio threshold for cell yield over global average. crops in cells below threshold will be ignored
#' @param selectyears defaults to all years available
#'
#' @return magpie object in cellular resolution
#' @author Kristine Karstens
#'
#' @examples
#' \dontrun{ calcOutput("GrowingPeriod", aggregate = FALSE) }
#'
#' @importFrom madrat toolGetMapping toolAggregate
#' @importFrom magclass collapseNames getNames new.magpie getYears dimSums
#' @importFrom mrcommons toolHarmonize2Baseline
#'
#' @export

calcGrowingPeriod <- function(version="LPJmL5", climatetype="CRU_4", time="raw", averaging_range=NULL, dof=NULL,
                              harmonize_baseline=FALSE, ref_year="y2015", yield_ratio=0.1, selectyears="all") {

  if(harmonize_baseline==FALSE){

    if(time=="raw"){

      ####################################################################################
      #Goal: calculate mean sowing date and growing period
      #Step 1 Take care of inconsistencies (harvest date or sowing date ==0 etc) and convert to magpie crop functional types
      #Step 2 remove wintercrops from both calculations for the northern hemisphere: sowd>180, hard<sowd
      #Step 3 remove crops that have an irrigated yield below 10% of global average (total cell area as aggregation weight for global value)
      #Step 4 Calculate growing period with the remaining crops
      #Step 5 remove sowd1 for sowing date calculation
      #Step 6 Calculate mean sowing date
      #Step 7 Set the sowd to 1 and growing period to 365 where dams are present and where they are NA (reflecting that all crops have been eliminated)
      #Step 8 Calculate the growing days per month for each cell and each year

      ####################################################################################

      ####################################################################################
      #Read sowing and harvest date input (new for LPJmL5)
      ####################################################################################

      LPJ2MAG      <- toolGetMapping("MAgPIE_LPJmL.csv", type = "sectoral", where = "mappingfolder")

      # Read yields first
      yields       <- collapseNames(calcOutput("Yields", version=version, climatetype=climatetype, time=time, averaging_range=averaging_range, dof=dof,
                                         harmonize_baseline=FALSE, calib_proxy=FALSE, split_cropcalc=FALSE, aggregate = FALSE)[,,"irrigated"])

      # Load Sowing dates from LPJmL (use just rainfed dates since they do not differ for irrigated and rainfed)
      sowd         <- collapseNames(calcOutput("LPJmL", version=version, climatetype=climatetype, subtype="sdate", time="raw",
                                               harmonize_baseline=FALSE, aggregate=FALSE)[,,"rainfed"])
      hard         <- collapseNames(calcOutput("LPJmL", version=version, climatetype=climatetype, subtype="hdate", time="raw",
                                               harmonize_baseline=FALSE, aggregate=FALSE)[,,"rainfed"])

      good_crops   <- LPJ2MAG$MAgPIE[which(LPJ2MAG$LPJmL%in%getNames(sowd))]
      bad_crops    <- LPJ2MAG$MAgPIE[which(!LPJ2MAG$LPJmL%in%getNames(sowd))]

      sowd         <- toolAggregate(sowd, rel=LPJ2MAG, from = "LPJmL", to = "MAgPIE", dim=3.1, partrel = TRUE)
      hard         <- toolAggregate(hard, rel=LPJ2MAG, from = "LPJmL", to = "MAgPIE", dim=3.1, partrel = TRUE)

      if(length(bad_crops)>0) vcat(2,"No information on the growing period found for those crops: ",paste(unique(bad_crops),collapse=", "))

      #####################################################################################

      ####################################################################################
      #Step 1 Take care of inconsistencies (hard==0 etc)
      ####################################################################################

      #Set sowd to 1 and hard to 365 where either sowd or hard are 0
      sowd[which(hard==0|sowd==0)] <- 1
      hard[which(hard==0|sowd==0)] <- 365

      #Set hard to sowd-1 where sowd and hard are equal
      hard[which(hard==sowd & sowd>1)]  <- sowd[which(hard==sowd)]-1 ### WHY???
      hard[which(hard==sowd & sowd==1)] <- 365

      ####################################################################################

      ####################################################################################
      #Step 2 remove crops that have an irrigated yield below 10% of global average (total cell area as aggregation weight)
      ####################################################################################

      area   <- dimSums(calcOutput("LUH2v2", cellular=TRUE, aggregate=FALSE, years="y1995"), dim = 3)
      yields <- collapseNames(yields[,,good_crops])

      cell2GLO     <- array(c(getCells(yields), rep("GLO",59199)), dim=c(59199,2))
      glo_yields   <- toolAggregate(yields, cell2GLO, weight = setYears(area, NULL))
      ratio_yields <- yields/glo_yields

      rm_lowyield   <- yields
      rm_lowyield[] <- 1
      rm_lowyield[ratio_yields < 0.1] <- NA

      rm(ratio_yields,yields,area,glo_yields)

      ####################################################################################

      ####################################################################################
      #Step 3 remove wintercrops from both calculations for the northern hemisphere: sowd>180, hard>365
      ####################################################################################

      cells_northern_hemisphere <- which(magpie_coord[,2] > 0)
      rm_wintercrops            <- new.magpie(getCells(sowd), getYears(sowd), names = getNames(sowd), sets = c("region", "year", "crop"), fill=1)

      #define all crops sowed after 180 days and where sowing date is after harvest date as wintercrops
      rm_wintercrops[cells_northern_hemisphere,,] <-
        ifelse(sowd[cells_northern_hemisphere,,]>180 &
                 hard[cells_northern_hemisphere,,]<sowd[cells_northern_hemisphere,,],          ### WHY???
               NA,1)

      ####################################################################################

      ####################################################################################
      #Step 4 Calculate mean growing period with the remaining crops
      ####################################################################################

      #calculate growing period as difference of sowing date to harvest date
      #sowd <- sowd[,years,]
      #hard <- hard[,years,]
      grow_period <- hard - sowd
      grow_period[which(hard < sowd)] <- 365 + grow_period[which(hard < sowd)]

      #calculate the mean after removing the before determined winter- and low yielding crops
      #rm_wintercrops <- rm_wintercrops[,years,]
      #rm_lowyield    <- rm_lowyield[,years,]

      n_crops          <- dimSums(rm_wintercrops * rm_lowyield, dim=3, na.rm=T)
      mean_grow_period <- dimSums(grow_period * rm_wintercrops * rm_lowyield, dim=3, na.rm=T)/n_crops
      mean_grow_period[is.infinite(mean_grow_period)] <- NA

      #############################################################################

      ####################################################################################
      #Step 5 remove sowd1 for sowing date calculation
      ####################################################################################

      rm_sowd1          <- sowd
      rm_sowd1[]        <- 1
      rm_sowd1[sowd==1] <- NA

      ####################################################################################

      ####################################################################################
      #Step 6 Calculate mean sowing date
      ####################################################################################

      n_crops   <- dimSums(rm_wintercrops * rm_lowyield * rm_sowd1, dim=3, na.rm=T)
      mean_sowd <- dimSums(grow_period * rm_wintercrops * rm_lowyield * rm_sowd1, dim=3, na.rm=T)/n_crops
      mean_sowd[is.infinite(mean_sowd)] <- NA

      rm(rm_wintercrops, rm_lowyield, rm_sowd1)
      rm(sowd, hard, grow_period)
      ####################################################################################

      ####################################################################################
      #Step 7 Set the sowd to 1 and growing period to 365 where dams are present and where they are NA (reflecting that all crops have been eliminated).
      ####################################################################################

      dams <- readSource("Dams", convert="onlycorrect")

      for(t in getYears(mean_sowd)){
        mean_sowd[which(dams==1),t]        <- 1
        mean_grow_period[which(dams==1),t] <- 365
      }

      mean_sowd[is.na(mean_sowd)]               <- 1
      mean_grow_period[is.na(mean_grow_period)] <- 365
      mean_sowd         <- round(mean_sowd)
      mean_grow_period  <- round(mean_grow_period)

      ####################################################################################

      ####################################################################################
      #Step 8 Calculate the growing days per month for each cell and each year
      ####################################################################################

      month        <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
      month_length <- c(   31,   28,   31,   30,   31,   30,   31,   31,   30,   31,   30,   31)
      names(month_length) <- month

      #Determine which day belongs to which month
      days_months        <- 1:365
      names(days_months) <- 1:365

      before <- 0
      for(i in 1:length(month_length)){
        days_months[(before+1):(before+month_length[i])]        <- i
        names(days_months)[(before+1):(before+month_length[i])] <- names(month_length)[i]
        before <- before+month_length[i]
      }

      #mag object for the growing days per month
      grow_days_per_month <- new.magpie(getCells(mean_sowd), getYears(mean_sowd), month, fill=0)

      #determine the harvest day, take care if it is greater than 365
      mean_hard <- (mean_sowd + mean_grow_period - 1)%%365
      mean_hard[mean_hard==0] <- 365

      mean_hard<-as.array(mean_hard)
      mean_sowd<-as.array(mean_sowd)
      grow_days_per_month <- as.array(grow_days_per_month)

      #Loop over the months to set the number of days that the growing period lasts in each month
      for(t in getYears(mean_sowd)){

        #goodcells are cells in which harvest date is after sowing date,
        #i.e. the cropping period does not cross the beginning of the year
        goodcells  <- ifelse(mean_hard[,t,] >= mean_sowd[,t,], 1, 0)
        badcells   <- ifelse(mean_hard[,t,] >= mean_sowd[,t,], 0, 1)

        for(month in 1:12){
          last_monthday  <- which(days_months==month)[length(which(days_months==month))]
          first_monthday <- which(days_months==month)[1]
          test_harvest_goodcells <- as.array(mean_hard[,t,]-first_monthday+1)
          days_in_this_month_goodcells<-as.array(last_monthday-mean_sowd[,t,]+1)
          days_in_this_month_goodcells[days_in_this_month_goodcells<0]<-0 #Month before sowing date
          days_in_this_month_goodcells[days_in_this_month_goodcells>month_length[month]]<-month_length[month] # Month is completely after sowing date
          days_in_this_month_goodcells[test_harvest_goodcells<0]<-0 #Month lies after harvest date
          days_in_this_month_goodcells[test_harvest_goodcells>0 & test_harvest_goodcells<month_length[month]]<-days_in_this_month_goodcells[test_harvest_goodcells>0 & test_harvest_goodcells<month_length[month]]-(last_monthday-mean_hard[test_harvest_goodcells>0 & test_harvest_goodcells<month_length[month],t,]) # Harvest date lies in the month. cut off the end of the month after harvest date
          days_in_this_month_goodcells<-days_in_this_month_goodcells<-days_in_this_month_goodcells*goodcells
          days_in_this_month_badcells_firstyear<-as.array(last_monthday-mean_sowd[,t,]+1)
          days_in_this_month_badcells_firstyear[days_in_this_month_badcells_firstyear<0]<-0 #Month before sowing date
          days_in_this_month_badcells_firstyear[days_in_this_month_badcells_firstyear>month_length[month]]<-month_length[month] # Month is completely after sowing date
          days_in_this_month_badcells_secondyear<-as.array(mean_hard[,t,]-first_monthday+1)
          days_in_this_month_badcells_secondyear[days_in_this_month_badcells_secondyear<0]<-0 #Month lies completely after harvest day
          days_in_this_month_badcells_secondyear[days_in_this_month_badcells_secondyear>month_length[month]]<-month_length[month] #Month lies completely before harvest day
          days_in_this_month_badcells<-(days_in_this_month_badcells_firstyear+days_in_this_month_badcells_secondyear)*badcells

          grow_days_per_month[,t,month]<-days_in_this_month_goodcells+days_in_this_month_badcells
        }
      }


      out <- as.magpie(grow_days_per_month)
      if(any(is.na(out))) stop("Found NAs.")


    } else {

      # Time smoothing:
      x <- calcOutput("GrowingPeriod", version=version, climatetype=climatetype, harmonize_baseline=FALSE,
                      time="raw", yield_ratio=yield_ratio, aggregate = FALSE)

      if(time=="average"){

        # Time smoothing through average:
        out <- toolTimeAverage(x, averaging_range=averaging_range)

      } else if(time=="spline"){

        # Time smoothing with spline method:
        out <- toolTimeSpline(x, dof=dof)
        # Replace value in 2100 with value from 2099 (LPJmL output ends in 2099)
        if ("y2099" %in% getYears(out)) out <- toolFillYears(out, c(getYears(out, as.integer=TRUE)[1]:2100))

        # replace values above days of a month with days of the month & negative values with 0
        month        <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
        month_length <- c(   31,   28,   31,   30,   31,   30,   31,   31,   30,   31,   30,   31)
        names(month_length) <- month
        out[out > as.magpie(month_length)] <- magpie_expand(as.magpie(month_length),out)[out > as.magpie(month_length)]
        out[out < 0] <- 0

      } else if(time!="raw"){
        stop("Time argument not supported!")
      }
    }

  } else {

    if(time=="raw") {
      stop("Harmonization with raw data not possible. Select time='spline' when applying harmonize_baseline=TRUE")
    } else {
      # load smoothed data
      baseline <- calcOutput("GrowingPeriod", version=version, climatetype=harmonize_baseline,
                             harmonize_baseline=FALSE, time=time, dof=dof, averaging_range=averaging_range,
                             yield_ratio=yield_ratio, aggregate = FALSE)
      x        <- calcOutput("GrowingPeriod", version=version, climatetype=climatetype,
                             harmonize_baseline=FALSE, time=time,dof=dof,averaging_range=averaging_range,
                             yield_ratio=yield_ratio, aggregate = FALSE)
      # Harmonize to baseline
      out <- toolHarmonize2Baseline(x=x, base=baseline, ref_year=ref_year, hard_cut=FALSE)

      # replace values above days of a month with days of the month & negative values with 0
      month        <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
      month_length <- c(   31,   28,   31,   30,   31,   30,   31,   31,   30,   31,   30,   31)
      names(month_length) <- month
      out[out > as.magpie(month_length)] <- magpie_expand(as.magpie(month_length),out)[out > as.magpie(month_length)]
      out[out < 0] <- 0
    }
  }

  if(selectyears!="all"){
    years       <- sort(findset(selectyears,noset = "original"))
    out         <- out[,years,]
  }

  return(list(
    x=out,
    weight=NULL,
    unit="d",
    description="Growing days per month in days",
    isocountries=FALSE))

}
