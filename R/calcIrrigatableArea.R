#' @title calcIrrigatableArea
#' @description This function calculates the area that can be irrigated given available area and available water for all magpie crops
#'
#' @return magpie object in cellular resolution
#' @author Felicitas Beier
#'
#' @param cells Switch between "lpjcell" (67420) and "magpiecell" (59199)
#'
#' @examples
#' \dontrun{ calcOutput("IrrigatableArea",aggregate=FALSE) }
#'
#' @import madrat
#' @import magclass

calcIrrigatableArea <- function(cells="lpjcell"){

  if (cells=="lpjcell") {
    NCELLS <- 67420
  } else if (cells=="magpiecell") {
    NCELLS <- 59199
  } else {
    stop("Cells argument not supported. Please select lpjcell for 67420 cells or magpiecell for 59199 cells")
  }

  systemnames <- c("drip","sprinkler","surface")
  cropnames   <- c("tece","maiz","trce","rice_pro","soybean","rapeseed","groundnut","sunflower","oilpalm","puls_pro","potato","cassav_sp","sugr_cane","sugr_beet","others","cottn_pro","foddr","begr","betr")
  type        <- c("consumption","withdrawal")
  systemnames <- paste(type,rep(systemnames,2),sep=".")

  # empty magpie object
  irrig_area <- new.magpie(1:NCELLS,"y1995",sort(paste(systemnames, rep(cropnames,6), sep=".")),sets=c("iso.cell","year","type.system.crop"))



  # read in land available for agricultural use (in mio. ha)
  land      <- collapseNames(calcOutput("AvlLandSi", aggregate=FALSE)[,,"si0"])

  for (type in c("consumption","withdrawal")) {
    # water requirement for full irrigation per crop per cell (in m^3)
    wat_req_full <- calcOutput("FullIrrigationRequirement", version=version, cells=cells, selectyears=y, climatetype=climatetype, harmonize_baseline=harmonize_baseline, time=time, dof=dof, irrig_requirement=type, aggregate=FALSE)

    # available water per cell (in mio. m^3)
    # read in available water for agricultural use (with argument: consumption vs. withdrawal; irrig_requirement="consumption" vs. "withdrawal)
    # transform to m^3: /1000000
    # [currently: placeholder]
    wat_avl <- new.magpie(1:NCELLS,"y1995",c("consumption","withdrawal"),sets=c("iso.cell","year","type"),fill=1000)


    # Calculations
    wat_req = wat_req_full/land

    if (wat_avl >= wat_req_full) {
      irrig_area <- land
    } else {
      # irrigatable area in ha given available water
      irrig_area <- wat_avl/wat_req
      # transform to mio. ha
      irrig_area <- irrig_area/1000000
    }

    irrig_area[,,type] <- irrig_area
  }

  x <- pmin(irrig_area[,,"consumption"],irrig_area[,,"withdrawal"])

  # Checks
  if(any(is.na(x))){
    stop("produced NA irrigatable area")
  }

  return(list(
    x=x,
    weight=NULL,
    unit="mio. ha",
    description="irrigatable area",
    isocountries=FALSE))
}
