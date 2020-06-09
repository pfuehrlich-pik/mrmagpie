#' @title readWATERGAP
#' @description Read in non-agricultural water demand data from WATERGAP model
#' @param subtype Data source to be read from
#' @return MAgPIE object of non-agricultural water demand at 0.5 cellular level in mio. m^3
#' @author Felicitas Beier, Abhijeet Mishra
#'
#' @examples
#' \dontrun{ readSource("WATERGAP", convert="onlycorrect") }
#'
#' @import madrat
#' @import magclass

readWATERGAP <- function(subtype="WATCH_ISIMIP_WATERGAP"){

  # Note to "WATCH_IMAGE_WATERGAP": old non-agricultural waterdemand data (will be deleted soon!)
  if(subtype=="WATCH_ISIMIP_WATERGAP"){
    # Non-agricultural water demands (withdrawals) from WATCH, ISIMIP and WATERGAP
      # A2: WATERGAP WATCH project; B1: WATERGAP WATCH project; SSP2: WATERGAP ISIMIP project
    x <- read.magpie(paste0(subtype,"/watdem_nonagr_0.5.mz"))

  } else if(subtype=="WATERGAP2020"){
    ### List of input files
    input <- list()
    input[["ssp1"]] <- list()
    input[["ssp1"]][["wc"]] <- list()
    input[["ssp1"]][["wc"]] <- c("watergap_ssp1_rcp4p5_elec_wc_annual_2005_2100.nc","watergap_ssp1_rcp4p5_man_wc_annual_2005_2100.nc","watergap_ssp1_rcp4p5_dom_wc_annual_2005_2100.nc")
    input[["ssp1"]][["ww"]] <- gsub("wc","ww",input[["ssp1"]][["wc"]])
    input[["ssp2"]][["wc"]] <- gsub("ssp1_rcp4p5","ssp2_rcp6p0",input[["ssp1"]][["wc"]])
    input[["ssp2"]][["ww"]] <- gsub("ssp1_rcp4p5","ssp2_rcp6p0",input[["ssp1"]][["ww"]])
    input[["ssp3"]][["wc"]] <- gsub("ssp1_rcp4p5","ssp3_rcp6p0",input[["ssp1"]][["wc"]])
    input[["ssp3"]][["ww"]] <- gsub("ssp1_rcp4p5","ssp3_rcp6p0",input[["ssp1"]][["ww"]])

    ### Reading in files and combining to one magpie object:
    x           <- read.magpie(paste0(subtype,"/",input[[1]][["wc"]][1]))
    getNames(x) <- paste0("sspX.",getNames(x))
    # Different SSPs:
    for (i in (1:length(input))){
      # Different water use types (withdrawal, consumption)
      for (j in (1:length(input[["ssp1"]]))){
        # Different industries (manufacturing, electricity, domestic)
        for (k in (1:length(input[["ssp1"]][["wc"]]))){
          tmp <- read.magpie(paste0(subtype,"/",input[[i]][[j]][k]))
          getNames(tmp) <- paste0("ssp",i,".",getNames(tmp))
          x <- mbind(x,tmp)
        }
      }
    }
    # Remove redundant scenario (was for temporary use in loop only)
    x <- x[,,"sspX",invert=T]

    ### Correct years dimension:
    years       <- as.integer(gsub("y","",getYears(x)))
    # Provided WATERGAP data starts with 2005:
    start_year  <- 2005
    years       <- years + start_year
    years       <- paste0("y",years)
    getYears(x) <- years

    # Unit transformation (from m3/yr to mio. m3/yr):
    x <- x/1000000

    ### Sum up over all non-agricultural water uses (domestic, industry)
    # water withdrawal:
    ww           <- dimSums(mbind(x[,,"elecww"],x[,,"domww"],x[,,"manww"]),dim=3.2)
    getNames(ww) <- paste0(getNames(ww),".withdrawal")
    # water consumption:
    wc           <- dimSums(mbind(x[,,"elecuse"],x[,,"domuse"],x[,,"manuse"]),dim=3.2)
    getNames(wc) <- paste0(getNames(wc),".consumption")
    x            <- mbind(ww,wc)
  }

  return(x)
}
