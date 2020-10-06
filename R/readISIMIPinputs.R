#' @title readISIMIPinputs
#' @description Read in non-agricultural water demand data from ISIMIP inputs
#' @param subtype Data source to be read from including path separated by ":", subtype separated by "."
#' @return MAgPIE object of non-agricultural water demand at 0.5 cellular level in mio. m^3
#' @author Felicitas Beier
#'
#' @examples
#' \dontrun{
#' readSource("ISIMIPinputs",
#' subtype="ISIMIP3b:water:histsoc.waterabstraction",
#' convert="onlycorrect")
#' }
#'
#' @import madrat
#' @import magclass

readISIMIPinputs <- function(subtype="ISIMIP3b:water:histsoc.waterabstraction"){

  if(grepl("\\.",subtype)){
    subtype     <- strsplit(gsub(":", "/", subtype), split="\\.")
    folder      <- unlist(subtype)[1]
    time        <- unlist(strsplit(folder, split="\\/"))[3]
    subtype     <- unlist(subtype)[2]
  }

  ### ISIMIP non-agricultural water withdrawals
  if(subtype=="waterabstraction"){

    # Time frame to be read in
    if(time=="histsoc"){
      readyears  <- paste0(time,"_annual_1901_2014")
    } else if(time=="2015soc"){
      readyears  <- paste0(time,"_annual_2015_2100")
    }

    # List of input files
    input <- list()
    input <- c(paste0(folder,"/domww_",readyears,".nc"),
               paste0(folder,"/indwc_",readyears,".nc"),
               paste0(folder,"/indww_",readyears,".nc"))

    # Read in data
    x <- read.magpie(paste0(folder,"/domwc_",readyears,".nc"))
    for (i in (1:length(input))){
      tmp <- read.magpie(input[i])
      x   <- mbind(x,tmp)
    }

    ### Correct years dimension:
    years       <- round(as.numeric(gsub("y","",getYears(x))),0)
    # Provided ISIMIP data starts in year 1901:
    start_year  <- 1901
    years       <- years + start_year
    years       <- paste0("y",years)
    getYears(x) <- years

    # Unit transformation (from m3/yr to mio. m3/yr):
    x <- x/1000000

    ### Sum up over all non-agricultural water uses (domestic, industry)
    # water withdrawal:
    ww           <- dimSums(mbind(x[,,"domww"],x[,,"indww"]),dim=3)
    getNames(ww) <- "withdrawal"
    # water consumption:
    wc           <- dimSums(mbind(x[,,"domwc"],x[,,"indwc"]),dim=3)
    getNames(wc) <- "consumption"
    x            <- mbind(ww,wc)

  }

  return(x)
}
