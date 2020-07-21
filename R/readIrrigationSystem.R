#' @title readIrrigationSystem
#' @description Read in irrigation system type for initialization
#' @param subtype Data source to be read from: LPJmL for LPJmL-dominant-irrigation-system-type or Jaegermeyr for country-specific shares and lpjcell (67420) or magpiecell (59199)
#' @return MAgPIE object of at country-level
#' @author Felicitas Beier
#'
#' @examples
#' \dontrun{ readSource("IrrigationSystem", convert="onlycorrect") }
#'
#' @import madrat
#' @import magclass

readIrrigationSystem <- function(subtype="Jaegermeyr_magpiecell"){

  # Mapping
  map <- toolGetMapping("CountryToCellMapping.csv",type="cell")
  map <- data.frame(celliso=map$celliso, iso=map$iso)

  # Irrigation system share by country
  if (grepl("Jaegermeyr", subtype)){

    # Read in and transform data
    x         <- read.csv("Jaegermeyr-2015-supplement_shr.csv")
    x$Country <- toolCountry2isocode(x$Country,mapping=c("bahamas, the"="BHS", "congo-brazzaville"="COG", "dr congo, former zaire"="COD",
                                            "gambia, the"="GMB", "myanmar (burma)"="MMR", "west bank"="PSE"))
    # Transform to MAgPIE object
    x <- as.magpie(x)
    x <- toolCountryFill(x)

    # Where system share is not given, it is assumed that 100% of irrigated land is surface irrigation
    x[,,"shr_AEI_surface"][is.na(x[,,"shr_AEI_surface"])] <- 1
    x[,,"shr_AEI_sprinkler"][is.na(x[,,"shr_AEI_sprinkler"])] <- 0
    x[,,"shr_AEI_drip"][is.na(x[,,"shr_AEI_drip"])] <- 0

    # Expand to cellular level
    x   <- x[map$iso,,]

  } else if(grepl("LPJmL", subtype)){

    # Read in and transform data
    x         <- read.csv("LPJmL_manage_irrig_IFT.csv")
    x         <- subset(x, x$country!="No_Land")
    x$country <- toolCountry2isocode(x$country,mapping=c("bahamas, the"="BHS", "congo-brazzaville"="COG", "dr congo, former zaire"="COD", "falkland islands (islas malvinas)"="FLK", "western samoa"="WSM",
                                                         "french southern and antarctica lands"="ATF", "gambia, the"="GMB", "myanmar (burma)"="MMR", "west bank"="PSE", "cocos keeling islands"="CCK",
                                                         "pitcairn islands"="PCN", "saint helena ascension and tristan da cunha"="SHN", "st. vincent and the grenadines"="VCT", "virgin islands"="VIR"))
    x         <- subset(x, !duplicated(x))

    # Transform to MAgPIE object
    x <- as.magpie(x)
    x   <- toolCountryFill(x)

    # Expand to cellular level
    x   <- x[map$iso,,]

    # Replace NAs with 1 (surface irrigation)
    # (Note: only affects country SJM: Svalbard and Jan Mayen)
    x <- toolConditionalReplace(x, conditions=c("is.na()","<0"), replaceby = 1)

  }

  # Object dimensions
  if (grepl("lpjcell", subtype)){
    lpj_map <- toolGetMapping("LPJ_CellBelongingsToCountries.csv",type="cell")
    tmp <- x
    getCells(tmp) <- paste(getCells(tmp),1:length(getCells(tmp)),sep=".")
    x <- new.magpie(cells_and_regions=paste(lpj_map$ISO,1:67420,sep="."),years=NULL,getNames(x),fill=0)
    x[magclassdata$cellbelongings$LPJ.Index,,] <- tmp[,,]
    x <- toolCell2isoCell(x,cells="lpjcell")
  } else if (grepl("magpiecell", subtype)){
    x <- toolCell2isoCell(x,cells="magpiecell")
  }

  return(x)
}
