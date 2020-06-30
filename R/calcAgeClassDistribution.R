#' @title calcAgeClassDistribution
#' @description This function calculates the share of each age class in secondary forests in each MAgPIE simulation cluster based on Global Forest Age Dataset from Poulter et al. 2019
#'
#' @return magpie object in cluster resolution
#' @author Abhijeet Mishra, Felicitas Beier
#'
#' @examples
#' \dontrun{ calcOutput("AgeClassDistribution", aggregate = FALSE) }
#'

calcAgeClassDistribution <- function(){

  poulter_dataset <- readSource("GFAD", convert="onlycorrect")  ## Poulter data is fraction of each cell

  #Area of cells
  mapping   <- toolGetMapping(type="cell",name="CountryToCellMapping.csv")
  lon <- mapping$lon
  lat <- mapping$lat

  cb <- as.data.frame(magpie_coord)
  cell_area  <- (111e3*0.5)*(111e3*0.5)*cos(cb$lat/180*pi)

  cell_area <- as.data.frame(cell_area)
  cell_area$cell <- mapping$cell
  cell_area_magpie <- as.magpie(cell_area[,c(2,1)])
  getNames(cell_area_magpie) <- NULL

  ######################

  getCells(poulter_dataset) <- mapping$cell

  forest_area <- poulter_dataset*cell_area_magpie

  forest_area <- dimSums(forest_area,dim=3.1)

  getNames(forest_area) <- paste0("ac",1:15*10)

  secdf_ac <- getNames(forest_area)[c(-1,-15)] ## Last age-class in primary forest. First age class is very young

  ac_distribution <- forest_area[,,secdf_ac]/dimSums(forest_area[,,secdf_ac],dim=3)
  ac_distribution <- add_columns(x = ac_distribution,addnm = "acx",dim = 3.1) ## Add back 15th age class which is primf. Keep at 0
  ac_distribution[,,"acx"] <- 0
  ac_distribution <- add_columns(x = ac_distribution,addnm = "ac10",dim = 3.1) ## Add back 1st age class which is young. Keep at 0
  ac_distribution[,,"ac10"] <- 0
  ac_distribution <- ac_distribution[,,c(paste0("ac",seq(10,140,10)),"acx")]
  ac_distribution[is.nan(ac_distribution)] <- 1/(length(getNames(ac_distribution))-2) ## Where no forest we don't want to create a mask with 0 so we allow uniform distribtuion

  out <- ac_distribution

  getCells(out) <- mapping$celliso

  names(dimnames(out))[1] <- "ISO.cell"

  weight <- setCells(cell_area_magpie,mapping$celliso)

  return(list(
    x=out,
    weight=weight,
    unit="1",
    description="Fraction of each age class in secondary forest from each spatially explicit cell",
    isocountries=FALSE))
}
