#' @title readGFAD
#' @description Read GLobal Forest Age Dataset derived from MODIS and COPENICUS satellite data
#'
#' @return magpie object in cellular resolution
#' @author Abhijeet Mishra, Felicitas Beier
#' @examples
#'
#' \dontrun{
#'   readSource("GFAD", convert="onlycorrect")
#' }
#' @importFrom raster brick subset as.matrix t
#' @importFrom ncdf4 nc_open ncvar_get
#' @importFrom madrat toolGetMapping
#' @importFrom magclass clean_magpie add_dimension setYears getSets
#' @export

readGFAD <- function(){
  poulter_data <- "GFAD_V1-1.nc"
  #poulter_data <- "C:/PIK/data_processing/inputdata/sources/GFAD/GFAD_V1-1.nc"
  #poulter_data <- "/p/projects/rd3mod/inputdata/sources/GFAD/GFAD_V1-1.nc"

  gfad <- nc_open(poulter_data)

  mapping   <- toolGetMapping(type="cell",name="CountryToCellMapping.csv")

  lon <- ncvar_get(gfad,"lon")
  lat <- ncvar_get(gfad,"lat")

  age_class <- paste0("X",1:15)
  forest_poulter <- c("NEEV","NEDE","BREV","BRDC")

  temp <- NULL

  for(forest_type in (1:length(forest_poulter))){
    for (i in age_class) {
      b <- brick(poulter_data,level=forest_type)
      b <- subset(b,i)

      cellNames <- mapping$celliso

      #Change longitude and latitude
      matrix <- t(as.matrix(b))

      #Create array for 59199 isocells, 1 year and 1 data dimension
      mag   <- array(NA,dim=c(59199,1,1),dimnames=list(cellNames,"y2010",i))

      #Fill array with data from raster object (note magpie_coord are

      for (j in 1:59199) {
        mag[j,,] <- matrix[which(magpie_coord[j, 1]==lon),
                           which(magpie_coord[j,2]==lat)]
      }

      #Convert array to magpie object and rename set
      x          <- clean_magpie(as.magpie(mag))
      getSets(x) <- c("cell","t","ac")
      x <- add_dimension(x = x,dim = 3.1,add = "type",nm = forest_poulter[forest_type])

      cat("Forest",getNames(x,dim=1),"ac",getNames(x,dim=2),"\n")

      temp <- mbind(temp,x)
    }
  }

  out <- setYears(temp,NULL)


  return(out)
}
