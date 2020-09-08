#' @title readLeifeld2018
#' @description read potential peatland area from Leifeld2018
#' @return List of magpie objects with results on cellular level, weight, unit and description.
#' @author Florian Humpenoeder
#' @examples
#'
#' \dontrun{
#'   readSource("Leifeld2018", convert="onlycorrect")
#' }
#' @importFrom raster raster projectRaster area writeRaster
#' @importFrom magclass read.magpie

readLeifeld2018 <- function(){
  
  x <- raster("Degradation_raster_homolosine_hires_rev4.tif",crs="+proj=igh")
  #crs(x) <- "+proj=igh"
  r <- raster(res=0.5)
  rp2 <- suppressWarnings(projectRaster(x,r,over=TRUE)) #re-project to regular grid
  aa <- area(rp2,na.rm=TRUE) #get area
  writeRaster(aa/10000,"Potential_Peatland_Area_0.5_Leifeld_et_al_2018.nc",overwrite=TRUE,varname="PotentialPeatlandArea", varunit="Mha",longname="Potential Peatland Area in 0.5-degree resolution based on Leifeld_et_al_2018")
  
  x <- read.magpie("Potential_Peatland_Area_0.5_Leifeld_et_al_2018.nc")[,,"PotentialPeatlandArea"]
  
  return(x)
}
