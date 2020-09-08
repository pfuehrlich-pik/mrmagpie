#' @title readGLW3
#' @description Read the gridded livestock of the world 3 dataset.
#' @param subtype Subtype of file to be opened (either Da or Aw)
#' @return Magpie objects
#' @author Marcos Alves
#' @examples
#'
#' \dontrun{
#' readSource("GLW3", subtype = "DA", convert="onlycorrect")
#' }
#'
#' @import madrat
#' @importFrom raster aggregate raster rasterToPoints
#'
# setConfig(forcecache=T)
# setConfig(globalenv = T)

readGLW3 <-
  function(subtype = "Da") {


    if (subtype == "Da") {
      str_name <- "5_Ct_2010_Da.tif"
      x <- raster(str_name)
      x <- raster::aggregate(x, fact=6, fun= sum)
      x <- rasterToPoints(x)
      mapping <- toolGetMapping(name = "CountryToCellMapping.csv", type = "cell")
      colnames(x) <- c("lon", "lat", "X5_Ct_2010_Da")
      x <- left_join(mapping,x, by = c("lat", "lon"), copy = TRUE)
      x <- as.magpie(x[,c(2,7)])
    }

    if (subtype == "Aw") {
      str_name <- "6_Ct_2010_Aw.tif"
      x <- raster(str_name)
      x <- raster::aggregate(x, fact=6, fun= sum)
      x <- rasterToPoints(x)
      mapping <- toolGetMapping(name = "CountryToCellMapping.csv", type = "cell")
      colnames(x) <- c("lon", "lat", "X5_Ct_2010_Aw")
      x <- left_join(mapping,x, by = c("lat", "lon"), copy = TRUE)
      x <- as.magpie(x[,c(2,7)])
    }

    return(x)

  }

