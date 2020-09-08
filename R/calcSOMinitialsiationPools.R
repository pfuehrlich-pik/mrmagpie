#' @title calcSOMinitialsiationPools
#' @description calculates Soil Organic Matter Pool, accounting for the management history as initialisation to magpie
#'
#' @return List of magpie object with results on country or cellular level, weight on cellular level, unit and description.
#' @author Benjamin Leon Bodirsky, Kristine Karstens
#'
#' @seealso
#' \code{\link[mrcommons]{calcSOM}}
#'
#' @examples
#' \dontrun{
#' calcOutput("SOMinitialsiationPools")
#' }
#'

calcSOMinitialsiationPools<-function(){

  past <- findset("past")
  som  <- calcOutput("SOM", aggregate=FALSE)
  som  <- collapseNames(som[,past,c("soilc")])

  return(list(
    x=som,
    weight=NULL,
    unit="Mt C",
    description="Soil carbon in cropland and non-cropland soils.",
    isocountries=FALSE,
    min=0,
    max=1000
  ))
}
