#' @title readCO2Atmosphere
#' @description Read CO2 global atmospheric concentration
#' @param subtype Switch between different inputs
#' @return Magpie objects with results on global level
#' @author Marcos Alves, Kristine Karstens
#' @seealso
#' \code{\link{readCO2Atmosphere}}
#' @examples
#'
#' \dontrun{
#' readSource("CO2Atmosphere_new", subtype = "ISIMIP3bv2:ssp126", convert = FALSE)
#' }
#'
#' @import madrat
#' @importFrom utils read.table

readCO2Atmosphere_new <-
  function(subtype = "ISIMIP3bv2:ssp126") {

    s <- toolSplitSubtype(subtype, list(version=NULL, scenario=NULL))

    file <- Sys.glob("*.dat")

    if (length(file) > 1) {
      stop("More than one file was found, please, check the source folder. readCO2Atmosphere")
    }

    y <- read.table(file)
    years <- y[, 1]
    x <- array(NA, dim = c(1, length(years), 1), dimnames = list(1, years, "co2"))
    for (i in 1:length(years)) {
      x[, i, ] <- y[i,2]
    }
    x <- clean_magpie(collapseNames(as.magpie(x, spatial = 1)))
    getNames(x) <- paste("CO2", s$version, s$scenario, sep = "_")
    getCells(x) <- "GLO"

    return(x)
  }
