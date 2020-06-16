#' @title calcSoilCharacteristics
#' @description Calculate Soil Characteristics based on a HWDS soil classification map
#'
#' @return Magpie objects with results on cellular level.
#' @author Marcos Alves
#' @seealso
#' \code{\link{readSoilClassification}},
#' @examples
#'
#' \dontrun{
#'   readSource("SoilClassification", subtype="HWSD.soil", convert="onlycorrect")
#' }
#'
#' @import madrat
#' @importFrom dplyr left_join

calcSoilCharacteristics <- function() {
  x  <- readSource("SoilClassification", subtype = "HWSD.soil", convert="onlycorrect")
  years <- seq(1995, 2100, 1)
  z <- array(NA, dim = c(dim(x)[1], length(years), 1),
             dimnames = list(1:dim(x)[1], years, "soil"))
  for (i in 1:length(years)) {
    z[, i, ] <- x@.Data[, 1,1]
  }

  soil_char <- toolGetMapping(name = "mappingSoil.csv", type = "sectoral")
  w = array(NA, dim = c(dim(x)[1], length(years), dim(soil_char[, -1])[2]),
            dimnames = list(1:dim(x)[1], years, dimnames(soil_char[, -1])[[2]]))

  for (i in 1:length(years)) {
    y           <- as.data.frame(z@.Data[, i, ])
    colnames(y) <- "soil"
    y           <- left_join(y, soil_char, by = "soil", keep = TRUE)
    y           <- data.matrix(y[, -1])
    w[, i, ]    <- y
  }

  w <- as.magpie(w, spatial = 1)
  x <- toolCell2isoCell(w)

  return(list(
    x=x,
    weight=NULL,
    unit=
"Ks: mm/h, Sf: mm ,
 w_pwp: % ,
 w_fc: % ,
 w_sat: % ,
 tdiff0: mm^2/s ,
 tdiff15: tmm^2/s ,
 tdiff100: mm^2/s ,
 cond_pwp:W/m^2/K) ,
 cond_100: W/m^2/K) ,
 cond_100_ice: W/m^2/K)",
    description =
"Ks: saturated hydraulic conductivity (mm/h) following Cosby (1984) ,
 Sf: Suction head (mm) in Green-Ampt equation following Rawls, Brakensiek and Miller (1982) ,
 w_pwp: water content at wilting point following Cosby (1984) ,
 w_fc: water content at field cpacity following Cosby (1984) ,
 w_sat: water content at saturation following Cosby (1984) ,
 tdiff0: thermal diffusivity (mm^2/s) at wilting point (0% whc) following Lawrence and Slater (2008) ,
 tdiff15: thermal diffusivity (mm^2/s) at 15% whc following Lawrence and Slater (2008) ,
 tdiff100: thermal diffusivity (mm^2/s) at field capacity (100% whc) following Lawrence and Slater (2008) ,
 cond_pwp: thermal conductivity (W/m^2/K) at wilting point following Lawrence and Slater (2008) ,
 cond_100: thermal conductivity (W/m^2/K) at saturation (all water) following Lawrence and Slater (2008) ,
 cond_100_ice: thermal conductivity (W/m^2/K) at saturation (all ice) Lawrence and Slater (2008)",
    isocountries=FALSE))
}

