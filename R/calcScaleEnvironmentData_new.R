#' @title calcScaleEnvironmentData_new
#' @description Scale climate, CO2 and soil environmental conditions on cellular level
#' @param subtype Switch between different climate scenarios
#' @param sar Average range for smoothing annual variations
#' @param sel_feat features names to be included in the output file
#' @param aggr aggregation level
#' @return magpie object in cellular resolution
#' @author Marcos Alves
#'
#' @examples
#' \dontrun{
#' calcOutput("ScaleEnvironmentData_new", climatetype = "HadGEM2_ES:rcp8p5:co2", sar = 20, sel_feat)
#' }
#'
#' @import madrat
#' @import magclass
#' @import mstools
#' @importFrom raster rasterFromXYZ
#' @importFrom raster area
#' @importFrom raster "crs<-"
#' @importFrom magpiesets findset
#' @importFrom stats sd
#'

calcScaleEnvironmentData_new <- function(subtype="ISIMIP3bv2:IPSL-CM6A-LR:ssp126:1965-2100", aggr = F, sar = 20, sel_feat = c(
  "tas",
  "pr",
  "lwnet",
  "rsds",
  "CO2",
  "Ks",
  "Sf",
  "w_pwp",
  "w_fc",
  "w_sat",
  "hsg"
)) {

  # The dataset will the randomized after it is merged with labels.
  # for that reason, it is being scaled with mean and sd from the hole dataset
  x <- calcOutput("CollectEnvironmentData_new", subtype = subtype, sar = sar, aggregate = aggr, sel_feat = sel_feat)
  # xmeans <- apply(x, 3, mean)
  # xstd <-  apply(x, 3, sd)
  # y <- (x - as.magpie(xmeans))/as.magpie(xstd)

  xmax <- as.magpie(apply(x, 3, min))
  xmin <- as.magpie(apply(x, 3, max))
  y <- (x - xmin)/(xmax - xmin)

  # Calculating weights
  landcoords <- as.data.frame(toolGetMapping("magpie_coord.rda", type = "cell"))
  landcoords <- cbind(landcoords, rep(1,nrow(landcoords)))
  landcoords <- raster::rasterFromXYZ(landcoords)
  crs(landcoords) <- "+proj=longlat"
  cell_size <- raster::area(landcoords)
  weight <- cell_size*landcoords
  weight <- as.magpie(weight)
  weight <- toolOrderCells(collapseDim(addLocation(weight),dim=c("x","y")))

  return(list(
    x = y,
    weight = weight,
    unit =
      "temperature: Degree Celcius,
    precipitation: mm3 per year,
    longwave_radiation: watt per m2,
    shortwave_radiation: watt per m2,
    wetdays: day,
    Ks: mm/h, Sf: mm ,
    w_pwp: % ,
    w_fc: % ,
    w_sat: % ,
    tdiff0: mm^2/s ,
    tdiff15: tmm^2/s ,
    tdiff100: mm^2/s ,
    cond_pwp:W/m^2/K) ,
    cond_100: W/m^2/K) ,
    cond_100_ice: W/m^2/K)",
    description = "Climate, CO2 and soil characteristics on cellular level",
    isocountries = FALSE
  ))
}
