#' @title calcTopsoilCarbon
#' @description This function extracts topsoil carbon densities from LPJ to MAgPIE
#'
#' @param version Switch between LPJmL4 and LPJmL4 (or LPJmL4+5)
#' @param climatetype Switch between different climate scenarios (default: "CRU_4")
#' @param time average, spline or raw (default)
#' @param averaging_range just specify for time=="average": number of time steps to average
#' @param dof             just specify for time=="spline": degrees of freedom
#' @param harmonize_baseline FALSE (default) nothing happens, if a baseline is specified here data is harmonized to that baseline (from ref_year on)
#' @param ref_year just specify for harmonize_baseline != FALSE : Reference year
#'
#' @return magpie object in cellular resolution
#' @author Kristine Karstens
#'
#' @examples
#' \dontrun{ calcOutput("TopsoilCarbon", aggregate = FALSE) }
#'
#' @importFrom magpiesets findset

calcTopsoilCarbon <- function(version="LPJmL4", climatetype="CRU_4", time="raw", averaging_range=NULL, dof=NULL,
                       harmonize_baseline=FALSE, ref_year="y2015"){

  if(version == "LPJmL4+5") version <- "LPJmL4"

  soilc_layer_natveg <-  calcOutput("LPJmL", version=version, climatetype=climatetype, subtype="soilc_layer",
                              time=time, dof=dof, harmonize_baseline=harmonize_baseline, ref_year=ref_year,
                              averaging_range=averaging_range,aggregate=FALSE)
  topsoilc           <- soilc_layer_natveg[,,1] + 1/3 * soilc_layer_natveg[,,2]
  getNames(topsoilc) <- "topsoilc"

  # Check for NAs
  if(any(is.na(topsoilc))){
    stop("produced NA Carbon")
  }

  weight <- calcOutput("LanduseInitialisation", aggregate=FALSE, cellular=TRUE, land="fao", input_magpie=TRUE, years="y1995", round=6)

  return(list(
    x=topsoilc,
    weight=weight,
    unit="t per ha",
    description="Topsoil carbon in tons per hectar for natural vegetation.",
    isocountries=FALSE))
}
