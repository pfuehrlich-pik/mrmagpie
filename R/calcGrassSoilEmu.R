#' @title calcGrassSoilEmu
#' @description Read files related to the training and optimization of the LPJml emulators.
#' @param subtype Subtype of file to be opened. Subtypes available:
#' 'weights', 'inputs',  'stddevs' and 'means'.
#' @param model trained model ID
#' @param mfile model file name
#' @return Magpie objects with a diverse inforamtion
#' @author Marcos Alves
#' @examples
#' \dontrun{
#' readSource("GrassSoilEmu",
#'   subtype = "ISIMIP3bv2:IPSL_CM6A_LR:ssp126:1965_2100",
#'   model = "5f5fa2", mfile = "weights"
#' )
#' }
#'
#' @import madrat
#' @importFrom stats runif
#'

calcGrassSoilEmu <-
  function(subtype = "ISIMIP3bv2:IPSL_CM6A_LR:ssp126:1965_2100", model = "5f5fa2", mfile = "weights") {
    subtype <- paste(subtype, model, mfile, sep = ":")

    x <- readSource("GrassSoilEmu", subtype = subtype, convert = FALSE)

    return(
      list(
        x = x,
        weight = NULL,
        unit = "none",
        description = "soil C ML emulator file",
        isocountries = FALSE
      )
    )
  }
