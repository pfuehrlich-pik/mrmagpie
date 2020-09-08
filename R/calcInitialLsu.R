#' @title calcInitialLsu
#' @description Loads the LSU that provides the maximum grass harvest as a initial values for MAgPIE
#' @return MAgPIE objects with optimal lsu on a cellular level.
#' @param model Grass harvest machine learning model ID
#' @author Marcos Alves
#' @examples
#'
#' \dontrun{
#' calOutput("InitialLsu", model = "f41f19be67")
#' }
#'
#' @import madrat
#' @import magpiesets
#' @export


calcInitialLsu <-
  function(model = "f41f19be67") {

    past <- findset("past")
    past <- past[7:length(past)]

    max_lsu           <- toolCell2isoCell(readSource("GrassYldEmu", subtype = paste(model,"max_lsu", sep = "."), convert="onlycorrect"))[,1:length(past),]
    getYears(max_lsu) <- past

    return(list(
      x = max_lsu,
      weight = NULL,
      unit = "lsu/ha",
      description = "Initial livestock densities that yields the maximum harvest per cell",
      isocountries = FALSE
    ))
  }
