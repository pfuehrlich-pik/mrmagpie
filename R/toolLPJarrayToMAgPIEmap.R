#' @title toolLPJarrayToMAgPIEmap
#' @description transforms object with 67420 cells to be able to be used in plotmap function
#'
#' This tool function selects the required cells from object with
#' 67420 cells to be able to plot using plotmap or plotmap2 that
#' are applied to magpie object of size 59199 cells
#'
#' @return magpie object with correct dimension for use in plotmap
#' @author Felicitas Beier

toolLPJarrayToMAgPIEmap <- function(x){
  # select correct cells and transform to magpie object
  out <- as.magpie(x[magclassdata$cellbelongings$LPJ_input.Index],spatial=1)

  return(out)
}
