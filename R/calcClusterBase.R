#' @title calcClusterBase
#' @description Reads a series of MAgPIE files and combines them to a matrix
#' which is then used for calculating a clustering.
#' @param years2use A vector with years with should be taken into account for
#' the clustering
#' @return A matrix containing the data
#' @author Jan Philipp Dietrich
#' @seealso \code{\link{calcCluster}}
#' @importFrom magclass wrap read.magpie
#' @importFrom madrat toolMappingFile
calcClusterBase <- function(years2use=1995) {

  d <- list()
  # read in data which should be used to determine cluster
  d$yld    <- calcOutput("Yields", selectyears=years2use, aggregate=FALSE)
  d$airrig <- calcOutput("Irrigation", selectyears=years2use, aggregate=FALSE)
  d$td     <- calcOutput("TransportDistance", aggregate=FALSE)[,,rep(1,16)]

  cdata <- do.call(cbind,lapply(d,wrap,list(1,c(2,3))))
  cdata <- scale(cdata)
  colnames(cdata) <- paste0("i",1:ncol(cdata))

  #make sure that spatial dimension contains country information
  iso <- toolMappingFile("cell", "CountryToCellMapping.csv", readcsv = TRUE)$iso
  dimnames(cdata)[[1]] <- paste(iso, 1:length(iso), sep = ".")

  return(list(
    x = as.magpie(cdata),
    weight = NULL,
    unit = "1",
    description = "Similarity matrix as basis for clustering",
    isocountries = FALSE))
}
