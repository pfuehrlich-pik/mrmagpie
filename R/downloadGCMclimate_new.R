#' @title downloadGCMClimate_new
#' @description Download GCM climate input used for Lpjml runs
#' @param subtype Switch between different inputs (eg. "ISIMIP3bv2:IPSL-CM6A-LR:historical:1850-2014:tas")
#' It consists of GCM version, climate model, scenario and variable.
#' @return metadata entry
#' @author  Marcos Alves
#' @examples
#' \dontrun{
#' readSource("GCMClimate_new",  convert = "onlycorrect")
#' }
#'
downloadGCMClimate_new <- function(subtype = "ISIMIP3bv2:IPSL-CM6A-LR:ssp126:2015-2100:tas") {

  x           <- toolSplitSubtype(subtype, list(version = NULL,  climatemodel = NULL,
                                                scenario = NULL, period = NULL,
                                                variable = NULL))

  if (x$climatemodel == "GSWP3-W5E5") {
         storage     <- "/p/projects/lpjml/input/historical/"
  } else storage     <- "/p/projects/lpjml/input/scenarios"

  path        <- file.path(storage,                     # historical or scenarios
                           x$version,                   # ISIMIP3[a|b](v2)
                           gsub("_", "/", x$scenario),  # obsclim e.g. ssp119
                           x$climatemodel)              # GCMs or GSWP3-W5E5

  if (!dir.exists(path)) {
    path      <- file.path(storage, x$version, gsub("_", "/", x$scenario),
                           gsub("_", "-", x$climatemodel))
  }

  list_files  <- list.files(path)
  file        <- grep(paste0(x$variable, "_"), list_files, value = TRUE)
  file        <- grep(x$period, file, value = TRUE)
  file_path   <- file.path(path, file)

  if (file.exists(file_path)) {
    file.copy(file_path, file)
  } else {
    stop("Data is not available so far!")
  }

  .getMetadata <- function(dataset, version) {
    out <- list(doi = NULL, version = NULL, title = NULL, description = NULL)
    return(out)
  }

  meta <- .getMetadata(x$dataset, x$version)

  # Compose meta data
  return(list(url           = paste0(storage, file_path),
              doi           = NULL,
              title         = x$version,
              author        = NULL,
              version       = x$version,
              release_date  = NULL,
              description   = NULL,
              license       = NULL,
              reference     = NULL)
  )
}
