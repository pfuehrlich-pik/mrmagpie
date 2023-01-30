#' @title downloadCO2Atmosphere_new
#' @description Download CO2 atm. inputs used for Lpjml runs
#' @param subtype Switch between different inputs (eg. "ISIMIP3bv2:IPSL-CM6A-LR:historical:1850-2014:tas")
#' It consists of GCM version, climate model, scenario and variable.
#' @return metadata entry
#' @author  Marcos Alves
#' @examples
#' \dontrun{
#' readSource("CO2Atmosphere_new", convert = "onlycorrect")
#' }
#'
downloadCO2Atmosphere_new <- function(subtype = "ISIMIP3bv2:ssp126") { # nolint: object_name_linter.
  ##### CONGIF #######
  substrRight <- function(x, n) {
    substr(x, nchar(x) - n + 1, nchar(x))
  }
  ##### CONFIG #######

  x           <- toolSplitSubtype(subtype, list(version = NULL, scenario = NULL))
  storage     <- "/p/projects/lpjml/input/scenarios/ISIMIP3b" # nolint: absolute_path_linter.
  files       <- grep(".dat", list.files(storage), value = TRUE)
  file        <- grep(substrRight(x$scenario, 2), files, value = TRUE)
  filePath   <- file.path(storage, file)

  if (length(file) > 1) {
    stop("More than one file was found, please, check the source folder")
  }

  if (file.exists(filePath)) {
    file.copy(filePath, file)
  } else {
    stop("Data is not available so far!")
  }

  # Compose meta data
  return(list(url           = paste0(storage, filePath),
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
