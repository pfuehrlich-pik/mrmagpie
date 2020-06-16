#' @title readGCMClimate
#' @description Read Climate data used as LPJmL inputs into MAgPIE objects
#' @param subtype Switch between different inputs
#' @return MAgPIE objects with results on cellular level.
#' @author Marcos Alves, Kristine Karstens
#' @seealso
#' \code{\link{readGCMClimate}}
#' @examples
#'
#' \dontrun{
#' readSource("GCMClimate", subtype="HadGEM2:rcp85.temperature", convert="onlycorrect")
#' }
#'
#' @importFrom lpjclass read.LPJ_input
#' @importFrom madrat toolSubtypeSelect
#' @export

readGCMClimate <-
  function(subtype = "GCMClimate:rcp85:HadGEM2.temperature") {
    if (grepl("\\.", subtype) & grepl("\\:", subtype)) {
      type     <- strsplit(gsub(":", "/" , subtype), split = "\\.")
      folder      <- unlist(type)[1]
      subtype     <- unlist(type)[2]

    } else {

      stop("Subtype did not define either: RCP, climate model or file. Please check the defintion of your subtype parameter.")

    }

    if (dir.exists(folder)) {
      files_list <- list.files(folder)
      files <-
        c(temperature           = files_list[grep("tas", files_list)],
          precipitation         = files_list[grep("pr", files_list)],
          longwave_radiation    = files_list[grep("lwnet", files_list)],
          shortwave_radiation   = files_list[grep("rsds", files_list)],
          wetdays               = files_list[grep("wet", files_list)])
    } else {

      stop(paste("Path", folder, "does not exist. Check the defition of your subtype or the folder structure you are trying to access."))

    }

    file_name <- toolSubtypeSelect(subtype, files)

    years <- as.numeric(unlist(regmatches(file_name, gregexpr("\\d{4}",file_name))))
    years <- seq(years[1],years[2],1)

    if (subtype %in% c("temperature", "precipitation", "wetdays")) {
      x <- read.LPJ_input(file_name = paste0(folder, "/", file_name),
                          out_years = paste0("y", years),
                          namesum = TRUE)

      x <- collapseNames(as.magpie(x))
      x <- x / 12
      col_name <-  unlist(strsplit(unlist(type)[1], split = "\\/"))[c(1,2)]
      getNames(x) <- paste(subtype,col_name[1],col_name[2], sep = "_")

    } else if (subtype %in% c("longwave_radiation", "shortwave_radiation")) {
      x <- read.LPJ_input(file_name = paste0(folder, "/", file_name),
                          out_years = paste0("y", years),
                          namesum = TRUE)
      x <- collapseNames(as.magpie(x))
      x <- x / 365
      col_name <-  unlist(strsplit(unlist(type)[1], split = "\\/"))[c(1,2)]
      getNames(x) <- paste(subtype,col_name[1],col_name[2], sep = "_")

    } else {
      stop(paste0("subtype ", subtype, " is not existing"))
    }

    return(x)

  }
