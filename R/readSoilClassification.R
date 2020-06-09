#' @title readSoilClassification
#' @description Read soil classification data used as input for lpjml
#' @param subtype Switch between different inputs
#' @return Magpie object with results on cellular level for soil types
#' @author Marcos Alves, Kristine Karstens
#' @examples
#'
#' \dontrun{
#' readSource("SoilClassification", subtype="HWSD.soil", convert="onlycorrect")
#' }
#'
#' @import madrat
#' @import magclass
#' @importFrom dplyr left_join

readSoilClassification <-
  function(subtype = "HWSD.soil") {
    if (grepl("\\.", subtype) & !grepl("\\:", subtype)) {
      subtype     <- strsplit(subtype, split = "\\.")
      folder      <- unlist(subtype)[1]
      subtype     <- unlist(subtype)[2]

    } else {
      stop(
        "Subtype did not define either: soil clasification source or file.
            Please check the defintion of your subtype parameter."
      )

    }

    if (dir.exists(folder)) {
      files_list <- list.files(folder)
      files <- c(soil = files_list[grep("soil", files_list)])
    } else {
      stop(paste("Path", folder, "does not exist. Check the defition of your",
          "subtype or the folder structure you are trying to access.")
      )
    }

    file_name <- toolSubtypeSelect(subtype, files)

    sk <- file(file.path(folder, file_name), "rb")
    y <- readBin(sk, integer(), n = 67460, size = 1)
    close(sk)
    y = y[which(magclassdata$grid_67420_59199 == 1)]
    y = y[magclassdata$cellbelongings$LPJ.Index]

    x <-  collapseNames(as.magpie(y, spatial = 1))
    getNames(x) <- subtype

    return(x)

  }
