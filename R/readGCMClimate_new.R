#' @title readGCMClimate_new
#' @description Read Climate data used as LPJmL inputs into MAgPIE objects
#' @param subtype Switch between different inputs,
#'                e.g. "ISIMIP3bv2:IPSL-CM6A-LR:historical:1850-2014:tas"
#'                Available variables are: * tas -
#'                                         * wet -
#'                                         * per -
#' @param subset Switch between different subsets of the same subtype
#' @return MAgPIE objects with results on cellular level.
#' @author Marcos Alves, Kristine Karstens, Felicitas Beier
#' @seealso
#' \code{\link{readGCMClimate_new}}
#' @examples
#' \dontrun{
#' readSource("GCMClimate_new", subtype, convert = "onlycorrect")
#' }
#'
#' @importFrom lpjclass read.LPJ_input
#' @importFrom madrat toolSplitSubtype
#' @importFrom magpiesets findset addLocation
#' @importFrom magclass collapseNames collapseDim as.magpie clean_magpie
#' @export

readGCMClimate_new <- function(subtype = "ISIMIP3bv2:IPSL-CM6A-LR:historical:1850-2014:tas",
                               subset  = "annual_mean") {

  subtype         <- toolSplitSubtype(subtype,
                                      list(version      = NULL,
                                           climatemodel = NULL,
                                           scenario     = NULL,
                                           period       = NULL,
                                           variable     = NULL))

  .prepareLPJ_input <- function(subset            = NULL) {

    file_name   <- Sys.glob(c("*.bin", "*.clm"))
    file_type   <- tail(unlist(strsplit(file_name, "\\.")), 1)

    if (file_type == "clm") {

      filedata     <- file(description = file_name, open = "rb",
                           blocking = TRUE, encoding = getOption("encoding"))
      seek(filedata, where = 15, origin = "start")
      in_header    <- as.numeric(readBin(filedata, what = integer(), size = 4,
                                         n = 5, endian = .Platform$endian))
      start_year   <- in_header[1]
      nyear        <- in_header[2]
      number_annual_predictions <- in_header[5]
      years        <- seq(start_year, start_year + nyear - 1, 1)
      close(filedata)

    } else {

      stop("File format of LPJmL input data unknown. Please provide .clm file format.")
    }

    if (subset == "wet") {

      x <- read.LPJ_input(file_name   = file_name,
                          out_years   = paste0("y", years),
                          namesum     = TRUE,
                          ncells      = 67420,
                          rule4binary = ">0") / number_annual_predictions

      class(x) <- "array"
      x        <- collapseNames(as.magpie(x, spatial = 1))


    } else if (subset == "annual_mean") {

      x <- read.LPJ_input(file_name   = file_name,
                          out_years   = paste0("y", years),
                          namesum     = TRUE,
                          ncells      = 67420) / number_annual_predictions

      class(x) <- "array"
      x        <- collapseNames(as.magpie(x, spatial = 1))


    } else if (subset == "annual_sum") {

      x <- read.LPJ_input(file_name   = file_name,
                          out_years   = paste0("y", years),
                          namesum     = TRUE,
                          ncells      = 67420)

      class(x) <- "array"
      x        <- collapseNames(as.magpie(x, spatial = 1))


    } else if (subset %in% c("monthly_mean", "monthly_sum")) {
      # define year sets (cut it in bunches)
      bunchLength <- 1
      yearsets    <- split(years, ceiling(seq_along(years) / bunchLength))

      # define month mapping
      monthLength  <- c(jan = 31, feb = 28, mar = 31, apr = 30,
                        may = 31, jun = 30, jul = 31, aug = 31,
                        sep = 30, oct = 31, nov = 30, dec = 31)
      daysMonth    <- NULL
      for (m in 1:12) daysMonth <- c(daysMonth, rep(names(monthLength[m]),
                                                    monthLength[m]))
      month2day    <- cbind(day   = 1:sum(monthLength),
                            month = daysMonth)

      # create output object
      x <- NULL

      # loop over bunches
      for (b in 1:length(yearsets)) {
        # read in a bunch of years
        tmp <- lpjclass::read.LPJ_input(file_name = file_name,
                                        out_years = paste0("y", yearsets[[b]]),
                                        namesum   = FALSE,
                                        ncells    = 67420)
        class(tmp) <- "array"
        tmp        <- as.magpie(tmp)

        # aggregate days to month
        tmp      <- toolAggregate(tmp,
                                  rel  = month2day,
                                  from = "day",
                                  to   = "month",
                                  dim  = 3)

        if (subset == "monthly_mean") tmp / as.magpie(monthLength)

        x      <- mbind(x, tmp)
      }

    } else if (grepl("\\d{4}:\\d{4}", subset)) {

      subYears <- eval(parse(text = subset))
      years    <- intersect(years, subYears)
      if (any(!(subYears %in% years))) {
        warning(paste0("Some subsetted years (subset = ", subset,
                       ") are not availabl\n in the original data.
                       Years set to:", years))
      }

      x <- read.LPJ_input(file_name   = file_name,
                          out_years   = paste0("y", years),
                          namesum     = FALSE,
                          ncells      = 67420)

      class(x) <- "array"
      x        <- collapseNames(as.magpie(x, spatial = 1))


    } else {
      stop("Subset argument unknown. Please check function help.")
    }

    return(x)
  }

  x        <- .prepareLPJ_input(subset) # maybe add conditionals on
  # which subtype subset combinations
  # should be allowed
  getNames(x) <- paste(subtype$variable, subset, sep = "_")
  x        <- collapseDim(addLocation(x), dim = c("N", "region"))
  x        <- clean_magpie(x)

  return(x)
}
