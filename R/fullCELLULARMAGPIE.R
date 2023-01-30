#' @title fullCELLULARMAGPIE
#' @description Function that produces the complete cellular data set required for running the
#' MAgPIE model.
#'
#' @param rev data revision which should be used as input (positive numeric).
#' @param ctype aggregation clustering type, which is a combination of a single letter, indicating the cluster
#' methodology, and a number, indicating the number of resulting clusters. Available methodologies are hierarchical
#' clustering (h), normalized k-means clustering (n) and combined hierarchical/normalized k-means clustering (c).
#' In the latter hierarchical clustering is used to determine the cluster distribution among regions whereas
#' normalized k-means is used for the clustering within a region.
#' @param dev development suffix to distinguish development versions for the same data revision.
#' This can be useful to distinguish parallel lines of development.
#' @param climatetype Global Circulation Model to be used
#' @param lpjml Defines LPJmL version for crop/grass and natveg specific inputs
#' @param isimip Defines isimip crop model input which replace maiz, tece, rice_pro and soybean
#' @param emu_id Pasture Soil carbon emulator ID
#' @param clusterweight Should specific regions be resolved with more or less detail? Values > 1 mean higher share,
#' < 1 lower share e.g. cfg$clusterweight <- c(LAM=2) means that a higher level of detail for region LAM if set to NULL
#' all weights will be assumed to be 1.
#' examples:
#' c(LAM=1.5,SSA=1.5,OAS=1.5)
#' c(LAM=2,SSA=2,OAS=2)
#' \code{\link{setConfig}} (e.g. for setting the mainfolder if not already set
#' properly).
#' @author Kristine Karstens, Jan Philipp Dietrich
#' @seealso
#' \code{\link{readSource}},\code{\link{getCalculations}},\code{\link{calcOutput}},\code{\link{setConfig}}
#' @examples
#' \dontrun{
#' retrieveData("CELLULARMAGPIE", revision = 12, mainfolder = "pathtowhereallfilesarestored")
#' }
#' @importFrom madrat setConfig getConfig
#' @importFrom magpiesets findset
#' @importFrom digest digest
#' @importFrom stringr str_split
#' @importFrom luplot plotregionscluster
#' @importFrom ggplot2 ggsave
#' @importFrom withr local_options

fullCELLULARMAGPIE <- function(rev = 0.1, dev = "",
                               ctype = "c200",
                               climatetype = "MRI-ESM2-0:ssp370",
                               lpjml = c(natveg = "LPJmL4_for_MAgPIE_44ac93de",
                                         crop = "ggcmi_phase3_nchecks_9ca735cb",
                                         grass = "lpjml5p2_pasture"),
                               isimip = NULL,
                               clusterweight = NULL,
                               emu_id = NULL) { # nolint

  "!# @pucArguments ctype clusterweight"

  withr::local_options(magclass_sizeLimit = 1e+12)

  ### Version settings ###
  if (rev < 4.66) stop("mrmagpie(>= 1.19.0) does not support revision below 4.66 anymore.
                       Please use a older snapshot/version of the library, if you need older revisions.")

  climatescen <- str_split(climatetype, ":")[[1]][2]

  message(paste0("Start preprocessing for \n climatescenario: ", climatetype,
    "\n LPJmL-Versions: ", paste(names(lpjml), lpjml, sep = "->", collapse = ", "),
    "\n clusterweight: ", paste(names(clusterweight), clusterweight, sep = ":", collapse = ", "),
    "\n isimip yield subtype: ", paste(names(isimip), isimip, sep = ":", collapse = ", ")))

  # Create version tag (will be returned at the very end of this function)
  versionTag <- paste(ctype,
                       gsub(":", "-", climatetype),
                       paste0("lpjml-", digest::digest(lpjml, algo = getConfig("hash"))),
                       sep = "_")
  versionTag <- ifelse(is.null(isimip), versionTag,
                        paste0(versionTag, "_isimip-",
                               digest::digest(isimip, algo = getConfig("hash"))))
  versionTag <- ifelse(is.null(clusterweight), versionTag,
                        paste0(versionTag, "_clusterweight-",
                               digest::digest(clusterweight, algo = getConfig("hash"))))
  versionTag <- ifelse(is.null(emu_id), versionTag,
                        paste0(versionTag, "_gsoilc-", emu_id))


  magYearsPastLong  <- c("y1995", "y2000", "y2005", "y2010", "y2015")
  magYears <- findset("time")
  shortYears <- findset("t_all")
  lpjYears <- seq(1995, 2100, by = 5)

  map      <- calcOutput("Cluster", ctype = ctype, weight = clusterweight, lpjml = lpjml,
    clusterdata = "yield_airrig", aggregate = FALSE)
  weightID <- ifelse(is.null(clusterweight), "", paste0("_", names(clusterweight), clusterweight, collapse = ""))
  clustermapname <- sub("\\.[^.]*$", ".rds", paste0("clustermap_rev", rev, dev, "_", ctype,
    weightID, "_", getConfig("regionmapping")))
  addMapping(clustermapname, map)

  # plot map with regions and clusters
  clustermap <- readRDS(clustermapname)
  p <- plotregionscluster(clustermap$cluster)
  ggsave(sub(".rds", ".pdf", sub("clustermap", "spamplot", clustermapname)), p, height = 6, width = 10, scale = 1)

  # 14 yields

  if (dev == "+calibYield") {

    calcOutput("YieldsCalibrated", aggregate = "cluster", source = c(lpjml = lpjml[["crop"]], isimip = isimip),
               climatetype = climatetype, round = 2, years = lpjYears, file = paste0("lpj_yields_", ctype, ".mz"))

  } else if (grepl("india", dev)) {

    calcOutput("Yields", aggregate = FALSE, source = c(lpjml = lpjml[["crop"]], isimip = isimip),
               climatetype = climatetype, round = 2, years = lpjYears, file = paste0("lpj_yields_0.5.mz"),
               weighting = "crop+irrigSpecific", indiaYields = TRUE, scaleFactor = 0.5)


    calcOutput("Yields", aggregate = "cluster", source = c(lpjml = lpjml[["crop"]], isimip = isimip),
               climatetype = climatetype, round = 2, years = lpjYears, file = paste0("lpj_yields_", ctype, ".mz"),
               weighting = "crop+irrigSpecific", indiaYields = TRUE, scaleFactor = 0.5)

  } else {

    calcOutput("Yields", aggregate = FALSE, source = c(lpjml = lpjml[["crop"]], isimip = isimip),
               climatetype = climatetype, round = 2, years = lpjYears, file = paste0("lpj_yields_0.5.mz"),
               weighting = ifelse(grepl("YieldWeights_", dev), gsub("YieldWeights_", "", dev), "totalCrop"))

    calcOutput("Yields", aggregate = "cluster", source = c(lpjml = lpjml[["crop"]], isimip = isimip),
               climatetype = climatetype, round = 2, years = lpjYears, file = paste0("lpj_yields_", ctype, ".mz"),
               weighting = ifelse(grepl("YieldWeights_", dev), gsub("YieldWeights_", "", dev), "totalCrop"))

  }


  # distinguish between region and superregion if mapping provides this distinction
  mapReg <- toolGetMapping(getConfig("regionmapping"), type = "regional")
  superregion <- ifelse("superregion" %in% colnames(mapReg), "superregion", "region")

  # 13 TC
  calcOutput("PastrTauHist", round = 2, past_mngmt = "mdef", file = "f13_pastr_tau_hist.csv", aggregate = superregion)

  # 09 drivers
  calcOutput("GridPop", source = "Gao", subtype = "all", cellular = TRUE, harmonize_until = 2015, urban = FALSE,
             aggregate = "cluster", years = magYears, round = 6, file = "f09_pop_grid.cs3")

  calcOutput("GridPop", source = "Gao", subtype = "all", cellular = TRUE, harmonize_until = 2015, urban = TRUE,
             aggregate = "cluster", years = magYears, round = 6, file = "f09_urbanpop_grid.cs3")

  # 10 land
  # seven land classes
  calcOutput("LanduseInitialisation", aggregate = FALSE, cellular = TRUE, cells = "magpiecell", nclasses = "seven",
             input_magpie = TRUE, selectyears = magYearsPastLong, round = 6, file = "avl_land_t_0.5.mz")
  calcOutput("LanduseInitialisation", aggregate = "cluster", cellular = TRUE, nclasses = "seven", input_magpie = TRUE,
             selectyears = magYearsPastLong, round = 6, file = paste0("avl_land_t_", ctype, ".mz"))
  calcOutput("LanduseInitialisation", aggregate = FALSE, cellular = FALSE, nclasses = "seven",
             input_magpie = TRUE, selectyears = magYearsPastLong, round = 6, file = paste0("avl_land_t_iso.cs3"))

  # nine land classes
  calcOutput("LanduseInitialisation", aggregate = FALSE, cellular = TRUE, cells = "magpiecell", nclasses = "nine",
            input_magpie = TRUE, selectyears = magYearsPastLong, round = 6, file = "avl_land_full_t_0.5.mz")
  calcOutput("LanduseInitialisation", aggregate = "cluster", cellular = TRUE, nclasses = "nine", input_magpie = TRUE,
             selectyears = magYearsPastLong, round = 6, file = paste0("avl_land_full_t_", ctype, ".mz"))
  calcOutput("LanduseInitialisation", aggregate = FALSE, cellular = FALSE, nclasses = "nine", input_magpie = TRUE,
             selectyears = magYearsPastLong, round = 6, file = paste0("avl_land_full_t_iso.cs3"))

  calcOutput("AvlLandSi", aggregate = FALSE, round = 6, file = "avl_land_si_0.5.mz")
  calcOutput("AvlLandSi", aggregate = "cluster", round = 6, file = paste0("avl_land_si_", ctype, ".mz"))

  # 22 land conservation
  calcOutput("ProtectedAreaBaseline", nclasses = "seven", cells = "magpiecell", magpie_input = TRUE,
    aggregate = FALSE, round = 6, file = "wdpa_baseline_0.5.mz")
  calcOutput("ProtectedAreaBaseline", nclasses = "seven", cells = "magpiecell", magpie_input = TRUE,
    aggregate = "cluster", round = 6, file = paste0("wdpa_baseline_", ctype, ".mz"))
  calcOutput("ConservationPriority", nclasses = "seven", cells = "magpiecell",
             aggregate = FALSE, round = 6, file = "consv_prio_areas_0.5.mz")
  calcOutput("ConservationPriority", nclasses = "seven", cells = "magpiecell",
             aggregate = "cluster", round = 6, file = paste0("consv_prio_areas_", ctype, ".mz"))

  calcOutput("ProtectArea", bhifl = ifelse(rev > 4.66, TRUE, FALSE), aggregate = "cluster", round = 6,
             file = paste0("protect_area_", ctype, ".mz"))

  # 30 crop
  calcOutput("Croparea", sectoral = "kcr", physical = TRUE, cellular = TRUE, irrigation = FALSE,
    aggregate = "cluster", file = paste0("f30_croparea_initialisation_", ctype, ".mz"))
  calcOutput("Croparea", sectoral = "kcr", physical = TRUE, cellular = TRUE, irrigation = TRUE,
    aggregate = "cluster", file = paste0("f30_croparea_w_initialisation_", ctype, ".mz"))

  calcOutput("AvlCropland", marginal_land = "magpie", cell_upper_bound = 0.9,  aggregate = FALSE, round = 6,
    file = "avl_cropland_0.5.mz")
  calcOutput("AvlCropland", marginal_land = "magpie", cell_upper_bound = 0.9, aggregate = "cluster", round = 6,
    file = paste0("avl_cropland_", ctype, ".mz"))
  calcOutput("AvlCropland", marginal_land = "magpie", cell_upper_bound = 0.9, aggregate = FALSE,
    country_level = TRUE, round = 6, file = paste0("avl_cropland_iso.cs3"))

  # 31 past

  calcOutput("GrasslandBiomass",  round = 3, file = "f31_grass_bio_hist.cs3", aggregate = "region")
  calcOutput("LUH2v2", aggregate = "cluster", landuse_types = "LUH2v2", cellular = TRUE,
             file = paste0("f31_LUH2v2_", ctype, ".mz"))
  # hard coded climate scenario for harmonization of data
  calcOutput("GrasslandsYields", lpjml = lpjml[["grass"]], climatetype = "MRI-ESM2-0:ssp126",
             subtype = "/co2/Nreturn0p5", # nolint
             lsu_levels = c(seq(0, 2.2, 0.2), 2.5), past_mngmt = "mdef",
             file = paste0("f31_grassl_yld_", ctype, ".mz"), years = magYears, aggregate = "cluster")
  calcOutput("GrasslandsYields", lpjml = lpjml[["grass"]], climatetype = "MRI-ESM2-0:ssp126",
             subtype = "/co2/Nreturn0p5", # nolint
             lsu_levels = c(seq(0, 2.2, 0.2), 2.5), past_mngmt = "mdef",
             file = paste0("f31_grassl_yld.mz"), years = magYears, aggregate = FALSE)
  calcOutput("PastureSuit",  subtype = paste("ISIMIP3bv2", "MRI-ESM2-0", "1850_2100", sep = ":"),
             file = paste0("f31_pastr_suitability_", ctype, ".mz"), years = magYears, aggregate = "cluster")
  calcOutput("PastureSuit",  subtype = paste("ISIMIP3bv2", "MRI-ESM2-0", "1850_2100", sep = ":"),
             file = "f31_pastr_suitability.mz", years = magYears, aggregate = FALSE)


  if (grepl("+PastrMngtLevels", dev)) {
    calcOutput("PastrMngtLevels", climatetype = paste0("MRI-ESM2-0", ":", climatescen),
               options = c("brazil_1", "brazil_2", "brazil_4"), cost_level = c(1, 2, 3),
               file = "PastrMngtLevels.mz", aggregate = FALSE)
  }

  calcOutput("ClimateClass", aggregate = "cluster", years = "y2015", file = paste0("koeppen_geiger_", ctype, ".mz"))
  calcOutput("ClimateClass", aggregate = "cluster", file = paste0("ipcc_climate_zones_", ctype, ".mz"))
  calcOutput("CellCountryFraction", aggregate = "cluster", file = paste0("cell_country_fraction_", ctype, ".mz"))

  # 32 forestry
  calcOutput("AfforestationMask", subtype = "noboreal",     aggregate = "cluster", round = 6,
    file = paste0("aff_noboreal_", ctype, ".mz"))
  calcOutput("AfforestationMask", subtype = "onlytropical", aggregate = "cluster", round = 6,
    file = paste0("aff_onlytropical_", ctype, ".mz"))
  calcOutput("AfforestationMask", subtype = "unrestricted", aggregate = "cluster", round = 6,
    file = paste0("aff_unrestricted_", ctype, ".mz"))

  calcOutput("NpiNdcAdAolcPol", aggregate = "cluster", round = 6, file = paste0("npi_ndc_ad_aolc_pol_", ctype, ".mz"))
  calcOutput("NpiNdcAffPol",    aggregate = "cluster", round = 6, file = paste0("npi_ndc_aff_pol_", ctype, ".mz"))

  calcOutput("BphEffect", aggregate = "cluster", file = paste0("f32_bph_effect_noTCRE_", ctype, ".mz"))
  calcOutput("BphTCRE",   aggregate = "cluster", file = paste0("f32_localTCRE_", ctype, ".mz"))
  calcOutput("BphMask",   aggregate = "cluster", file = paste0("f32_bph_mask_", ctype, ".mz"))

  # 34 urban land
  if (dev == "+GaoUrbanLand") {

    calcOutput("UrbanLandFuture", subtype = "Gao", aggregate = FALSE, round = 6, years = shortYears,
               file = "f34_urbanland_0.5.mz")
    calcOutput("UrbanLandFuture", subtype = "Gao", aggregate = "cluster", round = 6, years = shortYears,
               file = paste0("f34_urbanland_", ctype, ".mz"))
  } else {

  calcOutput("UrbanLandFuture", subtype = "LUH2v2", aggregate = FALSE, round = 6, years = shortYears,
             file = "f34_urbanland_0.5.mz")
  calcOutput("UrbanLandFuture", subtype = "LUH2v2", aggregate = "cluster", round = 6, years = shortYears,
    file = paste0("f34_urbanland_", ctype, ".mz"))
  }

  # 35 natveg
  calcOutput("AgeClassDistribution", aggregate = "cluster", round = 6, file = paste0("forestageclasses_", ctype, ".mz"))

  # 37 labour prod
  calcOutput("LabourProdImpactEmu", aggregate = "cluster", round = 6, subtype = "impact",
             file = paste0("f37_labourprodimpact_", ctype, ".mz"))
  calcOutput("LabourProdImpactEmu", aggregate = "cluster", round = 6, subtype = "relief",
             file = paste0("f37_labourprodrelief_", ctype, ".mz"))

  # 40
  calcOutput("TransportDistance", aggregate = "cluster", round = 6, file = paste0("transport_distance_", ctype, ".mz"))
  calcOutput("TransportDistance", aggregate = FALSE, round = 6, file = "transport_distance.mz")

  # 41 area equipped for irrigation
  calcOutput("AreaEquippedForIrrigation", aggregate = "cluster", cellular = TRUE, source = "Siebert",
    round = 6, file = paste0("avl_irrig_", ctype, ".mz"))
  calcOutput("AreaEquippedForIrrigation", aggregate = "cluster", cellular = TRUE, source = "LUH2v2",
    selectyears = magYearsPastLong, round = 6, file = paste0("avl_irrig_luh_t_", ctype, ".mz"))


  # 42 water demand
  calcOutput("Irrigation", lpjml = lpjml, years = lpjYears, climatetype = climatetype, aggregate = "cluster",
    round = 6, file = paste0("lpj_airrig_", ctype, ".mz"))

  # dummy Growing Period calc
  calcOutput("GrowingPeriod", lpjml = lpjml, years = lpjYears, climatetype = climatetype, yield_ratio = 0.1,
    aggregate = FALSE, round = 2, file = "lpj_grper_0.5.mz")

  # 43 water availability
  calcOutput("AvlWater", lpjml = lpjml, years = lpjYears, climatetype = climatetype, seasonality = "grper",
    aggregate = "cluster", round = 6, file = paste0("lpj_watavail_grper_", ctype, ".mz"))
  calcOutput("AvlWater", lpjml = lpjml, years = lpjYears, climatetype = climatetype, seasonality = "total",
    aggregate = "cluster", round = 6, file = paste0("lpj_watavail_total_", ctype, ".mz"))


  calcOutput("EnvmtlFlow", lpjml = lpjml, years = lpjYears, climatetype = climatetype, aggregate = "cluster",
    round = 6, seasonality = "grper", file = paste0("lpj_envflow_grper_", ctype, ".mz"))

  if (rev < 4.67) {
    calcOutput("WaterUseNonAg", datasource = "WATCH_ISIMIP_WATERGAP", years = lpjYears, seasonality = "grper",
               lpjml = lpjml, climatetype = climatetype, aggregate = "cluster",
               file = paste0("watdem_nonagr_grper_", ctype, ".mz"))
  } else {
    calcOutput("WaterUseNonAg", datasource = "WATERGAP_ISIMIP", usetype = "all:withdrawal",
               selectyears = lpjYears, seasonality = "grper", lpjml = lpjml, climatetype = climatetype,
               aggregate = "cluster", file = paste0("watdem_nonagr_grper_", ctype, ".mz"))
  }

  # 44 biodiversity
  calcOutput("BiomeType", aggregate = "cluster", cells = "magpiecell", round = 6,
             file = paste0("biorealm_biome_", ctype, ".mz"))
  calcOutput("Luh2SideLayers", aggregate = "cluster", round = 6, file = paste0("luh2_side_layers_", ctype, ".mz"))
  calcOutput("Luh2SideLayers", aggregate = FALSE, round = 6, file = "luh2_side_layers_0.5.mz")
  calcOutput("RRLayer", aggregate = "cluster", round = 6, file = paste0("rr_layer_", ctype, ".mz"))

  # 50 nitrogen
  calcOutput("AtmosphericDepositionRates", cellular = TRUE, aggregate = FALSE, round = 6,
    file = "f50_AtmosphericDepositionRates_0.5.mz")
  calcOutput("NitrogenFixationRateNatural",               aggregate = FALSE, round = 6,
    file = "f50_NitrogenFixationRateNatural_0.5.mz")

  calcOutput("AtmosphericDepositionRates", cellular = TRUE, aggregate = "cluster", round = 6,
    file = paste0("f50_AtmosphericDepositionRates_", ctype, ".mz"))
  calcOutput("NitrogenFixationRateNatural",               aggregate = "cluster", round = 6,
    file = paste0("f50_NitrogenFixationRateNatural_", ctype, ".mz"))

  calcOutput("Carbon", aggregate = FALSE, lpjml = lpjml, climatetype = climatetype,
    round = 6, years = "y1995", file = "lpj_carbon_stocks_0.5.mz")
  calcOutput("TopsoilCarbon", aggregate = FALSE, lpjml = lpjml, climatetype = climatetype,
    round = 6, years = "y1995", file = "lpj_carbon_topsoil_0.5.mz")

  calcOutput("Carbon", aggregate = "cluster", lpjml = lpjml, climatetype = climatetype,
    round = 6, years = lpjYears, file = paste0("lpj_carbon_stocks_", ctype, ".mz"))
  calcOutput("TopsoilCarbon", aggregate = "cluster", lpjml = lpjml, climatetype = climatetype,
    round = 6, years = lpjYears, file = paste0("lpj_carbon_topsoil_", ctype, ".mz"))

  # 58 peatland
  calcOutput("Peatland", subtype = "degraded", aggregate = FALSE, round = 6, file = "f58_peatland_degrad_0.5.mz")
  calcOutput("Peatland", subtype = "intact",   aggregate = FALSE, round = 6, file = "f58_peatland_intact_0.5.mz")
  calcOutput("Peatland", subtype = "degraded", aggregate = "cluster", round = 6,
    file = paste0("f58_peatland_degrad_", ctype, ".mz"))
  calcOutput("Peatland", subtype = "intact",   aggregate = "cluster", round = 6,
    file = paste0("f58_peatland_intact_", ctype, ".mz"))


  # 59 som
  calcOutput("SOMinitialsiationPools", aggregate = "cluster", round = 6,
    file = paste0("f59_som_initialisation_pools_", ctype, ".mz"))
  calcOutput("SOCLossShare",           aggregate = "cluster", rate = "loss", round = 6, years = "y1995",
    file = paste0("cshare_released_", ctype, ".mz"))

  ##### AGGREGATION ######

  # create info file
  writeInfo <- function(file, lpjmlData, resHigh, resOut, rev, cluster) {
    functioncall <- paste(deparse(sys.call(-3)), collapse = "")

    map <- toolGetMapping(type = "regional", name = getConfig("regionmapping"))
    regionscode <- regionscode(map)

    info <- c("lpj2magpie settings:",
      paste("* LPJmL data:", lpjmlData),
      paste("* Revision:", rev),
      "", "aggregation settings:",
      paste("* Input resolution:", resHigh),
      paste("* Output resolution:", resOut),
      paste("* Regionscode:", regionscode),
      "* Number of clusters per region:",
      paste(format(names(cluster), width = 5, justify = "right"), collapse = ""),
      paste(format(cluster, width = 5, justify = "right"), collapse = ""),
      paste("* Call:", functioncall)
      )
    base::cat(info, file = file, sep = "\n")
  }
  nrClusterPerRegion <- substr(attributes(p$data)$legend_text, 6,
                                  nchar(attributes(p$data)$legend_text) - 1)
  writeInfo(file = "info.txt", lpjmlData = climatetype,
    resHigh = "0.5", resOut = ctype, rev = rev, cluster = nrClusterPerRegion)

  return(list(tag = versionTag,
              pucTag = sub("^[^_]*_", "", versionTag)))

}
