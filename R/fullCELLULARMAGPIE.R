#' fullCELLULARMAGPIE
#'
#' Function that produces the complete cellular data set required for running the
#' MAgPIE model.
#'
#' @param rev data revision which should be used as input (positive numeric).
#' @param ctype aggregation clustering type, which is a combination of a single letter, indicating the cluster methodology, and a number,
#' indicating the number of resulting clusters. Available methodologies are hierarchical clustering (h), normalized k-means clustering
#' (n) and combined hierarchical/normalized k-means clustering (c). In the latter hierarchical clustering is used to determine the
#' cluster distribution among regions whereas normalized k-means is used for the clustering within a region.
#' @param dev development suffix to distinguish development versions for the same data revision. This can be useful to distinguish
#' parallel lines of development.
#' @param climatetype climate change scenario to be used
#' @param clusterweight Should specific regions be resolved with more or less detail? Values > 1 mean higher share, < 1 lower share
#' e.g. cfg$clusterweight <- c(LAM=2) means that a higher level of detail for region LAM if set to NULL all weights will be assumed to be 1.
#' examples:
#' c(LAM=1.5,SSA=1.5,OAS=1.5)
#' c(LAM=2,SSA=2,OAS=2)
#' \code{\link{setConfig}} (e.g. for setting the mainfolder if not already set
#' properly).
#' @author Kristine Karstens, Jan Philipp Dietrich
#' @seealso
#' \code{\link{readSource}},\code{\link{getCalculations}},\code{\link{calcOutput}},\code{\link{setConfig}}
#' @examples
#'
#' \dontrun{
#' retrieveData("CELLULARMAGPIE", revision=12, mainfolder="pathtowhereallfilesarestored")
#' }
#' @importFrom madrat setConfig getConfig
#' @importFrom magpiesets findset

fullCELLULARMAGPIE <- function(rev=0.1, dev="", ctype="c200", climatetype="HadGEM2_ES:rcp2p6:co2", clusterweight=NULL) {

  sizelimit <- getOption("magclass_sizeLimit")
  options(magclass_sizeLimit=1e+10)
  on.exit(options(magclass_sizeLimit=sizelimit))

  setConfig(debug=TRUE)

  mag_years_past_short <- c("y1995","y2000","y2005","y2010")
  mag_years_past_long  <- c("y1995","y2000","y2005","y2010","y2015")
  mag_years <- findset("time")
  short_years <- findset("t_all")
  lpj_years <- seq(1995, 2100,by=5)

  ### test settings (will be loaded from config in fina version)
  climatetype="HadGEM2_ES:rcp2p6:co2"
  harmonize_baseline="CRU_4"
  ref_year="y2015"

  map <- calcOutput("Cluster", ctype=ctype, weight=clusterweight, aggregate=FALSE)
  weightID <- ifelse(is.null(clusterweight),"",paste0("_",names(clusterweight),clusterweight,collapse=""))
  clustermapname <- sub("\\.[^.]*$",".rds",paste0("clustermap_rev",rev,dev,"_",ctype,weightID,"_",getConfig("regionmapping")))
  toolStoreMapping(map,clustermapname,type="regional",where=c("mappingfolder","outputfolder"),error.existing = FALSE)
  setConfig(extramappings = clustermapname)

  # 09 drivers
  ### gridded pop?

  # 14 yields
  calcOutput("Yields", version="LPJmL5", climatetype=climatetype, time="spline", dof=4,
             harmonize_baseline=harmonize_baseline, ref_year=ref_year, aggregate = FALSE,
             years="y1995", file=paste0("lpj_yields_0.5.mz"))
  calcOutput("Yields", version="LPJmL5", climatetype=climatetype, time="spline", dof=4,
             harmonize_baseline=harmonize_baseline, ref_year=ref_year, aggregate = "cluster",
             years=lpj_years, file=paste0("lpj_yields_", ctype, ".mz"))


  # These outputs need to be aggregated using weighted area mean
  calcOutput("GCMClimate", aggregate="cluster", file=paste0("rcp85.HadGEM2.temperature_", ctype, ".mz"), GCMModel = "HadGEM2", ClimateVariable = "temperature", rcp = "rcp85")
  calcOutput("GCMClimate", aggregate="cluster", file=paste0("rcp85.HadGEM2.precipitation_", ctype, ".mz"), GCMModel = "HadGEM2", ClimateVariable = "precipitation", rcp = "rcp85")
  calcOutput("GCMClimate", aggregate="cluster", file=paste0("rcp85.HadGEM2.longwave_radiation_", ctype, ".mz"), GCMModel = "HadGEM2", ClimateVariable = "longwave_radiation", rcp = "rcp85")
  calcOutput("GCMClimate", aggregate="cluster", file=paste0("rcp85.HadGEM2.shortwave_radiation_", ctype,  ".mz"), GCMModel = "HadGEM2", ClimateVariable = "shortwave_radiation", rcp = "rcp85")
  calcOutput("GCMClimate", aggregate="cluster", file=paste0("rcp85.HadGEM2.wetdays_", ctype, ".mz"), GCMModel = "HadGEM2", ClimateVariable = "wetdays", rcp = "rcp85")
  calcOutput("CO2Atmosphere", aggregate="cluster", file=paste0("calcCO2Atmosphere_", ctype, ".mz"), rcp="rcp85", level="cellular")
  calcOutput("SoilCharacteristics", aggregate="cluster", file=paste0("SoilCharacteristics_", ctype, ".mz"))
  calcOutput("ClimateClass", aggregate="cluster", years="y2015", file=paste0("koeppen_geiger_", ctype, ".mz"))

  #10 land
  calcOutput("LanduseInitialisation", aggregate=FALSE, cellular=TRUE, land="fao", input_magpie=TRUE, selectyears=mag_years_past_long, round=6, file="avl_land_t_0.5.mz")
  calcOutput("LanduseInitialisation", aggregate="cluster", cellular=TRUE, land="fao", input_magpie=TRUE, selectyears=mag_years_past_long, round=6, file=paste0("avl_land_t_", ctype, ".mz"))
  calcOutput("SeaLevelRise", aggregate=FALSE, cellular=TRUE, years=mag_years, round=6, file="f10_SeaLevelRise_0.5.mz")
  calcOutput("AvlLandSi", aggregate=FALSE, round=6, file="avl_land_si_0.5.mz")
  calcOutput("AvlLandSi", aggregate="cluster", round=6, file=paste0("avl_land_si_", ctype, ".mz"))


  #30 crop
  #calcOutput("Croparea", sectoral="kcr", physical=TRUE, cellular=TRUE, irrigation=FALSE, aggregate = FALSE,file="f30_croparea_initialisation_0.5.mz")
  #calcOutput("Croparea", sectoral="kcr", physical=TRUE, cellular=TRUE, irrigation=TRUE, aggregate = FALSE,file="f30_croparea_w_initialisation_0.5.mz")
  calcOutput("Croparea", sectoral="kcr", physical=TRUE, cellular=TRUE, irrigation=FALSE, aggregate = "cluster", file=paste0("f30_croparea_initialisation_", ctype, ".mz"))
  calcOutput("Croparea", sectoral="kcr", physical=TRUE, cellular=TRUE, irrigation=TRUE, aggregate = "cluster", file=paste0("f30_croparea_w_initialisation_", ctype, ".mz"))

  #32 forestry
  calcOutput("AfforestationMask", subtype="noboreal",     aggregate="cluster", round=6, file=paste0("aff_noboreal_", ctype, ".mz"))
  calcOutput("AfforestationMask", subtype="onlytropical", aggregate="cluster", round=6, file=paste0("aff_onlytropical_", ctype, ".mz"))
  calcOutput("AfforestationMask", subtype="unrestricted", aggregate="cluster", round=6, file=paste0("aff_unrestricted_", ctype, ".mz"))

  calcOutput("NpiNdcAdAolcPol", aggregate="cluster", round=6, file=paste0("npi_ndc_ad_aolc_pol_", ctype, ".mz"))
  calcOutput("NpiNdcAffPol",    aggregate="cluster", round=6, file=paste0("npi_ndc_aff_pol_", ctype, ".mz"))

  #35 natveg
  calcOutput("AgeClassDistribution", aggregate="cluster", round=6, file=paste0("secdf_age_class_dist_", ctype, ".mz"))
  calcOutput("ProtectArea",          aggregate="cluster", round=6, file=paste0("protect_area_", ctype, ".mz") )

  #34
  calcOutput("UrbanLandFuture", aggregate=FALSE, round=6, years=short_years, file="f34_UrbanLand_0.5.mz")
  calcOutput("UrbanLandFuture", aggregate="cluster", round=6, years=short_years, file=paste0("f34_UrbanLand_", ctype, ".mz"))

  #40
  calcOutput("TransportDistance", aggregate="cluster", round=6, file=paste0("transport_distance_", ctype, ".mz"))

  #41 area equipped for irrigation
  calcOutput("AreaEquippedForIrrigation", aggregate="cluster", cellular=TRUE, source="Siebert", round=6, file=paste0("avl_irrig_", ctype, ".mz"))
  calcOutput("AreaEquippedForIrrigation", aggregate="cluster", cellular=TRUE, source="LUH2v2",  selectyears=mag_years_past_long, round=6, file=paste0("avl_irrig_luh_t_", ctype, ".mz"))

  #42 water demand
  calcOutput("Irrigation", version="LPJmL5", years=lpj_years, climatetype=climatetype, harmonize_baseline=harmonize_baseline, ref_year=ref_year, time="spline", dof=4, aggregate="cluster", round=6, file=paste0("lpj_airrig_", ctype, ".mz"))
  calcOutput("EnvmtlFlow", version="LPJmL4", years=lpj_years, climatetype=climatetype, harmonize_baseline=harmonize_baseline, ref_year=ref_year, time="spline", dof=4, aggregate="cluster", round=6, seasonality="grper", file=paste0("lpj_envflow_grper_", ctype, ".mz"))
  calcOutput("NonAgWaterDemand", source="WATCH_ISIMIP_WATERGAP", years=lpj_years, seasonality="grper", aggregate="cluster", file="watdem_nonagr_grper_c200.mz")
  calcOutput("NonAgWaterDemand", source="WATERGAP2020", years=lpj_years, seasonality="grper", waterusetype="withdrawal", aggregate="cluster", file=paste0("watdem_nonagr_ww_grper_", ctype, ".mz"))
  calcOutput("NonAgWaterDemand", source="WATERGAP2020", years=lpj_years, seasonality="grper", waterusetype="consumption", aggregate="cluster", file=paste0("watdem_nonagr_wc_grper_", ctype, ".mz"))

  #43 water availability
  calcOutput("AvlWater", version="LPJmL4", years=lpj_years, climatetype=climatetype, harmonize_baseline=harmonize_baseline, ref_year=ref_year, time="spline", dof=4, seasonality="grper", aggregate="cluster", round=6, file=paste0("lpj_watavail_grper_", ctype, ".mz"))

  #44 biodiversity
  calcOutput("Luh2SideLayers", aggregate="cluster", round=6, file=paste0("luh2_side_layers_", ctype, ".mz"))
  calcOutput("RRLayer",        aggregate="cluster", round=6, file=paste0("rr_layer_", ctype, ".mz"))

  #50 nitrogen
  calcOutput("AtmosphericDepositionRates", cellular=TRUE, aggregate=FALSE, round=6, file="f50_AtmosphericDepositionRates_0.5.mz")
  calcOutput("NitrogenFixationRateNatural",               aggregate=FALSE, round=6, file="f50_NitrogenFixationRateNatural_0.5.mz")

  calcOutput("AtmosphericDepositionRates", cellular=TRUE, aggregate="cluster", round=6, file=paste0("f50_AtmosphericDepositionRates_", ctype, ".mz"))
  calcOutput("NitrogenFixationRateNatural",               aggregate="cluster", round=6, file=paste0("f50_NitrogenFixationRateNatural_", ctype, ".mz"))



  #52 carbon
  calcOutput("Carbon", aggregate = FALSE, version="LPJmL4+5", climatetype=climatetype,
             harmonize_baseline=harmonize_baseline, ref_year=ref_year,
             time="spline", dof=4, round=6, years="y1995", file="lpj_carbon_stocks_0.5.mz")
  calcOutput("TopsoilCarbon", aggregate = FALSE, version="LPJmL4", climatetype=climatetype,
             harmonize_baseline=harmonize_baseline, ref_year=ref_year,
             time="spline", dof=4, round=6, years="y1995", file="lpj_carbon_topsoil_0.5.mz")

  calcOutput("Carbon", aggregate = "cluster", version="LPJmL4+5", climatetype=climatetype,
             harmonize_baseline=harmonize_baseline, ref_year=ref_year,
             time="spline", dof=4, round=6, years=lpj_years, file=paste0("lpj_carbon_stocks_", ctype, ".mz"))
  calcOutput("TopsoilCarbon", aggregate = "cluster", version="LPJmL4", climatetype=climatetype,
             harmonize_baseline=harmonize_baseline, ref_year=ref_year,
             time="spline", dof=4, round=6, years=lpj_years, file=paste0("lpj_carbon_topsoil_", ctype, ".mz"))


  #58 peatland
  calcOutput("Peatland", subtype="degraded", aggregate = FALSE, round=6, file="f58_peatland_degrad_0.5.mz")
  calcOutput("Peatland", subtype="intact",   aggregate = FALSE, round=6, file="f58_peatland_intact_0.5.mz")
  calcOutput("Peatland", subtype="degraded", aggregate = "cluster", round=6, file=paste0("f58_peatland_degrad_", ctype, ".mz"))
  calcOutput("Peatland", subtype="intact",   aggregate = "cluster", round=6, file=paste0("f58_peatland_intact_", ctype, ".mz"))


  #59 som
  calcOutput("SOMinitialsiationPools", aggregate="cluster", round=6, file=paste0("f59_som_initialisation_pools_", ctype, ".mz"))
  calcOutput("SOCLossShare",           aggregate="cluster", rate="loss", round=6, years="y1995", file=paste0("cshare_released_", ctype, ".mz"))

  ##### AGGREGATION ######

  # create info file
  writeInfo <- function(file,lpjml_data, res_high, res_out, rev) {
    functioncall <- paste(deparse(sys.call(-2)), collapse = "")

    map <- toolMappingFile("regional", getConfig("regionmapping"),
                           readcsv = TRUE)
    regionscode <- regionscode(map)

    info <- c('lpj2magpie settings:',
              paste('* LPJmL data:',lpjml_data),
              paste('* Revision:', rev),
              '','aggregation settings:',
              paste('* Input resolution:',res_high),
              paste('* Output resolution:',res_out),
              paste('* Regionscode:',regionscode),
              paste('* Call:', functioncall))
    base::cat(info,file=file,sep='\n')
  }
  writeInfo(file=paste0(getConfig("outputfolder"),'/info.txt'), lpjml_data=climatetype, res_high="0.5", res_out=ctype, rev=rev)

}
