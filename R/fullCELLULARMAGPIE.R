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

fullCELLULARMAGPIE <- function(rev=0.1, ctype="c200") {

    mag_years_past_short <- c("y1995","y2000","y2005","y2010")
    mag_years_past_long  <- c("y1995","y2000","y2005","y2010","y2015")
    mag_years <- findset("time")
    short_years <- findset("t_all")

    ### test settings (will be loaded from config in fina version)
    climatetype="HadGEM2_ES:rcp2p6:co2"
    harmonize_baseline="CRU_4"
    ref_year="y2015"

    map <- calcOutput("Cluster", ctype=ctype, weight=NULL, aggregate=FALSE)
    toolStoreMapping(map,"clustermapping.csv",type="regional",error.existing = FALSE)
    setConfig(regionmapping = "clustermapping.csv")

    # 09 drivers
    ### gridded pop?

    # 14 yields
    # calcOutput("Yields")
    # These outputs need to be aggregated using weighted area mean
    calcOutput("GCMClimate", aggregate="cluster", file = "rcp85.HadGEM2.temperature_c200.mz", GCMModel = "HadGEM2", ClimateVariable = "temperature", rcp = "rcp85")
    calcOutput("GCMClimate", aggregate="cluster", file = "rcp85.HadGEM2.precipitation_c200.mz", GCMModel = "HadGEM2", ClimateVariable = "precipitation", rcp = "rcp85")
    calcOutput("GCMClimate", aggregate="cluster", file = "rcp85.HadGEM2.longwave_radiation_c200.mz", GCMModel = "HadGEM2", ClimateVariable = "longwave_radiation", rcp = "rcp85")
    calcOutput("GCMClimate", aggregate="cluster", file = "rcp85.HadGEM2.shortwave_radiation_c200.mz", GCMModel = "HadGEM2", ClimateVariable = "shortwave_radiation", rcp = "rcp85")
    calcOutput("GCMClimate", aggregate="cluster", file = "rcp85.HadGEM2.wetdays_c200.mz", GCMModel = "HadGEM2", ClimateVariable = "wetdays", rcp = "rcp85")
    calcOutput("CO2Atmosphere", aggregate="cluster", file="calcCO2Atmosphere_c200.mz", rcp="rcp85", level="cellular")
    calcOutput("SoilCharacteristics", aggregate="cluster", file="SoilCharacteristics_c200.mz")

    #10 land
    calcOutput("LanduseInitialisation", aggregate=FALSE, cellular=TRUE, land="fao", input_magpie=TRUE, years=mag_years_past_short, round=6, file="avl_land_t_0.5.mz")
    calcOutput("LanduseInitialisation", aggregate=TRUE, cellular=TRUE, land="fao", input_magpie=TRUE, years=mag_years_past_short, round=6, file="avl_land_t_c200.mz")
    calcOutput("SeaLevelRise", aggregate=FALSE, round=6, file="f10_SeaLevelRise_0.5.mz")
    calcOutput("AvlLandSi", aggregate=FALSE, round=6, file="avl_land_si_0.5.mz")



    #30 crop
    calcOutput("Croparea", sectoral="kcr", physical=TRUE, cellular=TRUE, irrigation=FALSE, aggregate = FALSE,file="f30_croparea_initialisation_0.5.mz")
    calcOutput("Croparea", sectoral="kcr", physical=TRUE, cellular=TRUE, irrigation=TRUE, aggregate = FALSE,file="f30_croparea_w_initialisation_0.5.mz")
    calcOutput("Croparea", sectoral="kcr", physical=TRUE, cellular=TRUE, irrigation=FALSE, aggregate = "cluster", file="f30_croparea_initialisation_c200.mz")
    calcOutput("Croparea", sectoral="kcr", physical=TRUE, cellular=TRUE, irrigation=TRUE, aggregate = "cluster",file="f30_croparea_w_initialisation_c200.mz")

    #32 forestry
    calcOutput("AfforestationMask",subtype="noboreal",aggregate=FALSE,round=6, file="aff_noboreal_0.5.mz")
    calcOutput("AfforestationMask",subtype="onlytropical",aggregate=FALSE,round=6, file="aff_onlytropical_0.5.mz")
    calcOutput("AfforestationMask",subtype="unrestricted",aggregate=FALSE,round=6, file="aff_unrestricted_0.5.mz")

    calcOutput("NpiNdcAdAolcPol", aggregate=FALSE, round=6, file="npi_ndc_ad_aolc_pol_0.5.mz")
    calcOutput("NpiNdcAffPol", aggregate=FALSE, round=6, file="npi_ndc_aff_pol_0.5.mz")

    #35 natveg
    calcOutput("AgeClassDistribution", aggregate="cluster", round=6, file="secdf_age_class_dist_c200.mz")
    calcOutput("ProtectArea", aggregate="cluster", round=6, file="protect_area_c200.mz" )

    #34
    calcOutput("UrbanLandFuture", aggregate=FALSE, round=6, file="f34_UrbanLand_0.5.mz")

    #40
    calcOutput("TransportDistance", aggregate=FALSE, round=6, file="transport_distance_0.5.mz")

    #41 area equipped for irrigation
    calcOutput("AreaEquippedForIrrigation", aggregate="cluster", cellular=TRUE, source="Siebert", round=6, file="avl_irrig_c200.mz")
    calcOutput("AreaEquippedForIrrigation", aggregate="cluster", cellular=TRUE, source="LUH2v2",  years=mag_years_past_short, round=6, file="avl_irrig_luh_t_c200.mz")

    #42 water demand
    calcOutput("Irrigation", version="LPJmL5", climatetype=climatetype, harmonize_baseline=harmonize_baseline, ref_year=ref_year, time="spline", dof=4, aggregate="cluster", round=6, file="lpj_airrig_c200.mz")
    calcOutput("EnvmtlFlow", version="LPJmL4", climatetype=climatetype, harmonize_baseline=harmonize_baseline, ref_year=ref_year, time="spline", dof=4, aggregate="cluster", round=6, seasonality="grper", file="lpj_envflow_grper_c200.mz")
    calcOutput("NonAgWaterDemand", source="WATCH_ISIMIP_WATERGAP", seasonality="grper", aggregate=FALSE, file="watdem_nonagr_grper_0.5.mz")
    calcOutput("NonAgWaterDemand", source="WATERGAP2020", seasonality="grper", waterusetype="withdrawal", aggregate="cluster", file="watdem_nonagr_ww_grper_c200.mz")
    calcOutput("NonAgWaterDemand", source="WATERGAP2020", seasonality="grper", waterusetype="consumption", aggregate="cluster", file="watdem_nonagr_wc_grper_c200.mz")

    #43 water availability
    calcOutput("AvlWater", version="LPJmL4", climatetype=climatetype, harmonize_baseline=harmonize_baseline, ref_year=ref_year, time="spline", dof=4, seasonality="grper", aggregate="cluster", round=6, file="lpj_watavail_grper_c200.mz")

    #44 biodiversity
    calcOutput("Luh2SideLayers", aggregate=FALSE, round=6, file="luh2_side_layers_0.5.mz")
    calcOutput("RRLayer", aggregate=FALSE, round=6, file="rr_layer_0.5.mz")

    #50 nitrogen
    calcOutput("AtmosphericDepositionRates", cellular=TRUE, aggregate=FALSE, round=6, file="f50_AtmosphericDepositionRates_0.5.mz")
    calcOutput("NitrogenFixationRateNatural", aggregate=FALSE, round=6, file="f50_NitrogenFixationRateNatural_0.5.mz")

    #52 carbon
    calcOutput("Carbon", aggregate = FALSE, round=6, file="lpj_carbon_stocks_0.5.mz")
    calcOutput("TopsoilCarbon", aggregate = FALSE, round=6, file="lpj_carbon_topsoil_0.5.mz")


    #59 som
    calcOutput("SOMinitialsiationPools", aggregate = FALSE, round=6, file="f59_som_initialisation_pools_0.5.mz")


    ## OTHER ##
    calcOutput("CalibratedArea", aggregate=FALSE, round=6, file="calibrated_area_0.5.mz" )
    calcOutput("CshareReleased", aggregate=FALSE, round=6, file="cshare_released_0.5.mz")
    calcOutput("Koeppen_geiger", aggregate=FALSE,round=6,file="koeppen_geiger_0.5.mz")


    ##### AGGREGATION ######

    # create info file
    writeInfo <- function(file,lpjml_data, res_high, res_out, rev) {
      functioncall <- paste(deparse(sys.call(-1)), collapse = "")

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
      cat(info,file=file,sep='\n')
    }
    writeInfo(file='info.txt', lpjml_data="default", res_high="0.5", res_out="IDK", rev=rev)


}
