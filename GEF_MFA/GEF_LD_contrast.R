source("/home/aiddata/Desktop/Github/geoML/geoML.R")

full.dta <- read.csv("/home/aiddata/Desktop/Github/geoML/GEF_MFA/GEF_LD_contrast.csv")

#==================================================================================
#==================================================================================
#Summarizing Financials and timestamps (this should be done in the data prep scripts
#in the future).
total <- as.numeric(gsub(",","",as.character(full.dta$Cofinance.CEO.endorse.stage))) +  
  as.numeric(gsub(",","",as.character(full.dta$GEF.Project.Grant.CEO.endorse.stage)))
cofinancing <- as.numeric(gsub(",","",as.character(full.dta$Cofinance.CEO.endorse.stage))) / total
full.dta$cofinance.ratio <- cofinancing 
full.dta$total.funding <- total
full.dta$GEF.funding <- as.numeric(gsub(",","",as.character(full.dta$GEF.Project.Grant.CEO.endorse.stage)))
full.dta$cofinance.funding <- as.numeric(gsub(",","",as.character(full.dta$Cofinance.CEO.endorse.stage)))

full.dta <- full.dta[!is.na(full.dta$Actual.date.of.implementation.start),]
full.dta <- full.dta[!(as.character(full.dta$Actual.date.of.implementation.start) == ""),]
full.dta$start.date <- as.Date(full.dta$Actual.date.of.implementation.start, format="%d-%b-%y")
full.dta$year <- as.numeric( format( full.dta$start.date, '%Y'))
full.dta$wdpa_5km.na.sum <- full.dta$wdpa_5km.na.sum / 4

#==================================================================================
#==================================================================================
#Define control variables
Vars <-  c("dist_to_all_rivers.na.mean", "dist_to_roads.na.mean",
           "srtm_elevation_500m.na.mean", "srtm_slope_500m.na.mean",
           "accessibility_map.na.mean", "gpw_v3_density.2000.mean",
           "wdpa_5km.na.sum", "treecover2000.na.mean", "latitude",
           "longitude", "udel_precip_v4_01_yearly_max.2002.mean", 
           "udel_precip_v4_01_yearly_min.2002.mean", 
           "udel_precip_v4_01_yearly_mean.2002.mean",
           "udel_air_temp_v4_01_yearly_max.2002.mean",   
           "udel_air_temp_v4_01_yearly_min.2002.mean",   
           "udel_air_temp_v4_01_yearly_mean.2002.mean",
           "v4composites_calibrated.2002.mean",
           "ltdr_yearly_ndvi_mean.2002.mean", 
           "Region",  "GEF.replenishment.phase",
           "cofinance.ratio", "total.funding", "GEF.funding")

VarNames <- c("Dist. to Rivers (m)", "Dist. to Roads (m)",
              "Elevation (m)", "Slope (degrees)",
              "Urb. Dist. (rel)", "Pop. Density (2000)",
              "Protected Area %", "Treecover (2000, %)", "Latitude",
              "Longitude", "Max Precip. (2002, mm)", 
              "Min Precip (2002, mm)", 
              "Mean Precip (2002, mm)",
              "Max Temp (2002, C)", 
              "Min Temp (2002, C)", 
              "Mean Temp (2002, C)",
              "Nightime Lights (2002, Relative)",
              "NDVI (2002, Unitless)", "Region", "Replenishment Phase",
              "Cofinance Ratio (% Total)", "Total Funding", "GEF funding"
)

out_path = "/home/aiddata/Desktop/Github/geoML/GEF_MFA/LD_MFA_Contrast/NDVI/"

t <- geoML(dta=full.dta, 
         trt=c("treatment", "GEF MFA Projects"), 
         ctrl=c(Vars, VarNames), 
         outcome=c("ltdr_yearly_ndvi_mean.2013.mean", "2013 NDVI"), 
         pth=out_path, 
         file.prefix="NDVI_max", 
         kvar=c("GEF.funding","cofinance.ratio",
                "treecover2000.na.mean","total.funding"),
         geog.fields = c("latitude", "longitude"),
         caliper=1.0,
         counterfactual.name = "Land Degradation",
         tree.ctrl = c(20,10),
         col.invert = FALSE,
         tree.cnt = 10000
)
