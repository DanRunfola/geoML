source("/home/aiddata/Desktop/Github/geoML/geoML.R")

full.dta <- read.csv("/home/aiddata/Desktop/Github/geoML/data/GEF_LD_contrast.csv")

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
           "ltdr_yearly_ndvi_mean.2002.mean")

VarNames <- c("Distance to Rivers (m)", "Distance to Roads (m)",
              "Elevation(m)", "Slope (°)",
              "Urban Travel Distance (relative)", "Population Density (2000)",
              "Protected Area %", "Treecover (2000, %)", "Latitude",
              "Longitude", "Maximum Precipitation (2002, mm)", 
              "Minimum Precipitation (2002, mm)", 
              "Mean Precipitation (2002, mm)",
              "Maximum Temperature (2002, C)", 
              "Minimum Temperature (2002, C)", 
              "Mean Temperature (2002, C)",
              "Nightime Lights (2002, Relative)",
              "NDVI (2002, Unitless)"                       
)

out_path = "/home/aiddata/Desktop/Github/geoML/output/"

source("/home/aiddata/Desktop/Github/geoML/geoML.R")
t <- geoML(dta=full.dta, 
         trt=c("treatment", "GEF MFA Projects"), 
         ctrl=c(Vars, VarNames), 
         outcome=c("ltdr_yearly_ndvi_mean.2013.mean", "2013 NDVI"), 
         pth=out_path, 
         file.prefix="NDVI_max", 
         kvar=c("v4composites_calibrated.2002.mean","dist_to_roads.na.mean",
                "accessibility_map.na.mean","srtm_slope_500m.na.mean"),
         geog.fields = c("latitude", "longitude"),
         caliper=0.5,
         counterfactual.name = "Land Degradation"
)
