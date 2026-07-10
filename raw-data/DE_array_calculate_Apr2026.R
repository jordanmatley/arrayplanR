###Evaluating effectiveness of array with pre-designated receiver locations###

#Assumes a user has some detection range data already conducted - see input files for comparison.

#Default approach uses 'receiver-centric' detection ranges - each receiver has its own detection range profile and cells around it have same profile, as opposed to 'cell-centric' detection ranges - each cell has its own detection range profile independent of the receivers around it. The main reason is because many studies will not have representative cell-specific detection range profiles throughout study area. Note, 'cell-centric' detection range capacity is still possible but needs to be manually adjusted (see optional code chunks).

#In example code, 'Habitat' is main factor distinguishing different detection range profiles across receivers. Other variables can easily be built in.

#Evaluation of array efficiency centered around prioritizing the 'effective area covered'  or 'effective sum of probabilities' above a defined probability threshold.

#The other main decision is whether the user wants cell probabilities based on being detected by at least 1 receiver or by more - for example an overlap of 3+ receivers is useful for fine-scale positioning objectives

#To incorporate temporal variability, it is suggested that separate evaluations are used.

#If processing is slow, it is likely because the selected 'cellsize' is very small relative to the study area. Consider increasing the 'cellsize' at least initially.

#It is suggested the user runs example code prior to using their own to understand the process.


library(dplyr)
library(ggplot2)
library(sf)
library(tidyr)
library(progressr)
library(data.table)
library(digest)
library(proxy)
library(purrr)
library(terra)
library(exactextractr)
library(ggspatial)

theme_set(theme_bw())

setwd("C:\\Users\\jorda\\OneDrive\\Research\\Flinders\\Research\\Other\\Array_DE\\R\\for_github")

#######SET NECESSARY PARAMETERS#########
cellsize=50 #Options: Numeric (in metres) (Default = 100). e.g., 100 = 100 x 100 m. Size of grids/cells used to calculate distances between receivers and cells (a smaller cell better for resolution but can slow computation)
max_distance=500 #Options: Numeric (in metres) (Default = 500). Distance from a receiver which you would expect ~0 (or low) detectability. Should be chosen based on range tests. Used to save computing time by ignoring areas not in range.
cell_prob_cutoff=0.1 #Options: Proportion (Default = 0.1). Equivalent to a x% chance of a transmission from a grid being heard (by n receivers). Should be chosen based on how stringent a user wants detection probabilities to be (e.g., only interested in grids with low vs high probability of detection).
n.receiver_overlap=1 #Options: Numeric (Default = 1). Set how many receivers to detect transmission from a grid (e.g., broad-scale tracking = 1, fine-scale positioning = 3+).
distance_approach="horizontal" #Options: "horizontal","relative_static","relative_dynamic" (Default = "horizontal"). Enables user to account for distances between an animal and a receiver that may occur at different depths in relation to each other. "horizontal" means no depth differences are accounted for (animal-receiver assumed to be aligned horizontally). "relative_static" allows for a singular relative depth difference (positive value) between animal and receiver (e.g., if study area and expected animal depth use are consistent). No bathymetry data are needed for these first two options. However, "relative_static" requires a numeric value (i.e., absolute relative depth difference)(see rel_depth_diff option). "relative_dynamic" incorporates bottom depth and the relative position an animal is expected to be in - "surface", "midwater", "bottom" (see water_column_position option), providing 3-dimensional estimates of detection probabilities.
rel_depth_diff<-NULL #Options: Numeric (Default = NULL). Only used when distance_approach="relative_static". The relative depth difference (static - i.e., only one value) in metres between animals and receiver locations. This option allows for non-horizontal differences in relatively homogeneous environments.
water_column_position=NULL #Options: "surface", "midwater", "bottom" (Default = NULL). Only used when distance_approach="relative_dynamic". The expected position of an animal in the water column to calculate non-horizontal distances between receiver and animal (where relevant).
#######SET NECESSARY PARAMETERS#########


#######LOAD AND PREPARE NECESSARY FILES#########
#1) study area file (gis file (e.g., .shp) or can be added manually from Google Earth)
study_area = st_read("BB_shp\\eaarlb_STT_250m_WaterMask_bbaypbay.shp") #using an existing shp file for study area (with land barriers)
#study_area2 = st_read("Study_area_cut.kml") #using a kml created directly in Google Earth using polygon function (this kml consisted of a smaller area of 'study_area')
#st_write(study_area, "Brewers_JM.geojson") #example of how to export a GeoJSON file to web application if desired (transform to correct crs first).

study_area<-st_transform(study_area, crs = 3857) #convert to utm (easier to work with) - adjust crs as needed.
#study_area2<-st_transform(study_area2, crs = 3857)

study_area_grid <- study_area %>%  #turn study area into gridded area based on 'cell size'
  st_make_grid(cellsize = cellsize, what = "polygons",square=TRUE)
study_area_grid <- st_as_sf(study_area_grid) %>%
  mutate(cell_id = 1:nrow(.)) %>%
  mutate(centroid=st_centroid(study_area_grid)) %>%
  st_intersection(study_area) %>% #remove grids outside study area
  rename(geometry=x)

#2) station file or automate based on receiver spacing for exploratory purposes
station_df<-read.csv("BB_array_df.csv")
station_df_dd <- st_as_sf(station_df, coords = c("lon_nad83", "lat_nad83"),crs = 4326) #crs of input locations
station_df_utm<-st_transform(station_df_dd, crs = 3857) #convert to metric crs
#in the above code, habitat is only associated via the station file, as this is often the approach taken and more simple. See other workflows "DE_array_optimize_1overlap" and DE_array_optimize_3overlap" to associate each cell with a habitat (either the same as the receiver, i.e. receiver-centric or based on cell habitat, i.e., cell-centric). For example, running affiliated code for 'study_area_grid' from the other workflows will create a more complete dataframe with those habitat values.

#######OPTIONAL#########
#instead of existing station locations, create station locations at regular spacing intervals to explore
#grid_points <- st_make_grid(study_area2,cellsize=200, what = "centers")
#grid_points <- st_sf(grid_points, crs = 3857) %>% st_intersection(study_area2)
#station_df_utm <- st_join(grid_points, station_df_utm %>%#allocate 'Habitat' of new station based on nearest 'Habitat'.
#                            select(Habitat,geometry),join = st_nearest_feature) %>%
#                            select(!c(Description,Name)) %>%
#                            rename(geometry=grid_points) %>%
#                            mutate(station=row_number()) %>%
#                            st_zm(drop = TRUE, what = "ZM")
#######OPTIONAL#########

#3) range testing file (or can be loaded further below as a specific model if already known)
range_df<-read.csv("de_range_df.csv") #be sure to include variables that affect detection range as desired (e.g., Habitat)

#######LOAD AND PREPARE NECESSARY FILES#########


#######OPTIONAL#########
#Below adds bathymetry info to each background cell - this is needed if distance_approach=="relative_dynamic"
#library(terra)
#library(exactextractr)
#grid_var<-rast("Depth_tif\\BB_depth.tif")
#grid_var <- project(grid_var, "EPSG:3857")
#grid_var <- aggregate(grid_var, fact = 5, fun = mean) #increase grid size accordingly if too slow (e.g., 5x coarser) or ignore
#plot(grid_var)
#study_area_grid$bottom_depth_cell <- abs(exact_extract(grid_var, study_area_grid, 'mean'))
#######OPTIONAL#########


#######VIEW STUDY AREA AND PARAMETERS#########
ggplot() +
  geom_sf(data = study_area_grid) +
  geom_sf(data=station_df_utm,aes(colour=Habitat),size=3)+
  geom_sf_text(data=station_df_utm,aes(label=station))
#ggsave("study_area.tiff",width = 12, height = 9, units = "in", compression = "lzw", bg = "white",dpi = 300)
#######VIEW STUDY AREA AND PARAMETERS#########


#######DETERMINE DETECTION RANGES AND RESULTANT CELL PROBABILITIES#########
#Calculated using range_df (i.e., range testing).  Assumes a horizontal distance between cell and receiver - which can be further adjusted (i.e., to account for different animal-receiver depths) with the 'distance_approach' option (and is applied in the compute_station_distances function below).
range_model<-glm(prop_detecs~dist_to_receiver+Habitat,data=range_df,family=binomial(link="logit"))

#plot range curves - adjust as needed
ggplot()+
  geom_smooth(data=range_df,aes(x=dist_to_receiver,y=prop_detecs,fill=Habitat,colour=Habitat),method = "glm", se = T,method.args = list(family = "binomial"),linewidth=2,alpha=0.1)+
  labs(x="Distance (m)",y="Detection Probability")+
  scale_y_continuous(breaks=seq(0, 1, 0.25),limits=c(0,1))+
  scale_fill_manual(values = c("Artificial" = "red", "Mixed" = "purple", "Seagrass" = "black"))+
  scale_color_manual(values = c("Artificial" = "red", "Mixed" = "purple", "Seagrass" = "black"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.text = element_text(size=22),legend.title = element_text(size=26),legend.key.size = unit(2,"line"),axis.title=element_text(size=28),axis.text=element_text(size=24))
#ggsave("range_curve.tiff",width = 12, height = 9, units = "in", compression = "lzw", bg = "white",dpi = 300)

#######DETERMINE DETECTION RANGES AND RESULTANT CELL PROBABILITIES#########


#######RECEIVER-GRID DISTANCE CALCULATIONS#########
#calculate distance of every grid from each designated receiver location
handlers(global = TRUE)  # needed once in session
dist_df <- compute_station_distances(
  station_df = station_df_utm,
  study_area_grid,
  study_area,
  max_distance,
  cellsize,
  distance_approach,
  rel_depth_diff,
  water_column_position,
  range_model)
head(dist_df)
#Habitat is locked to the receiver station (as this will be most common approach used - referred as 'receiver-centric detection range'). It can be adjusted to be relative to the habitat of grids (i.e., 'cell-centric detection range), simply assign 'Habitat' relative to cell_id (by merging an existing file or see optimize workflows for example code) and recalculate 'prop_detecs':
#dist_df$prop_detecs<-predict(range_model, newdata = dist_df, type="response")
#dist_df$prop_detecs <- ifelse(dist_df$landinway, 0, dist_df$prop_detecs)
#note that 'distance_approach' parameters have already been accounted for and do not need to be applied again.
#Name is left empty in case user wants to use it for other qualifications
#######RECEIVER-GRID DISTANCE CALCULATIONS#########


#######CALCULATE OVERALL CELL PROBABILITIES#########
dist_df_prob <- compute_cell_prob(
  dist_df,
  study_area_grid,
  n.receiver_overlap
)
head(dist_df_prob)
#######CALCULATE OVERALL CELL PROBABILITIES#########


#######SUMMARY TABLE#########
if(nrow(dist_df_prob %>% filter(prob_cell_all>cell_prob_cutoff))==0) {stop("Error: all grids failed to meet set conditions - cannot go forward")
} else{
dist_df_prob_summ<-dist_df_prob %>% filter(prob_cell_all>cell_prob_cutoff) %>% summarise(n_receiver_overlap=n.receiver_overlap,cellsize=cellsize,max_distance=max_distance,cell_prob_cutoff=cell_prob_cutoff,mean_prob=mean(prob_cell_all),sum_prob_cells=sum(prob_cell_all),ncells=length(cell_id),total_area_m2=length(cell_id)*cellsize^2) %>% st_drop_geometry()
dist_df_prob_summ
#saveRDS(dist_df_prob_summ,"dist_df_prob_summary.rds") #save if desired.
}

#######SUMMARY TABLE#########


#######MAIN PLOT#########
#adjust as necessary
dist_df_prob %>% filter(prob_cell_all>cell_prob_cutoff) %>% #toggle to show cut-off cells
  ggplot() +
  geom_sf(data = study_area,fill="white") +
  geom_sf(aes(fill=prob_cell_all),colour=NA)+
  geom_sf(data=station_df_utm,aes(colour=Habitat))+
  geom_sf(data = study_area,fill=NA) +
  scale_color_manual(values = c("Artificial" = "red", "Mixed" = "purple", "Seagrass" = "black"))+
  labs(fill="Probability\nof detection")+
  ggtitle(paste("n.receiver_overlap = ",n.receiver_overlap))+
  scale_fill_viridis_c(breaks=seq(0, 1, 0.25),limits=c(0,1),option = "viridis",alpha=1)+
  theme(panel.background = element_rect(fill = "grey", colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.text = element_text(size=16),legend.title = element_text(size=18),legend.key.size = unit(1.25,"line"))+
  annotation_scale(line_width = 2, height = unit(0.5, "cm"),location="tr")
#ggsave("Receiver_array_probabilities.tiff",width = 12, height = 9, units = "in", compression = "lzw", bg = "white",dpi = 200)
#######MAIN PLOT#########


