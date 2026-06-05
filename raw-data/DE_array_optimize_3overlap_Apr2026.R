###Optimizing arrays without pre-designated locations by estimating the best receiver locations for >1-receiver overlaps (only) based on input parameters###

#This code is only optimized for when the input n.receiver_overlap>1 (e.g., fine-scale tracking). It is suggested to use 'DE_array_optimize_1overlap' workflow if n.receiver_overlap=1.

#Assumes a user has some detection range data already conducted - see input files for comparison.

#Default approach uses 'receiver-centric' detection ranges - each receiver has its own detection range profile and cells around it have same profile, as opposed to 'cell-centric' detection ranges - each cell has its own detection range profile independent of the receivers around it. The main reason is because many studies will not have representative cell-specific detection range profiles throughout study area. Note, 'cell-centric' detection range capacity is still possible but needs to be manually adjusted (see optional code chunks). 

#In example code, 'Habitat' is main factor distinguishing different detection range profiles across receivers. Other variables can easily be built in, but the term 'Habitat' is used to control specific parameters within functions (e.g., optim_method="habitat_specific") and is a needed entry for those to work.

#Evaluation of array efficiency centered around prioritizing the 'effective area covered'  or 'effective sum of probabilities' above a defined probability threshold.  

#To incorporate temporal variability, it is suggested that separate evaluations are used.

#It is suggested the user runs example code prior to using their own to understand the process.

#Optimizing gridded array positions can be computationally cumbersome. There are typically hundreds of millions of possible receiver-cell combinations when cell sizes are small (and study areas are large). Therefore, an iterative stepwise estimation process is used to minimize the need for all combinations. Specifically, this approach selects individual cells (square or hexagonal and that overlay a square background grid layer) one at a time and tests performance for nrecs-1 locations that are nearest to the focal receiver. This ensures that spread angles are high (e.g., triangular arrangement of receivers). Setting 'rec_dist_iterations' allows for different receiver spacings to optimize performance (across all tested cell locations). Note, further investigations can be performed using the web application. It is suggested that users start broad to evaluate processing times (e.g., large 'cellsize', low 'subsampl_perc', low 'max_i_iterations') before progressing further.


library(dplyr)
library(ggplot2)
library(sf)
library(tidyr)
library(progressr)
library(data.table)
library(digest)
library(proxy) 
library(purrr)
theme_set(theme_bw())

setwd("C:\\Users\\jorda\\OneDrive\\Research\\Flinders\\Research\\Other\\Array_DE\\R\\for_github")

#######SET NECESSARY PARAMETERS#########
seed<-99 #Options: Numeric (Default = 99). Controls reproducibility since optimal solutions may change 
cellsize=100 #Options: Numeric (in metres) (Default = 100). e.g., 50 = 50 x 50 m. Size of grids/cells used to calculate distances between receivers and cells (a smaller cell better for resolution but can slow computation)
max_distance=500 #Options: Numeric (in metres) (Default = 500). Distance from a receiver which you would expect ~0 (or low) detectability. Should be chosen based on range tests. Used to save computing time by ignoring areas not in range. 
cell_prob_cutoff=0.1 #Options: Proportion (Default = 0.1). Equivalent to a x% chance of a transmission from a grid being heard (by n receivers). Should be chosen based on how stringent a user wants detection probabilities to be (e.g., only interested in grids with low vs high probability of detection).
n.receiver_overlap=3 #Options: Numeric (Default = 3). Set how many receivers to detect transmission from a grid (e.g., broad-scale tracking = 1, fine-scale positioning = 3+). Note, this code is only optimized for n.receiver_overlap>1. If n.receiver_overlap=1, it is suggested to use 'DE_array_optimize_1overlap' workflow.
nrecs=9 #Options: Numeric (Default = 9). The number of receivers to be deployed.
optim_method="all_possible" #Options: "all_possible","min_distance","habitat_specific" (Default = "all_possible"). The method used to optimize receivre array. Options are based on whether all cells are considered ("all_possible") OR whether a specific subset of receivers by habitat ("habitat_specific") OR whether stations are to be set at a minimum distance from each other ("min_distance"). Additional input parameters needed if "habitat_specific" is called (see 'CALCULATE OVERALL CELL PROBABILITIES' section below). Similarly, additional input parameter ('min_rec_distance') needed for "min_distance" approach.
optim_metric="sum_prob_cells" #Options: "sum_prob_cells", "total_area_m2" (Default = "sum_prob_cells"). Approach used to optimize array performance. Either based on the sum of all cell probability values above cell_prob_cutoff ("sum_prob_cells") or total area of cells with probability above cell_prob_cutoff ("total_area_m2")
min_rec_distance=100 #Options: Numeric (Default = 100). Only used when optim_method="min_distance". The minimum distance that receivers should be apart.
starting_loc=NULL #Options: NULL or c(lon,lat) (Default = NULL). Setting a forced receiver location to optimize the array around. When not NULL, assumes input is decimal degrees from Google Earth: EPSG:4326). Example location for this tutorial is c(-64.980, 18.340).
rec_dist_iterations=round(seq(from = 100, to = 250, length.out = 4)) #Options: vector of receiver spacing (Default: c(100,150,200,250)). Different receiver spacing distances (h iterations) to test for optimization estimates across adjacent n receivers (i iterations). To clarify, an h iteration consists of a pre-defined receiver spacing distance (e.g., all receivers are 100m apart). An i iteration loops through each receiver in the h iteration and combines with n of the nearest receivers to it to estimate detection probabilities. 
square_array_config=FALSE #Options: TRUE or FALSE (Default = FALSE). FALSE provides a hexagonal/triangular receiver array configuration (suggested at least for n.receiver_overlap>1), whereas TRUE provides a square receiver array configuration. Note, there is another (presumably smaller) square cell/grid layer as background to this that habitat/distances are drawn from.
rec_dist_variation=FALSE #Options: TRUE or FALSE (Default = FALSE). TRUE allows for receivers to be spaced at different distances (e.g., combinations of rec_dist_iterations) in addition to the static spacing designated. FALSE only incorporates iterations with the same receiver spacing. If habitat (or similar) is expected to drastically change detection range, then TRUE may be a more suitable option; however, it will increase processing time and only a subset of additional options are tested. Note, there is no guarantee 'rec_dist_variation=TRUE' will result in well-aligned configurations (high spread angles) and should be verified for feasibility. Specifically, this approach consists of re-analyzing k iterations - a k iteration is a single row from the main summary of all configuration performances (loop_hi_to_subet). A k iteration goes through the output of each row and works through j iterations, which are different mixed receiver distance combinations to be tested. For example, if there are 9 receivers across 3 habitats (3 receivers in each), each j iteration will freeze locations of receivers in 2 habitats and cycle through new receiver locations based on 'rec_dist_iterations' (selecting new locations in different sized grids that are closest to the receivers they are replacing - to ensure they remain in the general area). To ensure stability as much as possible (i.e., keeping optimal spread angles), only one habitat is cycled at a time (instead of two in the above example) - additional investigations are suggested with the interactive web application.  
distance_approach="horizontal" #Options: "horizontal","relative_static","relative_dynamic" (Default = "horizontal"). Enables user to account for distances between an animal and a receiver that may occur at different depths in relation to each other. "horizontal" means no depth differences are accounted for (animal-receiver assumed to be aligned horizontally). "relative_static" allows for a singular relative depth difference (positive value) between animal and receiver (e.g., if study area and expected animal depth use consistent). No bathymetry data are needed for these first two options. However, "relative_static" requires a numeric value (i.e., absolute relative depth difference)(see rel_depth_diff option). "relative_dynamic" incorporates bottom depth and the relative position an animal is expected to be in - "surface", "midwater", "bottom" (see water_column_position option), providing 3-dimensional estimates of detection probabilities.
rel_depth_diff<-NULL #Options: Numeric (Default = NULL). Only used when distance_approach="relative_static". The relative depth difference (static - i.e., only one value) in metres between animals and receiver locations. This option allows for non-horizontal differences in relatively homogeneous environments.
water_column_position=NULL #Options: "surface", "midwater", "bottom" (Default = NULL). Only used when distance_approach="relative_dynamic". The expected position of an animal in the water column to calculate non-horizontal distances between receiver and animal (where relevant).

#######SET NECESSARY PARAMETERS#########


#######LOAD AND PREPARE NECESSARY FILES#########
#1) study area file (gis file (e.g., .shp) or can be added manually from Google Earth)
study_area = st_read("BB_shp\\eaarlb_STT_250m_WaterMask_bbaypbay.shp") #using existing shp file for study area
study_area2 = st_read("Study_area_cut.kml") #created a smaller polygon for receivers while allowing detection range to spill over into study_area. Kml object created directly in Google Earth using polygon function. If not desired, simply use study_area2<-study_area.
#st_write(study_area, "Brewers_JM.geojson") #to export as a GeoJSON file to web application if desired (transform to correct crs first).
study_area<-st_transform(study_area, crs = 3857) #convert to utm (much easier to work with and understand)
study_area2<-st_transform(study_area2, crs = 3857)

#creating the larger background square grid to allow for detection range spill over (from smaller receiver polygon)
study_area_grid_pre <- study_area %>% 
  st_make_grid(cellsize = cellsize, what = "polygons",square=TRUE) 
study_area_grid_pre <- st_as_sf(study_area_grid_pre) %>% 
  mutate(cell_id = 1:nrow(.)) %>%
  mutate(centroid=st_centroid(study_area_grid_pre)) %>%
  st_intersection(study_area) %>%
  mutate(Name="study_spillover")

#layout smaller study area grid where receivers will go
study_area_grid_pre<-study_area_grid_pre %>% dplyr::select(x,cell_id,Name,centroid) %>% rename(geometry=x)
study_area_grid_pre<-st_join(study_area_grid_pre,st_as_sf(study_area2 %>% dplyr::select(!Name) %>% mutate(Description="study_area")),join = st_intersects,left = TRUE, largest = TRUE)

#2) Incorporate habitat/environmental variables that will impact detection range (will be needed for every cell of study_area2 with receiver-centric approach)

#Note, if no 'Habitat' or equivalent is being used, simply ignore section #2) below and run the following: 
#study_area_grid_pre$Habitat="no_habitat"
#study_area_grid<-study_area_grid_pre

#Habitat shp/kml (here example of loading 3 separate polygons made in Google Earth)
Hab1 = st_read("Artificial_Habitat.kml")
Hab2 = st_read("Seagrass_Habitat.kml")
Hab3 = st_read("Mixed_Habitat.kml")
Habitat_areas<-rbind(Hab1,Hab2,Hab3)
Habitat_areas<-st_transform(Habitat_areas, crs = 3857)
plot(Habitat_areas)

#add habitat areas per cell (within study area only)
study_area_grid<-st_join(study_area_grid_pre,st_as_sf(Habitat_areas %>% rename(Habitat=Name) %>% dplyr::select(!Description)),join = st_intersects,left = TRUE, largest = TRUE)
anyNA(study_area_grid$Habitat)  #does every cell has a habitat designation - FALSE means it does (TRUE is fine for us here but if using cell-centric approach it should be FALSE)
study_area_grid$Description<-ifelse(!is.na(study_area_grid$Habitat),"study_area","study_spillover") #differentiating study area where receivers can go
points_within_study<-study_area_grid$centroid %>% st_within(study_area,sparse = FALSE) #clip cells with centroids specifically that are not on water
study_area_grid <- study_area_grid[points_within_study, ] 
head(study_area_grid) #Description is NA if not in main receiver study area (i.e., study_area2). Habitats will be NA except when in main study area

#######OPTIONAL#########
#Below adds bathymetry info to each background cell - this is needed if distance_approach=="relative_dynamic"
#library(terra)
#library(exactextractr)
#grid_var<-rast("Depth_tif\\BB_depth.tif")
#grid_var <- project(grid_var, "EPSG:3857")
#grid_var <- aggregate(grid_var, fact = 5, fun = mean) #increase grid size accordingly if too slow (e.g., 5x coarser) or ignore
#plot(grid_var)
#study_area_grid$bottom_depth_cell <- abs(exact_extract(grid_var, study_area_grid, 'mean'))
#Above is also example of how a cell-centric approach (as opposed to receiver-centric) could be used - i.e., associate habitat based on study_area_grid for each cell_id instead of the habitat associated with each receiver station. Adjust code as needed.
#######OPTIONAL#########

#3) range testing file (or can be loaded further below as a specific model if already known)
range_df<-read.csv("de_range_df.csv")

#######LOAD AND PREPARE NECESSARY FILES#########


#######VIEW STUDY AREA AND PARAMETERS#########
ggplot() + 
  geom_sf(data=study_area,fill="white")+
  geom_sf(data = study_area_grid %>% filter(!is.na(Habitat)),aes(fill=Habitat),colour="white")+
  scale_fill_manual(values = c("Artificial" = "red", "Mixed" = "purple", "Seagrass" = "black"))+
  geom_sf(data = study_area,fill=NA)+
  labs(fill="Habitat")+
  theme(panel.background = element_rect(fill = "grey", colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text=element_text(size=15),legend.text = element_text(size=16),legend.title = element_text(size=18),legend.key.size = unit(1.25,"line"))+
  ggspatial::annotation_scale(line_width = 2, height = unit(1, "cm"),location="tr",text_cex = 1.5)
#ggsave("study_area_3optimizing.tiff",width = 12, height = 9, units = "in", compression = "lzw", bg = "white",dpi = 200)

#######VIEW STUDY AREA AND PARAMETERS#########


#######DETERMINE DETECTION RANGES AND RESULTANT CELL PROBABILITIES#########
#Calculated using range_df (i.e., range testing).  Assumes a horizontal distance between cell and receiver - which can be further adjusted (i.e., to account for different animal-receiver depths) with the 'distance_approach' option (and is applied in the optimize_3receivers function).
range_model<-glm(prop_detecs~dist_to_receiver+Habitat,data=range_df,family=binomial(link="logit"))

ggplot()+
  geom_smooth(data=range_df,aes(x=dist_to_receiver,y=prop_detecs,fill=Habitat,colour=Habitat),method = "glm", se = T,method.args = list(family = "binomial"),linewidth=2,alpha=0.1)+
  labs(x="Distance (m)",y="Detection Probability")+
  scale_y_continuous(breaks=seq(0, 1, 0.25),limits=c(0,1))+
  scale_fill_manual(values = c("Artificial" = "red", "Mixed" = "purple", "Seagrass" = "black"))+
  scale_color_manual(values = c("Artificial" = "red", "Mixed" = "purple", "Seagrass" = "black"))+ 
  theme(panel.background = element_rect(fill = "white", colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.text = element_text(size=22),legend.title = element_text(size=26),legend.key.size = unit(2,"line"),axis.title=element_text(size=28),axis.text=element_text(size=24))
#ggsave("range_curve_optimizing.tiff",width = 12, height = 9, units = "in", compression = "lzw", bg = "white",dpi = 300)

#######DETERMINE DETECTION RANGES AND RESULTANT CELL PROBABILITIES#########


#######CALCULATE OVERALL CELL PROBABILITIES#########
#if using "habitat_specific" optim_method, habitat-receiver allocation needs to be set here 
if(optim_method=="habitat_specific"){
  habitat_levels=unique(sort(range_df$Habitat))
  habitat_levels #order for line below
  hab_nrec_alloc=c(4,0,5) #adjust to desired habitat areas (e.g., 1 receiver in "Artificial", 3 in "Mixed", 1 in "Seagrass)
  hab_nrec_alloc <- rep(habitat_levels, times = hab_nrec_alloc)
  hab_nrec_alloc
}else{hab_nrec_alloc=NULL}

#main optimization function
loop_hi_to_subset_list <- optimize_3receivers(study_area_grid,hab_nrec_alloc,starting_loc)

loop_hi_to_subset<-loop_hi_to_subset_list$loop_hi_to_subset #access summary info (all h and i outputs provided)
head(loop_hi_to_subset)
#saveRDS(loop_hi_to_subset,"loop_hi_output_optim3.rds") #save if desired.

#store iterations for future plotting/output
if(rec_dist_variation == "TRUE"){array_size_configs_all<-loop_hi_to_subset_list$array_size_configs_all}
rec_distance_tests<-loop_hi_to_subset_list$rec_distance_tests
dist_df_iteration<-loop_hi_to_subset_list$dist_df_iteration
station_possible_iteration<-loop_hi_to_subset_list$station_possible_iteration

#######CALCULATE OVERALL CELL PROBABILITIES#########


#######SUMMARY TABLE/INFO#########
#loop_hi_to_subset<-loop_hi_to_subset[sample(1:nrow(loop_hi_to_subset)),] #randomize to alter max chosen if ties (if desired)

#enforcing suggested use of convex hulls of receivers (with small buffer of the minimum of h-iteration distances) as areas to optimize instead of absolute total area or sum prob cells because the latter do not necessarily provide realistic results (i.e., they tend to group receivers close together to increase area covered around them)
if(optim_metric=="total_area_m2"){optim_metric2<-"convex_total_area_m2"} else if(optim_metric=="sum_prob_cells") {optim_metric2<-"convex_sum_prob_cells"} 
loop_hi_subset<-loop_hi_to_subset %>% filter(status=="mixed_distances_iteration"|is.na(status)) %>% slice_max(order_by = .data[[optim_metric2]], n = 1, with_ties = FALSE) #select top iteration

#to access/plot other results manually
#note duplicated receiver stations (with same performance) may occur when using changing receiver spacing option (they have not been removed for posterity)
#loop_hi_subset<-loop_hi_to_subset %>% filter(hi_iteration=="1-51" & is.na(config_id)) #for static receiver distances
#loop_hi_subset<-loop_hi_to_subset %>% filter(hi_iteration=="3-14" & config_id==9) # for mixed distances iterations

loop_hi_subset
#######SUMMARY TABLE/INFO#########

#######RE-ACCESS BEST ITERATION AND PREPARE FOR PLOTTING#########
if (any(is.na(loop_hi_subset$status) | loop_hi_subset$status != "mixed_distances_iteration")) { #different ways needed to re-access data depending on whether it is static rec_dist_iterations or mixed distances
  
  stations_to_use<-unlist(loop_hi_subset$stations) #access those best stations
  #select best h iteration array layout and display to user coordinates of best array configuration
  dist_df_best<-dist_df_iteration[[as.numeric(loop_hi_subset$h_iteration)]]
  
  if (!is.null(starting_loc)){ 
    starting_loc_plot <- st_sfc(st_point(starting_loc), crs = 4326)
  }
  
} else{
  message("Warning: since rec_dist_variation=TRUE, the provided configuration may not have optimal spread angles (compared to rec_dist_variation=FALSE). If produced layout is not functional, it is suggested to use the next best configuration (which may be a loop_hi_to_subset row where status is not derived from a mixed_distances_iteration).")
  best_mixed_array<-array_size_configs_all[[loop_hi_subset$k_iteration]] %>% filter(config_id %in% loop_hi_subset$config_id)
  best_mixed_array
  stations_to_use<-unique(best_mixed_array$station)
  best_mixed_array<-best_mixed_array %>% rename(hex_id=station) #converting so it matches input for first time function run....
  study_area_buffer2<-st_buffer(best_mixed_array,max_distance)
  
  dist_df_best <- compute_grid_distances3(
    station_df = best_mixed_array,
    study_area_buffer = study_area_buffer2,
    study_area_grid = study_area_grid,
    study_area = study_area
  )
  dist_df_best$prop_detecs <- predict(range_model, newdata = dist_df_best, type = "response")
  if (is.null(dist_df_best$Habitat)) dist_df_best$Habitat <- "no_habitat"
  dist_df_best$prop_detecs[dist_df_best$landinway == TRUE] <- 0
  dist_df_best<-dist_df_best %>% rename(station=hex_id) #converting so it matches input...
  
  if (!is.null(starting_loc)){ 
    starting_loc_plot <- st_sfc(st_point(starting_loc), crs = 4326)  
  }
}

prep <- prepare_dist_df_best(loop_hi_subset,dist_df_iteration,array_size_configs_all,study_area_grid,study_area,max_distance,range_model,starting_loc)

# unpack into standalone objects
dist_df_best      <- prep$dist_df_best
stations_to_use   <- prep$stations_to_use
starting_loc_plot <- prep$starting_loc_plot

#######RE-ACCESS BEST ITERATION AND PREPARE FOR PLOTTING#########

#######PLOT (with FUNCTION)#########
#to prepare for plotting (e.g., detection probabilities for receiver locations)
to_plot_sf<-plot_output_prep(dist_df_best, study_area_grid, stations_to_use, n.receiver_overlap)

#adjust as necessary
to_plot_sf %>% filter(prob_cell_all>=cell_prob_cutoff) %>% #toggle to show cut-off cells
  ggplot() + 
  geom_sf(data = study_area,fill="white") +
  geom_sf(aes(fill=prob_cell_all),colour=NA,show.legend=TRUE)+
  #geom_sf(data=station_possible_iteration[[loop_hi_subset$h_iteration]],aes(geometry=station_loc),colour="pink")+
  #geom_sf_text(data=dist_df_best %>% filter(station %in% stations_to_use) %>% distinct(station,station_loc),aes(geometry=station_loc, label=station),colour="red")+
  geom_sf(data=dist_df_best %>% filter(station %in% stations_to_use) %>% distinct(station,station_loc,.keep_all=TRUE),aes(geometry=station_loc,colour=Habitat),show.legend=TRUE)+
  #geom_sf(data=starting_loc_plot,colour="Yellow")+ #do not use if starting_loc=NULL
  geom_sf(data = study_area,fill=NA) +
  labs(fill="Probability\nof detection",title=paste("Optimisation:",optim_method,"| Metric:",optim_metric,"| Overlap:",n.receiver_overlap))+
  scale_colour_manual(values = c("Artificial" = "red", "Mixed" = "purple", "Seagrass" = "black"))+
  scale_fill_viridis_c(breaks=seq(0, 1, 0.25),limits=c(0,1),option = "viridis",alpha=1)+
  theme(panel.background = element_rect(fill = "grey", colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.ticks=element_blank(),axis.text=element_blank(),axis.title=element_blank(),legend.text = element_text(size=16),legend.title = element_text(size=18),legend.key.size = unit(1.25,"line"))+
  ggspatial::annotation_scale(line_width = 2, height = unit(1, "cm"),location="tr",text_cex = 1.25)
#ggsave("Receiver_array_probabilities_3optimized.tiff",width = 12, height = 9, units = "in", compression = "lzw", bg = "white",dpi = 200)

#to access array info
stations_locs_best1<-dist_df_best %>% filter(station %in% stations_to_use) %>% distinct(station,station_loc) %>% rename(station_loc_proj=station_loc)
stations_locs_best2 <- st_as_sf(stations_locs_best1, crs = st_crs(study_area))
stations_locs_best2 <- st_transform(stations_locs_best2, crs=4326) %>% dplyr::select(-station) %>% rename(station_loc=station_loc_proj) #or adjust crs as needed to original input
stations_locs_best<-cbind(stations_locs_best1,stations_locs_best2)
stations_locs_best
#######PLOT (with FUNCTION)#########