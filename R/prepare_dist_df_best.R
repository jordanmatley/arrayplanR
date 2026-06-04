#' @import data.table
NULL

#' Preparing stored data for plotting
#'
#' organizing stored data for later plotting
#'
#' @param loop_hi_subset Dataframe of best performing configuration
#' @param dist_df_iteration Stored info of best performing configuration
#' @param array_size_configs_all Stored info of all configurations
#' @param study_area_grid An sf object containing the gridded study area 
#' @param study_area An sf object containing the study area.
#' @param max_distance Maximum distance from a receiver to include in calculations
#' @param range_model Detection range model based on range testing
#' @param starting_loc An optional identifier for a forced receiver location
#'
#' @return A list of various performance metrics for optimal configuration
#' @export 



#######DE ARRAY OPTIMIZE 3-OVERLAP FUNCTIONS#########

###FUNCTION 4/5
prepare_dist_df_best <- function(loop_hi_subset, dist_df_iteration, array_size_configs_all, study_area_grid, study_area, max_distance, range_model, starting_loc = NULL) {
  
  if (any(is.na(loop_hi_subset$status) | 
          loop_hi_subset$status != "mixed_distances_iteration")) {
    
    stations_to_use <- unlist(loop_hi_subset$stations)
    
    dist_df_best <- dist_df_iteration[[as.numeric(loop_hi_subset$h_iteration)]]
    
  } else {
    message("Warning: since rec_dist_variation=TRUE, the provided configuration may not have optimal spread angles (compared to rec_dist_variation=FALSE). If produced layout is not functional, it is suggested to use the next best configuration (which may be a loop_hi_to_subset row where status is not derived from a mixed_distances_iteration).")
    
    best_mixed_array <- array_size_configs_all[[loop_hi_subset$k_iteration]] %>%
      dplyr::filter(config_id %in% loop_hi_subset$config_id)
    
    stations_to_use <- unique(best_mixed_array$station)
    best_mixed_array<-best_mixed_array %>% rename(hex_id=station) #converting so it matches input for first time function run....
    study_area_buffer2 <- sf::st_buffer(best_mixed_array, max_distance)
    
    dist_df_best <- compute_grid_distances3(
      station_df = best_mixed_array,
      study_area_buffer = study_area_buffer2,
      study_area_grid = study_area_grid,
      study_area = study_area
    )
    
    dist_df_best$prop_detecs <- predict(range_model, newdata = dist_df_best, type = "response")
    dist_df_best<-dist_df_best %>% rename(station=hex_id) #converting so it matches input for first time function run....
    
    if (is.null(dist_df_best$Habitat)) {
      dist_df_best$Habitat <- "no_habitat"
    }
    
    dist_df_best$prop_detecs[dist_df_best$landinway == TRUE] <- 0
  }
  
  # Optional starting location
  starting_loc_plot <- NULL
  if (!is.null(starting_loc)) {
    starting_loc_plot <- sf::st_sfc(sf::st_point(starting_loc), crs = 4326)
  }
  
  list(
    dist_df_best = dist_df_best,
    stations_to_use = stations_to_use,
    starting_loc_plot = starting_loc_plot
  )
}
