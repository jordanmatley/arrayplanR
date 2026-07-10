#' @importFrom sf st_as_sf
#' @importFrom sf st_buffer
#' @importFrom sf st_point
#' @importFrom sf st_sfc
#' @importFrom dplyr filter
#' @importFrom dplyr rename
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
prepare_dist_df_best <- function(loop_hi_subset,
                                 dist_df_iteration,
                                 array_size_configs_all,
                                 study_area_grid,
                                 study_area,
                                 max_distance,
                                 range_model,
                                 starting_loc = NULL) {

  if (any(is.na(loop_hi_subset$status) |
          loop_hi_subset$status != "mixed_distances_iteration")) {

    stations_to_use <- unlist(loop_hi_subset$stations)

    dist_df_best <- dist_df_iteration[[as.numeric(loop_hi_subset$h_iteration)]]

  } else {
    message("Warning: since rec_dist_variation=variable, the provided configuration may not have optimal spread angles (compared to rec_dist_variation=static). If produced layout is not functional, it is suggested to use the next best configuration (which may be a loop_hi_to_subset row where status is not derived from a mixed_distances_iteration).")

    best_mixed_array <- array_size_configs_all[[loop_hi_subset$k_iteration]] %>%
      filter(config_id %in% loop_hi_subset$config_id)

    stations_to_use <- unique(best_mixed_array$station)
    best_mixed_array<-best_mixed_array %>% rename(hex_id=station) #converting so it matches input for first time function run....
    study_area_buffer2 <- st_buffer(best_mixed_array, max_distance)

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
    starting_loc_plot <- st_sfc(st_point(starting_loc), crs = 4326)
  }

  list(
    dist_df_best = dist_df_best,
    stations_to_use = stations_to_use,
    starting_loc_plot = starting_loc_plot
  )
}

###FUNCTION 5/5
#Plotting function
#duplicate
plot_output_prep <- function(dist_df_best,
                             study_area_grid,
                             stations_to_use,
                             n.receiver_overlap=3) {

  # Assumes compute_cell_prob_optim(p_vec, thr) exists in the environment
  dist.df <- as.data.table(dist_df_best)

  # Subset to relevant stations first
  dist_sub <- dist.df[station %in% stations_to_use]

  if (n.receiver_overlap == 1) {
    probs_combined <- dist_sub[, .(opp_prob = 1 - prop_detecs), by = cell_id
    ][, .(prob_cell_all = 1 - prod(opp_prob)), by = cell_id]
  } else {
    # Build one list of prop_detecs per cell and compute tail via FFT function
    dt_prop_list <- dist_sub[, .(prop_vec = list(prop_detecs)), by = cell_id]

    # If some cells have fewer contributors than the required overlap, return 0 for them
    dt_prop_list[, prob_cell_all := vapply(prop_vec, function(v) {
      if (length(v) < n.receiver_overlap) return(0)
      compute_cell_prob_optim(v, n.receiver_overlap)
    }, numeric(1))]
    probs_combined <- dt_prop_list[, .(cell_id, prob_cell_all)]
  }

  # Merge with study_area_grid (keep all grid cells), set NA -> 0
  to_plot <- merge(
    unique(study_area_grid[, c("cell_id", "geometry", "centroid")]),
    probs_combined,
    by = "cell_id",
    all.x = TRUE
  )
  to_plot$prob_cell_all[is.na(to_plot$prob_cell_all)] <- 0
  st_as_sf(to_plot)
}

