#' @import data.table
NULL

#' Preparing stored data for plotting
#'
#' organizing stored data for later plotting
#'
#' @param station_df An sf receiver-cell distance dataframe (from compute_grid_distances1 function).
#' @param study_area_grid An sf object containing the gridded study area 
#' @param stations_to_use Vector of station names to be used in final array configuration
#' @param n.receiver_overlap  The number of receivers to detect a transmission across grids (e.g., fine-scale vs broad-scale) 
#'
#' @return A data frame with detection probabilities across study area.
#' @export 


#######DE ARRAY OPTIMIZE 1-OVERLAP FUNCTIONS#########
###FUNCTION 4/4

#######DE ARRAY OPTIMIZE 3-OVERLAP FUNCTIONS#########
###FUNCTION 5/5

#Plotting function

plot_output_prep <- function(dist_df, study_area_grid, stations_to_use, n.receiver_overlap) {
  # Assumes compute_cell_prob_optim(p_vec, thr) exists in the environment
  dist.df <- as.data.table(dist_df)
  
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
