#' @importFrom sf st_cast
#' @importFrom sf st_coordinates
#' @importFrom sf st_crs
#' @importFrom sf st_distance
#' @importFrom sf st_geometry
#' @importFrom sf st_intersection
#' @importFrom sf st_intersects
#' @importFrom sf st_linestring
#' @importFrom sf st_sfc
#' @importFrom sf st_transform
#' @importFrom dplyr bind_rows
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom progressr handlers
#' @importFrom progressr progressor
#' @importFrom progressr with_progress
NULL

#' Compute distances between stations (>1-receiver overlap)
#'
#' calculate distance of every grid from every other possible receiver location
#'
#' @param station_df #station information for configuration iteration
#' @param study_area_buffer #distance buffer around all receivers from station_df
#' @param study_area_grid An sf object containing the gridded study area
#' @param study_area An sf object containing the study area.
#'
#' @return A data frame of receiver-grid distances.
#' @export

#######DE ARRAY OPTIMIZE 3-OVERLAP FUNCTIONS#########

###FUNCTION 1/5
#calculate distance of every grid from each designated receiver location
compute_grid_distances3 <- function(station_df,
                                    study_area_buffer,
                                    study_area_grid,
                                    study_area) {

  #initial checks
  target_crs <- st_crs(study_area)
  if (st_crs(station_df) != target_crs) station_df <- st_transform(station_df, target_crs)
  if (st_crs(study_area_buffer) != target_crs) study_area_buffer <- st_transform(study_area_buffer, target_crs)
  if (st_crs(study_area_grid) != target_crs) study_area_grid <- st_transform(study_area_grid, target_crs)

  # For land barrier
  study_area_line <- st_cast(study_area, "MULTILINESTRING")

  station_split <- split(station_df, station_df$hex_id) # Split hex_id stations [not background cells because they are just reference habitats and distances - including potentially duplicate stations since hex_ids may fall in duplicate station/cell_ids]

  result_list <- with_progress({
    handlers("cli")
    p <- progressor(along = station_split)

    lapply(station_split, function(station_df_loop) {
      p()
      station_i <- station_df_loop$hex_id[1] #is hex_id/proper station id
      station_geom <- station_df_loop$geometry[1] #is geometry of hex_id


      # Subset buffer & intersect with grid
      suppressWarnings(
        study_buffer_loop <- study_area_buffer %>%
          filter(hex_id == station_i) %>% #with new hex_id fix
          st_intersection(study_area_grid %>% select(-Name) %>% rename(background_grid_habitat=Habitat))
      )
      #centroid = background cell_id centroid
      #geometry = geometry of used cell_id (may be partially used)

      if (nrow(study_buffer_loop) == 0L) return(NULL)

      study_buffer_loop$centroid <- study_area_grid$centroid[
        match(study_buffer_loop$cell_id, study_area_grid$cell_id)
      ]
      centroids_sfc <- st_geometry(study_buffer_loop$centroid)

      # Compute distances
      dist_to_receiver <- as.numeric(st_distance(centroids_sfc, station_geom)) #background cell centroids and hex_id station loc

      # Barrier check
      station_coords <- matrix(st_coordinates(st_geometry(station_geom)), ncol = 2)
      centroid_coords <- st_coordinates(centroids_sfc)[,1:2, drop = FALSE]
      # vectorized LINESTRING creation using lapply
      station_cell_line <- st_sfc(
        lapply(seq_len(nrow(centroid_coords)), function(j) {
          st_linestring(rbind(station_coords, centroid_coords[j, ]))
        }),
        crs = target_crs
      )
      landinway <- lengths(st_intersects(station_cell_line, study_area_line)) > 0

      # Combine results
      study_buffer_loop %>%
        mutate(dist_to_receiver = dist_to_receiver,
               station_loc = station_geom,
               landinway = landinway)

    })
  })

  # Drop empty entries and merge
  result_list <- Filter(Negate(is.null), result_list)
  bind_rows(result_list)
}

