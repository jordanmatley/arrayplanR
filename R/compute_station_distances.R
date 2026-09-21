#' @importFrom sf st_buffer
#' @importFrom sf st_cast
#' @importFrom sf st_coordinates
#' @importFrom sf st_crs
#' @importFrom sf st_distance
#' @importFrom sf st_drop_geometry
#' @importFrom sf st_intersects
#' @importFrom sf st_linestring
#' @importFrom sf st_sfc
#' @importFrom dplyr bind_rows
#' @importFrom dplyr mutate
#' @importFrom purrr map
#' @importFrom progressr handlers
#' @importFrom progressr progressor
#' @importFrom progressr with_progress
NULL

#' Compute distances between stations
#'
#' calculate distance of every grid from each designated receiver location
#'
#' @param station_df A .csv object containing station locations.
#' @param study_area_grid An sf object containing the gridded study area
#' @param study_area An sf object containing the study area.
#' @param max_distance Distance from a receiver which you would expect ~0 (or low) detectability
#' @param cellsize Size of grids/cells used to calculate distances between receivers and cells
#' @param distance_approach An optional identifier that enables three-dimensional distances between an animal and a receiver to be included
#' @param rel_depth_diff The relative depth difference (static - i.e., only one value) in metres between animals and receiver locations. Only used when distance_approach="relative_static"
#' @param water_column_position The expected position of an animal in the water column to calculate non-horizontal distances between receiver and animal. Only used when distance_approach="relative_dynamic"
#' @param range_model Detection range model based on range testing
#'
#' @return A data frame of receiver-grid distances.
#' @export


#######DE ARRAY QUANTIFY FUNCTION#########

###FUNCTION 1/2
#calculate distance of every grid from each designated receiver location
compute_station_distances <- function(station_df,
                                      study_area_grid,
                                      study_area,
                                      max_distance,
                                      cellsize,
                                      distance_approach="horizontal",
                                      rel_depth_diff=NULL,
                                      water_column_position=NULL,
                                      range_model)   {

  if (max_distance <= cellsize) {
    stop("Error: your grid size (cellsize) is larger than max_distance.")
  }

  # Pre-cast study area to multilinestring for land barrier tests
  study_area_line <- st_cast(study_area, "MULTILINESTRING")

  # Pre-buffer all stations once
  station_buffers <- st_buffer(station_df, max_distance)

  # Split by station
  station_split <- split(station_df, station_df$station)

  handlers("txtprogressbar")
  result_list <- with_progress({
    p <- progressor(along = station_split)

    map(names(station_split), function(st_id) {
      p()

      # station_df row for this station
      station_loop <- station_split[[st_id]]

      # station attributes (all non-geometry columns)
      station_atts <- station_loop %>% st_drop_geometry()

      # buffer for this station (precomputed)
      station_buf <- station_buffers[station_buffers$station == st_id, ]

      # FAST grid filtering (bounding intersects)
      idx <- st_intersects(study_area_grid, station_buf, sparse = FALSE)[, 1]
      study_buffer_loop <- study_area_grid[idx, ]

      n_grid <- nrow(study_buffer_loop)
      if (n_grid == 0) return(NULL)

      # ---- FAST VECTORISED DISTANCE ----
      dist_to_receiver <- as.numeric(
        st_distance(study_buffer_loop$centroid, station_loop$geometry)
      )

      # ---- FAST PURE-SF LINES ----
      station_cell_line <- st_sfc(
        mapply(
          function(p1, p2)
            st_linestring(rbind(st_coordinates(p1), st_coordinates(p2))),
          station_loop$geometry,
          study_buffer_loop$centroid,
          SIMPLIFY = FALSE
        ),
        crs = st_crs(station_loop)
      )

      # ---- Land barrier test ----
      landinway <- st_intersects(station_cell_line, study_area_line, sparse = FALSE)[, 1]

      # ---- Combine results (preserve ALL station + grid attributes) ----
      out <- study_buffer_loop %>%
        mutate(
          dist_to_receiver = dist_to_receiver,
          landinway = landinway
        )

      # attach station attributes to every row
      for (col in names(station_atts)) {
        out[[col]] <- station_atts[[col]]
      }

      # attach station geometry explicitly
      out$station_loc <- station_loop$geometry

      ######
      #if using bathymetry - need to add depth of station (from background study_area grids)
      if ("bottom_depth_cell" %in% names(study_area_grid)) {

        station_depth_idx <- st_intersects(
          station_loop$geometry,
          study_area_grid
        )

        station_depth <- if (length(station_depth_idx[[1]]) == 0) {
          NA_real_
        } else {
          study_area_grid$bottom_depth_cell[
            station_depth_idx[[1]][1]
          ]
        }

        out$station_depth <- station_depth
      }

      # ==========================================================
      # FOR DISTANCE APPROACHES
      # ==========================================================
      #if using a relative depth approach (because tag-receiver is not necessarily horizontal) - assumes receiver deployed on bottom

      if (distance_approach == "relative_static" &&
          !is.null(rel_depth_diff)) {

        out$dist_to_receiver <- sqrt(
          out$dist_to_receiver^2 +
            abs(rel_depth_diff)^2
        )

      } else if (
        distance_approach == "relative_dynamic" &&
        "bottom_depth_cell" %in% names(out) &&
        !is.null(water_column_position)
      ) {

        if (water_column_position == "surface") {

          out$dist_to_receiver <- sqrt(
            out$dist_to_receiver^2 +
              abs(out$station_depth)^2
          )

        } else if (water_column_position == "midwater") {

          out$dist_to_receiver <- sqrt(
            out$dist_to_receiver^2 +
              abs(
                out$station_depth -
                  (out$bottom_depth_cell / 2)
              )^2
          )

        } else if (water_column_position == "bottom") {

          out$dist_to_receiver <- sqrt(
            out$dist_to_receiver^2 +
              abs(
                out$station_depth -
                  out$bottom_depth_cell
              )^2
          )

        } else {
          stop("Error - water_column_position not appropriately organized")
        }

      } else if (
        is.na(distance_approach) |
        is.null(distance_approach) |
        distance_approach == "horizontal"
      ) {

        # keep horizontal distance only

      } else {
        stop("Error - distance_approach not appropriately organized")
      }

      # ==========================================================
      # Detection probability
      # ==========================================================

      out$prop_detecs <- predict(
        range_model,
        newdata = out,
        type = "response"
      )

      out$prop_detecs <- ifelse(
        out$landinway == TRUE,
        0,
        out$prop_detecs
      )
      ######
      return(out)

    })
  })

  bind_rows(result_list)
}
