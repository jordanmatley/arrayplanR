#' @importFrom sf st_cast
#' @importFrom sf st_centroid
#' @importFrom sf st_coordinates
#' @importFrom sf st_crs
#' @importFrom sf st_drop_geometry
#' @importFrom sf st_geometry
#' @importFrom sf st_intersects
#' @importFrom sf st_linestring
#' @importFrom sf st_sfc
#' @importFrom dplyr left_join
#' @importFrom dplyr select
#' @importFrom progressr handlers
#' @importFrom progressr progressor
#' @importFrom proxy dist
#' @import data.table
NULL

#' Compute distances between stations (1-receiver overlap)
#'
#' calculate distance of every grid from every other possible receiver location
#'
#' @param study_area_grid An sf object containing the gridded study area
#' @param study_area An sf object containing the study area
#' @param max_distance Distance from a receiver which you would expect ~0 (or low) detectability
#' @param cellsize Size of grids/cells used to calculate distances between receivers and cells
#' @param distance_approach An optional identifier that enables three-dimensional distances between an animal and a receiver to be included
#' @param rel_depth_diff The relative depth difference (static - i.e., only one value) in metres between animals and receiver locations. Only used when distance_approach="relative_static"
#' @param water_column_position The expected position of an animal in the water column to calculate non-horizontal distances between receiver and animal. Only used when distance_approach="relative_dynamic"
#'
#' @return A data frame of receiver-grid distances
#' @export


#######DE ARRAY OPTIMIZE 1-OVERLAP FUNCTIONS#########

###FUNCTION 1/4
#calculate distance of every grid from every other possible receiver location
compute_grid_distances1 <- function(study_area_grid,
                                    study_area,
                                    max_distance,
                                    cellsize,
                                    distance_approach="horizontal",
                                    rel_depth_diff=NULL,
                                    water_column_position=NULL) {
  handlers("txtprogressbar")

  #initial checks
  if(max_distance<=cellsize){stop("Error: your grid size (cellsize) is larger than the expected max range (max_distance) - not useful to continue")}
  if (!st_crs(study_area_grid) == st_crs(study_area)) {
    st_crs(study_area) <- st_crs(study_area_grid)
  }

  #Subset station (filtered study area)
  station <- study_area_grid[study_area_grid$Description == "study_area", ]
  station$centroid <- st_centroid(st_geometry(station))

  if (nrow(station) == 0) {
    return(data.table(station = character(),
                      cell_id = character(),
                      dist_to_receiver = numeric(),
                      landinway = logical()))
  }

  cell_id <- study_area_grid
  cell_id$centroid <- st_centroid(st_geometry(cell_id))

  #Compute cross-distances: station -> cell_id
  from_coords <- st_coordinates(station$centroid)
  to_coords   <- st_coordinates(cell_id$centroid)

  message("Calculating all relevant receiver-cell distances - this may take a while depending on cellsize")
  distmat <- as.matrix(dist(from_coords, to_coords))

  # Apply max_distance filter
  if (!is.null(max_distance)) {
    distmat[distmat > max_distance] <- NA
  }

  #Flatten distance matrix into data.table
  pairs_idx <- which(!is.na(distmat), arr.ind = TRUE)

  #######
  if("bottom_depth_cell" %in% names(study_area_grid)){
    dt <- data.table(
      station = station$cell_id[pairs_idx[,1]],
      cell_id   = cell_id$cell_id[pairs_idx[,2]],
      dist_to_receiver    = distmat[pairs_idx],
      bottom_depth_cell=cell_id$bottom_depth_cell[pairs_idx[,2]]
    )
  } else {
    dt <- data.table(
      station = station$cell_id[pairs_idx[,1]],
      cell_id   = cell_id$cell_id[pairs_idx[,2]],
      dist_to_receiver    = distmat[pairs_idx]
    )
  }
  #######

  #land barrier
  coast_line <- tryCatch({
    st_cast(study_area, "MULTILINESTRING")
  }, error = function(e) NULL)

  dt[, landinway := FALSE]

  if (!is.null(coast_line) && length(coast_line) > 0 && nrow(dt) > 0) {
    a_list <- station$centroid[pairs_idx[,1]]
    b_list <- cell_id$centroid[pairs_idx[,2]]

    n_total <- length(a_list)
    chunk <- 100000
    p <- progressor(steps = n_total %/% chunk + 1)  # total steps for chunk

    line_list <- vector("list", n_total)

    for (i in seq_len(n_total)) {
      coords_pair <- rbind(st_coordinates(a_list[[i]]), st_coordinates(b_list[[i]]))
      line_list[[i]] <- st_linestring(coords_pair)

      # update progress bar every 'chunk' iterations
      if (i %% chunk == 0 || i == n_total) {
        p(1)  # increment one step
      }
    }

    lines_sfc <- st_sfc(line_list, crs = st_crs(study_area_grid))
    inter_mat <- st_intersects(lines_sfc, coast_line, sparse = FALSE)

    if (is.matrix(inter_mat)) {
      dt[, landinway := apply(inter_mat, 1, any)]
    } else {
      dt[, landinway := as.logical(inter_mat)]
    }
  }
  dist_df_pre<-as.data.frame(dt)
  #Merge relevant study area habitat and grid information. In this example ('receiver-centric'), receivers and their buffer have the same habitat). Below can be altered for 'cell-centric' as needed
  dist_df<-dist_df_pre %>% left_join(study_area_grid %>% select(cell_id,Habitat,centroid),by=c("station"="cell_id")) %>% rename(station_loc=centroid) %>% select(-geometry) #centroid/station_loc relative to station
  dist_df<-dist_df %>% left_join(study_area_grid %>% select(cell_id,centroid,geometry),by=c("cell_id"="cell_id")) #centroid/geometry relative to cell_id grid

  ######
  #if using bathymetry - need to add depth of station (from background study_area grids)
  if("bottom_depth_cell" %in% names(study_area_grid)){
    station_lookup <- unique(st_drop_geometry(dist_df)[, c("station", "station_loc")])
    station_depth_idx <- st_intersects(station_lookup$station_loc,study_area_grid)
    station_lookup$station_depth <- vapply(station_depth_idx,function(x) {
      if(length(x) == 0) NA_real_
      else study_area_grid$bottom_depth_cell[x[1]]
    },numeric(1))

    dist_df <- merge(dist_df,station_lookup[, c("station", "station_depth")],
                     by = "station",
                     all.x = TRUE)
  }

  # ==========================================================
  # FOR DISTANCE APPROACHES
  # ==========================================================

  #if using a relative depth approach (because tag-receiver is not necessarily horizontal) - assumes receiver deployed on bottom
  if (is.null(distance_approach) ||
      is.na(distance_approach) ||
      distance_approach == "horizontal") {
    dist_df$prop_detecs <- predict(range_model,
                                   newdata = dist_df,
                                   type = "response")
    dist_df$prop_detecs <- ifelse(dist_df$landinway, 0, dist_df$prop_detecs)
  } else if (distance_approach == "relative_static") {
    if (is.null(rel_depth_diff)) {
      stop("'rel_depth_diff' must be supplied when distance_approach = 'relative_static'.")
    }
    dist_df$dist_to_receiver <-
      sqrt(dist_df$dist_to_receiver^2 + abs(rel_depth_diff)^2)
    dist_df$prop_detecs <- predict(range_model,
                                   newdata = dist_df,
                                   type = "response")
    dist_df$prop_detecs <- ifelse(dist_df$landinway, 0, dist_df$prop_detecs)
  } else if (distance_approach == "relative_dynamic") {
    if (!("bottom_depth_cell" %in% names(dist_df))) {
      stop("'bottom_depth_cell' is required for distance_approach = 'relative_dynamic'.")
    }
    if (is.null(water_column_position)) {
      stop("'water_column_position' must be supplied when distance_approach = 'relative_dynamic'.")
    }

    dist_df$prop_detecs <- predict(range_model, newdata = dist_df, type = "response")
    dist_df$prop_detecs <- ifelse(dist_df$landinway, 0, dist_df$prop_detecs)
  }
  #########

  return(dist_df)
}
