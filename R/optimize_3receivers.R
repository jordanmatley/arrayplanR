#' @importFrom sf st_area
#' @importFrom sf st_as_sf
#' @importFrom sf st_buffer
#' @importFrom sf st_cast
#' @importFrom sf st_centroid
#' @importFrom sf st_convex_hull
#' @importFrom sf st_coordinates
#' @importFrom sf st_crs
#' @importFrom sf st_distance
#' @importFrom sf st_drop_geometry
#' @importFrom sf st_geometry
#' @importFrom sf st_intersection
#' @importFrom sf st_intersects
#' @importFrom sf st_linestring
#' @importFrom sf st_make_grid
#' @importFrom sf st_point
#' @importFrom sf st_sfc
#' @importFrom sf st_sf
#' @importFrom sf st_transform
#' @importFrom sf st_union
#' @importFrom dplyr filter
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom progressr handlers
#' @importFrom progressr progressor
#' @importFrom progressr with_progress
#' @importFrom progressr without_progress
#' @importFrom proxy dist
#' @import data.table
NULL

#' Optimize receivers (>1-receiver overlap)
#'
#' estimate optimized receivers based on >1 n.receiver_overlap
#'
#' @param study_area_grid An sf object containing the gridded study area
#' @param hab_nrec_alloc An optional identifier for habitat-specific receiver locations
#' @param starting_loc An optional identifier for a forced receiver location
#' @param seed Iteration seed
#' @param rec_dist_iterations Different receiver spacing distances (h iterations) to test for optimization estimates across adjacent n receivers (i iterations)
#' @param cellsize Size of grids/cells used to calculate distances between receivers and cells
#' @param min_rec_distance The minimum distance that receivers should be apart. Only used when optim_method="min_distance"
#' @param max_distance Distance from a receiver which you would expect ~0 (or low) detectability
#' @param study_area2 An sf object containing a reduced area within the study area where receivers can be positioned
#' @param optim_method The method used to optimize receiver array
#' @param square_array_config Selection between square or triangular receiver array configuration
#' @param study_area An sf object containing the study area
#' @param distance_approach An optional identifier that enables three-dimensional distances between an animal and a receiver to be included
#' @param rel_depth_diff The relative depth difference (static - i.e., only one value) in metres between animals and receiver locations. Only used when distance_approach="relative_static"
#' @param water_column_position The expected position of an animal in the water column to calculate non-horizontal distances between receiver and animal. Only used when distance_approach="relative_dynamic"
#' @param range_model Detection range model based on range testing
#' @param cell_prob_cutoff Equivalent to a x% chance of a transmission from a grid being heard (by n receivers)
#' @param n.receiver_overlap Set how many receivers to detect transmission from a grid
#' @param nrecs The number of receivers to be deployed
#' @param rec_dist_variation Selection whether receivers are to be spaced at different distances (e.g., combinations of rec_dist_iterations) in addition to the static spacing designated, or consistent spacing
#'
#' @return A data frame of performance metrics across iterations.
#' @export


#######DE ARRAY OPTIMIZE 3-OVERLAP FUNCTIONS#########

###FUNCTION 3/5
#Estimate optimized receiver based on >1 n.receiver_overlap
optimize_3receivers <- function(study_area_grid,
                                hab_nrec_alloc = NULL,
                                starting_loc = NULL,
                                seed,
                                rec_dist_iterations,
                                cellsize,
                                min_rec_distance,
                                max_distance,
                                study_area2,
                                optim_method="all_possible",
                                square_array_config=FALSE,
                                study_area,
                                distance_approach="horizontal",
                                rel_depth_diff=NULL,
                                water_column_position=NULL,
                                range_model,
                                cell_prob_cutoff,
                                n.receiver_overlap=3,
                                nrecs,
                                rec_dist_variation="static") {

  message("First, calculating receiver-cell distances...")

  set.seed(seed)
  options(datatable.verbose = FALSE)
  options(datatable.print.progress = FALSE)

  if(is.null(study_area_grid$Habitat)){study_area_grid$Habitat="no_habitat"}

  dt_unique_grid <-study_area_grid %>% select(cell_id,geometry,Habitat,centroid) %>% rename(station=cell_id,station_loc=centroid)
  m2_possible_val <- sum(as.numeric(st_area(dt_unique_grid %>% filter(!is.na(Habitat))))) #should be static so precompute

  if (!is.null(starting_loc)){
    starting_loc_pre <- st_sfc(st_point(starting_loc), crs = 4326)
    starting_loc_proj <- st_transform(starting_loc_pre, st_crs(study_area))
    starting_xy <- st_coordinates(starting_loc_proj)
  }

  if(is.null(starting_loc)){
    starting_station_number<-NA
    starting_station<-NA
  }


  suppressWarnings({

    rec_distance_tests<-list()
    Hab_Locs_iterations<-list()
    dist_df_iteration<-list()
    station_possible_iteration<-list()
    if(is.null(rec_dist_iterations)){rec_dist_iterations<-c(100,150,200,250)}

    for (h in 1:length(rec_dist_iterations)) { #loop for each rec distance testing

      hexsize=rec_dist_iterations[h] #hex hereafter refers to the overlying iterated receiver grids (as hexagons or squares)

      if(hexsize<cellsize) {message("cell size is too large for this receiver spacing")
        next}
      if(optim_method == "min_distance" & hexsize>min_rec_distance) {message("rec_dist_iterations larger than min_rec_distance")
        next}
      if(max_distance<=cellsize) stop("your grid size (cellsize) is larger than the expected max range (max_distance) - not useful to continue")

      if (!is.null(starting_loc)){
        # Create initial hex grid over bounding box
        hex_locs <- st_make_grid(
          study_area2,
          cellsize = hexsize,
          what = "polygons",
          square = square_array_config
        )
        hex_locs <- st_sf(geometry = hex_locs) %>%
          mutate(hex_centroid=st_centroid(hex_locs))

        # Find the hex whose centroid is closest to desired starting location
        dist_to_start <- st_distance(hex_locs$hex_centroid, starting_loc_proj)
        closest_idx <- which.min(dist_to_start)

        # Compute the offset needed to shift the grid so the closest centroid matches starting_loc
        shift_vector <- starting_xy - st_coordinates(hex_locs$hex_centroid[closest_idx,])

        # Apply shift to all centroids and geometries
        hex_locs$geometry <- st_geometry(hex_locs) + shift_vector
        hex_locs$hex_centroid <- st_centroid(hex_locs)
        hex_locs<-hex_locs$hex_centroid$geometry
        hex_locs <- st_sf(geometry = st_set_crs(hex_locs, st_crs(study_area_grid)))

        # Clip hex grid to original study area
        hex_locs <- hex_locs %>% st_intersection(study_area2)
        hex_locs<-hex_locs %>% st_intersection(study_area) #in case additional points still exist
        hex_locs <- hex_locs %>% select(geometry)
        hex_locs$hex_id <- seq_len(nrow(hex_locs))

        #associate station ID
        distances <- st_distance(hex_locs$geometry, starting_loc_proj)
        closest_idx2 <- which.min(distances)
        starting_station <- hex_locs[closest_idx2, ]
        starting_station_number_pre <- starting_station$hex_id

      } else{

        #create hexagon grid to overlap
        hex_locs <- study_area2 %>%
          st_make_grid(cellsize = hexsize, what = "polygons",square=FALSE) #%>% square (TRUE) or hexagon (FALSE) arrangement?
        hex_locs <- st_as_sf(hex_locs) %>%
          mutate(hex_centroid=st_centroid(hex_locs))
        hex_locs<-hex_locs$hex_centroid %>%
          st_intersection(study_area2)
        hex_locs<-hex_locs %>%
          st_intersection(study_area) #in case additional points still exist
        hex_locs <- st_sf(geometry = hex_locs) #%>% rename(hex_geometry=geometry)
        hex_locs$hex_id <- seq_len(nrow(hex_locs))
      }

      hex_stations_all <- unique(hex_locs$hex_id)
      if (nrecs > length(hex_stations_all)) stop("nrecs greater than available stations")

      hex_locs_i<-hex_locs %>% st_intersection(dt_unique_grid %>% select(!c(station_loc,station)))
      hex_locs_i<-hex_locs_i %>% filter(!is.na(Habitat)) #be aware with this approach (receiver-centric) can have diff hex_id's in same station ID (which is from background grid cell_id)

      if (!is.null(starting_loc)) {
        starting_station_number_pre<-hex_locs_i %>% filter(hex_id %in% starting_station_number_pre)
        starting_station_number<-starting_station_number_pre$hex_id #hex_id fix

      }

      study_area_buffer<-st_buffer(hex_locs_i,max_distance)

      #use distance function here for hex_locs only...
      dist_df <- compute_grid_distances3(
        station_df = hex_locs_i,
        study_area_buffer = study_area_buffer,
        study_area_grid = study_area_grid,
        study_area = study_area
      )
      #hex_id = station ID (as hex)
      #station = background grid cell id/station (habitat identifier)
      #cell_id = background grid cell id (inlcudin spillover)
      #distance_to_receiver = distance between hex_id and cell_id (?)
      #centroid = cell_id centroid
      #geometry = hex_id
      #station_loc = hex_id location

      dist_df<-dist_df %>% rename(station=hex_id) #to match up with later and because these are stations really
      message("h = ", h, ": computing station distances for ", nrow(hex_locs_i), " possible stations...please be patient (or reduce study_area/cellsize)")


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

      if (is.null(distance_approach) || is.na(distance_approach)) {
        distance_approach <- "horizontal"
      }
      if (distance_approach == "horizontal") {
        # No adjustment required
      } else if (distance_approach == "relative_static") {
        if (is.null(rel_depth_diff)) {
          stop("'rel_depth_diff' must be supplied when distance_approach = 'relative_static'.")
        }
        dist_df$dist_to_receiver <-
          sqrt(dist_df$dist_to_receiver^2 + abs(rel_depth_diff)^2)
      } else if (distance_approach == "relative_dynamic") {
        if (!("bottom_depth_cell" %in% names(dist_df))) {
          stop("'bottom_depth_cell' is required when distance_approach = 'relative_dynamic'.")
        }
        if (is.null(water_column_position)) {
          stop("'water_column_position' must be supplied when distance_approach = 'relative_dynamic'.")
        }
        if (water_column_position == "surface") {
          dist_df$dist_to_receiver <-
            sqrt(dist_df$dist_to_receiver^2 +
                   abs(dist_df$station_depth)^2)
        } else if (water_column_position == "midwater") {
          dist_df$dist_to_receiver <-
            sqrt(dist_df$dist_to_receiver^2 +
                   abs(dist_df$station_depth -
                         (dist_df$bottom_depth_cell / 2))^2)
        } else if (water_column_position == "bottom") {
          dist_df$dist_to_receiver <-
            sqrt(dist_df$dist_to_receiver^2 +
                   abs(dist_df$station_depth -
                         dist_df$bottom_depth_cell)^2)
        } else {
          stop("'water_column_position' must be one of 'surface', 'midwater', or 'bottom'.")
        }
      } else {
        stop("'distance_approach' must be one of 'horizontal', 'relative_static', or 'relative_dynamic'.")
      }

      #########

      #calculate whole array for each iteration to estimate average coverage to compare across
      dist_df$prop_detecs <- predict(range_model, newdata = dist_df, type = "response")

      if(is.null(dist_df$Habitat)){dist_df$Habitat="no_habitat"}

      #set probability to 0 for receiver-grid combinations that have a land barrier
      dist_df$prop_detecs[dist_df$landinway == TRUE] <- 0

      dist_df_iteration[[h]]<-dist_df
      dt <- as.data.table(dist_df)
      dt_unique <- unique(dt[, .(station, station_loc, Habitat)], by = "station")

      stations_all<-dt_unique$station
      dt_unique_sf<-st_sf(dt_unique)
      station_possible_iteration[[h]]<-dt_unique_sf

      setkey(dt, station, cell_id)

      #just taking one receiver at a time, then selecting nearest n receivers to do the calcs
      loop_i_blank <- list()
      Hab_loop_i_blank<-list()
      any_valid_i <- FALSE  # track if any i iteration succeeded

      for (i in seq_along(stations_all)) {
        station_this_i <- stations_all[[i]] #i=101

        # Compute distances from station i to all other stations
        distances <- st_distance(dt_unique_sf[i, ], dt_unique_sf[-i, ])

        # Convert distances to numeric vector
        distances <- as.numeric(distances)
        dist_ranks <- rank(distances, ties.method = "random")

        # All candidate stations (excluding itself)
        candidate_stations <- dt_unique_sf$station[-i]

        # Barrier check
        station_coords1 <- matrix(st_coordinates(st_geometry(dt_unique_sf[i, ])), ncol = 2)
        centroid_coords1 <- st_coordinates(dt_unique_sf[-i, ])[,1:2, drop = FALSE]
        target_crs <- st_crs(study_area)

        # vectorized LINESTRING creation using lapply
        station_cell_line1 <- st_sfc(
          lapply(seq_len(nrow(centroid_coords1)), function(j) {
            st_linestring(rbind(station_coords1, centroid_coords1[j, ]))
          }),
          crs = target_crs
        )
        # Pre-cast study area to multilinestring for land barrier tests
        study_area_line <- st_cast(study_area, "MULTILINESTRING")

        landinway1 <- lengths(st_intersects(station_cell_line1, study_area_line)) > 0

        #Keep only distance-ranked stations that also have a clear line of sight
        nearest_idx <- which( dist_ranks <= (nrecs - 1) & landinway1==FALSE)


        # If all blocked, skip this iteration
        if (length(nearest_idx) == 0) next

        nearest_points <- dt_unique_sf[-i, ][nearest_idx, ]
        near_stations_i <- rbind(dt_unique_sf[i, ], nearest_points)
        new_stations <- unique(near_stations_i$station)

        #Habitat constraint
        if (optim_method == "habitat_specific") {
          if (!identical(sort(near_stations_i$Habitat), sort(hab_nrec_alloc))) {
            next  # skip this i
          }
        }

        #Compute detection probabilities
        dt_sub <- dt[, .(station, cell_id, prop_detecs, Habitat, centroid)]
        dt_sel <- dt_sub[station %in% new_stations]
        #head(dt_sel)

        dt_prop_list <- dt_sel[, .(prop_vec = list(prop_detecs)), by = cell_id]
        dt_prop_list[, prob_cell_all := vapply(prop_vec, function(v)
          compute_cell_prob_optim(v, n.receiver_overlap), numeric(1))]
        dt_prop_list[is.na(prob_cell_all), prob_cell_all := 0]
        dt_prob <- dt_prop_list[, .(cell_id, prob_cell_all)]

        #get convex hull of receivers
        rec_convexhull <- st_convex_hull(st_union(near_stations_i))
        rec_convexhull_buf <- st_buffer(rec_convexhull, dist = min(rec_dist_iterations))
        # Get only the grid cells in the convex hull
        cells_in_hull <- dt_unique_grid[st_intersects(dt_unique_grid, rec_convexhull_buf, sparse = FALSE), ]
        # Filter dt_prob to keep only those cell_ids because we want to optimize such that receivers are pushed outwards (and not clumped together which technically has large coverage - but is not valuable for triangulation)
        dt_prob_convex <- dt_prob[cell_id %in% cells_in_hull$station]

        #to skip this iteration if starting location is set and the station of interest is not being used in this iteration
        if (!is.null(starting_loc) && !(starting_station_number %in% new_stations)) {
          next  # skip this i
        }

        dt_prob_summ <- dt_prob[prob_cell_all >= cell_prob_cutoff,
                                .(ncells = .N,
                                  mean_prob = mean(prob_cell_all),
                                  sum_prob_cells = sum(prob_cell_all),
                                  total_area_m2 = .N * cellsize^2,
                                  m2_possible = m2_possible_val
                                )]

        dt_prob_convex_summ<-dt_prob_convex[prob_cell_all >= cell_prob_cutoff,
                                            .(convex_sum_prob_cells = sum(prob_cell_all),
                                              convex_ncells = .N,
                                              convex_total_area_m2 = .N * cellsize^2)]

        dt_prob_summ <- cbind(dt_prob_summ, dt_prob_convex_summ)

        if (nrow(dt_prob_summ) == 0) next

        #Add metadata
        dt_prob_summ[, `:=`(
          nrecs = nrecs,
          cell_prob_cutoff = cell_prob_cutoff,
          cellsize = cellsize,
          max_distance = max_distance,
          n_receiver_overlap = n.receiver_overlap,
          optim_method = optim_method,
          min_receiver_spacing = ifelse(optim_method == "min_distance",
                                        min_rec_distance, NA),
          rec_distance_iteration = hexsize,
          h_iteration = h,
          i_iteration = i,
          starting_loc_station=starting_station_number,
          n_habitats=n_distinct(near_stations_i$Habitat),
          focal_station=station_this_i,
          stations = list(new_stations),
          config_id=NA,
          status=ifelse(length(new_stations)!=nrecs,"nrecs not matched",NA),
          message=NA
        )]
        loop_i_blank[[i]] <- dt_prob_summ

        Hab_keep_track<-near_stations_i %>% mutate(rec_distance_iteration = hexsize,focal_station=station_this_i)
        Hab_keep_track$h_iteration<-h
        Hab_keep_track$i_iteration<-i
        Hab_loop_i_blank[[i]]<-Hab_keep_track

        any_valid_i <- TRUE  # mark at least one i iteration successful
      } # end i loop

      valid_i <- which(!vapply(loop_i_blank, is.null, logical(1)))

      if (length(valid_i) > 0) {
        loop_i_to_subset <- rbindlist(loop_i_blank[valid_i])
        Hab_i_to_subset <- rbindlist(Hab_loop_i_blank[valid_i])
        rec_distance_tests[[h]] <- loop_i_to_subset
        Hab_Locs_iterations[[h]] <- Hab_i_to_subset
        message(paste0("Completed h = ", h, " with ", length(valid_i), " valid i iterations."))
      } else {
        message(paste0("All i iterations failed for h = ", h, ". Are your requirements realistic? Moving to next h."))
        rec_distance_tests[[h]] <- data.table(
          h_iteration = h,
          status = "no_valid_i",
          message = paste("All i iterations invalid for h =", h)
        )
        Hab_Locs_iterations[[h]] <- data.table(
          h_iteration = h,
          status = "no_valid_i",
          message = paste("All i iterations invalid for h =", h)
        )
        next  # continue to next h iteration (NOT break)
      }

      used_iter <- which(!sapply(loop_i_blank, is.null))
      loop_i_to_subset <- rbindlist(loop_i_blank[used_iter])
      Hab_i_to_subset <- rbindlist(Hab_loop_i_blank[used_iter])
      Hab_i_to_subset$hi_iteration<-paste(Hab_i_to_subset$h_iteration,Hab_i_to_subset$i_iteration,sep="-")

      message(ifelse(length(used_iter) < length(stations_all),
                     paste0("Last possible i solution reached in ", length(used_iter), " iterations (max ", length(stations_all), ")."),
                     paste0("Completed ", length(stations_all), " i iterations.")))
      rec_distance_tests[[h]]<-loop_i_to_subset
      Hab_Locs_iterations[[h]]<-Hab_i_to_subset

    } #end of h loop
  }) #warning suppression

  used_iter2 <- which(!sapply(rec_distance_tests, is.null))

  if (length(used_iter2) > 0) {
    loop_hi_to_subset <- rbindlist(rec_distance_tests[used_iter2], fill = TRUE)
    Hab_hi_to_subset <- rbindlist(Hab_Locs_iterations[used_iter2], fill = TRUE)

    message(ifelse(length(used_iter2) < length(rec_dist_iterations),
                   paste0("Last possible h solution reached in ", length(used_iter2),
                          " iterations (max ", length(rec_dist_iterations), ")."),
                   paste0("Completed ", length(rec_dist_iterations), " h iterations.")))
  } else {
    loop_hi_to_subset <- data.table(
      h_iteration = NA,
      status = "no_valid_h",
      message = "All h iterations failed"
    )
    message("No valid h iterations overall.")
    Hab_hi_to_subset <- data.table(
      h_iteration = NA,
      status = "no_valid_h",
      message = "All h iterations failed"
    )
  }
  loop_hi_to_subset$hi_iteration<-paste(loop_hi_to_subset$h_iteration,loop_hi_to_subset$i_iteration,sep="-")
  loop_hi_to_subset$k_iteration<-NA

  ###################################################################
  ###main block if multiple receiver distances option (rec_dist_variation="variable") is selected
  if (rec_dist_variation == "variable") {

    #A k iteration is a single row from the main summary of all configuration performances (loop_hi_to_subet). A k iteration goes through each one and works through j iterations, which are different mixed receiver distance combinations to be tested. For example, if there are 9 receivers across 3 habitats (3 receivers in each), each j iteration will freeze locations of receivers in 2 habitats and cycle through new receiver distance combinations (selecting new locations in different sized grids that are closest to the ones they are replacing - to ensure they remain in the general area). To ensure stability as much as possible (i.e., keeping optimal spread angles), only one habitat is cycled at a time (instead of two in the above example) - additional investigations are suggested with the interactive web application.

    #loop each loop_hi_subset
    loop_hi_mainHAB<-list()
    array_size_configs_blank<-list()

    loop_hi_mainHAB_presplit <- loop_hi_to_subset[!(status %in% c("no_valid_i", "no_valid_h", "nrecs not matched")) ]
    if(nrow(loop_hi_mainHAB_presplit)>0){ #statement in case there are no valid hi_iterations
      loop_hi_mainHAB_split<-split(loop_hi_mainHAB_presplit,loop_hi_mainHAB_presplit$hi_iteration)
      # iterate configurations

      with_progress({
        p <- progressor(steps = length(loop_hi_mainHAB_split))

        message("Now running through mixed rec_dist_iterations combinations where different habitats may warrant mixed spacing...processing time will increase if this is set to TRUE. Only rec_dist_iterations that are valid/processed (e.g., based on habitat specification, minimum distance cutoffs, etc.) are considered further")

        for (k in 1:length(loop_hi_mainHAB_split)) {

          p(sprintf("k = %d / %d", k, length(loop_hi_mainHAB_split)))  # update FIRST

          loop_k_array <- loop_hi_mainHAB_split[[k]] #k=4

          if (is.na(loop_k_array$n_habitats[1]) || loop_k_array$n_habitats[1] <= 1) next #not needed if only 1 habitat or NA

          #access rec locations

          main_hab <- Hab_hi_to_subset[
            hi_iteration == loop_k_array$hi_iteration[1],
            !c("h_iteration", "i_iteration"),
            with = FALSE
          ]


          #find nearest focal receiver in other iterations
          tomatch_habs <- Hab_hi_to_subset[
            rec_distance_iteration != main_hab$rec_distance_iteration[1] &
              station == focal_station
          ]
          if(nrow(tomatch_habs)==0) next

          main_hab_focal<-main_hab[station==focal_station]

          tomatch_habs$dist_to_focal_m <- as.numeric(st_distance(tomatch_habs$station_loc, main_hab_focal$station_loc))

          tomatch_habs <- tomatch_habs[
            , .SD[which.min(dist_to_focal_m)], by = rec_distance_iteration
          ][
            , ite_focal := paste(rec_distance_iteration, focal_station, sep = "-")
          ]

          secondary_hab <- Hab_hi_to_subset[
            paste(rec_distance_iteration, focal_station, sep = "-") %in% tomatch_habs$ite_focal,
            !c("h_iteration", "i_iteration"),
            with = FALSE
          ]

          # Convert geometry column to sf if needed
          main_sf <- st_as_sf(main_hab, sf_column_name = "station_loc")
          sec_sf  <- st_as_sf(secondary_hab, sf_column_name = "station_loc")

          nearest_by_distance <- rbindlist(
            lapply(unique(secondary_hab$rec_distance_iteration), function(d) {

              sec_sub <- sec_sf[sec_sf$rec_distance_iteration == d, ]

              nearest_idx <- integer(nrow(main_sf))

              for (i in seq_len(nrow(main_sf))) {

                dists <- st_distance(main_sf[i, ], sec_sub)[1, ]

                if (optim_method == "min_distance") {
                  # apply minimum distance constraint
                  dists <- as.numeric(dists)
                  min_rec_distance <- as.numeric(min_rec_distance)
                  valid <- which(dists >= min_rec_distance)

                  if (length(valid) == 0) {
                    nearest_idx[i] <- NA_integer_
                  } else {
                    nearest_idx[i] <- valid[which.min(dists[valid])]
                  }

                } else {
                  # standard nearest neighbour (no constraint)
                  nearest_idx[i] <- which.min(dists)
                }
              }

              data.table(
                station = main_sf$station,
                main_habitat = main_sf$Habitat,
                main_station_loc = st_as_text(main_sf$station_loc),

                matched_distance_iteration = d,

                nearest_secondary_station = sec_sub$station[nearest_idx],
                nearest_secondary_habitat = sec_sub$Habitat[nearest_idx],
                secondary_station_loc = st_as_text(sec_sub$station_loc[nearest_idx]),

                nearest_distance = as.numeric(
                  mapply(function(i, j) {
                    if (is.na(j)) return(NA_real_)
                    st_distance(main_sf[i, ], sec_sub[j, ], by_element = TRUE)
                  }, seq_len(nrow(main_sf)), nearest_idx)
                )
              )
            })
          )

          #in case there are NAs in the nearest_distance - e.g., due to min_rec_dist
          nearest_by_distance <- nearest_by_distance[
            , if (!any(is.na(nearest_distance))) .SD,
            by = matched_distance_iteration
          ]
          if(nrow(nearest_by_distance)==0) next

          # Base habitat configuration (fixed receiver locations)
          # base configuration (fixed original array)
          base_config <- unique(
            nearest_by_distance[, .(
              station,
              main_station_loc,
              main_habitat
            )]
          )

          hab_types <- unique(base_config$main_habitat)

          array_size_configs <- rbindlist(
            lapply(hab_types, function(focal_hab) {

              rbindlist(
                lapply(unique(nearest_by_distance$matched_distance_iteration), function(d) {

                  cfg <- copy(base_config)

                  # pull full replacement table INCLUDING location
                  replacements <- nearest_by_distance[
                    matched_distance_iteration == d &
                      main_habitat == focal_hab,
                    .(
                      station,
                      new_habitat = nearest_secondary_habitat,
                      new_loc     = secondary_station_loc,
                      new_station = nearest_secondary_station
                    )
                  ]

                  # apply habitat + location swap
                  cfg[
                    replacements,
                    `:=`(
                      main_habitat = i.new_habitat,
                      main_station_loc = i.new_loc,
                      station = i.new_station
                    ),
                    on = "station"
                  ]

                  cfg[, focal_habitat := focal_hab]
                  cfg[, matched_distance_iteration := d]

                  cfg
                })
              )
            }),
            fill = TRUE
          )

          array_size_configs[, config_id := rep(
            1:(.N / nrow(base_config)),
            each = nrow(base_config)
          )]

          #head(nearest_by_distance) #contains already quantified main receiver config
          #head(array_size_configs) #new configurations attempts

          #need to filter out options if n_distinct locations are less than nrecs wanted (happens because may duplicate closest receiver to the main ones)
          array_size_configs <- array_size_configs[
            , if (!anyDuplicated(station)) .SD,
            by = config_id
          ]
          if (nrow(array_size_configs) == 0) next

          # ---- enforce nrecs ----
          array_size_configs <- array_size_configs[
            , if (uniqueN(station) == nrecs) .SD,
            by = config_id
          ]
          if (nrow(array_size_configs) == 0) next

          # ---- starting loc constraint ----
          if (!is.null(starting_loc)) {
            array_size_configs <- array_size_configs[
              , if (any(station == starting_station_number)) .SD,
              by = config_id
            ]
            if (nrow(array_size_configs) == 0) next
          }

          # ---- habitat constraint ----
          if (optim_method == "habitat_specific") {
            array_size_configs <- array_size_configs[
              , if (identical(sort(main_habitat), sort(hab_nrec_alloc))) .SD,
              by = config_id
            ]
            if (nrow(array_size_configs) == 0) next
          }

          array_size_configs <- st_as_sf(
            array_size_configs[
              , !c("focal_habitat", "matched_distance_iteration")
            ],
            wkt = "main_station_loc",
            crs = st_crs(study_area)
          )

          setnames(
            array_size_configs,
            old = c("main_habitat","main_station_loc"),
            new = c("Habitat","geometry")
          )
          st_geometry(array_size_configs) <- "geometry"


          # ---- min distance constraint ----
          if (optim_method == "min_distance") {
            array_size_configs <- as.data.table(array_size_configs)[
              , {
                d <- as.matrix(dist(st_coordinates(main_station_loc)))
                if (min(d[upper.tri(d)]) >= min_rec_distance) .SD else NULL
              },
              by = config_id
            ]
            if (nrow(array_size_configs) == 0) next
          }

          #calculating main metrics for later add-on to main output
          config_split <- split(array_size_configs,array_size_configs$config_id) #slowest bit

          process_config <- function(loop_j_array) { #loop_j_array is config_split

            loop_j_array <- copy(loop_j_array) #to avoid misnaming across iterations

            #Habitat constraint
            if (optim_method == "habitat_specific") {
              if (!identical(sort(loop_j_array$Habitat), sort(hab_nrec_alloc))) {
                return(NULL)
              }
            }
            setnames(loop_j_array, "station", "hex_id") #converting so it matches input for first time function run....

            study_area_buffer1<-st_buffer(loop_j_array,max_distance)

            loop_j_array_df <- st_drop_geometry(loop_j_array)
            #dist_df_hab <- progressr::without_progress({compute_grid_distances3(
            dist_df_hab <- local({
              old <- getOption("progressr.enable")
              options(progressr.enable = FALSE)
              on.exit(options(progressr.enable = old), add = TRUE)

              compute_grid_distances3(
                station_df = loop_j_array,
                study_area_buffer = study_area_buffer1,
                study_area_grid,
                study_area)
            })

            setnames(dist_df_hab, "hex_id", "station") #to match up with later and because these are stations really


            ##############
            #if using bathymetry - need to add depth of station (from background study_area grids)
            if("bottom_depth_cell" %in% names(study_area_grid)){
              station_lookup_hab <- unique(st_drop_geometry(dist_df_hab)[, c("station", "station_loc")])
              station_depth_idx_hab <- st_intersects(station_lookup_hab$station_loc,study_area_grid)
              station_lookup_hab$station_depth <- vapply(station_depth_idx_hab,function(x) {
                if(length(x) == 0) NA_real_
                else study_area_grid$bottom_depth_cell[x[1]]
              },numeric(1))
              dist_df_hab <- merge(dist_df_hab,station_lookup_hab[, c("station", "station_depth")],
                                   by = "station",
                                   all.x = TRUE)
            }


            # ==========================================================
            # FOR DISTANCE APPROACHES
            # ==========================================================
            #if using a relative depth approach (because tag-receiver is not necessarily horizontal) - assumes receiver deployed on bottom
            #add TO OTHER FUNCTIONS AND DEFINITIONS
            if (is.null(distance_approach) || is.na(distance_approach)) {
              distance_approach <- "horizontal"
            }
            if (distance_approach == "horizontal") {
              # No adjustment required
            } else if (distance_approach == "relative_static") {
              if (is.null(rel_depth_diff)) {
                stop("'rel_depth_diff' must be supplied when distance_approach = 'relative_static'.")
              }
              dist_df$dist_to_receiver <-
                sqrt(dist_df$dist_to_receiver^2 + abs(rel_depth_diff)^2)
            } else if (distance_approach == "relative_dynamic") {
              if (!("bottom_depth_cell" %in% names(dist_df))) {
                stop("'bottom_depth_cell' is required when distance_approach = 'relative_dynamic'.")
              }
              if (is.null(water_column_position)) {
                stop("'water_column_position' must be supplied when distance_approach = 'relative_dynamic'.")
              }
              if (water_column_position == "surface") {
                dist_df$dist_to_receiver <-
                  sqrt(dist_df$dist_to_receiver^2 +
                         abs(dist_df$station_depth)^2)
              } else if (water_column_position == "midwater") {
                dist_df$dist_to_receiver <-
                  sqrt(dist_df$dist_to_receiver^2 +
                         abs(dist_df$station_depth -
                               (dist_df$bottom_depth_cell / 2))^2)
              } else if (water_column_position == "bottom") {
                dist_df$dist_to_receiver <-
                  sqrt(dist_df$dist_to_receiver^2 +
                         abs(dist_df$station_depth -
                               dist_df$bottom_depth_cell)^2)
              } else {
                stop("'water_column_position' must be one of 'surface', 'midwater', or 'bottom'.")
              }
            } else {
              stop("'distance_approach' must be one of 'horizontal', 'relative_static', or 'relative_dynamic'.")
            }

            ####################


            dist_df_hab$prop_detecs <- predict(range_model, newdata = dist_df_hab, type = "response")

            if (is.null(dist_df_hab$Habitat)) dist_df_hab$Habitat <- "no_habitat"
            dist_df_hab$prop_detecs[dist_df_hab$landinway == TRUE] <- 0
            dt_hab <- as.data.table(dist_df_hab)
            dt_hab_sel <- dt_hab[
              station %in% loop_j_array$hex_id,
              .(station, cell_id, prop_detecs, Habitat, centroid)
            ]
            dt_hab_prop_list <- dt_hab_sel[
              ,
              .(prop_vec = list(prop_detecs)),
              by = cell_id
            ]
            dt_hab_prop_list[
              ,
              prob_cell_all := vapply(
                prop_vec,
                compute_cell_prob_optim,
                numeric(1),
                n.receiver_overlap
              )
            ]
            dt_hab_prop_list[is.na(prob_cell_all), prob_cell_all := 0]
            dt_hab_prob <- dt_hab_prop_list[, .(cell_id, prob_cell_all)]

            # convex hull
            rec_hab_convexhull <- st_convex_hull(st_union(loop_j_array))
            rec_hab_convexhull_buf <- st_buffer(rec_hab_convexhull, dist = min(rec_dist_iterations))
            cells_hab_in_hull <- dt_unique_grid[
              st_intersects(dt_unique_grid, rec_hab_convexhull_buf, sparse = FALSE),
            ]
            dt_hab_prob_convex <- dt_hab_prob[
              cell_id %in% cells_hab_in_hull$station
            ]
            dt_hab_prob_summ <- dt_hab_prob[
              prob_cell_all >= cell_prob_cutoff,
              .(
                ncells = .N,
                mean_prob = mean(prob_cell_all),
                sum_prob_cells = sum(prob_cell_all),
                total_area_m2 = .N * cellsize^2,
                m2_possible = m2_possible_val
              )
            ]
            dt_hab_prob_convex_summ <- dt_hab_prob_convex[
              prob_cell_all >= cell_prob_cutoff,
              .(
                convex_sum_prob_cells = sum(prob_cell_all),
                convex_ncells = .N,
                convex_total_area_m2 = .N * cellsize^2
              )
            ]
            if (nrow(dt_hab_prob_summ) == 0) return(NULL)
            dt_hab_prob_summ <- cbind(dt_hab_prob_summ, dt_hab_prob_convex_summ)
            dt_hab_prob_summ[
              ,
              `:=`(
                nrecs = nrecs,
                cell_prob_cutoff = cell_prob_cutoff,
                cellsize = cellsize,
                max_distance = max_distance,
                n_receiver_overlap = n.receiver_overlap,
                optim_method = optim_method,
                min_receiver_spacing = ifelse(optim_method == "min_distance", min_rec_distance, NA),
                rec_distance_iteration = hexsize,
                h_iteration = loop_k_array$h_iteration,
                i_iteration = loop_k_array$i_iteration,
                starting_loc_station = starting_station_number,
                n_habitats = uniqueN(loop_j_array$Habitat),
                focal_station = unique(main_hab_focal$station),
                stations = list(unique(loop_j_array$hex_id)),
                config_id=unique(loop_j_array_df$config_id),
                status="mixed_distances_iteration",
                message=NA,
                hi_iteration=unique(loop_k_array$hi_iteration),
                k_iteration=k
              )
            ]
            return(dt_hab_prob_summ)
          }

          array_configs_list <- without_progress(lapply(config_split, process_config))

          array_configs_df <- rbindlist(array_configs_list, fill = TRUE)

          #Need to merge all this with main loop file
          loop_hi_mainHAB[[k]]<-array_configs_df
          array_size_configs_blank[[k]] <- array_size_configs #store for later plotting

        } #end k loop

      }) #progress bar

      loop_hi_HAB_tojoin<-do.call(rbind,loop_hi_mainHAB)

      if(is.null(loop_hi_HAB_tojoin)) {
        loop_hi_to_subset<-loop_hi_to_subset} else{
          loop_hi_to_subset<-rbind(loop_hi_to_subset,loop_hi_HAB_tojoin)
        }

      #Return outputs
      list(
        loop_hi_to_subset = loop_hi_to_subset,
        array_size_configs_all = array_size_configs_blank, #will be empty if loop_hi_HAB_tojoin is null
        rec_distance_tests = rec_distance_tests,
        dist_df_iteration = dist_df_iteration,
        station_possible_iteration = station_possible_iteration
      )


      ###################################################################
      ###main block if multiple receiver distances option (rec_dist_variation="variable") is selected (above)
    }
    else{ #statement in case there are no valid hi_iterations
      list(
        loop_hi_to_subset = loop_hi_to_subset,
        rec_distance_tests = rec_distance_tests,
        dist_df_iteration = dist_df_iteration,
        station_possible_iteration = station_possible_iteration
      )
      if(all(as.character(loop_hi_to_subset$status)=="no_valid_i")) message("No valid i iterations across all h iterations. Are your parameter settings reasonable?")
    }
  } else {
    #Return outputs
    list(
      loop_hi_to_subset = loop_hi_to_subset,
      rec_distance_tests = rec_distance_tests,
      dist_df_iteration = dist_df_iteration,
      station_possible_iteration = station_possible_iteration
    )
  }
} #end of function
