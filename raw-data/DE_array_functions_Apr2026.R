
#######DE ARRAY CALCULATE FUNCTIONS#########

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

###FUNCTION 2/2
#calculate receiver-cell probabilities based on designated n.receiver_overlap
compute_cell_prob <- function(dist_df, 
                              study_area_grid, 
                              n.receiver_overlap, 
                              progress = TRUE) {

  #initial checks
  stopifnot("cell_id" %in% names(dist_df))
  stopifnot("prop_detecs" %in% names(dist_df))
  if (!inherits(study_area_grid, "sf")) stop("study_area_grid must be an sf object with geometry and cell_id.")
  stopifnot("cell_id" %in% names(study_area_grid))
  if (!is.numeric(n.receiver_overlap) || n.receiver_overlap < 1) stop("n.receiver_overlap must be integer >= 1")
  
  cell_split <- split(dist_df, dist_df$cell_id) #split cells
  cell_ids <- names(cell_split)
  
  if (progress) {
    handlers("txtprogressbar")
    p <- progressor(along = cell_ids)
  } else {
    p <- function(...) NULL
  }
  
  #internal probability calculation function
  compute_for_one <- function(cell_df, thr) {
    cid <- unique(cell_df$cell_id)
    n_cell <- nrow(cell_df)
    if (n_cell == 0 || n_cell < thr) {
      return(tibble(cell_id = cid, prob_cell_all = 0))
    }
    p_vec <- as.numeric(cell_df$prop_detecs)
    if (thr == 1) {
      prob_cell_all <- 1 - prod(1 - p_vec)
      return(tibble(cell_id = cid, prob_cell_all = prob_cell_all))
    }
    dp <- numeric(thr)
    dp[1] <- 1
    for (p_i in p_vec) {
      new_dp <- numeric(thr)
      new_dp[1] <- dp[1] * (1 - p_i)
      if (thr > 1) {
        for (k in 2:thr) {
          new_dp[k] <- dp[k] * (1 - p_i) + dp[k - 1] * p_i
        }
      }
      dp <- new_dp
    }
    tibble(cell_id = cid, prob_cell_all = 1 - sum(dp))
  }
  
  #apply internal function across cells
  out_list <- vector("list", length(cell_ids))
  for (i in seq_along(cell_ids)) {
    if (progress) p()
    out_list[[i]] <- compute_for_one(cell_split[[i]], n.receiver_overlap)
  }
  dist_df_prob <- bind_rows(out_list)
  
  #join back cells
  study_grid_unique <- study_area_grid %>%
    distinct(cell_id, .keep_all = TRUE)
  dist_df_prob_joined <- study_grid_unique %>%
    left_join(dist_df_prob, by = "cell_id") %>%
    mutate(prob_cell_all = replace_na(prob_cell_all, 0))
}


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

  #data.table::setDT(dt) #add during package creation
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
    #dist_df <- merge(dist_df,"bottom_depth_cell" %in% names(study_area_grid)[, c("station", "station_depth")],
    #                 by = "station",
    #                 all.x = TRUE)
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

###FUNCTION 2/4
#calculate optimized receiver-cell probabilities based on designated n.receiver_overlap
#duplicate
compute_cell_prob_optim <- function(p_vec, thr) {
  if (length(p_vec) == 0) return(0)
  if (thr <= 0) return(1)
  if (thr == 1) return(1 - prod(1 - p_vec))
  if (length(p_vec) < thr) return(0)
  
  N <- length(p_vec)
  L <- 1L
  while (L < (N + 1L)) L <- L * 2L
  
  j <- 0:(L - 1L)
  omega <- exp(-2i * pi * j / L)
  prod_vals <- rep(1 + 0i, L)
  for (p in p_vec) {
    prod_vals <- prod_vals * ((1 - p) + p * omega)
  }
  
  coeffs <- Re(fft(prod_vals, inverse = TRUE) / L)[1:(N + 1L)]
  coeffs[coeffs < 0] <- 0
  sum(coeffs[(thr + 1L):(N + 1L)])
}

###FUNCTION 3/4
#Estimate optimized receiver based on 1 n.receiver_overlap
optimize_1receivers <- function(dist_df, 
                                study_area_grid, 
                                hab_nrec_alloc = NULL,
                                seed,
                                starting_loc=NULL,
                                study_area,
                                optim_method="all_possible",
                                min_rec_distance,
                                max_i_iterations,
                                subsampl_perc_i=100,
                                n.receiver_overlap=1,
                                cell_prob_cutoff,
                                cellsize,
                                max_distance,
                                optim_metric="sum_prob_cells",
                                nrecs) {
  
  #creating temporary folder for potential memory issues (not tested for MacOS)
  #newtemp <- "C:/temp_R"
  #dir.create(newtemp, showWarnings = FALSE)
  #Sys.setenv(TMP = newtemp, TEMP = newtemp)
  
  message(paste("Identifying estimated optimized receiver locations based on: ", max_i_iterations,
                " i iterations. Each i iteration will include ~",
                round(choose(n_distinct(dist_df$station), n.receiver_overlap) * (subsampl_perc_i / 100), 0),
                "possible combinations (at", subsampl_perc_i, "% subsampling rate) from", nrecs, "receiver locations. This may take a while depending on cellsize..."))
  
  set.seed(seed) #iterative so results may change
  options(datatable.verbose = FALSE)
  options(datatable.print.progress = FALSE)
  
  if (!is.null(starting_loc)){ 
    # Project starting location to study area CRS
    starting_loc_pre <- st_sfc(st_point(starting_loc), crs = 4326)
    starting_loc_proj <- st_transform(starting_loc_pre, st_crs(study_area))
    starting_xy <- st_coordinates(starting_loc_proj)
    
    # Find the grid whose centroid is closest to desired starting location
    dist_to_start <- st_distance(study_area_grid[study_area_grid$Description == "study_area", ]$centroid,starting_loc_proj)    
    closest_idx <- which.min(dist_to_start)
    starting_station <- study_area_grid[closest_idx, ]
    starting_station_number <- starting_station$cell_id

    #re-calculate distances for starting station only
    setDT(dist_df)    
    dist_df1 <- dist_df[station == starting_station_number]
    dist_df2 <- dist_df[station != starting_station_number]
    dist_df1$station_loc<-starting_loc_proj
    dist_df1$dist_to_receiver <- round(mapply(function(a, b) st_distance(a, b),dist_df1$station_loc,dist_df1$centroid),1)
    dist_df1$prop_detecs <- predict(range_model, newdata = dist_df1, type = "response")
    dist_df1$prop_detecs <- ifelse(dist_df1$landinway, 0, dist_df1$prop_detecs)
    dist_df<-rbind(dist_df1,dist_df2)
    
  }
  if(is.null(starting_loc)){
    setDT(dist_df)
    starting_station_number<-NA
    starting_station<-NA
  }
  
  if(is.null(dist_df$Habitat)){dist_df$Habitat="no_habitat"}
  
  #subset
  dt<-dist_df
  dt_unique <- unique(dt[, .(station, station_loc, Habitat)], by = "station")
  dt_sub <- dt[, .(station,cell_id, prop_detecs, Habitat, centroid)]
  dt <- NULL
  stations_all <- unique(dt_unique$station)
  if (nrecs > length(stations_all)) stop("nrecs greater than available stations")
  
  #Initialize receiver locations (prior to first iteration)
  if (optim_method == "habitat_specific") {
    if (is.null(hab_nrec_alloc)) stop("Provide hab_nrec_alloc")
    habitat_levels <- unique(dt_unique$Habitat)
    if (length(hab_nrec_alloc) != length(habitat_levels)) stop("Length mismatch")
    names(hab_nrec_alloc) <- habitat_levels
    unique_stations_loop <- dt_unique[, .SD[sample(.N, hab_nrec_alloc[.BY$Habitat])], by = Habitat]$station
    
  } else if (optim_method == "min_distance") {
    find_recs <- st_sf(dt_unique[0, ])
    while (nrow(find_recs) < nrecs && nrow(dt_unique) > 0) {
      find_it <- dt_unique[sample(.N, 1)]
      find_it <- st_as_sf(find_it, sf_column_name = "station_loc", crs = st_crs(find_recs))
      find_recs <- rbind(find_recs, find_it)
      dists <- as.numeric(st_distance(find_it$station_loc, dt_unique$station_loc))
      locs_possible <- dt_unique[dists >= min_rec_distance, ]
    }
    unique_stations_loop <- unique(find_recs$station)
    if (length(unique_stations_loop) != nrecs)
      stop("Desired receiver number not met – likely min_rec_distance too large")
    
  } else {
    unique_stations_loop <- sample(stations_all, nrecs)
  }
  
  # Pre-index data for fast lookups
  setkey(dt_sub, station)
  setindex(dt_sub, cell_id)
  
  # Pre-split dt_sub by station
  prop_dt_list <- split(dt_sub, by = "station", keep.by = TRUE)
  
  # Precompute station centroids table (for min_distance checks)
  station_centroid <- dt_unique[, .(station, station_loc)]
  setnames(station_centroid, "station_loc", "centroid")
  
  setkey(station_centroid, station)
  
  #precompute full station distance matrix (if min_distance used)
  dist_mat <- NULL
  if (optim_method == "min_distance") {
    # convert centroids to numeric matrix
    coords_mat <- st_coordinates(station_centroid$centroid)
    # as.matrix(dist(...)) calculates Euclidean distances; acceptable substitute for st_distance for speed
    dist_mat_full <- as.matrix(dist(coords_mat))
    rownames(dist_mat_full) <- station_centroid$station
    colnames(dist_mat_full) <- station_centroid$station
    dist_mat <- dist_mat_full
  }
  
  # Used station combinations to avoid duplicates
  used_keys <- new.env(parent = emptyenv())
  
  with_progress({
    handlers("cli")
    p <- progressor(steps = max_i_iterations)
    
    loop_i_blank <- vector("list", max_i_iterations)
    combos_i_used <- vector("list", max_i_iterations)
    
    # Main iteration loop (have kept option for n.receiver_overlap>1 but advise against use)
    for (i in seq_len(max_i_iterations)) {
      
      if (!is.null(starting_loc) && i==1){ 
        pre_combos <- unique_stations_loop[-(1:(n.receiver_overlap+1))]
      } else if (!is.null(starting_loc) && i>1){
        pre_combos <- unique_stations_loop[-(1:(n.receiver_overlap))]
        pre_combos <- pre_combos[pre_combos != starting_station_number]
      } else {
        pre_combos <- unique_stations_loop[-(1:n.receiver_overlap)]
      }
      
      # Efficient random combination sampling: keep as a list of sorted integer vectors
      sample_random_combos <- function(stations_all, k, nsamples) {
        replicate(nsamples, sort(sample(stations_all, k)), simplify = FALSE)
      }
      
      n_possible <- choose(length(stations_all), n.receiver_overlap)
      nsamples <- ceiling(n_possible * subsampl_perc_i / 100)
      if (nsamples < 1) nsamples <- 1
      
      combos_k_using_random <- sample_random_combos(stations_all, n.receiver_overlap, nsamples)
      loop_k_blank <- vector("list", length(combos_k_using_random))
      
      # iterate combos
      for (k_idx in seq_along(combos_k_using_random)) {
        loop_k <- combos_k_using_random[[k_idx]]
        
        if (!is.null(starting_loc)){ 
          new_stations <- c(pre_combos, loop_k,starting_station_number)
        } else {
          new_stations <- c(pre_combos, loop_k)
        }
        
        if (length(unique(new_stations)) != nrecs) next
        
        # Fast hashed key based on sorted station vector
        stations_key <- digest(
          paste0(new_stations, collapse = ","),
          algo = "xxhash64",
          serialize = FALSE
        )
        
        if (exists(stations_key, envir = used_keys, inherits = FALSE)) next
        used_keys[[stations_key]] <- TRUE
        
        #Habitat constraints
        if (optim_method == "habitat_specific") {
          hab_removed <- dt_unique[station %in% unique_stations_loop[1:n.receiver_overlap],
                                   .SD[!duplicated(station)], .SDcols = c("station", "Habitat")][, Habitat]
          hab_kept <- dt_unique[station %in% loop_k,
                                .SD[!duplicated(station)], .SDcols = c("station", "Habitat")][, Habitat]
          if (!identical(sort(hab_removed), sort(hab_kept))) next
        }
        
        #Min-distance constraints
        if (optim_method == "min_distance") {
          dmat_sub <- dist_mat[as.character(new_stations), as.character(new_stations), drop = FALSE]
          diag(dmat_sub) <- Inf
          if (as.numeric(min(dmat_sub, na.rm = TRUE)) < min_rec_distance) next
        }
        
        # Compute per-cell detection probabilities
        rows_list <- prop_dt_list[as.character(new_stations)]
        if (any(vapply(rows_list, is.null, logical(1)))) next
        dt_sel <- rbindlist(rows_list, use.names = TRUE, fill = TRUE)
        
        # aggregate prop_detecs per cell
        dt_prop_list <- dt_sel[, .(prop_vec = list(prop_detecs)), by = cell_id]
        
        # compute cell-level prob depending on overlap
        if (n.receiver_overlap == 1) {
          dt_prop_list[, prob_cell_all := vapply(prop_vec, function(v) 1 - prod(1 - v), numeric(1))]
        } else {
          dt_prop_list[, prob_cell_all := vapply(prop_vec, function(v) {
            compute_cell_prob_optim(v, n.receiver_overlap)
          }, numeric(1))]
        }
        
        dt_prop_list[is.na(prob_cell_all), prob_cell_all := 0]
        dt_prob <- dt_prop_list[, .(cell_id, prob_cell_all)]
        
        # Skip if no cells above cutoff
        dt_prob_summ <- dt_prob[prob_cell_all >= cell_prob_cutoff,
                                .(ncells = .N,
                                  mean_prob = mean(prob_cell_all),
                                  sum_prob_cells = sum(prob_cell_all),
                                  total_area_m2 = .N * cellsize^2,
                                  m2_possible = sum(as.numeric(
                                    st_area(study_area_grid[study_area_grid$Description == "study_area", ])
                                  )))]
        
        if (nrow(dt_prob_summ) == 0) next
        
        # Add metadata to table
        dt_prob_summ[, `:=`(
          nrecs = nrecs,
          cell_prob_cutoff = cell_prob_cutoff,
          cellsize = cellsize,
          max_distance = max_distance,
          n_receiver_overlap = n.receiver_overlap,
          optim_method = optim_method,
          min_receiver_spacing = ifelse(optim_method == "min_distance",
                                        min_rec_distance, NA),
          i_iteration = i,
          starting_loc_station=starting_station_number,
          stations = list(new_stations)
        )]
        
        loop_k_blank[[k_idx]] <- dt_prob_summ
      } # end k loop
      
      valid_k <- which(!vapply(loop_k_blank, is.null, logical(1)))
      
      if (length(valid_k) > 0) {
        loop_k_to_subset <- rbindlist(loop_k_blank[valid_k])
        loop_k_to_subset <- loop_k_to_subset[sample(.N)]
        loop_k_subset <- loop_k_to_subset[which.max(get(optim_metric))]
        loop_i_blank[[i]] <- loop_k_subset
        
        unique_stations_loop <- unlist(loop_k_subset$stations)
        combos_i_used[[i]] <- sort(unique_stations_loop)
        
        p(sprintf("Iteration %d / %d complete", i, max_i_iterations))
      } else {
        message("All k combinations invalid in i iteration ", i, ". Increase subsampl_perc_i or try another seed.")
        break
      }
    } # end i loop
  }) # end with_progress
  
  used_iter <- which(!sapply(loop_i_blank, is.null))
  loop_i_to_subset <- if (length(used_iter) > 0) rbindlist(loop_i_blank[used_iter]) else NULL
  
  message(ifelse(length(used_iter) < max_i_iterations,
                 paste0("Last possible i solution reached in ", length(used_iter), " iterations (max ", max_i_iterations, ")."),
                 paste0("Completed ", max_i_iterations, " i iterations.")))
  
  list(
    loop_i_to_subset = loop_i_to_subset,
    loop_i_blank = loop_i_blank,
    combos_i_used = combos_i_used,
    unique_stations_final = unique_stations_loop
  )
}

###FUNCTION 4/4
#Plotting function
#duplicate
plot_output_prep <- function(dist_df, 
                             study_area_grid, 
                             stations_to_use, 
                             n.receiver_overlap=1) {
  
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
  
  #station_split <- split(station_df, station_df$station) # Split stations (technically uses this for second ime function run - so will have to adjust so both attempts run)
  station_split <- split(station_df, station_df$hex_id) # Split hex_id stations [not background cells because they are just reference habitats and distances - including potentially duplicate stations since hex_ids may fall in duplicate station/cell_ids]
  
  #handlers("txtprogressbar")  # progress bar
  result_list <- with_progress({
    handlers("cli")
    p <- progressor(along = station_split)
    
    lapply(station_split, function(station_df_loop) {
      p() 
      #station_i <- station_df_loop$station[1] #is background cell id
      station_i <- station_df_loop$hex_id[1] #is hex_id/proper station id
      station_geom <- station_df_loop$geometry[1] #is geometry of hex_id 
      

      # Subset buffer & intersect with grid
      suppressWarnings(
        study_buffer_loop <- study_area_buffer %>%
          #dplyr::filter(station == station_i) %>%
          filter(hex_id == station_i) %>% #with new hex_id fix
          st_intersection(study_area_grid %>% select(-Name) %>% rename(background_grid_habitat=Habitat))
      )
      #centroid = background cell_id centroid
      #geometry = geometry of used cell_id (may be partially used)

      #ggplot()+
      #  geom_sf(data = study_area,fill="white")+
      #  geom_sf(data=study_buffer_loop,aes(geometry=geometry))+
      #  geom_sf(data=study_buffer_loop,aes(geometry=centroid))
      
      #Okay to use centroids where grid only partially used
      
      
      if (nrow(study_buffer_loop) == 0L) return(NULL)
    
      study_buffer_loop$centroid <- study_area_grid$centroid[
        match(study_buffer_loop$cell_id, study_area_grid$cell_id)
      ]
      centroids_sfc <- st_geometry(study_buffer_loop$centroid)
      
      #ggplot()+
      #  geom_sf(data = study_area,fill="white")+
      #  geom_sf(data=study_buffer_loop,aes(geometry=geometry))+ #cell_id geometry
      #  geom_sf_text(data=study_buffer_loop,aes(geometry=centroid,label=cell_id))+ #cell_id centroid
      #  geom_sf(data=station_geom,colour="red") #station
      
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

###FUNCTION 2/5
#calculate optimized receiver-cell probabilities based on designated n.receiver_overlap
#duplicate
compute_cell_prob_optim <- function(p_vec, thr) {
  if (length(p_vec) == 0) return(0)
  if (thr <= 0) return(1)
  if (thr == 1) return(1 - prod(1 - p_vec))
  if (length(p_vec) < thr) return(0)
  
  N <- length(p_vec)
  L <- 1L
  while (L < (N + 1L)) L <- L * 2L
  
  j <- 0:(L - 1L)
  omega <- exp(-2i * pi * j / L)
  prod_vals <- rep(1 + 0i, L)
  for (p in p_vec) {
    prod_vals <- prod_vals * ((1 - p) + p * omega)
  }
  
  coeffs <- Re(fft(prod_vals, inverse = TRUE) / L)[1:(N + 1L)]
  coeffs[coeffs < 0] <- 0
  sum(coeffs[(thr + 1L):(N + 1L)])
}

###FUNCTION 3/5
#Estimate optimized receiver based on >1 n.receiver_overlap
optimize_3receivers <- function(study_area_grid, 
                                hab_nrec_alloc = NULL,
                                starting_loc = NULL,
                                seed,
                                study_area_buffer,
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
  
  #For temporary storage (not tested for MacOS)
  #newtemp <- "C:/temp_R"
  #dir.create(newtemp, showWarnings = FALSE)
  #Sys.setenv(TMP = newtemp, TEMP = newtemp)
  
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
      
      #this moves point to grid centroid
      #hex_locs_i<-st_join(hex_locs,dt_unique_grid,join = st_intersects,left = TRUE, largest = TRUE) %>% st_set_geometry("centroid") %>% select(-geometry)
      #hex_locs_i<-hex_locs_i %>% filter(!is.na(Habitat))
      #this keeps at specific location (and may skews some distances between grids,especially with larger grid sizes - can also switch study_area_grid_pre in "DE_array_optimize_3overlap.r" square to TRUE/FALSE to improve)
      #hex_locs_i<-hex_locs %>% st_intersection(dt_unique_grid %>% select(-station_loc)) #%>% rename(Habitat_station=Habitat) 
      hex_locs_i<-hex_locs %>% st_intersection(dt_unique_grid %>% select(!c(station_loc,station))) #do not need station i don't think and it confuses further down
      hex_locs_i<-hex_locs_i %>% filter(!is.na(Habitat)) #can have diff hex_id's in same station ID (which is from background grid cell_id - potentially problematic...need to make sure hex_id is actual new stationID because it is unique...should only be using station here in relation to what habitat it is and where it is from...)
      #head(hex_locs_i) #geometry is hex location
      
      if (!is.null(starting_loc)) {
        starting_station_number_pre<-hex_locs_i %>% filter(hex_id %in% starting_station_number_pre)
        #starting_station_number<-starting_station_number_pre$station
        starting_station_number<-starting_station_number_pre$hex_id #hex_id fix
        
      }
      
      study_area_buffer<-st_buffer(hex_locs_i,max_distance)
      
      #use distance function here for hex_locs only...
      #handlers(global = TRUE)
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
        #head(dist_df)
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
        
        #filter out land-blocked stations
        #land_blocks <- dt[.(station_this_i, candidate_stations), landinway] #won't work on 3 version because land in ways are between cell_id not stations/hex_id - so need to redo for stations specifically
        # Barrier check
        station_coords1 <- matrix(st_coordinates(st_geometry(dt_unique_sf[i, ])), ncol = 2)
        #centroid_coords1 <- st_coordinates(centroids_sfc)[,1:2, drop = FALSE]
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
        #candidate_stations1<-data.frame(station=c(station_this_i,candidate_stations),landinway1=c("FALSE",landinway1))
        #head(candidate_stations1) #issue if order important
        #dt<-dt %>% left_join(candidate_stations1,by="station") %>% filter(landinway1==FALSE)
        #head(dt)
        #there is already a landinway from station to cell_id - may need later?
        #unique(dt$landinway1)
        #land_blocks[is.na(land_blocks)] <- TRUE 
        #valid_connection <- !land_blocks
        
        #Keep only distance-ranked stations that also have a clear line of sight
        #nearest_idx <- which( dist_ranks <= (nrecs - 1) & valid_connection)
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
        rec_convexhull_buf <- st_buffer(rec_convexhull, dist = min(rec_dist_iterations))#could use receiver distances for this h-iteration, but no matter what we want receivers to optimiz
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
        Hab_keep_track$h_iteration<-h # i think i need this
        Hab_keep_track$i_iteration<-i
        #head(Hab_keep_track)
        Hab_loop_i_blank[[i]]<-Hab_keep_track
        
        any_valid_i <- TRUE  # mark at least one i iteration successful
      } # end i loop
      
      valid_i <- which(!vapply(loop_i_blank, is.null, logical(1)))
      
      if (length(valid_i) > 0) {
        loop_i_to_subset <- rbindlist(loop_i_blank[valid_i])
        Hab_i_to_subset <- rbindlist(Hab_loop_i_blank[valid_i])
        rec_distance_tests[[h]] <- loop_i_to_subset
        Hab_Locs_iterations[[h]] <- Hab_i_to_subset
        message(paste0("✅ Completed h = ", h, " with ", length(valid_i), " valid i iterations."))
      } else {
        message(paste0("⚠️ All i iterations failed for h = ", h, ". Are your requirements realistic? Moving to next h."))
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
    message("⚠️ No valid h iterations overall.")
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
          #head(loop_k_array)
          
          if (is.na(loop_k_array$n_habitats[1]) || loop_k_array$n_habitats[1] <= 1) next #not needed if only 1 habitat or NA
          
          #access rec locations
          #main_hab<-Hab_hi_to_subset %>% filter(hi_iteration %in% unique(loop_k_array$hi_iteration)) %>% dplyr::select(!c(h_iteration,i_iteration))
          main_hab <- Hab_hi_to_subset[
            hi_iteration == loop_k_array$hi_iteration[1],
            !c("h_iteration", "i_iteration"),
            with = FALSE
          ]
          #main_hab
          
          #find nearest focal receiver in other iterations
          #tomatch_habs<-Hab_hi_to_subset  %>% filter(rec_distance_iteration!=unique(main_hab$rec_distance_iteration)) %>% filter(station==focal_station)
          tomatch_habs <- Hab_hi_to_subset[
            rec_distance_iteration != main_hab$rec_distance_iteration[1] &
              station == focal_station
          ]
          #head(tomatch_habs)
          if(nrow(tomatch_habs)==0) next 
          
          #main_hab_focal<-main_hab %>% filter(station==focal_station)
          main_hab_focal<-main_hab[station==focal_station]
          #main_hab_focal
          
          tomatch_habs$dist_to_focal_m <- as.numeric(st_distance(tomatch_habs$station_loc, main_hab_focal$station_loc))
          #head(tomatch_habs)
          
          #tomatch_habs<-tomatch_habs %>% group_by(rec_distance_iteration) %>% slice_min(dist_to_focal_m,n=1) %>% mutate(ite_focal=paste(rec_distance_iteration,focal_station,sep="-"))
          tomatch_habs <- tomatch_habs[
            , .SD[which.min(dist_to_focal_m)], by = rec_distance_iteration
          ][
            , ite_focal := paste(rec_distance_iteration, focal_station, sep = "-")
          ]
          
          
          #secondary_hab<-Hab_hi_to_subset %>% mutate(ite_focal=paste(rec_distance_iteration,focal_station,sep="-")) %>% filter(ite_focal %in% unique(tomatch_habs$ite_focal)) %>% dplyr::select(!c(h_iteration,i_iteration))
          secondary_hab <- Hab_hi_to_subset[
            paste(rec_distance_iteration, focal_station, sep = "-") %in% tomatch_habs$ite_focal,
            !c("h_iteration", "i_iteration"),
            with = FALSE
          ]
          #head(secondary_hab)
          
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
          
          #nearest_by_distance
          
          #in case there are NAs in the nearest_distance - e.g., due to min_rec_dist
          nearest_by_distance <- nearest_by_distance[
            , if (!any(is.na(nearest_distance))) .SD,
            by = matched_distance_iteration
          ]
          if(nrow(nearest_by_distance)==0) next
          
          #setDT(nearest_by_distance)
          
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
          
          # ---- sf conversion ----
          # array_size_configs <- st_as_sf(
          #    array_size_configs[, !c("focal_habitat", "matched_distance_iteration")],
          #    wkt = "main_station_loc",
          #    crs = st_crs(study_area)
          #  ) %>% rename(geometry=main_station_loc,Habitat=main_habitat) 
          #  #names(array_size_configs)[names(array_size_configs) == "main_habitat"] <- "Habitat"
          
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
          #array_size_configs
          
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
            ####MAY BE AN ISSUE
            #loop_j_array<-loop_j_array %>% rename(hex_id=station) #converting so it matches input for first time function run....
            setnames(loop_j_array, "station", "hex_id") #converting so it matches input for first time function run....
            ####MAY BE AN ISSUE
            
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
            #}) #silence progress bar within function
            #dist_df_hab<-dist_df_hab %>% rename(station=hex_id) #to match up with later and because these are stations really
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
              #head(dist_df_hab)
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
          
          #suppressWarnings({
          array_configs_list <- without_progress(lapply(config_split, process_config))
          #})
          #array_configs_list <- progressr::without_progress({
          #  options(progressr.enable = FALSE)
          #  on.exit(options(progressr.enable = TRUE), add = TRUE)
          #  
          #  lapply(config_split, process_config)
          #})
          
          array_configs_df <- rbindlist(array_configs_list, fill = TRUE)
          #head(array_configs_df)
          
          #Need to merge all this with main loop file
          
          loop_hi_mainHAB[[k]]<-array_configs_df
          #array_size_configs_blank[[k]]<-array_size_configs
          array_size_configs_blank[[k]] <- array_size_configs #store for later plotting
          
          #})
          
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
      if(all(as.character(loop_hi_to_subset$status)=="no_valid_i")) message("⚠️ No valid i iterations across all h iterations. Are your parameter settings reasonable?")  
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

