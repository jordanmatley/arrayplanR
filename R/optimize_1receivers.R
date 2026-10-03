#' @importFrom sf st_area
#' @importFrom sf st_as_sf
#' @importFrom sf st_coordinates
#' @importFrom sf st_crs
#' @importFrom sf st_distance
#' @importFrom sf st_point
#' @importFrom sf st_sfc
#' @importFrom sf st_sf
#' @importFrom sf st_transform
#' @importFrom progressr handlers
#' @importFrom progressr progressor
#' @importFrom progressr with_progress
#' @importFrom proxy dist
#' @importFrom digest digest
#' @import data.table
NULL

#' Optimize receivers (1-receiver overlap)
#'
#' estimate optimized receivers based on 1 n.receiver_overlap
#'
#' @param dist_df An sf receiver-cell distance dataframe (from compute_grid_distances1 function)
#' @param study_area_grid An sf object containing the gridded study area
#' @param hab_nrec_alloc An optional identifier for habitat-specific receiver locations
#' @param seed Iteration seed
#' @param starting_loc An optional identifier for a forced receiver location
#' @param study_area An sf object containing the study area.
#' @param optim_method The method used to optimize receiver array
#' @param min_rec_distance The minimum distance that receivers should be apart. Only used when optim_method="min_distance"
#' @param max_i_iterations Number of i iterations to process as a balance between computational costs and maximizing optimization solution
#' @param subsampl_perc_i The percent of all possible receiver-grid combinations taken (and randomized) in each i iteration
#' @param n.receiver_overlap Set how many receivers to detect transmission from a grid
#' @param cell_prob_cutoff Equivalent to a x% chance of a transmission from a grid being heard (by n receivers)
#' @param cellsize Size of grids/cells used to calculate distances between receivers and cells
#' @param max_distance Distance from a receiver which you would expect ~0 (or low) detectability
#' @param optim_metric Approach used to optimize array performance
#' @param nrecs The number of receivers to be deployed
#'
#' @return A data frame of performance metrics across iterations.
#' @export


#######DE ARRAY OPTIMIZE 1-OVERLAP FUNCTIONS#########

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

  message(paste("Identifying estimated optimized receiver locations based on: ", max_i_iterations,
                " i iterations. Each i iteration will include ~",
                round(choose(n_distinct(dist_df$station), n.receiver_overlap) * (subsampl_perc_i / 100), 0),
                "possible combinations (at", subsampl_perc_i, "% subsampling rate) from", nrecs, "receiver locations. This may take a while depending on cellsize..."))

  set.seed(seed) #iterative so results may change
  options(datatable.verbose = FALSE)
  options(datatable.print.progress = FALSE)

  #precompute for later (since static)
  #m2_possible_val = sum(as.numeric(st_area(study_area_grid[study_area_grid$Description == "study_area", ]))) #study area only
  m2_possible_val = sum(as.numeric(st_area(study_area_grid))) #study area + spillover

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
    select_stations <- function(dt_unique, nrecs, min_distance, n_attempts = 1000) {
      dmat <- units::drop_units(sf::st_distance(dt_unique$station_loc))
      for (attempt in seq_len(n_attempts)) {
        selected <- integer(0)
        for (i in sample.int(nrow(dt_unique))) {
          if (length(selected) == 0 ||
              all(dmat[i, selected] >= min_distance)) {
            selected <- c(selected, i)
          }
          if (length(selected) == nrecs) {
            return(dt_unique$station[selected])
          }}}
      stop(
        "Unable to generate a starting configuration of ", nrecs,
        " stations with a minimum separation of ", min_distance, " m. Likely min_rec_distance too large for the iteration.")}
    unique_stations_loop<-select_stations(dt_unique, nrecs, min_distance)
    if (length(unique_stations_loop) != nrecs)
      stop("Desired receiver number not met – likely min_rec_distance too large for the iteration")

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
          paste0(sort(new_stations), collapse = ","),
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
                                  m2_possible = m2_possible_val
                                )]

        #evenness calculation
        ncells_study_area <- as.data.table(st_drop_geometry(study_area_grid))
        ncells_study_area <- ncells_study_area[Description == "study_area",.(cell_id)] #if restricting to study area only - advised
        #ncells_study_area <- ncells_study_area[, .(cell_id)] #if including spillover area as well
        ncells_study_area <- merge(ncells_study_area,dt_prob,by = "cell_id",all.x = TRUE,sort = FALSE)
        ncells_study_area[is.na(prob_cell_all), prob_cell_all := 0]
        cv_prob <- (sd(ncells_study_area$prob_cell_all, na.rm = TRUE) /  mean(ncells_study_area$prob_cell_all, na.rm = TRUE))*-1
        dt_prob_summ[, evenness := cv_prob]

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
