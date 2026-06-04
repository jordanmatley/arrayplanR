#' @import data.table
NULL

#' Compute receiver-cell probabilities
#'
#' calculate receiver-cell probabilities based on designated n.receiver_overlap
#'
#' @param station_df An sf receiver-cell distance dataframe (from compute_station_distances function).
#' @param study_area_grid An sf object containing the gridded study area 
#' @param study_area The number of receivers to detect transmissions from each grid .
#'
#' @return A data frame of all cell detection probabilities.
#' @export 


#######DE ARRAY CALCULATE FUNCTION#########

###FUNCTION 2/2
#calculate receiver-cell probabilities based on designated n.receiver_overlap
compute_cell_prob <- function(dist_df, study_area_grid, n.receiver_overlap, progress = TRUE) {
  
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
    mutate(prob_cell_all = tidyr::replace_na(prob_cell_all, 0))
}

