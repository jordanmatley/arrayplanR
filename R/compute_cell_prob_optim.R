
#######DE ARRAY OPTIMIZE 1-OVERLAP FUNCTIONS#########
###FUNCTION 2/4

#######DE ARRAY OPTIMIZE 3-OVERLAP FUNCTIONS#########
###FUNCTION 2/5

#calculate optimized receiver-cell probabilities based on designated n.receiver_overlap

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

  coeffs <- Re(stats::fft(prod_vals, inverse = TRUE) / L)[1:(N + 1L)]
  coeffs[coeffs < 0] <- 0
  sum(coeffs[(thr + 1L):(N + 1L)])
}
