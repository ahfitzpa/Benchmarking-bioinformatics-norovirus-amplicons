#https://stackoverflow.com/questions/24845909/generate-n-random-integers-that-sum-to-m-in-r
rand_vect <- function(N, M, sd = 1, pos.only = TRUE) {
  vec <- rlnorm(N, M/N, sd)
  if (abs(sum(vec)) < 0.01) vec <- vec + 1
  vec <- round(vec / sum(vec) * M)
  deviation <- M - sum(vec)
  for (. in seq_len(abs(deviation))) {
    vec[i] <- vec[i <- sample(N, 1)] + sign(deviation)
  }
  if (pos.only) while (any(vec < 0)) {
    negs <- vec < 0
    pos  <- vec > 0
    vec[negs][i] <- vec[negs][i <- sample(sum(negs), 1)] + 1
    vec[pos][i]  <- vec[pos ][i <- sample(sum(pos ), 1)] - 1
  }
  vec
}

add_abundances <- function(df){
  updated_df <- df %>%  filter_all(any_vars(complete.cases(.)))
  n_value <- nrow(updated_df)
  values <- rand_vect(n_value,100)
  updated_df %>% dplyr::mutate(abundance = values) 
}
