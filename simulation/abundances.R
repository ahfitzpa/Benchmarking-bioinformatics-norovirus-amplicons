library(tidyr)
library(dplyr)
library(magrittr)
library(janitor)
library(tibble)


#https://stackoverflow.com/questions/24845909/generate-n-random-integers-that-sum-to-m-in-r
rand_vect <- function(N, M, sd = 1, pos.only = TRUE) {
  vec <- rlnorm(N, M/N, sd)
  if (abs(sum(vec)) < 0.001) vec <- vec + .001
  vec <-vec / sum(vec) * M
  deviation <- M - sum(vec)
  for (. in seq_len(abs(deviation))) {
    vec[i] <- vec[i <- sample(N, 1)] + sign(deviation)
  }
  if (pos.only) while (any(vec < 0)) {
    negs <- vec < 0.001
    pos  <- vec > 0.01
    vec[negs][i] <- vec[negs][i <- sample(sum(negs), 1)] + .01
    vec[pos][i]  <- vec[pos ][i <- sample(sum(pos ), 1)] - .01
  }
  vec
}

add_abundances <- function(df){
  updated_df <- df %>%  filter_all(any_vars(complete.cases(.)))
  n_value <- nrow(updated_df)
  values <- rand_vect(n_value,1)
  updated_df %>% dplyr::mutate(abundance = values) 
}

accessions <- read.table(file = "accession_ids.txt", sep ="", header=F)

nmax <- max(stringr::str_count(accessions$V2, "_")) + 1
check <- accessions %>% separate(V2, paste0("col", seq_len(nmax)), sep = "_", fill = "right")%>%
  pivot_longer(-V1) %>%
  pivot_wider(names_from=V1, values_from=value)  %>%
  select(-name)%>%
  as_tibble()

# Creates a new dataframe from each column, maintains the original column names in the new dataframes
my_list <- list()               # Create empty list
for(i in 1:ncol(check)){
  temp <- data.frame(check[,i])
  colnames(temp) <- colnames(check)[i]
  assign(colnames(check)[i], temp)
  my_list[[i]] <- temp
  names(my_list)[[i]] <-  colnames(temp)
}

updated_list <- lapply(my_list, add_abundances)

for(i in seq_along(updated_list)) {
  write.table(updated_list[[i]], paste0("abundances/",names(updated_list)[i], ".txt", sep = ""), 
              col.names = FALSE, row.names = FALSE,quote = FALSE)
}
