# set up a file with 40 sample
# between 2-10 genotypes per sample 
library(dplyr)
library(magrittr)

sample <- as.data.frame(matrix(paste0("sample", "", 01:40), ncol = 1, nrow = 40)) %>% 
  dplyr::rename(sample =V1) %>%
  dplyr::mutate(n_genotypes =sample(2:10,40, replace=T ))

list <- c(sample[,2])

write.csv(list, "genotypes.txt", row.names=F)
write.csv(sample[,1],"sample.txt", quote=FALSE,row.names=F)

# distribution based on input reads in statsdada2 files from various sequencing runs
#clean lognormal distribution reads accross samples and features
# use truncated log normal distribution to prevent most samples centering at 0 reads
library(EnvStats)

rand_vect_cont <- function(N, S, sdlog, min,max) {
  vec <-  data.frame(reads=abs(EnvStats::rlnormTrunc(N, meanlog=8, sdlog, min, max)))
  all <- (vec / sum(vec)) * S
  even <- all %>% 
    dplyr::mutate(reads=as.integer(reads)) %>%
    dplyr::mutate(reads=ifelse((reads%%2)==0,reads, reads+1)) 
  even
}

max_reads <-  50000000
vec <- rand_vect_cont(40, S=max_reads,sdlog=10, min=100, max=3000000)

write.csv(vec, "reads.txt", row.names=F)

## normal distribution of quality scores for MiSeq runs
library(EnvStats)
quality_scores <- data.frame(high_q = as.integer(rnormTrunc(40, mean = 36, sd = 2, min = 34, max = 40)), 
                             low_q =as.integer(rnormTrunc(40, mean = 33, sd = 2, min = 20, max = 33)))
write.table(quality_scores, "quality_scores.txt",quote=FALSE, col.names=F, row.names=F,sep = "\t")
