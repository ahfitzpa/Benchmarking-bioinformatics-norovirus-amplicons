if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("DECIPHER")
library(DECIPHER)
library(dplyr)
library(magrittr)
library(tidyverse)
library(stringr)

train_classifier <- function(fasta, tax) {
# specify the path to your file of training sequences:
seqs_path <- fasta
# read the sequences into memory
seqs <- readDNAStringSet(seqs_path)
# Alternatively use readAAStringSet or readRNAStringSet

# (optionally) specify a path to the taxid file:
rank_path <- tax
taxid <- read.table(rank_path,
                      header=FALSE,
                      col.names=c('Index', 'Name', 'Parent', 'Level', 'Rank'),
                      sep="*", # asterisks delimited
                      quote="", # preserve quotes
                      stringsAsFactors=FALSE)

# obtain the taxonomic assignments
groups <- names(seqs) # sequence names
# assume the taxonomy begins with 'Root;'
groups <- gsub("(.*)(Root;)", "\\2", groups) # extract the group label
groupCounts <- table(groups)
u_groups <- names(groupCounts) # unique groups
length(u_groups) # number of groups

maxGroupSize <- 10 # max sequences per label (>= 1)
remove <- logical(length(seqs))
for (i in which(groupCounts > maxGroupSize)) {
  index <- which(groups==u_groups[i])
  keep <- sample(length(index),
                 maxGroupSize)
  remove[index[-keep]] <- TRUE
}
sum(remove) # number of sequences eliminated

maxIterations <- 3 # must be >= 1
allowGroupRemoval <- FALSE
probSeqsPrev <- integer() # suspected problem sequences from prior iteration
for (i in seq_len(maxIterations)) {
  cat("Training iteration: ", i, "\n", sep="")
  # train the classifier
  trainingSet <- LearnTaxa(seqs[!remove],
                           names(seqs)[!remove],
                           taxid)
  # look for problem sequences
  probSeqs <- trainingSet$problemSequences$Index
  if (length(probSeqs)==0) {
    cat("No problem sequences remaining.\n")
    break
  } else if (length(probSeqs)==length(probSeqsPrev) &&
             all(probSeqsPrev==probSeqs)) {
    cat("Iterations converged.\n")
    break
  }
  if (i==maxIterations)
    break
  probSeqsPrev <- probSeqs
  # remove any problem sequences
  index <- which(!remove)[probSeqs]
  remove[index] <- TRUE # remove all problem sequences
  if (!allowGroupRemoval) {
    # replace any removed groups
    missing <- !(u_groups %in% groups[!remove])
    missing <- u_groups[missing]
    if (length(missing) > 0) {
      index <- index[groups[index] %in% missing]
      remove[index] <- FALSE # don't remove
    }
  }
}

sum(remove) # total number of sequences eliminated
length(probSeqs) # number of remaining problem sequences
return(trainingSet)
}
fasta_files <- Sys.glob("*.fasta")
tax_files <- Sys.glob("*.tax")

classifiers <- purrr::map2(fasta_files,tax_files, ~ train_classifier(fasta=.x, tax=.y))
names(classifiers) <- c("custom", "HucaT", "noronet")

vsearch_output <- Sys.glob("/data/Food/analysis/R6564_NGS/amy_fitzpatrick/denoise_comparison/classifier_comparison/vsearch/*.fasta")

classify_vsearch <- function(seq_output, classifier){
  test <-  readDNAStringSet(seq_output)
  trainingSet <- classifier
  filename <- names(classifier)
  
  ids <- IdTaxa(test,
                trainingSet,
                type="extended",
                strand="top",
                threshold=80,
                processors=1)
  
  assignment <- sapply(ids,
                       function(x)
                         paste(x$taxon,
                               collapse=";"))
  
  clean_assignment <- assignment %>% 
    as.data.frame() %>%
    rename(taxonomy= ".") %>%
    tidyr::separate(taxonomy, c("Class", "Order", "Family", "Genus", "genogroup", "genotype", "genotype_capsid"), ";") %>%
    mutate(across(Class:genotype_capsid, ~ str_remove(.x, "\\[.*")))
 return(clean_assignment)
}

custom <- list()
noronet <- list()
HuCaT <- list()
list_output <- list()
for (i in vsearch_output){
  custom[i] <- list(classify_vsearch(i, classifiers[[1]]))
  HuCaT[i] <- list(classify_vsearch(i, classifiers[[2]]))
  noronet[i] <- list(classify_vsearch(i, classifiers[[3]]))
}

names(custom) <- c("custom_001","custom_002","custom_003","custom_004","custom_005","custom_006","custom_007","custom_008","custom_009","custom_010")
names(HuCaT) <- c("HuCaT_001","HuCaT_002","HuCaT_003","HuCaT_004","HuCat_005","HuCaT_006","HuCaT_007","HuCaT_008","HuCaT_009","HuCaT_010")
names(noronet) <- c("noronet_001","noronet_002","noronet_003","noronet_004","noronet_005","noronet_006","noronet_007","noronet_008","noronet_009","noronet_010")

lapply(names(HuCaT), function(x) write.table(HuCaT[[x]], file=paste(x,"_idtaxa.tsv"), quote = FALSE, row.names = F, sep="\t"))
lapply(names(custom), function(x) write.table(custom[[x]], file=paste(x,"_idtaxa.tsv"), quote = FALSE, row.names = F, sep="\t"))
lapply(names(noronet), function(x) write.table(noronet[[x]], file=paste(x,"_idtaxa.tsv"), quote = FALSE, row.names = F, sep="\t"))
