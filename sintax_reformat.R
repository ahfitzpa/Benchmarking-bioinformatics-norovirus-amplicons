library(dplyr)
library(magrittr)
library(tidyverse)
library(stringr)

fasta_table <- function(fasta){
  fastaFile <- Biostrings::readDNAStringSet(fasta, format = "fasta",use.names = TRUE)
  OTU <- names(fastaFile)
  sequence <-  paste(fastaFile)
  df <- data.frame(OTU, sequence)
  return(df)
}

#functions
writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

fasta_list <- Sys.glob("*.fasta")
tax_files <- Sys.glob("*.tsv")

fasta_dfs <- lapply(fasta_list, fasta_table)
tax_dfs <- lapply(tax_files, read.table, sep="\t") %>%
  lapply(., dplyr::rename, OTU=V1, taxonomy=V2 )
names(fasta_dfs) <- c("custom", "HuCaT", "noronet")
names(tax_dfs) <- c("custom", "HuCaT", "noronet")

custom <- left_join(fasta_dfs[[1]], tax_dfs[[1]])
HuCaT <- left_join(fasta_dfs[[2]], tax_dfs[[2]])
noronet <- left_join(fasta_dfs[[3]], tax_dfs[[3]])

list_classifiers <- list(custom, HuCaT, noronet)

reformat_fasta <- function(classifier){
classifier %>%
    rename(seq=sequence) %>%
 tidyr::separate(taxonomy, c("d", "p", "c", "o", "f", "g", "s"), ";") %>%
  mutate(d= paste0("d:", d)) %>%
  mutate(p= paste0("p:", p)) %>%
  mutate(c= paste0("c:", c)) %>%
  mutate(o= paste0("o:", o)) %>%
  mutate(f= paste0("f:", f)) %>%
  mutate(g= paste0("g:", g)) %>%
  mutate(s= paste0("s:", s)) %>%
  mutate(name=paste0(OTU,";tax=",d,",",p,",",c,",",o,",",f,",",g,",",s))
}

updated_fasta <- lapply(list_classifiers, reformat_fasta)

writeFasta(data=updated_fasta[[1]], filename="custom_sintax.fasta")
writeFasta(data=updated_fasta[[3]], filename="HuCaT_sintax.fasta")
writeFasta(data=updated_fasta[[3]], filename="noronet_sintax.fasta")
