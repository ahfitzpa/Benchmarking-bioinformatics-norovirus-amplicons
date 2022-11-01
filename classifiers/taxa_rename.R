##############################################################################################################
################################# FIX TAXONOMY FILE FOR RDP/BLAST #################################################

library(dplyr)
library(tidyr)
library(stringr)
library(Biostrings)
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

noronet_seq <- fasta_table("noronet_capsid_db.fasta")
noronet_tax <- read.table("noronet_taxonomy.tsv", sep = "\t", header=F) %>%
  rename(OTU=V1, tax=V2)

noronet <- left_join(noronet_seq, noronet_tax) %>%
  separate(tax, c("Class","Order","Family","Genus","genogroup","genotype", "genotype_capsid"),";") %>%
  unite(name, c("OTU", "genotype_capsid"), remove=F) %>%
  unite(taxonomy, c("Class","Order","Family","Genus","genogroup","genotype", "genotype_capsid")) %>%
  select(-OTU) %>%
  rename(seq=sequence)

updated_taxonomy <- noronet %>%
  select(name, taxonomy) %>%
  write.table(., "noronet_taxonomy.tsv", quote=F, col.names=F, row.names=F, sep= "\t")

updated_fasta <- noronet %>%
  select(name,seq) %>%
  writeFasta(., "noronet_capsid_db.fasta")

###############################################################################################################################
HuCaT_seq <- fasta_table("calicinet_capsid_db.fasta")
HuCaT_tax <- read.table("HuCaT_taxonomy.tsv", sep = "\t", header=F) %>%
  dplyr::rename(OTU="V1", tax="V2")

custom_seq <- fasta_table("custom_capsid_db.fasta")
custom_tax <- read.table("custom_taxonomy.tsv", sep = "\t", header=F) %>%
  dplyr::rename(OTU=V1, tax=V2)

HuCaT <- left_join(HuCaT_seq, HuCaT_tax) %>%
  dplyr::rename(name=OTU) %>%
  mutate(name=str_replace_all(name, " ", "_"))%>%
  dplyr::rename(seq=sequence)

custom <- left_join(custom_seq, custom_tax) %>%
  dplyr::rename(name=OTU) %>%
  mutate(name=str_replace_all(name, " ", "_"))%>%
  dplyr::rename(seq=sequence)

HuCaT_taxonomy <- HuCaT %>%
  select(name, tax) %>%
  write.table(., "HuCaT_taxonomy.tsv", quote=F, col.names=F, row.names=F, sep= "\t")

custom_taxonomy <- custom %>%
  select(name, tax) %>%
  write.table(., "custom_taxonomy.tsv", quote=F, col.names=F, row.names=F, sep= "\t")

HuCaT_fasta <- HuCaT %>%
  select(name,seq) %>%
  writeFasta(., "calicinet_capsid_db.fasta")

custom_fasta <- custom %>%
  select(name,seq) %>%
  writeFasta(., "custom_capsid_db.fasta")
