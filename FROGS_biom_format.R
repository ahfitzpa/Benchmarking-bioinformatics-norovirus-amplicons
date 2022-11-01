library(dplyr)
library(biomformat)
library(stringr)
library(janitor)

otu_table <- read.delim(file = "otu_table.txt", sep = "\t", header =F)%>% 
	slice(-1)%>%
	janitor::row_to_names(row_number = 1) %>%
	dplyr::rename_all(~str_replace_all(.,"_","\\.")) 
write.table(otu_table, file='input_otu.tsv', quote=FALSE, sep='\t', row.names=F)
