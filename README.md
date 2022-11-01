# Benchmarking-bioinformatic-tools-for-amplicon-based-sequencing-of-norovirus-

Using a simulated sequencing dataset, denoising-based (DADA2, Deblur and USEARCH-UNOISE3) and clustering-based pipelines (VSEARCH and FROGS) were compared with respect to their ability to represent composition and sequence information. Open source classifiers (RDP, BLASTn, IDTAXA, QIIME2 naive Bayes and sintax) were trained using three different databases; a custom database, the NoroNet database and the Human Calicivirus database. Each classifier and database combination was compared from the perspective of their classification accuracy. 

This repository contains the final version of pipelines and classifier scripts used to process the simulated sequencing data. 
