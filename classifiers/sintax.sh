###################################################################################################################################################################################################################
# fetch noronet, HuCaT and custom tax and fasta files, generated for RDP
mkdir classifier_comparison/sintax
#noronet
cp classifier_comparison/RDP/clean_noronet.fasta classifier_comparison/sintax/noronet.fasta 
cp classifier_comparison/RDP/noronet_taxonomy.tsv classifier_comparison/sintax/noronet.tsv

# HuCaT
cp classifier_comparison/RDP/clean_HuCaT.fasta classifier_comparison/sintax/HuCaT.fasta 
cp classifier_comparison/RDP/HuCaT_taxonomy.tsv classifier_comparison/sintax/HuCaT.tsv

# custom
cp classifier_comparison/RDP/clean_custom.fasta classifier_comparison/sintax/custom.fasta 
cp classifier_comparison/RDP/custom_taxonomy.tsv classifier_comparison/sintax/custom.tsv
###################################################################################################################################################################################################################
cp sintax_reformat.R classifier_comparison/sintax/sintax_reformat.R
cd classifier_comparison/sintax

module load R/4.0.2
Rscript sintax_reformat.R
module unload R/4.0.2

############################################################################
usearch11="usearch11.0.667_i86linux32

$usearch11 -makeudb_usearch noronet_sintax.fasta -output noronet_sintax.udb
$usearch11 -makeudb_usearch HuCaT_sintax.fasta -output HuCaT_sintax.udb
$usearch11 -makeudb_usearch custom_sintax.fasta -output custom_sintax.udb

###################################################################################################################################################################################################################
mkdir classified_sintax
for i in classifier_comparison/vsearch/*.fasta;
do
name=$(basename "$i" .fasta)
$usearch11 -sintax $i -db custom_sintax.udb -tabbedout classified_sintax/custom_${name}.tsv -strand both -sintax_cutoff 0.8
$usearch11 -sintax $i -db noronet_sintax.udb -tabbedout classified_sintax/noronet_${name}.tsv -strand both -sintax_cutoff 0.8
$usearch11 -sintax $i -db HuCaT_sintax.udb -tabbedout classified_sintax/HuCaT_${name}.tsv -strand both -sintax_cutoff 0.8
done

###############################################################################################################################
for i in vsearch/*.fasta;
do
filename=$(basename "$i" .fasta)
name=$(basename "$i" .fasta |cut -f1 -d.)
awk -v var="$name" 'BEGIN{FS=OFS="\t"}{$1=$1 "_" var "_" "custom"}1' classified_sintax/custom_${filename}.tsv  > classified_sintax/custom_${filename}1.tsv 
awk -v var="$name" 'BEGIN{FS=OFS="\t"}{$1=$1 "_" var "_" "HuCaT"}1' classified_sintax/HuCaT_${filename}.tsv  > classified_sintax/HuCaT_${filename}1.tsv 
awk -v var="$name" 'BEGIN{FS=OFS="\t"}{$1=$1 "_" var "_" "noronet"}1' classified_sintax/noronet_${filename}.tsv  > classified_sintax/noronet_${filename}1.tsv 
done

#################################################################################################################################################################################################################### 
# merge taxonomy outputs for all databases
awk 'FNR>1 || NR==1 {print $0","FILENAME}'  classified_sintax/*1.tsv  >  classified_sintax/taxonomy_sintax.tsv
