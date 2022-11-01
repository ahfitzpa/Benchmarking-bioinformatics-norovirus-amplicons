###################################################################################################################################################################################################################
# fetch noronet, HuCaT and custom tax and fasta files, generated for RDP
mkdir classifier_comparison/IDTAXA 
#noronet
cp classifier_comparison/RDP/noronet_RDP.fasta classifier_comparison/IDTAXA/noronet_RDP.fasta 
cp classifier_comparison/RDP/noronet_RDP.tax classifier_comparison/IDTAXA/noronet_RDP.tax

# HuCaT
cp classifier_comparison/RDP/HuCaT_RDP.fasta classifier_comparison/IDTAXA/HuCaT_RDP.fasta 
cp classifier_comparison/RDP/HuCaT_RDP.tax classifier_comparison/IDTAXA/HuCaT_RDP.tax

# custom
cp classifier_comparison/RDP/custom_RDP.fasta classifier_comparison/IDTAXA/custom_RDP.fasta 
cp classifier_comparison/RDP/custom_RDP.tax classifier_comparison/IDTAXA/custom_RDP.tax

cp IDTAXA.R classifier_comparison/IDTAXA/IDTAXA.R
cd classifier_comparison/IDTAXA
###################################################################################################################################################################################################################
module load R/4.0.2
Rscript IDTAXA.R
module unload R/4.0.2
###############################################################################################################################
for i in classifier_comparison/vsearch/*.fasta;
do
filename=$(basename "$i" .fasta)
name=$(basename "$i" .fasta |cut -f1 -d.)
echo $name
awk -v var="$name" 'BEGIN{FS=OFS="\t"}{$1=$1 "_" var "_" "custom"}1' custom_${name}_idtaxa.tsv  > custom_${nname}_1idtaxa.tsv 
awk -v var="$name" 'BEGIN{FS=OFS="\t"}{$1=$1 "_" var "_" "HuCaT"}1' HuCaT_${name}_idtaxa.tsv  > HuCaT_${name}_1idtaxa.tsv 
awk -v var="$name" 'BEGIN{FS=OFS="\t"}{$1=$1 "_" var "_" "noronet"}1' noronet_${name}_idtaxa.tsv  > noronet_${name}_1idtaxa.tsv 
done

#################################################################################################################################################################################################################### 
# merge taxonomy outputs for all databases
awk 'FNR>1 || NR==1 {print $0 " " FILENAME}' *1idtaxa.tsv  >  taxonomy_IDTAXA.tsv
awk 'FNR>1 || NR==1 {print $0","FILENAME}'  *_1idtaxa.tsv   >   taxonomy_IDTAXA.tsv
