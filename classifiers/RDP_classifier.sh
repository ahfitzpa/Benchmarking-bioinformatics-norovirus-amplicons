
###################################################################################################################################################################################################################
# fetch noronet, HuCaT and custom tsv and fasta files
mkdir classifier_comparison/RDP
#noronet
cp norovirus_capsid_db_taxonomy.tsv classifier_comparison/RDP/noronet_taxonomy.tsv
cp noronet_capsid_db.fasta  classifier_comparison/RDP/noronet_capsid_db.fasta

# HuCaT
cp calicinet_capsid_db_taxonomy.tsv classifier_comparison/RDP/HuCaT_taxonomy.tsv
cp calicinet_capsid_db.fasta classifier_comparison/RDP/calicinet_capsid_db.fasta

# custom
cp custom_capsid_db_taxonomy.tsv classifier_comparison/RDP/custom_taxonomy.tsv
cp custom_capsid_db.fasta  classifier_comparison/RDP/custom_capsid_db.fasta

cd classifier_comparison/RDP

#################################################################################################################################################################################
cp taxa_rename.R taxa_rename.R
module load R/4.0.2
Rscript taxa_rename.R
module unload R/4.0.2
###################################################################################################################################################################################################################

# count samples with full taxonomy 
for i in *taxonomy.tsv; do 
cat $i | while read -r line; do 
    echo -n "${line} " 
    grep -o ";" <<<"$line"| wc -l 
  done > 1_${i}

# remove samples without full taxonomy and remove count column of levels of classification
awk '$3 >= 6' FS=' ' 1_${i} | awk '!($3="")' > 2_${i}

# send accession names to file
awk '{print $1}' 2_${i} > ids_${i}

done

rm 1*
rm 2*

#remove sequences without complete taxonomy 
module load seqkit/1.4  
seqkit grep -f ids_noronet_taxonomy.tsv noronet_capsid_db.fasta > clean_noronet.fasta
seqkit grep -f ids_HuCaT_taxonomy.tsv calicinet_capsid_db.fasta > clean_HuCaT.fasta
seqkit grep -f ids_custom_taxonomy.tsv custom_capsid_db.fasta > clean_custom.fasta
module unload seqkit/1.4  

###################################################################################################################################################################################################################
for i in clean*; do
  outdir=$(basename "$i" .fasta |cut -f2 -d"_")
  echo "$outdir"
 
module load python3.7/i3
python3.7 fasta2RDP.py -d $i -t ${outdir}_taxonomy.tsv -r Domain Phylum Class Order Family Genus Species --rdp-taxonomy ${outdir}_RDP.tax --rdp-fasta ${outdir}_RDP.fasta 
module unload python3.7/i3

### copy properties file ###
# look in preformated database we provided, and copy properties file and renamed it
cp silva_128_16S/silva_128_16S.fasta.properties ${outdir}_RDP.fasta.properties

### format for blast ###
module load blast/2.8.1
makeblastdb -in ${outdir}_RDP.fasta -dbtype nucl 
module unload blast/2.8.1

done

####################################################################################################################################################################################################################
# download and unzip RDP classifier
classifier_RDP="rdp_classifier_2.13/dist/classifier.jar"
java -Xmx30g -jar $classifier_RDP  train -o HuCaT -s HuCaT_RDP.fasta -t HuCaT_RDP.tax
java -Xmx30g -jar $classifier_RDP train -o noronet -s noronet_RDP.fasta -t noronet_RDP.tax
java -Xmx30g -jar $classifier_RDP  train -o custom -s custom_RDP.fasta -t custom_RDP.tax

cp noronet_RDP.fasta.properties noronet/noronet_RDP.fasta.properties
cp HuCaT_RDP.fasta.properties HuCaT/HuCaT_RDP.fasta.properties
cp custom_RDP.fasta.properties custom/custom_RDP.fasta.properties

#################################################################################################################################################################################################################### 
# classify 10 vsearch datasets
mkdir classified_RDP
for i in classifier_comparison/vsearch/*.fasta;
do
name=$(basename "$i" .fasta)
java -Xmx8g -jar $classifier_RDP -t HuCaT/HuCaT_RDP.fasta.properties -c 0.5 -o classified_RDP/HuCaT_${name}.classified.tsv $i
java -Xmx8g -jar $classifier_RDP -t noronet/noronet_RDP.fasta.properties -c 0.5 -o classified_RDP/noronet_${name}.classified.tsv $i
java -Xmx8g -jar $classifier_RDP -t custom/custom_RDP.fasta.properties -c 0.5 -o classified_RDP/custom_${name}.classified.tsv $i

done

#################################################################################################################################################################################################################### 
# merge taxonomy outputs for all databases
awk 'FNR>1 || NR==1 {print $0","FILENAME}'  classified_RDP/*.classified.tsv  >  classified_RDP/taxonomy_RDP.tsv
