
mkdir raw_sequences
cp simulated_sequencing/insilicoseq/*.fastq raw_sequences

result=${PWD##*/}
annotate=$(echo "$result"| awk -F'[_.]' '{print $4}')

cd raw_sequences
for i in *R1.fastq
  do
  name=$(basename "$i" | sed -e 's/.........$//');
  mv $i ${name}.${annotate}_L001_R1_001.fastq;
done

for i in *R2.fastq
  do
  name=$(basename "$i" | sed -e 's/.........$//');
  mv $i ${name}.${annotate}_L001_R2_001.fastq;
done

gzip *.fastq

cd ..
######################################################################################################################################################################################################################################
mkdir metadata

cd simulated_sequencing/insilicoseq
# combine expected abundances into one file and remove individual txt files
#awk 'BEGIN{print "filename","accession.id", "abundance"}FNR>1{print FILENAME" "$0}' *.txt > expected_composition1.txt

#rm *_abundance.txt

cd ../..

# set up updated metadata file based on expected composition and correct genotypic classification of included sequences

cp simulated_sequencing/insilicoseq/expected_composition1.txt metadata/expected_composition1.txt 
# cut accession ids out of expected composition file
awk '{print $2}' metadata/expected_composition1.txt | sed '1d' > metadata/accession.txt
awk '$3>0 {print $2}' metadata/expected_composition1.txt | sed '1d' > metadata/accession.txt

# fetch sequences included in simulation by accession id#
cp test_norovirus.part_*.fasta  metadata/test_norovirus.fasta
cat metadata/test_norovirus.fasta | tr -d "_" >  metadata/test_norovirus1.fasta
cat metadata/accession.txt| tr -d "_" >  metadata/accession1.txt

# change to single line fasta file and 
module load fastx_toolkit/0.0.14
fasta_formatter -i  metadata/test_norovirus1.fasta -o metadata/test_norovirus2.fasta -w 0
module unload fastx_toolkit/0.0.14

module load seqkit/1.4
grep -w -A 1 -f metadata/accession1.txt metadata/test_norovirus2.fasta --no-group-separator > metadata/expected_sequences.fasta
module unload seqkit/1.4

grep -Fwf metadata/accession1.txt metadata/expected_composition1.txt > metadata/expected_composition.txt

################################################################################################################################################################################################
cd metadata
module load seqkit/1.4 
#cut amplicons from full length genomes
seqkit amplicon -F CNTGGGAGGGCGATCGCAA -R CCRCCNGCATRNCCRTTRTACAT -m 7 expected_sequences.fasta > expected_GII.fasta
seqkit seq -M 400 expected_GII.fasta  > expected_GII_400.fasta
seqkit amplicon -F CTGCCCGAATTYGTAAATGA -R CCAACCCARCCATTRTACA -m 5 expected_sequences.fasta > expected_GI.fasta
seqkit seq -M 400 expected_GI.fasta > expected_GI_400.fasta

# check average length of sequences
awk '{/>/&&++a||b+=length()}END{print b/a}' expected_GII_400.fasta
awk '{/>/&&++a||b+=length()}END{print b/a}' expected_GI_400.fasta
module unload seqkit/1.4

# common non norovirus and norovirus sequences to be tested in classifiers
cat expected_GII_400.fasta expected_GI_400.fasta > expected_sequences1.fasta
# remove empty fasta lines
#https://www.biostars.org/p/78786/
awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' expected_sequences1.fasta > expected_sequences2.fasta

# classify sequences using noronet external tool
