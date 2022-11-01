
#DOWNLOAD ALL AVAILABLE NOROVRUS SEQUENCES FROM NCBI AND SAVE IN input_files folder
cd input_files
mkdir cv_input

# common non norovirus and norovirus sequences to be tested in classifiers
#cat *.fasta > test_all.fasta

# mask ambigious sequences
module load emboss/6.6.0
maskambignuc -sequence cv_input/test_all.fasta -outseq test_masked.fasta
module unload emboss/6.6.0

# Remove duplicates
module load bbmap/38.22
dedupe.sh in=cv_input/test_masked.fasta out=nodupes.fasta
module unload bbmap/38.22

# change to single line fasta file and 
module load fastx_toolkit/0.0.14
fasta_formatter -i nodupes.fasta -o test_norovirus2.fasta -w 0
module unload fastx_toolkit/0.0.14

# keep only accession_id
sed '/^>/s/^>\([^ ]*\) .*/>\1 /' test_norovirus2.fasta > test_norovirus.fasta

module load cdhit/4.7
cd-hit-est -i test_norovirus.fasta -o test_norovirus_99.fasta -c 0.99 -n 6 -d 0 -M 16000 -T 8
module unload cdhit/4.7

out=test_norovirus.fasta.split

module load seqkit
seqkit split2 --by-size 1 test_norovirus_99.fasta -O $out

find $out -name "*.fasta" \
    | while read f; do \
        mv $f $out/$(seqkit seq --name --only-id $f).fasta; \
    done

module unload seqkit
