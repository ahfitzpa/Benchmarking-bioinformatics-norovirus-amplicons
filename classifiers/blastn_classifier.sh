mkdir classifier_comparison/blast
cp classifier_comparison/RDP/custom_RDP.fasta classifier_comparison/blast/custom_RDP.fasta
cp classifier_comparison/RDP/custom_RDP.tax classifier_comparison/blast/custom_RDP.tax
cp classifier_comparison/RDP/custom_taxonomy.tsv classifier_comparison/blast/custom_taxonomy.tsv

cp classifier_comparison/RDP/HuCaT_RDP.fasta classifier_comparison/blast/HuCaT_RDP.fasta
cp classifier_comparison/RDP/HuCaT_RDP.tax classifier_comparison/blast/HuCaT_RDP.tax
cp classifier_comparison/RDP/HuCaT_taxonomy.tsv classifier_comparison/blast/HuCaT_taxonomy.tsv

cp classifier_comparison/RDP/noronet_RDP.fasta classifier_comparison/blast/noronet_RDP.fasta
cp classifier_comparison/RDP/noronet_RDP.tax classifier_comparison/blast/noronet_RDP.tax
cp classifier_comparison/RDP/noronet_taxonomy.tsv classifier_comparison/blast/noronet_taxonomy.tsv

cd classifier_comparison/blast
module load blast/2.8.1
for i in *.fasta;
do
name=$(basename "$i" .fasta)
makeblastdb -in $i -parse_seqids -taxid_map ${name}.tax -dbtype nucl -out ${name}
done

#######################################################################################################################################################################
mkdir classified_blast
module load python3.7/i3
for i in vsearch/*.fasta;
do
name=$(basename "$i" .fasta)
assign_taxonomy="classifier_comparison/blast/Assign-Taxonomy-with-BLAST-master/taxonomy_assignment_BLAST.py"

blastn -db HuCaT_RDP -query $i -out classified_blast/HuCaT_${name}.blastout.txt -num_threads 5 -outfmt "6 qseqid qlen sseqid pident length qstart qend sstart send evalue bitscore" 
blastn -db noronet_RDP -query $i -out classified_blast/noronet_${name}.blastout.txt -num_threads 5 -outfmt "6 qseqid qlen sseqid pident length qstart qend sstart send evalue bitscore" 
blastn -db custom_RDP -query $i -out classified_blast/custom_${name}.blastout.txt -num_threads 5 -outfmt "6 qseqid qlen sseqid pident length qstart qend sstart send evalue bitscore" 

python3 $assign_taxonomy --blast_flavor blastn \
--cutoff_species 90 --cutoff_family 70 --cutoff_phylum 60 --length_percentage 0.50 \
--blast_database HuCat_RDP \
--output_dir classified_blast/HuCaT_${name} \
--blast_file classified_blast/HuCaT_${name}.blastout.txt $i HuCaT_taxonomy.tsv

python3 $assign_taxonomy --blast_flavor blastn \
--cutoff_species 90 --cutoff_family 70 --cutoff_phylum 60 --length_percentage 0.50 \
--blast_database noronet_RDP \
--output_dir classified_blast/noronet_${name} \
--blast_file classified_blast/noronet_${name}.blastout.txt $i noronet_taxonomy.tsv

python3 $assign_taxonomy --blast_flavor blastn \
--cutoff_species 90 --cutoff_family 70 --cutoff_phylum 60 --length_percentage 0.50 \
--blast_database custom_RDP \
--output_dir classified_blast/custom_${name} \
--blast_file classified_blast/custom_${name}.blastout.txt $i custom_taxonomy.tsv

done

#######################################################################################################################################################################
#extract taxonomy files out of subdirectories
#rename sample.fasta after parent directory and move up directory
find "classified_blast/" -type f -iname 'taxonomy_assignment_per_sequence.tsv' -exec sh -c '
    path="${1%/*}"; filename="${1##*/}";
    cp -nv "${1}" "classified_blast/${path##*/}.tsv" ' sh_cp {} \;

# Remove subdirectories and content. 
rm -rf classified_blast/*/

# merge all taxonomy outputs but include file name

awk 'FNR>1 || NR==1 {print $0","FILENAME}' classified_blast/*.vsearch.tsv  > classified_blast/taxonomy_blast.tsv
