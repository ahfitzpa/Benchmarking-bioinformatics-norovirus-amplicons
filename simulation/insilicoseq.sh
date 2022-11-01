cp abundances.R simulated_sequencing
cp functions.R simulated_sequencing
cp tags.txt simulated_sequencing

cd simulated_sequencing 
##########################################################################################################################################################################################
# merge all genome files (contain accession numbers of genotypes in each fasta file) and include file name 
for y in $(ls sample*genotypes.txt);do a="$(cut -d 'g' -f 1 <<< "$y")";for x in $(cat "$y" );do echo "$a" "$x"  >> genome_list.txt ;done;done
rm sample_*genotypes.txt

#create file with sample name, assigned reads and list of accession numbers in each fasta
cut -f 1,2 genome_list.txt > accession_ids.txt
paste  genome_list.txt reads2.txt quality_scores.txt tags.txt| column -s $'\t' -t > outputfile.txt

################################################################################################################################################
# add indexes to sequences, based on Nextera XT kit
module load seqkit/1.4 
# use seqkit module to add indexes to sequences (each index is 8 bp long)
awk -F ' ' '{print"cat "$1".fasta | seqkit mutate -i 0:"$8" > "$1"_a.fasta"}' outputfile.txt > forward_indexes.sh
sh forward_indexes.sh
awk -F ' ' '{print"cat "$1"_a.fasta | seqkit mutate -i -1:"$10" > "$1"_b.fasta"}' outputfile.txt > reverse_indexes.sh
sh reverse_indexes.sh
module unload seqkit/1.4 

# remove intermediary files
awk -F ' ' '{print"rm "$1".fasta"}' outputfile.txt > remove_orignal_fasta.sh
awk -F ' ' '{print"rm "$1"_a.fasta"}' outputfile.txt > remove_intermediary_fasta.sh
sh remove_orignal_fasta.sh
sh remove_intermediary_fasta.sh

# check average length of sequences
awk '{/>/&&++a||b+=length()}END{print b/a}' *.fasta > length_indexed_fasta.txt

################################################################################################################################################
module load python3.7/i3
pip3 install --user InSilicoSeq
pip install biopython==1.78
#############################################################################################################################################################################################  
#create file
mkdir insilicoseq
awk -F ' ' '{print"iss generate --cpus 20 --genomes "$1"_b.fasta --abundance 'zero_inflated_lognormal' --n_reads "$3" --model miseq --output insilicoseq/"$1""}' outputfile.txt > insilico_simulations.sh 

sh insilico_simulations.sh ########## cuts up genomes into 300bp length

#############################################################################################################################################################################################    	

module unload python3.7/i3
