### IF SCRIPT NOT WORKING check FASTA headers (formatted correctly) and check that the capsid cut products are of the correct size #####
dir=$(pwd)
mkdir simulated_sequencing 
cd simulated_sequencing 
#############################################################################################################################################################################################   
# set up sample files for simulated data including, number of genotypes, reads and names of samples
module load R/4.0.2
Rscript randomdiversity2.R
module unload R/4.0.2

tail -n +2 sample.txt > sample2.txt
tail -n +2 reads.txt > reads2.txt ######## creates random reads in a truncated log normal distribution
tail -n +2 genotypes.txt > genotypes2.txt ###### names for each output genome

rm -f sample.txt
rm -f reads.txt
rm -f genotypes.txt
#############################################################################################################################################################################################  
# pick random files from available fasta files for norovirus

# create a folder for each sample
cat sample2.txt | xargs -L 1 mkdir

# 1. create commands for to pick n number of fasta files for each folder
awk '{print "find ../../capsid_sequences/capsid3.fasta.split/*.fasta -type f | shuf -n "$0" > genotypes_list.txt"}' genotypes2.txt >> genotypes_list.txt  

# 2. Pick one line from file and place in each subdirectory
for d in */; do
  (cd $d && cat ../genotypes_list.txt | shuf -n 1  > genotypes_list.sh)
done

# 3. For each subdirectory, carry out commands from part 1
for d in */;  
do 
	(cd $d
 	sh genotypes_list.sh

	awk -F '\t' '{print"head -n 1 "$1"  | cut -f 2 >> species.txt"}' genotypes_list.txt > species.sh
	# extract names of fasta sequences
	sh species.sh

	# combine all fasta files for input into random reads. 
	xargs -a genotypes_list.txt cp -t .
	cat *.fasta > sample.fasta 

	# extract only sequnece name of fasta file - accession number
	awk '{print $1}' species.txt > Gnomenames.txt

	# combine names into one row. 
	cat Gnomenames.txt | sed 's/>//' > genotypes2.txt
	cat genotypes2.txt | paste -s -d "_" > genotypes.txt)
done

# 4. rename sample.fasta after parent directory and move up directory
find "../simulated_sequencing/" -type f -iname 'sample.fasta' -exec sh -c '
    path="${1%/*}"; filename="${1##*/}";
    cp -nv "${1}" ../"simulated_sequencing/${path##*/}.fasta" ' sh_cp {} \;

# 5. rename genotypes.txt file after parent directory and move up directory
find "../simulated_sequencing/" -type f -iname '*genotypes.txt' -exec sh -c '
    path="${1%/*}"; filename="${1##*/}";
    cp -nv "${1}" ../"simulated_sequencing/${path##*/}${filename}" ' sh_cp {} \;

# 6. Remove subdirectories and content. 
rm -rf ../simulated_sequencing/*/
