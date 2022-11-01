mkdir vsearch 

cp metadata/updated_metadata.tsv vsearch/metadata1.tsv

cd vsearch
cp map.pl map.pl
cp FROGS_biom_format.R FROGS_biom_format.R

pipeline_name=${PWD##*/}
simulation=$(echo ${PWD} | sed -e "s/.*\/\([^/]*\)\/[^/]*/\1/" | awk -F'[_.]' '{print $4}')
awk -v var="$pipeline_name" 'BEGIN{FS=OFS="\t"}NR>1{$1=$1 "." var}1' metadata1.tsv > metadata.tsv
rm metadata1.tsv

date

mkdir raw_data
cp ../cutadapt/removed_primers/*.fastq.gz raw_data/
gunzip raw_data/*.gz

cd raw_data

for f in *_L001_R2_001.fastq; do
    mv -- "$f" "${f%.${simulation}_L001_R2_001.fastq}_${simulation}_${pipeline_name}_L001_R2_001.fastq"
done

for f in *_L001_R1_001.fastq; do
   mv -- "$f" "${f%.${simulation}_L001_R1_001.fastq}_${simulation}_${pipeline_name}_L001_R1_001.fastq"
done
cd ..

###################################################################################################################################################################################################
module load vsearch/2.3.4
raw_data="raw_data"
merged_data="2.merged_reads"
# Maximum no of mismatches in the alignment - Default 5. Consider increasing if you have long overlaps.
maxdiffs="20"
threads="10"
# Discard pair if alignment is shorter than this value - Default 16 bp
overlap="50"

mkdir ${merged_data}
mkdir working1

#*****************************************************************************************
# Step1: merge data with vsearch

for file1 in ${raw_data}/*R1_001.fastq
    do

        echo ""
        echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        echo Merging paired reads
        echo forward reads are:
        echo $(basename ${file1})
        echo reverse reads are:
        echo $(basename ${file1} R1_001.fastq)R2_001.fastq

    vsearch --fastq_mergepairs ${file1} --threads ${threads} --reverse "${raw_data}/$(basename -s R1_001.fastq ${file1})R2_001.fastq" --fastq_maxdiffs ${maxdiffs} --fastq_minovlen ${overlap} --fastqout "working1/$(basename "$file1")" --fastq_eeout 
done
#*****************************************************************************************
# Step 2: Remove "_L001_R1_001" from filenames

for file2 in working1/*.fastq
    do

        rename="$(basename ${file2} _L001_R1_001.fastq).fastq"

        mv ${file2} ${merged_data}/${rename}
done

#*****************************************************************************************
# Removing working directory

        rm -r working1
###################################################################################################################################################################################################
Q_stats="3.quality_stats"

    echo
    echo Calculate quality statistics

mkdir ${Q_stats}

for file3 in ${merged_data}/*.fastq
    do
        echo ""
        echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        echo Quality control and removing dimer seqs
        echo input is:
        echo ${file3}

    vsearch --fastq_eestats  ${file3} --output "${Q_stats}/$(basename "$file3" .fastq).stats"
done    

###################################################################################################################################################################################################
QF_reads="4.quality_filtered"
# Discard reads with > E total expected errors for all bases in the read after any truncation options have been applied. Usually set at 1, however may want to lower to 0.5 or 0.25 for more stringent filtering.
max_ee="1"
minlen="100"
maxlen="400"
maxns="0"

    echo Quality filtering

mkdir ${QF_reads}

for file4 in ${merged_data}/*.fastq
    do
        echo ""
        echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        echo Quality control and removing dimer seqs
        echo input is:
        echo ${file4}


    vsearch --fastq_filter ${file4} \
        --fastq_maxee ${max_ee} \
        --fastq_minlen ${minlen} \
        --fastq_maxlen ${maxlen} \
        --fastq_maxns ${maxns}\
        --fastaout "${QF_reads}/$(basename "$file4" .fastq).fasta" \
        --fasta_width 0

done
###################################################################################################################################################################################################
derep_dir="6.derep_data"
low_abund_seqs="7.low_abund_sequences"
SF="8.singleton_filtered"
# Enter max replicate cluster size (eg. to remove singletons enter 1, for duplicates enter 2)
maxsize="1"

echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo Removing low abundant seqs singletons per sample
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo ""

# Creating directories

mkdir ${derep_dir}

#*****************************************************************************************
# Step 1: Dereplicating

for file7 in ${QF_reads}/*.fasta
    do

        echo ""
        echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        echo Removing singletons step 1: derep_fulllength
        echo input is:
        echo ${file7}


 vsearch --derep_fulllength ${file7} \
        --strand plus \
        --output "${derep_dir}/$(basename "$file7" .fasta).fasta" \
        --sizeout \
        --uc  "${derep_dir}/$(basename "$file7" .fasta).uc" \
        --relabel $(basename "$file7" .fasta). \
        --fasta_width 0

echo Sum of unique sequences in each sample: $(cat "${derep_dir}/$(basename "$file7" .fasta).fasta" | grep -c "^>")
done

#*****************************************************************************************
# Step 2: Merge samples
# At this point there should be one fasta file for each sample                  
# It should be quality filtered and dereplicated.                               

echo
echo ====================================
echo Processing all samples together
echo ====================================

cat 6.derep_data/*.fasta > 6.derep_data/all.fasta
###################################################################################################################################################################################################
# Removal of low abundance and chimeric sequences

#*****************************************************************************************
# Step 1: Filtering low abundant seqs {maxsize}
mkdir 7.chimera_removal
echo
echo Dereplicate across samples and remove singletons
vsearch --derep_fulllength 6.derep_data/all.fasta \
    --minuniquesize 2 \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc 7.chimera_removal/all.derep.uc \
    --output 7.chimera_removal/all.derep.fasta

echo Unique non-singleton sequences: $(grep -c "^>" 7.chimera_removal/all.derep.fasta)

#*****************************************************************************************
# Step 2: Clustering prior to chimera removal
echo
echo Precluster at 98% before chimera detection

vsearch --cluster_size 7.chimera_removal/all.derep.fasta \
    --threads $threads \
    --id 0.99 \
    --strand plus \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc 7.chimera_removal/all.preclustered.uc \
    --centroids 7.chimera_removal/all.preclustered.fasta

echo Unique sequences after preclustering: $(grep -c "^>"  7.chimera_removal/all.preclustered.fasta)

#*****************************************************************************************
# Step 3: De novo chimera detection

echo
echo De novo chimera detection

vsearch --uchime_denovo 7.chimera_removal/all.preclustered.fasta \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --nonchimeras 7.chimera_removal/all.denovo.nonchimeras.fasta \

echo Unique sequences after de novo chimera detection: $(grep -c "^>" 7.chimera_removal/all.denovo.nonchimeras.fasta)

#*****************************************************************************************
# Step 4: Reference chimera detection

echo
echo Reference chimera detection
ref_db="/data/Food/analysis/R6564_NGS/amy_fitzpatrick/qiime_NoV_classifier/capsid_db85.fasta"

vsearch --uchime_ref 7.chimera_removal/all.denovo.nonchimeras.fasta \
    --threads $threads \
    --db $ref_db \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --nonchimeras 7.chimera_removal/all.ref.nonchimeras.fasta


echo
echo Extract all non-chimeric, non-singleton sequences, dereplicated

cd 7.chimera_removal
perl ../map.pl all.derep.fasta all.preclustered.uc all.ref.nonchimeras.fasta > all.nonchimeras.derep.fasta

echo Unique non-chimeric, non-singleton sequences: $(grep -c "^>" all.nonchimeras.derep.fasta)

echo
echo Extract all non-chimeric, non-singleton sequences in each sample

perl ../map.pl ../6.derep_data/all.fasta all.derep.uc all.nonchimeras.derep.fasta > all.nonchimeras.fasta

echo Sum of unique non-chimeric, non-singleton sequences in each sample: $(grep -c "^>" all.nonchimeras.fasta)

###################################################################################################################################################################################################
echo
echo Cluster at 97% and relabel with OTU_n, generate OTU table
cd ..

mkdir 8.cluster
vsearch --cluster_size 7.chimera_removal/all.nonchimeras.fasta \
    --threads $threads\
    --id 0.99 \
    --strand plus \
    --sizein \
    --fasta_width 0 \
    --relabel OTU_ \
    --uc 8.cluster/all.clustered.uc \
    --centroids 8.cluster/all.otus.fasta \
    --otutabout 8.cluster/all.otutab.txt \
    --biomout 8.cluster/all.otutab.biom


echo
echo Number of OTUs: $(grep -c "^>" 8.cluster/all.otus.fasta)

echo
echo Done

date

module unload vsearch/2.3.4 
###################################################################################################################################################################################################
##################################################################################################################################################################################################
# format biom file for input
biom convert -i 8.cluster/all.otutab.biom -o otu_table.txt --to-tsv
module load R/4.0.2
Rscript FROGS_biom_format.R
module unload R/4.0.2
biom convert -i input_otu.tsv -o input_otu.biom --table-type="OTU table" --to-json
