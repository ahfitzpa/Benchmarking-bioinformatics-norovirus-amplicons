mkdir unoise3
mkdir unoise3/output
mkdir unoise3/logs
usearch11="usearch11.0.667_i86linux32"

cp metadata/updated_metadata.tsv unoise3/metadata1.tsv

cd unoise3
pipeline_name=${PWD##*/}
simulation=$(echo ${PWD} | sed -e "s/.*\/\([^/]*\)\/[^/]*/\1/" | awk -F'[_.]' '{print $4}')
awk -v var="$pipeline_name" 'BEGIN{FS=OFS="\t"}NR>1{$1=$1 "." var}1' metadata1.tsv > metadata.tsv
rm metadata1.tsv

mkdir raw_data
cp ../cutadapt/removed_primers/*.fastq.gz raw_data/
gunzip raw_data/*.gz
cd raw_data

#*****************************************************************************************
for f in *_L001_R2_001.fastq; do
    mv -- "$f" "${f%_L001_R2_001.fastq}.${pipeline_name}_L001_R2_001.fastq"
done

for f in *_L001_R1_001.fastq; do
    mv -- "$f" "${f%_L001_R1_001.fastq}.${pipeline_name}_L001_R1_001.fastq"
done
cd ..

###################################################################################################################################################################################################
raw_data="raw_data"
merged_data="2.merged_reads"
# Maximum no of mismatches in the alignment - Default 5. Consider increasing if you have long overlaps.
maxdiffs="20"
# Discard pair if alignment is shorter than this value - Default 16 bp
overlap="50"

mkdir ${merged_data}
mkdir working1

#*****************************************************************************************
# Step1: merge data with usearch9 -fastq-filter

for file1 in ${raw_data}/*R1_001.fastq
    do

        echo ""
        echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        echo Merging paired reads
        echo forward reads are:
        echo $(basename ${file1})
        echo reverse reads are:
        echo $(basename ${file1} R1_001.fastq)R2_001.fastq

    ${usearch11} -fastq_mergepairs ${file1} -reverse "${raw_data}/$(basename -s R1_001.fastq ${file1})R2_001.fastq" -fastqout "working1/$(basename "$file1")" -fastq_maxdiffs ${maxdiffs} -fastq_minovlen ${overlap} -report ${merged_data}/2a_merging_seqs_report.txt -tabbedout ${merged_data}/2b_tabbedout.txt
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
QF_reads="3.quality_filtered"
# Discard reads with > E total expected errors for all bases in the read after any truncation options have been applied. Usually set at 1, however may want to lower to 0.5 or 0.25 for more stringent filtering.
max_ee="1"
# Enter min length of sequence for trimming in bp (eg. to keep all seqs above 200 bp enter "200")
minlen="100"

mkdir ${QF_reads}

for file4 in ${merged_data}/*.fastq
    do
        echo ""
        echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        echo Quality control and removing dimer seqs
        echo input is:
        echo ${file4}

    ${usearch11} -fastq_filter ${file4} -fastaout "${QF_reads}/$(basename "$file4" .fastq).fasta" -fastq_maxee ${max_ee} -fastq_minlen ${minlen}
done

###################################################################################################################################################################################################
labeled_data="4.labeled_data"

echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo Renameing sequences with ">barcodelabel=sample_id;sequence_id"
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mkdir ${labeled_data}
mkdir working2

#*****************************************************************************************
# Step 1: Remove ">" from start of sequence_ID

for file5 in ${QF_reads}/*.fasta
    do

    sed -e 's/>/>barcodelabel=;/g' ${file5} > working2/$(basename "$file5" .fasta).txt
done

#*****************************************************************************************
# Step 2: Add sample_ID (should be filename) to produce ">barcodelabel=sample_ID;sequence_ID"

for file6 in working2/*.txt
    do

    sample_id=$(basename ${file6} .txt)
    echo ${sample_id}

    sed -e "s/;/${sample_id};/g" ${file6} > "${labeled_data}/$(basename "$file6" .txt).fasta"
done

#*****************************************************************************************
# Remove working directories

rm -r working2

###################################################################################################################################################################################################
derep_dir="5.derep_data"
low_abund_seqs="6.low_abund_sequences"
SF="7.singleton_filtered"
# Enter max replicate cluster size (eg. to remove singletons enter 1, for duplicates enter 2)
maxsize="1"

echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo Removing low abundant seqs singletons per sample
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo ""

# Creating directories

mkdir ${derep_dir}
mkdir ${low_abund_seqs}
mkdir ${SF}

#*****************************************************************************************
# Step 1: Dereplicating

for file7 in ${labeled_data}/*.fasta
    do

        echo ""
        echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        echo Removing singletons step 1: derep_fulllength
        echo input is:
        echo ${file7}

     ${usearch11} -fastx_uniques ${file7} -fastaout "${derep_dir}/$(basename "$file7" .fasta).fasta" -sizeout
done


#*****************************************************************************************
# Step 2: Filtering low abundant seqs {maxsize}

for file8 in ${derep_dir}/*.fasta
    do

        echo ""
        echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        echo Removing singletons step 2: sorting uniques
        echo input is:
        echo ${file8}

     ${usearch11} -sortbysize ${file8} -fastaout "${low_abund_seqs}/$(basename "$file8" .fasta).fasta" -maxsize ${maxsize}
done

#*****************************************************************************************
# Step 3: Mapping reads

for file9 in ${labeled_data}/*.fasta
    do

        echo ""
        echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        echo Removing singletons step 3: mapping reads to low abundant uniques
        echo input is:
        echo ${file9}

     ${usearch11} -search_exact ${file9} -db "${low_abund_seqs}/$(basename "$file9" .fasta).fasta" -strand plus -notmatched "${SF}/$(basename "$file9" .fasta).fasta"
done


*****************************************************************************************
###########################################################################################################################################################
SF="7.singleton_filtered"
cluster="8.cluster"
uparse_otus="8a.uparse_otus"

echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo UPARSE on all
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mkdir ${cluster}

   cat ${SF}/*.fasta > ${cluster}/all_SF.fasta

cd ${cluster}

  ${usearch11} -fastx_uniques all_SF.fasta -fastaout all_SF_DR.fasta -sizeout

mkdir ${uparse_otus}
cd ${uparse_otus}

  ${usearch11} -cluster_otus ../all_SF_DR.fasta -otus uparse_otus.fasta -relabel OTU

  ${usearch11} -usearch_global ../all_SF.fasta -db uparse_otus.fasta -strand both -id 0.97 -otutabout uparse_otu_tab.txt -biomout uparse_otu_biom.biom

  ${usearch11} -calc_distmx uparse_otus.fasta -tabbedout uparse_otus_distmx.txt -maxdist 0.2 -termdist 0.3
  
  ${usearch11} -cluster_aggd uparse_otus_distmx.txt -treeout uparse_otus_clusters.tree -clusterout uparse_otus_clusters.txt \
      -id 0.80 -linkage min

cd ..
##################################################################################################################################################################################################
cluster="8.cluster"
unoise_zotus="8b.unoise_zotus"

echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo UNOISE on all
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mkdir ${unoise_zotus}
cd ${unoise_zotus}

  ${usearch11} -unoise3 ../all_SF_DR.fasta -zotus unoise_zotus.fasta -tabbedout unoise_tab.txt

  ${usearch11} -fastx_relabel unoise_zotus.fasta -prefix Otu -fastaout unoise_zotus_relabelled.fasta -keep_annots

  ${usearch11} -otutab ../all_SF.fasta -zotus unoise_zotus_relabelled.fasta -otutabout unoise_otu_tab.txt -biomout unoise_otu.biom -mapout unoise_map.txt -notmatched unoise_notmatched.fasta -dbmatched dbmatches.fasta -sizeout

  ${usearch11} -calc_distmx unoise_zotus.fasta -tabbedout unoise_zotus_distmx.txt -maxdist 0.2 -termdist 0.3
  
  ${usearch11} -cluster_aggd unoise_zotus_distmx.txt -treeout unoise_zotus_clusters.tree -clusterout unoise_zotus_clusters.txt \
      -id 0.80 -linkage min

cd ..
cd ..

##################################################################################################################################################################################################
#Load QIIME2
module load qiime2/2021.2
source activate qiime2-2021.2

#Convert zOTU table and rep-seqs to q2 artifacts
mkdir q2
mkdir q2/quality-control
mkdir q2/quality-control
cp ../metadata/expected_data.tsv q2/quality-control/expected_composition.tsv
cp ../metadata/expected_sequences3.fasta q2/quality-control/expected_sequences.fasta

qiime tools import \
--input-path 8.cluster/8b.unoise_zotus/unoise_zotus_relabelled.fasta \
--output-path q2/rep-seqs.qza \
--type "FeatureData[Sequence]"

biom convert -i 8.cluster/8b.unoise_zotus/unoise_otu_tab.txt -o q2/table.from_txt_json.biom --table-type="OTU table" --to-json

module unload qiime2/2021.2
conda deactivate
