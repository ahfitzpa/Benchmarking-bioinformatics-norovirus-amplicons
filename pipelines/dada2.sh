
################################################################################################################################
# create folder for data generated using a q20 phred threshold
mkdir dada2
cp metadata/updated_metadata.tsv dada2/metadata1.tsv
cd dada2

pipeline_name=${PWD##*/}
simulation=$(echo ${PWD} | sed -e "s/.*\/\([^/]*\)\/[^/]*/\1/" | awk -F'[_.]' '{print $4}')
awk -v var="$pipeline_name" 'BEGIN{FS=OFS="\t"}NR>1{$1=$1 "." var}1' metadata1.tsv > metadata.tsv
rm metadata1.tsv
date

mkdir raw_data
cp ../cutadapt/removed_primers/*.fastq.gz raw_data/

cd raw_data

#*****************************************************************************************
for f in *_L001_R2_001.fastq.gz; do
    mv -- "$f" "${f%_L001_R2_001.fastq.gz}.${pipeline_name}_L001_R2_001.fastq.gz"
done

for f in *_L001_R1_001.fastq.gz; do
    mv -- "$f" "${f%_L001_R1_001.fastq.gz}.${pipeline_name}_L001_R1_001.fastq.gz"
done
cd ..

# load qiime on HPC
module load qiime2/2021.2
source activate qiime2-2021.2

## import sequence reads. Sequences need to have been subjected to cutadapt prior to import into dada2 as they have been generated using two different primers ##
## CasavaOneEightSingleLanePerSampleDirFmt refers to naming convention. Check this if errors encountered ##
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path raw_data/ \
  --input-format CasavaOneEightLanelessPerSampleDirFmt \
  --output-path demultiplexed-seqs.qza

#rm -r -d raw_data

# make directory to store output file
mkdir output 

# load qiime on HPC
module load qiime2/2021.2
source activate qiime2-2021.2

## denoise using dada2. Q20: Reads truncated at the first instance of a quality score less than or equal to this value ##
## sequences stored in rep-seqs ##
## statistics of quality metrics (merged reads, denoised, non chimeric etc) stored in stats-dada2.qza ##
## dada2.qza file later converted to qzv, shows frequency of features ##
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demultiplexed-seqs.qza \
  --p-trunc-len-f 0 \
  --p-trunc-len-r 0 \
  --p-trim-left-f 3 \
  --p-trim-left-r 5 \
  --p-max-ee-f 2 \
  --p-max-ee-r 2 \
  --p-trunc-q 20 \
  --p-chimera-method none \
  --o-representative-sequences output/rep-seqs.qza \
  --o-table output/dada2_table.qza \
  --o-denoising-stats output/stats-dada2.qza \
  --verbose 

## next few lines convert to qza files to visual qzv files  ##

qiime metadata tabulate \
  --m-input-file output/stats-dada2.qza \
  --o-visualization output/stats-dada2.qzv

qiime feature-table tabulate-seqs \
  --i-data output/rep-seqs.qza \
  --o-visualization output/rep-seqs.qzv

qiime feature-table summarize \
  --i-table output/dada2_table.qza \
  --m-sample-metadata-file metadata.tsv \
  --o-visualization output/dada2_table.qzv 

module unload qiime2/2021.2
conda deactivate
