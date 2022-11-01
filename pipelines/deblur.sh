# create folder for data generated using a q15 phred threshold
mkdir deblur
cp metadata/updated_metadata.tsv deblur/metadata1.tsv
cd deblur

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
#*****************************************************************************************

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
################################################################################################################################
# make directory to store output file
mkdir output 

# load qiime on HPC

qiime vsearch join-pairs \
  --i-demultiplexed-seqs demultiplexed-seqs.qza \
  --p-minovlen 50 \
  --p-allowmergestagger \
  --p-minlen 100 \
  --o-joined-sequences output/demux-joined.qza

qiime quality-filter q-score \
  --i-demux output/demux-joined.qza \
  --p-min-quality 20 \
  --o-filtered-sequences output/demux-filtered.qza \
  --o-filter-stats output/demux-filter-stats.qza

qiime deblur denoise-other \
  --i-demultiplexed-seqs output/demux-filtered.qza \
  --i-reference-seqs NoV_seqs.qza \
  --p-trim-length 260 \
  --p-min-reads 1 \
  --p-mean-error 0.05 \
  --p-sample-stats \
  --o-representative-sequences output/rep-seqs.qza \
  --o-table output/table.qza \
  --o-stats output/deblur-stats.qza

################################################################################################################################
qiime deblur visualize-stats \
  --i-deblur-stats output/deblur-stats.qza \
  --o-visualization deblur-stats.qzv

qiime feature-table summarize \
  --i-table output/table.qza \
  --o-visualization output/table.qzv \
  --m-sample-metadata-file metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data output/rep-seqs.qza \
  --o-visualization output/rep-seqs.qzv

conda deactivate
module unload qiime2
