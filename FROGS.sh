mkdir FROGS
mkdir FROGS/output
mkdir FROGS/logs
mkdir FROGS/raw_data
###################################################################################################################################################################################################
# prepare directories and files
cp cutadapt/removed_primers/*.fastq.gz FROGS/raw_data/
cp metadata/updated_metadata.tsv FROGS/metadata1.tsv
cp metadata/expected_data.tsv FROGS/expected_composition.tsv
cp metdata/expected.fasta FROGS/expected_sequences.fasta
cp FROGS_biom_format.R FROGS/output
gunzip FROGS/raw_data/*.gz
cp frogs_db/capsid_db* FROGS/

cd FROGS
pipeline_name=${PWD##*/}
simulation=$(echo ${PWD} | sed -e "s/.*\/\([^/]*\)\/[^/]*/\1/" | awk -F'[_.]' '{print $4}')
awk -v var="$pipeline_name" 'BEGIN{FS=OFS="\t"}NR>1{$1=$1 "." var}1' metadata1.tsv > metadata.tsv
rm metadata1.tsv

cd raw_data

#*****************************************************************************************
for f in *_L001_R2_001.fastq; do
    mv -- "$f" "${f%.${simulation}_L001_R2_001.fastq}_${simulation}_${pipeline_name}_L001_R2_001.fastq"
done

for f in *_L001_R1_001.fastq; do
   mv -- "$f" "${f%.${simulation}_L001_R1_001.fastq}_${simulation}_${pipeline_name}_L001_R1_001.fastq"
done

for i in *_L001_R1_001.fastq
  do
  name=$(basename "$i" | sed -e 's/..................$//')
  mv $i ${name}_R1.fastq;
done

for i in *_L001_R2_001.fastq
  do
  name=$(basename "$i" | sed -e 's/..................$//')
  mv $i ${name}_R2.fastq;
done

(basename -a *R1.fastq | sed 's/.........$//') > filenames

samplenames="cat filenames"
$samplenames
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

for file1 in ${raw_data}/*R1.fastq
    do

        echo ""
        echo ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        echo Merging paired reads
        echo forward reads are:
        echo $(basename ${file1})
        echo reverse reads are:
        echo $(basename ${file1} R1.fastq)R2.fastq

    vsearch --fastq_mergepairs ${file1} --threads ${threads} --reverse "${raw_data}/$(basename -s R1.fastq ${file1})R2.fastq" --fastq_maxdiffs ${maxdiffs} --fastq_minovlen ${overlap} --fastqout "working1/$(basename "$file1")" --fastq_eeout 
done


#*****************************************************************************************
# Rename and move
for file2 in working1/*R1.fastq
  do
  mv ${file2} ${merged_data}/$(basename "${file2}" |  sed -e 's/.........$//').fastq
done


#*****************************************************************************************
# Removing working directory

        rm -r working1


module unload vsearch/2.3.4
#######################################################################################################################################################################################################################
echo "Step preprocess : Vsearch `date`"
R1_FILES="2.merged_reads/*.fastq"
$R1_FILES
cpu="5"

module load frogs/3.2.3
source activate frogs@3.2.3

preprocess.py illumina \
 --input-R1 ${R1_FILES} \
 --R1-size 300 \
 --already-contiged \
 --min-amplicon-size 100 --max-amplicon-size 344 \
 --without-primers \
 --nb-cpus ${cpu}  \
 --expected-amplicon-size 344 \
 --output-dereplicated output/01-prepro-vsearch.fasta \
 --output-count output/01-prepro-vsearch.tsv \
 --summary output/01-prepro-vsearch.html \
 --log-file logs/01-prepro-vsearch.log 
 
if [ $? -ne 0 ]
then
	echo "Error in preprocess : Vsearch" >&2
	exit 1;
fi

#######################################################################################################################################################################################################################
echo "Step clustering1 `date`"

clustering.py \
 --distance 1 --denoising \
 --input-fasta output/01-prepro-vsearch.fasta \
 --input-count output/01-prepro-vsearch.tsv \
 --output-b output/02-clustering1.biom \
 --output-fasta output/02-clustering1.fasta \
 --output-compo output/02-clustering1_compo.tsv \
 --log-file logs/02-clustering1.log \
 --nb-cpus 5
if [ $? -ne 0 ]
then
	echo "Error in clustering" >&2
	exit 1;
fi

echo "Step remove_chimera `date`"

remove_chimera.py \
 --input-fasta output/02-clustering1.fasta \
 --input-biom output/02-clustering1.biom \
 --non-chimera output/03-chimera.fasta \
 --out-abundance output/03-chimera.biom \
 --summary output/03-chimera.html \
 --log-file output/03-chimera.log \
 --nb-cpus 5
if [ $? -ne 0 ]
then
	echo "Error in remove_chimera" >&2
	exit 1;
fi


echo "Step otu filters `date`"

otu_filters.py \
 --min-abundance 0.00005 \
 --min-sample-presence 1 \
 --contaminant /data/Food/analysis/R6564_NGS/amy_fitzpatrick/denoise_comparison/phi.fa \
 --nb-cpus 5 \
 --input-biom output/03-chimera.biom \
 --input-fasta output/03-chimera.fasta \
 --output-fasta output/04-filters.fasta \
 --output-biom output/04-filters.biom \
 --excluded output/04-filters.excluded \
 --summary output/04-filters.html \
 --log-file logs/04-filters.log 

if [ $? -ne 0 ]
then
	echo "Error in filters" >&2
	exit 1;
fi

echo "Step affiliation_OTU `date`"

affiliation_OTU.py \
 --reference capsid_db.fasta \
 --input-fasta output/04-filters.fasta \
 --input-biom output/04-filters.biom \
 --output-biom output/06-affiliation.biom \
 --summary output/06-affiliation.html \
 --log-file logs/06-affiliation.log \
 --nb-cpus 5  

echo "Step affiliation_filter: masking mode `date`"

affiliation_filters.py \
--input-biom output/06-affiliation.biom \
--input-fasta output/04-filters.fasta \
--output-biom output/07-affiliation_masked.biom \
--summary output/07-affiliation_masked.html \
--impacted output/07-impacted_OTU_masked.tsv \
--impacted-multihit output/07-impacted_OTU_masked_multihit.tsv \
--log-file logs/07-affiliation_filter_maskMode.log \
--min-blast-length 100 \
--min-blast-identity 0.7 \
--min-blast-coverage 0.5 \
--max-blast-evalue 1e-150 \
--mask \
--taxonomic-ranks Domain Phylum Class Order Family Genus Species  # this is the default value of this option

echo "Step affiliation_postprocess `date`"

affiliation_postprocess.py \
 --input-biom output/06-affiliation.biom \
 --input-fasta output/04-filters.fasta \
 --reference capsid_db.fasta  \
 --output-biom output/08-affiliation_postprocessed.biom \
 --output-compo output/08-affiliation_postprocessed.compo.tsv \
 --output-fasta output/08-affiliation_postprocessed.fasta \
 --log-file logs/08-affiliation_postprocessed.log
 
echo "Step clusters_stat `date`"

clusters_stat.py \
 --input-biom output/08-affiliation_postprocessed.biom  \
 --output-file output/10-clustersStat.html \
 --log-file logs/10-clustersStat.log
					
echo "Step affiliations_stat `date`"

affiliations_stat.py \
 --input-biom output/08-affiliation_postprocessed.biom  \
 --output-file output/11-affiliationsStat.html \
 --log-file output/11-affiliationsStat.log \
 --tax-consensus-tag "blast_taxonomy" \
 --identity-tag "perc_identity" \
 --coverage-tag "perc_query_coverage" \
 --multiple-tag "blast_affiliations" \
 --rarefaction-ranks Family Genus Species \
 --taxonomic-ranks Domain Phylum Class Order Family Genus Species # this is the default value of this option

echo "Step biom_to_tsv `date`"

biom_to_tsv.py \
 --input-biom output/08-affiliation_postprocessed.biom\
 --input-fasta output/08-affiliation_postprocessed.fasta\
 --output-tsv output/12-biom2tsv.tsv \
 --output-multi-affi output/12-biom2tsv-affiliation_multihit.tsv \
 --log-file output/12-biom2tsv.log

echo "Step biom_to_stdBiom `date`"


biom_to_stdBiom.py \
 --input-biom output/08-affiliation_postprocessed.biom \
 --output-biom output/13-affiliation_std.biom \
 --output-metadata output/13-affiliation_multihit.tsv \
 --log-file output/13-biom2stdbiom.log

echo "Step tsv_to_biom `date`"


tsv_to_biom.py \
 --input-tsv output/12-biom2tsv.tsv \
 --input-multi-affi output/12-biom2tsv-affiliation_multihit.tsv \
 --output-biom output/14-tsv2biom.biom \
 --output-fasta output/14-tsv2biom.fasta \
 --log-file output/14-tsv2biom.log 

echo "Step tree `date`"

tree.py \
 --nb-cpus 3 \
 --input-sequences output/08-affiliation_postprocessed.fasta \
 --biom-file output/08-affiliation_postprocessed.biom \
 --out-tree output/15-tree-mafft.nwk \
 --html output/15-tree-mafft.html \
 --log-file output/15-tree-mafft.log

echo "Step phyloseq_import_data `date`"

phyloseq_import_data.py  \
 --biomfile output/08-affiliation_postprocessed.biom \
 --samplefile metadata.tsv \
 --treefile output/15-tree-mafft.nwk \
 --rdata output/16-phylo_import.Rdata \
 --html output/16-phylo_import.nb.html \
 --log-file output/16-phylo_import.log

echo "Step phyloseq_composition `date`"

phyloseq_composition.py  \
 --varExp n_genotype --taxaRank1 Class --taxaSet1 Norovirus --taxaRank2 Species --numberOfTaxa 47 \
 --rdata output/16-phylo_import.Rdata \
 --html output/17-phylo_composition.nb.html \
 --log-file output/17-phylo_composition.log

echo "Step phyloseq_alpha_diversity `date`"

phyloseq_alpha_diversity.py  \
 --varExp n_genotype \
 --rdata output/16-phylo_import.Rdata --alpha-measures Observed Chao1 Shannon \
 --alpha-out output/18-phylo_alpha_div.tsv \
 --html output/18-phylo_alpha_div.nb.html \
 --log-file output/18-phylo_alpha_div.log

echo "Step phyloseq_beta_diversity `date`"

phyloseq_beta_diversity.py  \
 --varExp n_genotype --distance-methods cc,unifrac \
 --rdata output/16-phylo_import.Rdata \
 --matrix-outdir output \
 --html output/19-phylo_beta_div.nb.html \
 --log-file output/19-phylo_beta_div.log

 
module unload frogs/3.2.3
conda deactivate
##################################################################################################################################################################################################
# format biom file for input
cd output
biom convert -i 08-affiliation_postprocessed.biom  -o otu_table.txt --to-tsv
module load R/4.0.2
Rscript FROGS_biom_format.R
module unload R/4.0.2
biom convert -i input_otu.tsv -o input_otu.biom --table-type="OTU table" --to-json
