mkdir cutadapt
module load cutadapt/2.6

#Directories
DATADIR='raw_sequences'
Output='cutadapt'
adapter_fw='TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'
adapter_rev='GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'
PrimerF_GI='CTGCCCGAATTNGTAAATGA'
PrimerF_GII='CNTGGGAGGGCGATCGCAA'
PrimerR_GI='CCAACCCANCCATTNTACA'
PrimerR_GII='CCNCCNGCATNNCCNTTNTACAT'
mkdir $Output/removed_primers

for i in $DATADIR/*.*_L001_R1_001.fastq.gz
do
R1=$i;

#replace "_R1_" with "_R2_" to get R2 file name
R2=${R1%_L001_R1_001.fastq.gz}_L001_R2_001.fastq.gz;

SAMPLE1=`basename ${R1%.fastq.gz}`;
SAMPLE2=`basename ${R2%.fastq.gz}`;

echo "processing $fqfileR1"
cutadapt \
-a $adapter_fw \
-A $adapter_rev \
-b $PrimerF_GI \
-b $PrimerF_GII \
-B $PrimerR_GI \
-B $PrimerR_GII \
-e 0.1 \
--minimum-length 100 \
--discard-untrimmed \
-o $Output/removed_primers/$SAMPLE1.fastq.gz --paired-output $Output/removed_primers/$SAMPLE2.fastq.gz $R1 $R2

done

module unload cutadapt/2.6

