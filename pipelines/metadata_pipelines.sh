cp metadata/updated_metadata.tsv unoise3/metadata1.tsv
cp metadata/updated_metadata.tsv vsearch/metadata1.tsv
cp metadata/updated_metadata.tsv deblur/metadata1.tsv
cp metadata/updated_metadata.tsv dada2/metadata1.tsv
cp metadata/updated_metadata.tsv FROGS/metadata1.tsv

for j in vsearch unoise3 FROGS deblur dada2; do 
 (cd $j &&
	pipeline_name=${PWD##*/}
	simulation=$(echo ${PWD} | sed -e "s/.*\/\([^/]*\)\/[^/]*/\1/" | awk -F'[_.]' '{print $4}')
	awk -v var="$pipeline_name" 'BEGIN{FS=OFS="\t"}NR>1{$1=$1 "." var}1' metadata1.tsv > metadata.tsv)
  done

