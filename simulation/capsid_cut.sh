mkdir capsid_sequences
cd capsid_sequences

module load seqkit/1.4 
#cut amplicons from full length genomes
seqkit amplicon -F CNTGGGAGGGCGATCGCAA -R CCRCCNGCATRNCCRTTRTACAT -m 7 ../test*.fasta > capsid_GII.fasta
seqkit seq -M 400 capsid_GII.fasta  > capsid_GII_400.fasta 
seqkit amplicon -F CTGCCCGAATTYGTAAATGA -R CCAACCCARCCATTRTACA -m 5 ../test*.fasta > capsid_GI.fasta
seqkit seq -M 400 capsid_GI.fasta > capsid_GI_400.fasta 

# check average length of sequences
awk '{/>/&&++a||b+=length()}END{print b/a}' capsid_GII_400.fasta
awk '{/>/&&++a||b+=length()}END{print b/a}' capsid_GI_400.fasta

# add adapter sequences to samples
cat capsid_GII_400.fasta | seqkit mutate -i 0:TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCNTGGGAGGGCGATCGCAA > capsid_GII_a.fasta
cat capsid_GII_a.fasta |seqkit mutate -i -1:TACATNTTNCCNNTACGNCCNCCGACAGAGAATATGTGTAGAGGCTCGGGTGCTCTG > capsid_GII_b.fasta
cat capsid_GI_400.fasta |seqkit mutate -i 0:TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCTGCCCGAATTNGTAAATGA > capsid_GI_a.fasta
cat capsid_GI_a.fasta | seqkit mutate -i -1:ACATNTTACCNACCCAACCGACAGAGAATATGTGTAGAGGCTCGGGTGCTCTG  > capsid_GI_b.fasta
module unload seqkit/1.4 

# check number of sequences
grep ">" capsid_GII_b.fasta | wc -l
grep ">" capsid_GI_b.fasta | wc -l

# common non norovirus and norovirus sequences to be tested in classifiers
cat capsid_GI_b.fasta capsid_GII_b.fasta > capsid.fasta

# check average length of sequences
awk '{/>/&&++a||b+=length()}END{print b/a}' capsid.fasta

# remove empty fasta lines
#https://www.biostars.org/p/78786/
awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' capsid.fasta > capsid1.fasta

rm capsid_GII_b.fasta
rm capsid_GI_b.fasta

module load seqkit/1.4 
seqkit seq -M 500 capsid1.fasta > capsid2.fasta
seqkit seq -m 400 capsid2.fasta > capsid3.fasta
seqkit split --by-id capsid3.fasta
module unload seqkit/1.4

rm capsid.fasta 
rm capsid2.fasta 
