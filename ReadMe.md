
Downloading human genome in colorspace fasta from Bowtie site:
wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/hg19_c.ebwt.zip
unzip hg19_c.ebwt.zip -d hg19

Downloading bacterial metagenomes in fasta from IGC site:
wget ftp://climb.genomics.cn/pub/10.5524/100001_101000/100064/1.GeneCatalogs/IGC.fa.gz
gunzip IGC.fa.gz

Downloading bacterial metagenomes annotation from IGC site:
ftp://climb.genomics.cn/pub/10.5524/100001_101000/100064/3.IGC.AnnotationInfo/IGC.annotation_OF.summary.gz
gunzip IGC.annotation_OF.summary.gz

Or make your own annotation with FASTA headers (if you're using GenBank-formatted reference): 
grep "^>" file.fasta | sed -e 's/|/\t/g' | sed -e 's/>gb\t//g' > annotation.txt

Add a header from file (note the delimiter and the line separator):
"Gene_ID	Gene_name	Gene_length	Gene_completeness_status	Cohort_origin	Phylum	Genus	KEGG	eggNOG	Sample_occurence_frequency	Individual_occurence_frequency	KEGG_function	eggNOG_function	Cohort_assembled
"
cat header.txt IGC.annotation.summary.v2 > IGC.annotation.summary.v3

FASTA chunking
I've used Fasta File Splitter tool under Windows:
https://github.com/PNNL-Comp-Mass-Spec/Fasta-File-Splitter

Its command line was like:
FastaFileSplitter.exe 760MetaHit_139HMP_368PKU_511Bac.fa.90_95 /N:2
pause

Make IGC genome colorspace indexes:

bowtie-build -C 760MetaHit_139HMP_368PKU_511Bac.fa_2x_1.fasta igc1 > bowtie-build_1.log
bowtie-build -C 760MetaHit_139HMP_368PKU_511Bac.fa_2x_2.fasta igc2 > bowtie-build_2.log

samtools faidx 760MetaHit_139HMP_368PKU_511Bac.fa_2x_1.fasta
samtools faidx 760MetaHit_139HMP_368PKU_511Bac.fa_2x_2.fasta

awk -v OFS='\t' {'print $1,$2'} 760MetaHit_139HMP_368PKU_511Bac.fa_2x_1.fasta.fai > 760MetaHit_139HMP_368PKU_511Bac.fa_2x_1.genome
awk -v OFS='\t' {'print $1,$2'} 760MetaHit_139HMP_368PKU_511Bac.fa_2x_2.fasta.fai > 760MetaHit_139HMP_368PKU_511Bac.fa_2x_2.genome
cat "760MetaHit_139HMP_368PKU_511Bac.fa_2x_1.genome" "760MetaHit_139HMP_368PKU_511Bac.fa_2x_2.genome" > "760MetaHit_139HMP_368PKU_511Bac.fa.90_95.genome"

You can use *.genome indexes as your annotation for the fail-safe processing. You can also use sed to add a header:
sed -i '1s/^/Gene_name\tGene_length\tGene_ID\n/' annotation.txt

If you want add a file as an extra column:
paste file1 file2 > file3
