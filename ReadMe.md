# ThreeBees

**ThreeBees** is an attempt to make a functional seamless pipeline for symbiotic metagenome data analysis originally developed for optimization of dataflow performed above sequenced human gut microbiomes. It uses 3 consecutive bowtie reads alignments (meant as "B"). 

>**WARNINGS/DISCLAIMERS**

> This pipeline is published for educational purposes only. Its code is presented "as is" and might not work in some cases. Make sure you have made all necessary changes in configuration parts!

### How does it work

 1. Pre-configured *ThreeBeesGenerator.py* makes *ThreeBees.sh* script in the folder specified after -f flag. This folder must contain *"LXX"* subdirectories, and **.xsq* file must have path like *"LXX/result/"*;
 2. *ThreeBees.sh* runs main data processing including the next step;
 3. *ThreeBeesCombiner.py* collects all required coverage data into *bp.txt* or *pos.txt* files created in *"ThreeBees/5_Statistics/Finalized/"* path. *"bp.txt"* contains total coverage and *"pos.txt"* contains numbers of covered nucleotides per gene (entry). 

### Data processing

 1. ABI SOLiD Reads preprocessing *(XSQ Tools Software, reads-filter_new.pl)*;
 2. Colorspace FASTA Reads processing *(SOLiD Accuracy Enhancer Tool (SAET), csfasta_quality_filter.pl)*;
 3. Negative genome mapping (*bowtie* with the *-un* flag);
 4. Unmapped reads positive mapping (*bowtie* without the *-un* flag);
 5. Statistics gathering *(SAMtools, BEDTools, get_cov_from_bedtools_hist.py)*

### A small preparation guide

 - You will need a Linux-based mainframe. We've used the computational cluster with 24 CPUs, 64 Gb RAM, 30 Tb RAID storage running on CentOS;
 - Install all mentioned software. You have to also check *perl* and *python3* with *pandas* package;
 - Download and build sequence data (see below);
 - Paste all paths into SoftwarePaths and SequencesPaths sections of *ThreeBeesGenerator.py* file.
 - Launch commands:
```
cd <full path to folder with "LXX" subdirectories>  # Navigation into main working  directory
python ~/scripts/ThreeBees/ThreeBeesGenerator.py -f <full path to folder with "LXX" subdirectories>  # Bash script generation
nohup bash ThreeBees.sh > /dev/null 2>&1 & echo $! > run.pid  # Launch bash script with free console and save master control process ID, killing it will stop pipeline
top  # Watch the progress, press "q" to exit
```

#### Making References: How We Did It

 - Download human *hg19* genome in colorspace fasta from Bowtie site:
```
wget ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/hg19_c.ebwt.zip
unzip hg19_c.ebwt.zip -d hg19
```
    
 - Download bacterial metagenomes in fasta from IGC site:
```
wget ftp://climb.genomics.cn/pub/10.5524/100001_101000/100064/1.GeneCatalogs/IGC.fa.gz
gunzip IGC.fa.gz
```

 - Download bacterial metagenomes annotation from IGC site:
```
ftp://climb.genomics.cn/pub/10.5524/100001_101000/100064/3.IGC.AnnotationInfo/IGC.annotation_OF.summary.gz
gunzip IGC.annotation_OF.summary.gz`
```

 - Or make own annotation with FASTA headers (if using GenBank-formatted reference): 

```
grep "^>" file.fasta | sed -e 's/|/\t/g' | sed -e 's/>gb\t//g' > annotation.txt
```

 - Add a header from file (note the delimiter and the line separator): 
 > Gene_ID	Gene_name	Gene_length	Gene_completeness_status	Cohort_origin	Phylum	Genus	KEGG	eggNOG	Sample_occurence_frequency	Individual_occurence_frequency	KEGG_function	eggNOG_function	Cohort_assembled
 
```
cat header.txt IGC.annotation.summary.v2 > IGC.annotation.summary.v3
```

 - For FASTA chunking we've used [Fasta File Splitter](https://github.com/PNNL-Comp-Mass-Spec/Fasta-File-Splitter) tool under Windows. Its command line was like:

```
FastaFileSplitter.exe 760MetaHit_139HMP_368PKU_511Bac.fa.90_95 /N:2
pause
```

 - Make IGC genome colorspace indexes:

```
bowtie-build -C 760MetaHit_139HMP_368PKU_511Bac.fa_2x_1.fasta igc1 > bowtie-build_1.log

bowtie-build -C 760MetaHit_139HMP_368PKU_511Bac.fa_2x_2.fasta igc2 > bowtie-build_2.log

samtools faidx 760MetaHit_139HMP_368PKU_511Bac.fa_2x_1.fasta

samtools faidx 760MetaHit_139HMP_368PKU_511Bac.fa_2x_2.fasta

awk -v OFS='\t' {'print $1,$2'} 760MetaHit_139HMP_368PKU_511Bac.fa_2x_1.fasta.fai > 760MetaHit_139HMP_368PKU_511Bac.fa_2x_1.genome

awk -v OFS='\t' {'print $1,$2'} 760MetaHit_139HMP_368PKU_511Bac.fa_2x_2.fasta.fai > 760MetaHit_139HMP_368PKU_511Bac.fa_2x_2.genome

cat "760MetaHit_139HMP_368PKU_511Bac.fa_2x_1.genome" "760MetaHit_139HMP_368PKU_511Bac.fa_2x_2.genome" > "760MetaHit_139HMP_368PKU_511Bac.fa.90_95.genome"

```

 - Add a header to reserve **.genome* indexes as annotation for the fail-safe processing (we did not use it):

```
sed -i '1s/^/Gene_name\tGene_length\tGene_ID\n/' annotation.txt
```

 - Add a file as an extra column:

```
paste file1 file2 > file3
```
