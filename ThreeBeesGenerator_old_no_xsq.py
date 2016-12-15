# Paths should be ending at "/". No spaces are allowed.
# Make sure the script is placed into root lanes directory containing "LXX" directories.

# [SoftwarePaths]
qualityFilterPath = ""
SAETPath = ""
qualityTrimmerPath = ""
BowtieFolder = ""

# [SequencesPaths]
humanGenomePath = "/data/reference/homo_sapiens/ucsc/hg19/sequence/ColorSpaseIndex/genome"

bacterialGenomeBwtPath1 = "/data/reference/IGC/fs/igc1"
bacterialGenomeBwtPath2 = "/data/reference/IGC/fs/igc2"

bacterialGenomeFai1 = "/data/reference/IGC/fs/760MetaHit_139HMP_368PKU_511Bac.fa_2x_1.fasta.fai"
bacterialGenomeFai2 = "/data/reference/IGC/fs/760MetaHit_139HMP_368PKU_511Bac.fa_2x_2.fasta.fai"

bacterialGenomeFasta1 = "/data/reference/IGC/fs/760MetaHit_139HMP_368PKU_511Bac.fa_2x_1.fasta"
bacterialGenomeFasta2 = "/data/reference/IGC/fs/760MetaHit_139HMP_368PKU_511Bac.fa_2x_2.fasta"

# Script begin
import inspect, os, re, csv

# root folder:  /data1/bio/kazan_solid_metagenome_data/kazan5500_2016_05_31_HP_A
# lane folder:  /data1/bio/kazan_solid_metagenome_data/kazan5500_2016_05_31_HP_A/L01_conversion/Libraries/

# demultiplexed samples folder:
#               /data1/bio/kazan_solid_metagenome_data/kazan5500_2016_05_31_HP_A/L01_conversion

# csfasta file: /data1/bio/kazan_solid_metagenome_data/kazan5500_2016_05_31_HP_A/L01_conversion/Libraries/130HP/F3/reads/
#                            kazan5500_2016_05_31_HP_A_L01_130HP_F3.csfasta

# Workflow folders creation
for newDir in ["1_Preprocessed_reads/filtered/fixed", "2_Processed_reads", "3_Non-mapped_reads", "4_Mapped_reads", "5_Statistics"]:
    os.makedirs(str("ThreeBees/" + newDir))

# Lane proessing
threeBees = open("ThreeBees.sh", 'w')
threeBees.write(str("LaneDir=$PWD\n" +
                    "ScriptDir=$LaneDir/ThreeBees\n"))

csfastaPaths = []

for lane in ["L01", "L02", "L03", "L04", "L05", "L06"]:
    multiplexFilenames = []

    for sampleName in os.listdir(str(lane + "_conversion/Libraries/")):
        if not any(badWord in sampleName for badWord in ["##Library", "Unassigned", "Unclassified", "All"]):
            multiplexFilenames.append(sampleName)

    for sampleName in multiplexFilenames:
        for file in os.listdir(str(lane + "_conversion/Libraries/" + sampleName + "/F3/reads/")):
            if file.endswith(".csfasta"):
                csfastaPaths.append(os.path.abspath(str(lane + "_conversion/Libraries/" + sampleName + "/F3/reads/" + file)))

for csfastaPath in csfastaPaths:
    csfastaFileName = str(csfastaPath.rsplit('/', 1)[-1]).replace(".csfasta", "")

    threeBees.write("# FOR SAMPLE FILE: " + csfastaFileName + "\n" +
                    "cd $ScriptDir/1_Preprocessed_reads/" + "\n" +
                    qualityFilterPath + "reads-filter_new.pl -f " + csfastaPath + " -q " + str(csfastaPath).replace(".csfasta", ".QV.qual") + " -o filtered -t 15 &>$ScriptDir/5_Statistics/" + csfastaFileName + "_reads-filter_new.log" + "\n" +
                    "\n" +
                    "cd filtered" + "\n" +
                    SAETPath + "saet " + csfastaFileName + ".15.filtered.csfasta 200000000 -qual " + csfastaFileName + ".QV.15.filtered.qual -qvupdate -trustprefix 25 -localrounds 3 -globalrounds 2 -numcores 20 -log $ScriptDir/5_Statistics/" + csfastaFileName + "_SAET.log" + "\n" +
                    "\n" +
                    "cd fixed" + "\n" +
                    "mv " + csfastaFileName + ".15.filtered.csfasta " + csfastaFileName + ".15.filtered.fixed.csfasta" + "\n" +
                    "mv none " + csfastaFileName + ".QV.15.filtered.fixed.qual" + "\n" +
                    "\n" +
                    qualityTrimmerPath + "csfasta_quality_filter.pl -f " + csfastaFileName + ".15.filtered.fixed.csfasta -q " + csfastaFileName + ".QV.15.filtered.fixed.qual -o $ScriptDir/2_Processed_reads/" + csfastaFileName + ".15.filtered.fixed.trimmed.csfasta -l 30 -s 30 &>$ScriptDir/5_Statistics/" + csfastaFileName + "_csfasta_quality_filter.log" + "\n" +
                    "\n" +
                    "cd $ScriptDir" +
                    "\n" +
                    "bowtie -f -C -S -t -v 3 -k 1 --threads 20 --un $ScriptDir/3_Non-mapped_reads/" + csfastaFileName + ".nohum.csfasta " + humanGenomePath + " $ScriptDir/2_Processed_reads/" + csfastaFileName + ".15.filtered.fixed.trimmed.csfasta $ScriptDir/4_Mapped_reads/" + csfastaFileName + ".hum.sam &>$ScriptDir/5_Statistics/" + csfastaFileName + "_bowtieA.log" + "\n" +
                    "\n" +
                    "bowtie -f -C -S -t -v 3 -k 1 --threads 20 --un $ScriptDir/3_Non-mapped_reads/" + csfastaFileName + ".nohum.noigc1.csfasta " + bacterialGenomeBwtPath1 + " $ScriptDir/3_Non-mapped_reads/" + csfastaFileName + ".nohum.csfasta $ScriptDir/4_Mapped_reads/" + csfastaFileName + ".nohum.igc1.sam &>$ScriptDir/5_Statistics/" + csfastaFileName + "_bowtieB1.log" + "\n" +
                    "bowtie -f -C -S -t -v 3 -k 1 --threads 20 --un $ScriptDir/3_Non-mapped_reads/" + csfastaFileName + ".nohum.noigc1.noigc2.csfasta " + bacterialGenomeBwtPath2 + " $ScriptDir/3_Non-mapped_reads/" + csfastaFileName + ".nohum.noigc1.csfasta $ScriptDir/4_Mapped_reads/" + csfastaFileName + ".nohum.noigc1.igc2.sam &>$ScriptDir/5_Statistics/" + csfastaFileName + "_bowtieB2.log" + "\n" +
                    "\n" +
                    "samtools import " + bacterialGenomeFai1 + " $ScriptDir/4_Mapped_reads/" + csfastaFileName + ".nohum.igc1.sam $ScriptDir/4_Mapped_reads/" + csfastaFileName + ".nohum.igc1.bam" + "\n" +
                    "samtools import " + bacterialGenomeFai2 + " $ScriptDir/4_Mapped_reads/" + csfastaFileName + ".nohum.noigc1.igc2.sam $ScriptDir/4_Mapped_reads/" + csfastaFileName + ".nohum.noigc1.igc2.bam" + "\n" +
                    "\n" +
                    "samtools sort $ScriptDir/4_Mapped_reads/" + csfastaFileName + ".nohum.igc1.bam $ScriptDir/4_Mapped_reads/" + csfastaFileName + ".nohum.igc1.sorted" + "\n" +
                    "samtools sort $ScriptDir/4_Mapped_reads/" + csfastaFileName + ".nohum.noigc1.igc2.bam $ScriptDir/4_Mapped_reads/" + csfastaFileName + ".nohum.noigc1.igc2.sorted" + "\n" +
                    "\n" +
                    "genomeCoverageBed -ibam $ScriptDir/4_Mapped_reads/" + csfastaFileName + ".nohum.igc1.sorted.bam -g " + bacterialGenomeFasta1 + " > $ScriptDir/5_Statistics/" + csfastaFileName + "_genomeCoverageBed_B1.txt" + "\n" +
                    "genomeCoverageBed -ibam $ScriptDir/4_Mapped_reads/" + csfastaFileName + ".nohum.noigc1.igc2.sorted.bam -g " + bacterialGenomeFasta2 + " > $ScriptDir/5_Statistics/" + csfastaFileName + "_genomeCoverageBed_B2.txt" + "\n" +
                    "\n" +
                    "cat $ScriptDir/5_Statistics/" + csfastaFileName + "_genomeCoverageBed_B1.txt  $ScriptDir/5_Statistics/" + csfastaFileName + "_genomeCoverageBed_B2.txt > $ScriptDir/5_Statistics/" + csfastaFileName + "_genomeCoverageBed_B_fused.txt" + "\n" +
                    "\n"
                    )

threeBees.close()
print("YAY")




