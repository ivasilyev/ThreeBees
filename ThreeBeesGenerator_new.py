# Paths should be ending at "/". No spaces are allowed.
# Make sure the script is placed into root lanes directory containing "LXX" directories.

# [SoftwarePaths]
XSQToolsPath = "/data/projects/lifescope_xsq/XSQ_Tools/"
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

# Get lanes
lanePaths = []

for lane in ["L01", "L02", "L03", "L04", "L05", "L06"]:
    for file in os.listdir(str(lane + "/result/")):
        if file.endswith(".xsq"):
            lanePaths.append(os.path.abspath(str(lane + "/result/" + file)))

# root folder:  /data2/bio/sandbox/kazan5500_2016_06_21_1/
# lane file:    /data2/bio/sandbox/kazan5500_2016_06_21_1/L01/result/
#                            kazan5500_2016_06_21_1_L01.xsq

# demultiplexed sample file: /data2/bio/sandbox/kazan5500_2016_06_21_1/ThreeBees/1_Preprocessed_reads/
#                            kazan5500_2016_06_21_1_L01_1_002_01.xsq

# csfasta file: /data2/bio/sandbox/kazan5500_2016_06_21_1/
#                            ThreeBees/1_Preprocessed_reads/Libraries/1_002_01/F3/reads/
#                            " + csfastaFileName + ".csfasta

# Workflow folders creation
for newDir in ["1_Preprocessed_reads", "2_Processed_reads", "3_Non-mapped_reads", "4_Mapped_reads", "5_Statistics"]:
    os.makedirs(str("ThreeBees/" + newDir))

# Lane proessing
threeBees = open("ThreeBees.sh", 'w')
threeBees.write(str("LaneDir=$PWD\n" +
                    "ScriptDir=$LaneDir/ThreeBees\n"))

for laneFile in lanePaths:
    threeBees.write(str("# FOR LANE FILE: " + str(laneFile.rsplit('/', 1)[-1]) + "\n" +
                        "cd " + XSQToolsPath + "\n" +
                        "bash convertFromXSQ.sh -s -c " + laneFile + " -o $ScriptDir/1_Preprocessed_reads &>$ScriptDir/5_Statistics/" + str(laneFile.rsplit('/', 1)[-1]).replace(".xsq", "") + "_convertFromXSQ.log\n" +
                        "mv $ScriptDir/1_Preprocessed_reads/Libraries $ScriptDir/1_Preprocessed_reads/Libraries_" + str(laneFile.rsplit('/', 1)[-1]).replace(".xsq", "") + "\n" +
                        "\n"
                        )
                    )

    # Predicting demultiplexed folders
    multiplexParser = open(str(laneFile).replace(".xsq", "_Multiplex_.txt"), 'rU')
    multiplexFilenames = []

    for row in multiplexParser:
        if not any(badWord in row for badWord in ["##Library", "Unassigned", "Unclassified", "All"]):
            multiplexFilenames.append(str(row).split('\t', 1)[0])

    # Main workflow
    for sampleName in multiplexFilenames:
        csfastaFileName = str(laneFile.rsplit('/', 1)[-1]).replace(".xsq", "") + "_" + sampleName + "_F3"

        threeBees.write("# FOR SAMPLE FILE: " + csfastaFileName + "\n" +
                        "cd $ScriptDir/1_Preprocessed_reads/Libraries_" + str(laneFile.rsplit('/', 1)[-1]).replace(".xsq", "") + "/" + sampleName + "/F3/reads" + "\n" +
                        "mkdir filtered" + "\n" +
                        qualityFilterPath + "reads-filter_new.pl -f " + csfastaFileName + ".csfasta -q " + csfastaFileName + ".QV.qual -o filtered -t 15 &>$ScriptDir/5_Statistics/" + csfastaFileName + "_reads-filter_new.log" + "\n" +
                        "\n" +
                        "cd filtered" + "\n" +
                        "mkdir fixed" + "\n" +
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
