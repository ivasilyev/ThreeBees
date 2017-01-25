# Paths should be ending at "/". No spaces are allowed.
# Make sure the script is placed into the root lanes directory containing "LXX" subdirectories.

# [SoftwarePaths]
selfPath = ""
qualityFilterPath = ""
saetPath = ""
qualityTrimmerPath = ""
bowtiePath = ""
samtoolsPath = ""
bedtoolsPath = ""
coverageExporterPath = "~/scripts/"

# [SequencesPaths]
humanGenomePath = "/data/reference/homo_sapiens/ucsc/hg19/sequence/ColorSpaseIndex/genome"

bacterialGenomeBwtMask1 = "/data/reference/IGC/fs/igc1"
bacterialGenomeBwtMask2 = "/data/reference/IGC/fs/igc2"

bacterialGenomeFai1 = "/data/reference/IGC/fs/760MetaHit_139HMP_368PKU_511Bac.fa_2x_1.fasta.fai"
bacterialGenomeFai2 = "/data/reference/IGC/fs/760MetaHit_139HMP_368PKU_511Bac.fa_2x_2.fasta.fai"

bacterialGenomeLengths1 = "/data/reference/IGC/fs/760MetaHit_139HMP_368PKU_511Bac.fa_2x_1.genome"
bacterialGenomeLengths2 = "/data/reference/IGC/fs/760MetaHit_139HMP_368PKU_511Bac.fa_2x_2.genome"

bacterialGenomeTags = "/data/reference/IGC/IGC.annotation.summary.v2"


# [Script begin]
import os, sys, getopt


def usage():
    print("Usage: " + sys.argv[0] + " -f/--folder <folder>\tWorking folder containing \"LXX\" subfolders")
    sys.exit(2)


def main():
    opts = ""
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hf:", ["help", "folder="])
    except getopt.GetoptError as arg_err:
        print(str(arg_err))
        usage()
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
        elif opt in ("-f", "--folder"):
            working_folder = os.path.abspath(arg)
            return working_folder
    print("The working folder is not specified!")
    usage()

workingFolder = main()
os.chdir(workingFolder)


# Workflow folders creation
for newDir in ["1_Preprocessed_reads/filtered/fixed", "2_Processed_reads", "3_Non-mapped_reads", "4_Mapped_reads", "5_Statistics"]:
    os.makedirs(str("ThreeBees/" + newDir))

# Lane processing
threeBees = open("ThreeBees.sh", 'w')
threeBees.write(str("# !/bin/sh" + "\n" +
                    "\n" +
                    "# FOR THE BEST USABILITY USE LAUNCH COMMAND LIKE: nohup bash ThreeBees.sh > /dev/null 2>&1 & echo $! > run.pid" + "\n" +
                    "\n" +
                    "LaneDir=$PWD" + "\n" +
                    "ScriptDir=$LaneDir/ThreeBees" + "\n" +
                    "\n"))

# Getting lanes
lanePaths = []

for lane in ["L01", "L02", "L03", "L04", "L05", "L06"]:
    for file in os.listdir(str(lane + "/result/")):
        if file.endswith(".xsq"):
            lanePaths.append(os.path.abspath(str(lane + "/result/" + file)))

for laneFile in lanePaths:
    threeBees.write(str("# FOR LANE FILE: " + str(laneFile.rsplit('/', 1)[-1]) + "\n" +
                        "cd " + XSQToolsPath + "\n" +
                        "bash convertFromXSQ.sh -s -c " + laneFile + " -o $ScriptDir/1_Preprocessed_reads &>$ScriptDir/5_Statistics/" + str(laneFile.rsplit('/', 1)[-1]).replace(".xsq", "") + "_convertFromXSQ.log\n" +
                        "mv $ScriptDir/1_Preprocessed_reads/Libraries $ScriptDir/1_Preprocessed_reads/Libraries_" + str(laneFile.rsplit('/', 1)[-1]).replace(".xsq", "") + "\n" +
                        "\n"
                        )
                    )
    multiplexParser = open(str(laneFile).replace(".xsq", "_Multiplex_.txt"), 'rU')

    # Predicting demultiplexed reads folders
    multiplexFilenames = []

    for row in multiplexParser:
        if not any(badWord in row for badWord in ["##Library", "Unassigned", "Unclassified", "All"]):
            multiplexFilenames.append(str(row).split('\t', 1)[0])

    # Main workflow
    for sampleName in multiplexFilenames:
        csfastaFileName = str(laneFile.rsplit('/', 1)[-1]).replace(".xsq", "") + "_" + sampleName + "_F3"

        threeBees.write("# FOR SAMPLE FILE: " + csfastaFileName + "\n" +
                        "\n" +
                        "cd $ScriptDir/1_Preprocessed_reads/" + "\n" +
                        qualityFilterPath + "reads-filter_new.pl -f " + csfastaPath + " -q " + str(csfastaPath).replace(".csfasta", ".QV.qual") + " -o filtered -t 15 &>$ScriptDir/5_Statistics/" + csfastaFileName + "_reads-filter_new.log" + "\n" +
                        "\n" +
                        "cd filtered" + "\n" +
                        saetPath + "saet " + csfastaFileName + ".15.filtered.csfasta 200000000 -qual " + csfastaFileName + ".QV.15.filtered.qual -qvupdate -trustprefix 25 -localrounds 3 -globalrounds 2 -numcores 20 -log $ScriptDir/5_Statistics/" + csfastaFileName + "_SAET.log" + "\n" +
                        "\n" +
                        "cd fixed" + "\n" +
                        "mv " + csfastaFileName + ".15.filtered.csfasta " + csfastaFileName + ".15.filtered.fixed.csfasta" + "\n" +
                        "mv none " + csfastaFileName + ".QV.15.filtered.fixed.qual" + "\n" +
                        "\n" +
                        qualityTrimmerPath + "csfasta_quality_filter.pl -f " + csfastaFileName + ".15.filtered.fixed.csfasta -q " + csfastaFileName + ".QV.15.filtered.fixed.qual -o $ScriptDir/2_Processed_reads/" + csfastaFileName + ".15.filtered.fixed.trimmed.csfasta -l 30 -s 30 &>$ScriptDir/5_Statistics/" + csfastaFileName + "_csfasta_quality_filter.log" + "\n" +
                        "\n" +
                        "cd $ScriptDir" +
                        "\n" +
                        bowtiePath + "bowtie -f -C -S -t -v 3 -k 1 --threads 20 --un $ScriptDir/3_Non-mapped_reads/" + csfastaFileName + ".nohum.csfasta " + humanGenomePath + " $ScriptDir/2_Processed_reads/" + csfastaFileName + ".15.filtered.fixed.trimmed.csfasta $ScriptDir/4_Mapped_reads/" + csfastaFileName + ".hum.sam &>$ScriptDir/5_Statistics/" + csfastaFileName + "_bowtieA.log" + "\n" +
                        "\n" +
                        bowtiePath + "bowtie -f -C -S -t -v 3 -k 1 --threads 20 --un $ScriptDir/3_Non-mapped_reads/" + csfastaFileName + ".nohum.noigc1.csfasta " + bacterialGenomeBwtMask1 + " $ScriptDir/3_Non-mapped_reads/" + csfastaFileName + ".nohum.csfasta $ScriptDir/4_Mapped_reads/" + csfastaFileName + ".nohum.igc1.sam &>$ScriptDir/5_Statistics/" + csfastaFileName + "_bowtieB1.log" + "\n" +
                        bowtiePath + "bowtie -f -C -S -t -v 3 -k 1 --threads 20 --un $ScriptDir/3_Non-mapped_reads/" + csfastaFileName + ".nohum.noigc1.noigc2.csfasta " + bacterialGenomeBwtMask2 + " $ScriptDir/3_Non-mapped_reads/" + csfastaFileName + ".nohum.noigc1.csfasta $ScriptDir/4_Mapped_reads/" + csfastaFileName + ".nohum.noigc1.igc2.sam &>$ScriptDir/5_Statistics/" + csfastaFileName + "_bowtieB2.log" + "\n" +
                        "\n" +
                        samtoolsPath + "samtools import " + bacterialGenomeFai1 + " $ScriptDir/4_Mapped_reads/" + csfastaFileName + ".nohum.igc1.sam $ScriptDir/4_Mapped_reads/" + csfastaFileName + ".nohum.igc1.bam" + "\n" +
                        samtoolsPath + "samtools import " + bacterialGenomeFai2 + " $ScriptDir/4_Mapped_reads/" + csfastaFileName + ".nohum.noigc1.igc2.sam $ScriptDir/4_Mapped_reads/" + csfastaFileName + ".nohum.noigc1.igc2.bam" + "\n" +
                        "\n" +
                        samtoolsPath + "samtools sort $ScriptDir/4_Mapped_reads/" + csfastaFileName + ".nohum.igc1.bam $ScriptDir/4_Mapped_reads/" + csfastaFileName + ".nohum.igc1.sorted" + "\n" +
                        samtoolsPath + "samtools sort $ScriptDir/4_Mapped_reads/" + csfastaFileName + ".nohum.noigc1.igc2.bam $ScriptDir/4_Mapped_reads/" + csfastaFileName + ".nohum.noigc1.igc2.sorted" + "\n" +
                        "\n" +
                        bedtoolsPath + "genomeCoverageBed -ibam $ScriptDir/4_Mapped_reads/" + csfastaFileName + ".nohum.igc1.sorted.bam -g " + bacterialGenomeLengths1 + " > $ScriptDir/5_Statistics/" + csfastaFileName + "_genomeCoverageBed_B1.txt" + "\n" +
                        bedtoolsPath + "genomeCoverageBed -ibam $ScriptDir/4_Mapped_reads/" + csfastaFileName + ".nohum.noigc1.igc2.sorted.bam -g " + bacterialGenomeLengths2 + " > $ScriptDir/5_Statistics/" + csfastaFileName + "_genomeCoverageBed_B2.txt" + "\n" +
                        "\n" +
                        "cd $ScriptDir/5_Statistics/"
                        "\n" +
                        "python " + coverageExporterPath + "get_cov_from_bedtools_hist.py " + csfastaFileName + "_genomeCoverageBed_B1.txt " + bacterialGenomeLengths1 + " " + csfastaFileName + "_bp_coverage_B1.txt " + csfastaFileName + "_pos_coverage_B1.txt" + "\n" +
                        "python " + coverageExporterPath + "get_cov_from_bedtools_hist.py " + csfastaFileName + "_genomeCoverageBed_B2.txt " + bacterialGenomeLengths2 + " " + csfastaFileName + "_bp_coverage_B2.txt " + csfastaFileName + "_pos_coverage_B2.txt" + "\n" +
                        "\n" +
                        "cat " + csfastaFileName + "_bp_coverage_B1.txt " + csfastaFileName + "_bp_coverage_B2.txt > " + csfastaFileName + "_bp_coverage_B_merged.txt" + "\n" +
                        "cat " + csfastaFileName + "_pos_coverage_B1.txt " + csfastaFileName + "_pos_coverage_B2.txt > " + csfastaFileName + "_pos_coverage_B_merged.txt" + "\n" +
                        "\n" +
                        "cd $ScriptDir" +
                        "\n\n")

threeBees.write("############################################################################################" + "\n" +
                "python " + selfPath + "ThreeBeesCombiner.py -f 5_Statistics -t " + bacterialGenomeTags + "\n" +
                "\n" +
                "##############" + "\n" +
                "# SCRIPT END #" + "\n" +
                "##############" + "\n")
threeBees.close()


print("The generation has been completed successfully. \n" +
      "Your working directory is: \"" + os.getcwd() + "\"\n" +
      "Please check your working directory and start \"ThreeBees.sh\", the main script. \n" +
      "The example command: \"nohup bash ThreeBees.sh > /dev/null 2>&1 & echo $! > run.pid\". \n" +
      "You have to repeat the generation for the each sample collection.")
