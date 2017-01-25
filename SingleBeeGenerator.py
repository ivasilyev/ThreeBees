# This script will generate the single alignment & postprocessing script for every *.csfasta file.
# Make sure the script is placed into the same folder!
# Paths should be ending at "/". No spaces are allowed.

# [SoftwarePaths]
bowtiePath = ""
samtoolsPath = ""
bedtoolsPath = ""
coverageExporterPath = "~/scripts/"
selfPath = "/data1/bio/kazan_solid_metagenome_data/Malanin/"

# [SequencesPaths]
bacterialGenomeBwtMask = "/data1/bio/kazan_solid_metagenome_data/Malanin/bwt"
bacterialGenomeFai = "/data1/bio/kazan_solid_metagenome_data/Malanin/nucleotide_fasta_protein_homolog_model.fasta.fai"
bacterialGenomeLengths = "/data1/bio/kazan_solid_metagenome_data/Malanin/nucleotide_fasta_protein_homolog_model.fasta.genome"
bacterialGenomeTags = "/data1/bio/kazan_solid_metagenome_data/Malanin/annotation.txt"

# [Script begin]

import os

csfastaPaths = []

for file in os.listdir(os.getcwd()):
    if file.endswith(".csfasta"):
        csfastaPaths.append(os.path.abspath(file))
        
for newDir in ["Non-mapped_reads", "Mapped_reads", "Statistics"]:
    os.makedirs(newDir)
        
singleBee = open("SingleBee.sh", 'w')

for csfastaPath in csfastaPaths:
    csfastaFileName = str(csfastaPath.rsplit('/', 1)[-1]).replace(".nohum.csfasta", "")

    singleBee.write("# FOR SAMPLE FILE: " + csfastaFileName + "\n" +
                    bowtiePath + "bowtie -f -C -S -t -v 3 -k 1 --threads 20 --un Non-mapped_reads/" + csfastaFileName + ".csfasta " + bacterialGenomeBwtMask + " " + csfastaFileName + ".nohum.csfasta Mapped_reads/" + csfastaFileName + ".sam &>Statistics/" + csfastaFileName + "_bowtie.log" + "\n" +
                    "\n" +
                    samtoolsPath + "samtools import " + bacterialGenomeFai + " Mapped_reads/" + csfastaFileName + ".sam Mapped_reads/" + csfastaFileName + ".bam" + "\n" +
                    "\n" +
                    samtoolsPath + "samtools sort Mapped_reads/" + csfastaFileName + ".bam Mapped_reads/" + csfastaFileName + ".sorted" + "\n" +
                    "\n" +
                    bedtoolsPath + "genomeCoverageBed -ibam Mapped_reads/" + csfastaFileName + ".sorted.bam -g " + bacterialGenomeLengths + " > Statistics/" + csfastaFileName + "_genomeCoverageBed.txt" + "\n" +
                    "\n" +
                    "cd Statistics/"
                    "\n" +
                    "python " + coverageExporterPath + "get_cov_from_bedtools_hist.py " + csfastaFileName + "_genomeCoverageBed.txt " + bacterialGenomeLengths + " " + csfastaFileName + "_bp_coverage_B_merged.txt " + csfastaFileName + "_pos_coverage_B_merged.txt" + "\n" +
                    "\n" +
                    "cd " + os.getcwd() + "\n" +
                    "\n\n")

singleBee.write("############################################################################################" + "\n" +
                "python " + selfPath + "Combiner.py -f Statistics -t " + bacterialGenomeTags + "\n" +
                "\n" +
                "##############" + "\n" +
                "# SCRIPT END #" + "\n" +
                "##############" + "\n")
singleBee.close()

print("The generation has been completed successfully. \n" +
      "Your working directory is: \"" + os.getcwd() + "\"\n" +
      "Please check your working directory and start \"SingleBee.sh\", the main script. \n" +
      "The example command: \"nohup bash SingleBee.sh > /dev/null 2>&1 & echo $! > run.pid\". \n" +
      "You have to repeat the generation for the each sample collection.")

