
#########################################################################
# An internal script. Use it when the statistics gathering is complete. #
#########################################################################

# The tag file must contain a header with "Gene_name" column.
# Other files must not contain a header and have only two columns.

import os, sys, getopt, pandas


def usage():
    print("Usage: " + sys.argv[0] + " -f <folder> -t <tags>" + "\n\n" +
          "-f/--folder <folder> \t\tWorking folder containing files ending with \"_coverage_B_merged.txt\" " + "\n" +
          "-t/--tags <file> \t\tFASTA headers annotation file")
    sys.exit(2)


def main():
    opts = ""
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hf:t:", ["help", "folder=", "tags="])
    except getopt.GetoptError as arg_err:
        print(str(arg_err))
        usage()

    results_folder = None
    tag_file = None

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
        elif opt in ("-f", "--folder"):
            results_folder = os.path.abspath(arg)
        elif opt in ("-t", "--tags"):
            tag_file = os.path.abspath(arg)

    if (results_folder is not None) and (tag_file is not None):
        return results_folder, tag_file

    print("The parameters are not yet specified!")
    usage()

resultsPath, tagFile = main()

os.makedirs(str(resultsPath + "/Finalized"))

mainDF = pandas.read_table(tagFile, sep='\t', header=0, engine='python')

for bpOrPos in ["bp", "pos"]:
    resultFileList = []

    for file in os.listdir(resultsPath):
        if file.endswith("_F3_" + bpOrPos + "_coverage_B_merged.txt"):
            resultFileList.append(os.path.abspath(str(resultsPath + "/" + file)))

    resultFileList.sort()
    resultCollectionTable = mainDF

    for resultFile in resultFileList:
        sampleName = str(resultFile.rsplit('/', 1)[-1]).replace(str("_F3_" + bpOrPos + "_coverage_B_merged.txt"), "")
        resultSampleTable = pandas.read_table(resultFile, sep='\t', header='infer', names=["Gene_name", sampleName], engine='python')
        resultCollectionTable = pandas.merge(resultCollectionTable, resultSampleTable, on="Gene_name", how='left')

    resultCollectionTable.to_csv(str(resultsPath + "/Finalized/" + bpOrPos + ".txt"), sep='\t', index=False)

