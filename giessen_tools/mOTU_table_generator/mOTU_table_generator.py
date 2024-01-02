#!/usr/bin/env python3

__author__ = "mkuemmel"
__projekt__ = "mOTU_table_generator_for_Galaxy"
__date__ = "26.09.18"
__version__ = "1.0"


"""
program to create a "mOTUR_table".
Needs a ".tabular" file of the program usearch with the "-cluster_fast -uc" settings.
And a blast ".tabular" file, with the following entries:
    "-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames'"
"""

import argparse
from pathlib import Path
from collections import Counter

def converter(args):
    """
    Method to convert an ".tabular" usearch file and a blast ".tabular" file in to a csv data.
    This csv shows which dataset is in which cluster.

    :param args: list of the method "argparse". contents the arguments for the pipeline.

    """
    # File paths
    inputFilePath = Path(args.input[0])
    outputFilePath = Path(args.output[0])
    blastinputFilePath = Path(args.input_Blast[0])
    blastoutputFilePath = Path(args.output_Blast[0])

    # Lists
    dictCluster = {}  # Lists all file names for each cluster. The cluster are the keys (dict)
    listDataName = []  # List wih all data names
    dictTaxa = {}  # Lists all taxa(sseqid) for each read in the blast file. The read names are the keys (dict)
    dict_Read_to_Cluster = {}  # Dict of the cluster numbers (keys). The read names are the value.
    dict_Cluster_to_Read = {}  # Dict of the read names (keys). The cluster are teh value.
    dict_for_best_hit = {}  # Read as Key and ssciname as best hit (if there are more than one, the entrys are separated with an ";")
    dict_for_best_hit_bitscore = {}  # Read as Key and tuple of bitscore and pident as entry. Same order as the entrys of dict_for_best_hit


    # opens the given .blast input file path
    blastinputFile = open(str(blastinputFilePath.absolute()))


    # Iterates over the .blast input file
    for line in blastinputFile:
        # header are not needed
        if line.startswith("#"):
            continue
        lineentrys = line.split("\t")  # gets all columns of the actual line
        readname = lineentrys[0].split(";")[0]  # saves the read name
        # sets the read name as key and the ssciname as entry
        # only the top match wil be saved
        if readname in dictTaxa:
            pass
        else:
            dictTaxa[readname] = lineentrys[13].strip()
        bitscore = lineentrys[11].strip()
        if readname in dict_for_best_hit:
            # if there is the same bitscore like the loop before the read will be appended to dict_for_best_hit
            if float(bitscore) == float(before_loop_bitscore):
                dict_for_best_hit[readname].append(lineentrys[13].strip())  # append the ssciname to list for the reads best hits
                dict_for_best_hit_bitscore[readname].append((lineentrys[13].strip(), bitscore, lineentrys[2]))  # append a list of ssciname, bitscore and pident fot the reads best hits
            # if the bitscore is bigger than the loop before the read will be set as new list
            elif float(bitscore) > float(before_loop_bitscore):
                dict_for_best_hit[readname] = [lineentrys[13].strip()]  # sets the ssciname as new list for the reads best hits in a list
                dict_for_best_hit_bitscore[readname] = [(lineentrys[13].strip(), bitscore, lineentrys[2])]  # sets a list of ssciname, bitscore and pident for the reads best hits in a list
                before_loop_bitscore = bitscore # if the bitscore is higher, a new reference bit score is needed
        # if the readname is not in the dict. a new entry will be generated
        else:
            dict_for_best_hit[readname] = [lineentrys[13].strip()]  # sets the ssciname as new list for the reads best hits in a list
            dict_for_best_hit_bitscore[readname] = [(lineentrys[13].strip(), bitscore, lineentrys[2])]  # sets a list of ssciname, bitscore and pident for the reads best hits in a list
            before_loop_bitscore = bitscore  # the new reference bit score is set

    # closes the .blast input file
    blastinputFile.close()

    # opens the given usearch .tabular input file path an greps the cluster numbers and the names of the files, where the reads comes from.
    inputFile = open(str(inputFilePath.absolute()))
    for line in inputFile:
        lineentrys = line.split()
        # names of the data files
        dataName = lineentrys[8]
        dataName = dataName.split("___")[0]
        # if the cluster number already exists as key append the data name to the entry
        if int(lineentrys[1]) in dictCluster:
            dictCluster[int(lineentrys[1])].append(dataName)
        # if the cluster number not exist as key generate a new entry as a list
        else:
            dictCluster[int(lineentrys[1])] = [dataName]

        # saves the data name into a list
        if dataName in listDataName:
            pass
        else:
            listDataName.append(dataName)

        # if a read name is in the dictTaxa, the read is in the blast data an will be needed
        if lineentrys[8].split(";", 1)[0] in dictTaxa:
            # in the dict_Read_to_Cluster the cluster number will be saved as key and the read name will be the value
            dict_Read_to_Cluster[int(lineentrys[1])] = lineentrys[8].split(";", 1)[0]
            # needed to convert the read name in the blast file in to a cluster number
            # in the dict_Cluster_to_Read the read name will be saved as key and the cluster number will be the value
            dict_Cluster_to_Read[lineentrys[8].split(";", 1)[0]] = int(lineentrys[1])

    # closes the input file
    inputFile.close()

    # writes the string in the new output file. Path given by "outputFilePath"
    newOutputFile = open(str(outputFilePath.absolute()), "w")

    # sets the first part of the header of the first column in the ".mOTU_table" file
    newText = "OUT_Nr\tssciname\tpident\tbitscore"

    # sets the header for the rest of the first column
    for name in sorted(listDataName):
        newText = newText + "\t" + name
    newText = newText + "\n"
    newOutputFile.write(newText)

    # writes the matrix in to the string
    for cluster in sorted(dictCluster):
        # only if the cluster also exists in the blast data, the line will be written
        if cluster in dict_Read_to_Cluster:
            # variables for the newText
            ssciname = ""
            pident_for_text = ""
            bitscore_for_text = ""
            # uses the Counter method, to count each identical elements in dict_for_best_hit
            countervar = Counter(dict_for_best_hit[dict_Read_to_Cluster[cluster]])
            # returns a list of the 10 most elements
            varvar = countervar.most_common(10)  # 10 is a random number

            # iterates over the list generated by the Counter method
            for i, e in enumerate(varvar):
                # in the first iteration the variable can just be set
                if i == 0:
                    ssciname = ssciname + str((varvar[i][0]))
                    # iterates over the dict for best hit bitscore and sets the best pident and bitscore for the ssciname of the readname (blast data)
                    for elements in dict_for_best_hit_bitscore[dict_Read_to_Cluster[cluster]]:
                        if str((varvar[i][0])) == elements[0]:
                            pident_for_text = elements[2]
                            bitscore_for_text = elements[1]
                            break
                else:
                    # the first element is already the highest bit score. All other elements are only same or worse.
                    if not (varvar[i-1][1] == varvar[i][1]):
                        break
                    else:
                        # if the elements are same the entrys will be written behind the variables for the text separated by ";"
                        ssciname = ssciname + ";" + str((varvar[i][0]))
                        # iterates over the dict for best hit bitscore and sets the best pident and bitscore for the ssciname of the readname (blast data)
                        for elements in dict_for_best_hit_bitscore[dict_Read_to_Cluster[cluster]]:
                            if str((varvar[i][0])) == elements[0]:
                                pident_for_text = elements[2]
                                bitscore_for_text = elements[1]
                                break

            # writes all variables into the newText variable
            newText = str(cluster) + "\t" + ssciname + "\t" + pident_for_text + "\t" + bitscore_for_text
            # write the number of names in the list in the dict.
            for name in sorted(listDataName):
                newText = newText + "\t" + str(dictCluster[cluster].count(name))
            newText = newText + "\n"
            newOutputFile.write(newText)
        # if the cluster don't exists in the blast file this loop will be passed
        else:
            pass

    newOutputFile.close()


    print("write blastout")
    newBlastoutCSVFile = open(str(blastoutputFilePath.absolute()), "w")

    # opens the given blastFile Path a second time
    blastinputFile = open(str(blastinputFilePath.absolute()))
    newText = "qsequi (OUT) (Clusternumber)\tsseqid (Genbank)\tpident (identity)\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tstaxids\tssciname\n"
    newBlastoutCSVFile.write(newText)

    print("start blastinput")
    for line in blastinputFile:
        if line.startswith("#"):
            continue
        # splits the line into the readname and the rest of the line
        lineentrys = line.split("\t", 1)

        # replace the readname with the clustername
        newText = str(dict_Cluster_to_Read[lineentrys[0].split(";", 1)[0]])

        # writes the rest of the line behind it
        newText = newText + "\t" + lineentrys[1]
        newBlastoutCSVFile.write(newText)

    # closes the input blast file
    blastinputFile.close()

    newBlastoutCSVFile.close()



if __name__ == "__main__":

    """
    main method.
    Implements the argparser and the help.
    If all important arguments are given, the pipe will be started and gets the args input.
    If some arguments are missing the method prints the help and ends.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", nargs=1, help="tabular file. USEARCH cluster format from the -cluster_fast -uc settings")
    parser.add_argument("-o", "--output", nargs=1, help="Path/Name of the new mOTUR table")
    parser.add_argument("-ib", "--input_Blast", nargs=1, help="Tab separated blast file. need following entrys:\n"
                                                              "-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore  staxids sscinames'")
    parser.add_argument("-ob", "--output_Blast", nargs=1, help="Path/Name of the new output blast data.")

    args = parser.parse_args()

    if args.input and args.output:
        converter(args)
    else:
        parser.print_help()
