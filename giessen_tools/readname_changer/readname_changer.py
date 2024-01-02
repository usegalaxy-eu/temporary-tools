#!/usr/bin/env python3

__author__ = "mkuemmel"
__projekt__ = "file name in read names integrator"
__date__ = "21.09.2018"
__version__ = "1.0"

"""
Programm, to integrate the name of the input file in to the the readnames of the fastq file itself.
"""

import argparse
import subprocess
from pathlib import Path

def readnameChanger(args):
    # saves the path of the given input file
    inputFilePath = Path(args.input[0])
    # saves the path of the given input file
    outputFileFilePath = Path(args.output[0])
    # saves the new name for adding bevore the read names
    nameString = str(args.name[0])
    nameString = nameString.split(".")[0]


    # start to change fastq header
    print("-------------change headers of the Fastq file for \"" + nameString + "\"-----------------------\n")
    # before: @NameRead
    # after: @Dataname___NameRead
    process = subprocess.Popen("sed -e '1~4s/^@\(.*\)/@" + nameString + "___\\1/' " + str(inputFilePath.absolute()) + " > " + str(outputFileFilePath.absolute()), shell=True,
                               stdout=subprocess.PIPE)
    process.wait()
    print("---------------------------------------------done-----------------------------------------------\n")


if __name__ == "__main__":

    """
    main method.
    Implements the argparser and the help.
    If all important arguments are given, the program will be started and gets the args input.
    If some arguments are missing the method prints the help and ends.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", nargs=1, help= "Path to the file where the readnames will be changed. The file must be in fastq!")
    parser.add_argument("-n", "--name", nargs=1, help= "Name to insert bevore the readnames.")
    parser.add_argument("-o", "--output", nargs=1, help= "Path to the new changed file. The file must be in fastq!")


    args = parser.parse_args()

    if args.input and args.output and args.name:
        path_new_File = readnameChanger(args)
    else:
        parser.print_help()