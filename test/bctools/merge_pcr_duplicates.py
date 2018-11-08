#!/usr/bin/env python

import argparse
import logging
from sys import stdout
from subprocess import check_call
from shutil import rmtree
from tempfile import mkdtemp
from os.path import isfile
# avoid ugly python IOError when stdout output is piped into another program
# and then truncated (such as piping to head)
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)

tool_description = """
Merge PCR duplicates according to random barcode library.

Barcodes containing uncalled base 'N' are removed.

Input:
* bed6 file containing alignments with fastq read-id in name field
* fastq library of random barcodes

Output:
* bed6 file with random barcode in name field and number of PCR duplicates as
  score, sorted by fields chrom, start, stop, strand, name

Example usage:
- read PCR duplicates from file duplicates.bed and write merged results to file
  merged.bed:
merge_pcr_duplicates.py duplicates.bed bclibrary.fa --outfile merged.bed
"""

epilog = """
Author: Daniel Maticzka
Copyright: 2015
License: Apache
Email: maticzkd@informatik.uni-freiburg.de
Status: Testing
"""

# parse command line arguments
parser = argparse.ArgumentParser(description=tool_description,
                                 epilog=epilog,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
# positional arguments
parser.add_argument(
    "alignments",
    help="Path to bed6 file containing alignments.")
parser.add_argument(
    "bclib",
    help="Path to fastq barcode library.")
# optional arguments
parser.add_argument(
    "-o", "--outfile",
    required=True,
    help="Write results to this file.")
# misc arguments
parser.add_argument(
    "-v", "--verbose",
    help="Be verbose.",
    action="store_true")
parser.add_argument(
    "-d", "--debug",
    help="Print lots of debugging information",
    action="store_true")
parser.add_argument(
    '--version',
    action='version',
    version='0.2.0')

args = parser.parse_args()

if args.debug:
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(filename)s - %(levelname)s - %(message)s")
elif args.verbose:
    logging.basicConfig(level=logging.INFO, format="%(filename)s - %(levelname)s - %(message)s")
else:
    logging.basicConfig(format="%(filename)s - %(levelname)s - %(message)s")
logging.info("Parsed arguments:")
logging.info("  alignments: '{}'".format(args.alignments))
logging.info("  bclib: '{}'".format(args.bclib))
if args.outfile:
    logging.info("  outfile: enabled writing to file")
    logging.info("  outfile: '{}'".format(args.outfile))
logging.info("")

# see if alignments are empty and the tool can quit
n_alns = sum(1 for line in open(args.alignments))
if n_alns == 0:
    logging.warning("WARNING: Working on empty set of alignments, writing empty output.")
    eventalnout = (open(args.outfile, "w") if args.outfile is not None else stdout)
    eventalnout.close()
    exit(0)

# check input filenames
if not isfile(args.bclib):
    raise Exception("ERROR: barcode library '{}' not found.")
if not isfile(args.alignments):
    raise Exception("ERROR: alignments '{}' not found.")

try:
    tmpdir = mkdtemp()
    logging.debug("tmpdir: " + tmpdir)

    # prepare alinments
    syscall2 = "cat " + args.alignments + " | awk -F \"\\t\" 'BEGIN{OFS=\"\\t\"}{split($4, a, \" \"); $4 = a[1]; print}'| sort --compress-program=gzip -k4,4 > " + tmpdir + "/alns.csv"
    check_call(syscall2, shell=True)

    # join barcode library and alignments
    # after join: id, bc, chr, start, stop, mapscore, strand
    # after datamash: bc, chr, start, stop, strand, ndupes, idrepresentative
    syscall3 = "cat " + \
        args.bclib + \
        " | awk 'BEGIN{OFS=\"\\t\"}NR%4==1{gsub(/^@/,\"\"); id=$1}NR%4==2{bc=$1}NR%4==3{print id,bc}' " + \
        " | sort --compress-program=gzip -k1,1 | join -1 1 -2 4 - " + tmpdir + "/alns.csv " + \
        " | awk 'BEGIN{OFS=\"\\t\"}$2!~/N/{print $1,$2,$3,$4,$5,$6,$7}' " + \
        " | datamash --sort -g 2,3,4,5,7 count 2 first 1 " + \
        " | awk 'BEGIN{OFS=\"\\t\"}{print $2,$3,$4,$7,$6,$5}' > " + args.outfile
    check_call(syscall3, shell=True)
finally:
    logging.debug("removed tmpdir: " + tmpdir)
    rmtree(tmpdir)
