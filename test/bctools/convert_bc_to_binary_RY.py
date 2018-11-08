#!/usr/bin/env python

import argparse
import logging
from string import maketrans
from sys import stdout
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

tool_description = """
Convert standard nucleotides to IUPAC nucleotide codes used for binary barcodes.

A and G are converted to nucleotide code R. T, U and C are converted to Y. By
default output is written to stdout.

Example usage:
- write converted sequences from file in.fa to file file out.fa:
convert_bc_to_binary_RY.py in.fa --outfile out.fa
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
    "infile",
    help="Path to fastq input file.")
# optional arguments
parser.add_argument(
    "-o", "--outfile",
    help="Write results to this file.")
parser.add_argument(
    "--fasta-format",
    dest="fasta_format",
    help="Read and write fasta instead of fastq format.",
    action="store_true")
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
    version='0.1.0')


def translate_nt_to_RY(seq):
    """Translates nucleotides to RY (A,G -> R; C,U,T -> Y).

    >>> translate_nt_to_RY("ACGUTACGUT")
    RYRYYRYRYY
    """
    trans_table = maketrans("AGCUT", "RRYYY")
    trans_seq = seq.translate(trans_table)
    logging.debug(seq + " -> " + trans_seq)
    return trans_seq


def translate_nt_to_RY_iterator(robj):
    """Translate SeqRecords sequences to RY alphabet."""
    for record in robj:
        if not args.fasta_format:
            saved_letter_annotations = record.letter_annotations
        record.letter_annotations = {}
        record.seq = Seq(translate_nt_to_RY(str(record.seq)),
                         IUPAC.unambiguous_dna)
        if not args.fasta_format:
            record.letter_annotations = saved_letter_annotations
        yield record

# handle arguments
args = parser.parse_args()
if args.debug:
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(filename)s - %(levelname)s - %(message)s")
elif args.verbose:
    logging.basicConfig(level=logging.INFO, format="%(filename)s - %(levelname)s - %(message)s")
else:
    logging.basicConfig(format="%(filename)s - %(levelname)s - %(message)s")
logging.info("Parsed arguments:")
if args.outfile:
    logging.info("  outfile: enabled writing to file")
    logging.info("  outfile: '{}'".format(args.outfile))
logging.info("  outfile: '{}'".format(args.outfile))
logging.info("")

# get input iterator
input_handle = open(args.infile, "rU")
if args.fasta_format:
    input_seq_iterator = SeqIO.parse(input_handle, "fasta")
else:
    input_seq_iterator = SeqIO.parse(input_handle, "fastq")
convert_seq_iterator = translate_nt_to_RY_iterator(input_seq_iterator)
output_handle = (open(args.outfile, "w") if args.outfile is not None else stdout)
if args.fasta_format:
    SeqIO.write(convert_seq_iterator, output_handle, "fasta")
else:
    SeqIO.write(convert_seq_iterator, output_handle, "fastq")
output_handle.close()
