#!/usr/bin/env python

import argparse
import logging
from subprocess import check_call
import os

tool_description = """
Remove spurious events originating from errors in random sequence tags.

This script compares all events sharing the same coordinates. Among each group
of events the maximum number of PCR duplicates is determined. All events that
are supported by less than 10 percent of this maximum count are removed.

Input:
* bed6 file containing crosslinking events with score field set to number of PCR
  duplicates

Output:
* bed6 file with spurious crosslinking events removed, sorted by fields chrom,
  start, stop, strand

Example usage:
- remove spurious events from spurious.bed and write results to file cleaned.bed
rm_spurious_events.py spurious.bed --oufile cleaned.bed
"""

epilog = """
Author: Daniel Maticzka
Copyright: 2015
License: Apache
Email: maticzkd@informatik.uni-freiburg.de
Status: Testing
"""


class DefaultsRawDescriptionHelpFormatter(argparse.ArgumentDefaultsHelpFormatter,
                                          argparse.RawDescriptionHelpFormatter):
    # To join the behaviour of RawDescriptionHelpFormatter with that of ArgumentDefaultsHelpFormatter
    pass


def main():
    # parse command line arguments
    parser = argparse.ArgumentParser(description=tool_description,
                                     epilog=epilog,
                                     formatter_class=DefaultsRawDescriptionHelpFormatter)
    # positional arguments
    parser.add_argument(
        "events",
        help="Path to bed6 file containing alignments.")
    # optional arguments
    parser.add_argument(
        "-o", "--outfile",
        required=True,
        help="Write results to this file.")
    parser.add_argument(
        "-t", "--threshold",
        type=float,
        default=0.1,
        help="Threshold for spurious event removal."
    )
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
        version='0.1.0')

    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(filename)s - %(levelname)s - %(message)s")
    elif args.verbose:
        logging.basicConfig(level=logging.INFO, format="%(filename)s - %(levelname)s - %(message)s")
    else:
        logging.basicConfig(format="%(filename)s - %(levelname)s - %(message)s")
    logging.info("Parsed arguments:")
    logging.info("  alignments: '{}'".format(args.events))
    logging.info("  threshold: '{}'".format(args.threshold))
    if args.outfile:
        logging.info("  outfile: enabled writing to file")
        logging.info("  outfile: '{}'".format(args.outfile))
    logging.info("")

    # check threshold parameter value
    if args.threshold < 0 or args.threshold > 1:
        raise ValueError("Threshold must be in [0,1].")

    if not os.path.isfile(args.events):
        raise Exception("ERROR: file '{}' not found.")

    # prepare barcode library
    syscall = "cat " + args.events + " | sort -k1,1V -k6,6 -k2,2n -k3,3 -k5,5nr | perl " + os.path.dirname(os.path.realpath(__file__)) + "/rm_spurious_events.pl --frac_max " + str(args.threshold) + "| sort -k1,1V -k2,2n -k3,3n -k6,6 -k4,4 -k5,5nr > " + args.outfile
    check_call(syscall, shell=True)


if __name__ == "__main__":
    main()
