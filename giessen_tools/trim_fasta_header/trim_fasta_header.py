#!/usr/bin/env python3

__author__ = "sgriep"
__projekt__ = "trim_fasta_header"
__date__ = "2021-02-15"
__version__ = "1.0"


import sys
import re

input_file = sys.argv[1]
output_file = sys.argv[2]

fh_out = open(output_file, 'w')
fh_in = open(input_file, 'r')

for line in fh_in:
    if line[0] != '>':
        fh_out.write(line)
    else:
        cols = line.split()
        fh_out.write('{}\n'.format(cols[0]))

fh_out.close
fh_in.close

