#!/usr/bin/env python3

__author__ = "sgriep"
__projekt__ = "demultiplex_by_primers"
__date__ = "2021-02-15"
__version__ = "1.0"


import sys
import os
import re
import argparse
import subprocess as sp

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help='input fasta file to demultiplex')
    parser.add_argument('-o', '--output', required=True, help='output fasta file with demultiplexed readname')
    parser.add_argument('-p', '--primers', required=True, help='primer table')
    args = parser.parse_args()
    return args

def process_primers(fasta_file, primer_file, outfile):
    sample = get_sample(fasta_file)
    print('sample: {}'.format(sample))

    deplexed_fastas = []
    fasta_for_next_primer = fasta_file
    with open(primer_file) as fh_primer:
        for line in fh_primer:
            cols = line.strip().split('\t')
            if len(cols) > 1 and cols[1] == sample:
                primer = cols[0]
                sample_primer_file = 'process.primer.txt'
                with open(sample_primer_file, 'w') as fh_sample_primer_file:
                    fh_sample_primer_file.write(line)
                demultiplexed_fasta, fasta_for_next_primer = ngsfilter(sample, primer, fasta_for_next_primer, sample_primer_file)
                deplexed_fasta = readname_changer(primer, demultiplexed_fasta)
                deplexed_fastas.append(deplexed_fasta)
    concat_fastas(deplexed_fastas, outfile)

    return

def get_sample(fasta_file):
    read_sample = ''
    with open(fasta_file) as fh:
        for line in fh:
            if line[0] == '>':
                read_id_line = line.split('\t')
                read_id = read_id_line[0][1:].split('___')
                read_sample = re.sub('_L001_R[12]_001', '', read_id[0])
                break
    return read_sample

def ngsfilter(sample, primer, fasta_file, primer_file):

    prefix = '{sample}_{primer}'.format(sample=sample, primer=primer)
    
    demultiplexed_fasta = '{prefix}.fasta'.format(prefix=prefix)
    fh_demultiplexed_fasta = open(demultiplexed_fasta, 'w')
    
    fasta_for_next_primer = '{prefix}.next.fasta'.format(prefix=prefix)
    fh_fasta_for_next_primer = open(fasta_for_next_primer, 'w')
    
    temp_fasta = 'process.fasta'

    cmd_ngsfilter = [
        'ngsfilter',
        fasta_file,
        '-t', primer_file,
        '-u', temp_fasta
    ]
    cmd_obigrep_1 = [
        'obigrep',
        '-a', 'error:Cannot assign sequence to a sample',
        temp_fasta
    ]
    cmd_obigrep_2 = [
        'obigrep',
        '-a', 'error:No reverse primer match',
        temp_fasta,
    ]
    cmd_obigrep_3 = [
        'obigrep',
        '-a', 'error:No primer match',
        temp_fasta,
    ]

    sp.check_call(cmd_ngsfilter)
    sp.call(cmd_obigrep_1,
            stdout=fh_demultiplexed_fasta,
           )
    sp.call(cmd_obigrep_2,
            stdout=fh_demultiplexed_fasta,
           )
    sp.call(cmd_obigrep_3,
            stdout=fh_fasta_for_next_primer,
           )

    fh_demultiplexed_fasta.close()
    fh_fasta_for_next_primer.close()

    return demultiplexed_fasta, fasta_for_next_primer

def readname_changer(primer, fasta_file):
    fasta_renamed_reads = '{prefix}.renamed.fasta'.format(prefix=os.path.splitext(fasta_file)[0])
    fh_fasta_renamed_reads = open(fasta_renamed_reads, 'w')

    cmd_prepend_readnames = [
        'sed',
        '-e', 's/^>\([^ ]\+\)\(___\)/>\\1_{primer}\\2/'.format(
            primer=primer,
        ),
        fasta_file
    ]

    sp.call(cmd_prepend_readnames,
            stdout=fh_fasta_renamed_reads,
           )
    
    return fasta_renamed_reads

def concat_fastas(fasta_files, concat_file):

    cmd_cat_fastas = ['cat']
    cmd_cat_fastas.extend(fasta_files)

    sp.call(cmd_cat_fastas,
            stdout=open(concat_file, 'w'),
           )

    return

def main():
    args = get_args()

    if os.path.exists(args.input) and os.path.exists(args.primers):
        process_primers(args.input, args.primers, args.output)

        return 0
    return 1

if __name__ == "__main__":
    sys.exit(main())

