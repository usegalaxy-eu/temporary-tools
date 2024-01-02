#!/usr/bin/env python
# coding: utf-8

__author__ = "lkoehl"
__projekt__ = "index_read_filter"
__date__ = "2020-04-16"
__version__ = "1.0"



import sys
import csv
import gzip
import math
import numpy as np
from Bio import SeqIO

# The program filters the reads for each quality in the range 0 <= q <= MAX_QUAL
# Possible BUG: might be necessary to expose this as an external parameter of the script and galaxy tool, but 36 may be sufficient
MAX_QUAL = 40

# For a given i7 index i, find all i5 indices that are not paired with it
# That's set(i5s) - set(partners)
# The sample sheet must contain all pairs (i7,i5) that identify a sample

def partners(i7):
    """Provide a list of i5 sequences the given i7 is partnered with"""
    return [i5s()[i] for i, x in enumerate(i7s()) if x == i7]


def control_pairs_single(i7):
    """Build a set of pairs {(i7, a), (i7, b), ...} that don't occur in the sample sheet
    The index sequences from the sample sheet are represented as internal indexes (int)
    They are the controls for estimating misassignments"""
    non_partners = set(i5s()) - set(partners(i7))
    return {(I7I[i7], I5I[np]) for np in non_partners}


from itertools import chain


def control_pairs():
    return set(chain(*map(control_pairs_single, i7s())))


i7s = lambda: [pair[0] for pair in SAMPLE_SHEET]
i5s = lambda: [pair[1] for pair in SAMPLE_SHEET]


def initialize_sample_sheet_globals(sample_sheet_fname):
    # sample sheet data are globals for now. They are immutable once initialized
    # XXX: if the programs grows, this could be refactored, eg. into a class sample_sheet
    global SAMPLE_SHEET, I7I, I5I, CONTROL_PAIRS, I7S, I5S
    with open(sample_sheet_fname) as f:
        r = csv.reader(f)
        # the sample sheet is just the true indexes as rows, eg. i7,i5 ...
        SAMPLE_SHEET = [line for line in r]

    # map sequence index (str) to internal indexes (int) for identifying the i7/i5 sequences from the sample sheet
    I7I = dict((v, k) for k,v in enumerate(set(i7s())))
    I5I = dict((v, k) for k,v in enumerate(set(i5s())))
    I7S = dict((I7I[k], k) for k in I7I)
    I5S = dict((I5I[k], k) for k in I5I)
    # The control pairs used to estimate the misassignment rate.
    CONTROL_PAIRS = control_pairs()


# Build a list of tuples (i7i, i5i, qavg) from the two given fastq files



def average_qscore_vec(fastq_qual):
    """Same as above but vectorized"""
    qual = np.array(fastq_qual, dtype=np.int8)
    prob = np.power(10, (-qual) / 10)
    mean = np.mean(prob)
    qavg = -10 * math.log10(mean)
    return int(qavg)


# def fuzzy_index_equality(genuine: str, aspires: str):
#     """Return genuine if the second string is equal to the first except for maybe one character position"""
#     if len(genuine) != len(aspires): return False
#     s = sum(int(i == j) for i, j in zip(genuine, aspires))
#     return genuine if 0 <= len(genuine) - s <= 1 else False

others = []

def build_index_tuples(i7fq, i5fq):
    """Yield the index pairs from the samples fastq files, include the average quality for each index read"""
    for (i7, i5) in zip(SeqIO.parse(gzip.open(i7fq, 'rt'), 'fastq'), SeqIO.parse(gzip.open(i5fq, 'rt'), 'fastq')):
        (i7_name, i7_seq, i7_qual) = (i7.id, i7.seq, i7.letter_annotations["phred_quality"])
        (i5_name, i5_seq, i5_qual) = (i5.id, i5.seq, i5.letter_annotations["phred_quality"])
        try:
            i7i = I7I[i7_seq]
            i5i = I5I[i5_seq]
            qavg = average_qscore_vec(i7_qual + i5_qual)
            #assert i7_name==i5_name
            yield (i7_name, i7i, i5i, qavg)
        except KeyError:
            # Maybe TODO: handle index reads with base errors using fuzzy_index_equality
            # This requires the program to know the expected indexes for the sample its currently handling
            # Galaxy provides just the file paths to the datasets. Can we assume this information is always in the fastq comment?
            others.append((str(i7_seq), str(i5_seq)))
            continue


# Filter the list by quality 0 <= q <= 36
# Increasing the quality score required for a reads to pass the cutoff filter mean the misassignment rate should go down


def filter_triples(index_triples, quality_cutoff):
    """Filter the given list of triples (i7, i5, q) to contain only the (i7, i5) where q >= to a minimal required quality"""
    return [(name, i7, i5) for name, i7, i5, q in index_triples if q >= quality_cutoff]


# Count the number of times the control pairs occur in the filtered list
# Do that for each quality score

from collections import Counter

def get_sample_index_filter_lists(index_triples, quality):
    """ Provide the quality and the two key numbers after filtering the sample at a q-score,
        which are the percent of reads remaining and number of control hits"""
    filtered_pairs = filter_triples(index_triples, quality)
    sample_index_counts = Counter([(i7s, i5s) for name, i7s, i5s in filtered_pairs])
    sample_index_counts2 = [(I7S[counts[0]], I5S[counts[1]], sample_index_counts[counts]) for counts in sample_index_counts]
    control_counts = [sample_index_counts[cp] for cp in CONTROL_PAIRS]
    return (filtered_pairs, sample_index_counts, control_counts)


def counts_at_quality(index_triples, quality):
    #""" Provide the quality and the two key numbers after filtering the sample at a q-score,
    #    which are the percent of reads remaining and number of control hits"""
    #filtered_pairs = filter_triples(index_triples, quality)
    #sample_index_counts = Counter([(i7s, i5s) for name, i7s, i5s in filtered_pairs])
    #control_counts = [sample_index_counts[cp] for cp in CONTROL_PAIRS]
    filtered_pairs, sample_index_counts, control_counts = get_sample_index_filter_lists(index_triples, quality)
    #sample_index_counts2 = [(I7S[counts[0]], I5S[counts[1]], sample_index_counts[counts]) for counts in sample_index_counts]
    #correct_pairs = [(I7I[pair[0]], I5I[pair[1]]) for pair in SAMPLE_SHEET]
    #correct_counts = [sample_index_counts[cp] for cp in correct_pairs]
    #print('{}\t{}\t{}'.format(quality, len(filtered_pairs), filtered_pairs))
    #print('{}'.format(sample_index_counts))
    #print('{}'.format(sample_index_counts2))
    #print('{}'.format(control_counts))
    #print('{}'.format(correct_counts))
    return {
        'phred_score': quality,
        'remaining_reads': len(filtered_pairs),
        'total_reads': len(index_triples),
        'control_hits': sum(control_counts),
        #'misassign_hits': sum(control_counts),
        #'correct_hits': sum(correct_counts),
    }

def get_filtered_sample_ids_correct(index_triples, quality, misassign):
    filtered_pairs, sample_index_counts, control_counts = get_sample_index_filter_lists(index_triples, quality)
    if misassign:
        #print('{}'.format(SAMPLE_SHEET))
        #print('{}'.format(CONTROL_PAIRS))
        #print('{}'.format(I7I))
        #print('{}'.format(I5I))
        correct_pairs = [(I7I[pair[0]], I5I[pair[1]]) for pair in SAMPLE_SHEET]
        correct_counts = [sample_index_counts[cp] for cp in correct_pairs]
        #print('correct_counts: {}'.format(correct_counts))
        #print('correct_pairs({}): {}'.format(len(correct_pairs), correct_pairs))
        correct_ids = [fp for fp in filtered_pairs if (fp[1], fp[2]) in correct_pairs]
        #print('correct_ids({}): {}'.format(len(correct_ids), correct_ids))
        #correct_counts = [sample_index_counts[cp] for cp in SAMPLE_SHEET]
        return correct_ids
    else:
        filtered_ids = [fp for fp in filtered_pairs]
        #print('filtered_ids({}): {}'.format(len(filtered_ids), filtered_ids))
        return filtered_ids

def count_sample(i7fq_fname, i5fq_fname, qual=None, misassign=False):
    """Produce a list of dicts with counts for each key metric"""
    index_triples = list(build_index_tuples(i7fq_fname, i5fq_fname))
    if qual:
        return get_filtered_sample_ids_correct(index_triples, qual, misassign)
        #return counts_at_quality(index_triples, qual)
    else:
        #return filter_triples(index_triples, qual)
        return [counts_at_quality(index_triples, qcut) for qcut in range(MAX_QUAL + 1)]

def filter_fwd_rev_reads(fwd_fq, rev_fq, id_list, out_fwd_fq, out_rev_fq):
    """"""
    with gzip.open(out_fwd_fq, 'wt') as out_fwd:
        with gzip.open(out_rev_fq, 'wt') as out_rev:
            for (fwd, rev) in zip(SeqIO.parse(gzip.open(fwd_fq, 'rt'), 'fastq'), SeqIO.parse(gzip.open(rev_fq, 'rt'), 'fastq')):
                #assert fwd.id == rev.id
                if fwd.id in id_list:
                    SeqIO.write(fwd, out_fwd, 'fastq')
                    SeqIO.write(rev, out_rev, 'fastq')


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='Estimate misassignment rates at different qualities in a dual-index sample from a multi-sample Illumina sequencing run',
    )
    parser.add_argument('-ss', action="store", required=True, dest="ss", help="Sample sheet CSV with true i7/i5 combinations")
    parser.add_argument('-i7', action="store", required=True, dest="i7", help="i7 index reads")
    parser.add_argument('-i5', action="store", required=True, dest="i5", help="i5 index reads")
    parser.add_argument('-q', action='store', required=False, dest='quality', type=int, help="quality filter threshold")
    parser.add_argument('-ma', action="store_true", required=False, dest='misassign', help="also discard misassigned reads")
    parser.add_argument('-fwd', action="store", required=False, dest="fwd", help="fwd reads")
    parser.add_argument('-rev', action="store", required=False, dest="rev", help="rev reads")
    parser.add_argument('-out_fwd', action="store", required=False, dest="out_fwd", help="fwd reads")
    parser.add_argument('-out_rev', action="store", required=False, dest="out_rev", help="rev reads")
    args = parser.parse_args()

    initialize_sample_sheet_globals(args.ss)
    if not args.quality:
        rows = count_sample(args.i7, args.i5)

        w = csv.DictWriter(sys.stdout, ['phred_score', 'remaining_reads', 'total_reads', 'control_hits'])
        #w = csv.DictWriter(sys.stdout, ['phred_score', 'remaining_reads', 'total_reads', 'misassign_hits', 'correct_hits'])
        w.writeheader()
        w.writerows(rows)
    else:
        filtered_list = count_sample(args.i7, args.i5, args.quality, args.misassign)
        #print('{}'.format(filtered_list))
        id_list = [entry[0] for entry in filtered_list]
        filter_fwd_rev_reads(args.fwd, args.rev, id_list, args.out_fwd, args.out_rev)

