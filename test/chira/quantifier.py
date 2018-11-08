import os
import sys
import re
from collections import defaultdict
import pkg_resources

import utilities


def merge_bed(bed, mergedbed, merge_overlap):
    fh_bed = open(bed, "r")
    d_desc = defaultdict(lambda: defaultdict(list))
    for line in fh_bed:
        f = line.rstrip('\n').split("\t")
        d_desc[f[0] + "\t" + f[5]][(int(f[1]), int(f[2]))].append(f[3])
    fh_bed.close()
    print("size new:", sys.getsizeof(d_desc))
    print("Read bed file")

    fh_mergedbed = open(mergedbed, "w")
    for chr in sorted(d_desc.keys()):
        startpos_sorted = sorted(d_desc[chr], key=lambda tup: tup[0])
        t_merged = []
        d_mergeddesc = defaultdict(list)
        for currentpos in startpos_sorted:
            if not t_merged:
                t_merged.append(currentpos)
                d_mergeddesc[currentpos[0]].extend(d_desc[chr][currentpos])
            else:
                prevpos = t_merged[-1]
                if currentpos[0] <= prevpos[1]:
                    overlap_len = min(currentpos[1], prevpos[1]) - currentpos[0] + 1
                    if overlap_len/float(currentpos[1]-currentpos[0]+1) >= merge_overlap \
                            or overlap_len/float(prevpos[1]-prevpos[0]+1) >= merge_overlap:
                        # if overlap_len > 15:
                        end = max(prevpos[1], currentpos[1])
                        t_merged[-1] = (prevpos[0], end)  # replace by merged interval
                        d_mergeddesc[(prevpos[0])].extend(d_desc[chr][currentpos])
                    else:
                        t_merged.append(currentpos)
                        d_mergeddesc[currentpos[0]].extend(d_desc[chr][currentpos])
                        continue
                else:
                    t_merged.append(currentpos)
                    d_mergeddesc[currentpos[0]].extend(d_desc[chr][currentpos])

        for pos in t_merged:
            [chrid, strand] = chr.split("\t")
            fh_mergedbed.write(chrid + "\t" +
                               str(pos[0]) + "\t" +
                               str(pos[1]) + "\t" +
                               strand + "\t" +
                               ";".join(set(d_mergeddesc[pos[0]])) + "\n")
    fh_mergedbed.close()


def maximize(d_tg, d_group_reads, d_read_groups, d_group_total_counts, em_threshold, i=1):
    d_tmp = defaultdict(float)
    print("iteration: "+str(i))
    # print("d_tg keys: " + str(len(d_tg)))
    for both in d_tg.keys():
        [readid, groupid] = both.split("\t")
        group_total_count = 0.0
        read_total_count = 0.0
        for read in set(d_group_reads[groupid]):
            group_total_count += d_tg[read + "\t" + groupid]
        for group in set(d_read_groups[readid]):
            read_total_count += d_group_total_counts[group]
        d_tmp[readid + "\t" + groupid] = group_total_count / read_total_count
    s = 0.0
    for groupid in d_group_reads.keys():
        c = 0.0
        for readid in set(d_group_reads[groupid]):
            c += d_tmp[readid+"\t"+groupid]
        d_group_total_counts[groupid] = c
        s += c

    equal = True
    for k in d_tg:
        if abs(d_tg[k] - d_tmp[k]) >= em_threshold:
            equal = False
    if equal:
        return d_tmp
    else:
        i += 1
        return maximize(d_tmp, d_group_reads, d_read_groups, d_group_total_counts, em_threshold, i)


def create_mmgs(merged_bed, groups_file, crg_share):
    # without d_mmgreads
    fh_merged_bed = open(merged_bed, "r")
    l_locireads = []
    d_readlocus_transcripts = defaultdict(lambda: defaultdict(list))
    l_locipos = []
    l_locilen = []
    print(merged_bed)
    for n_locus, line in enumerate(fh_merged_bed):
        # chr14     64814786    64814804    -   tag_1308593|1|r,ENSMUST00000176386;tag_1308594|2|r,ENSMUST00000176386
        f = line.rstrip('\n').split('\t')
        locuslen = float(f[2])-float(f[1])+1
        l_locilen.append(locuslen)
        pos = ':'.join([f[0], f[1], f[2], f[3]])
        l_locipos.append(pos)
        alignments = f[4].split(';')  # in the description column of the bed file,alignments are seperated by ';'
        l_locusreads = []
        for alignment in alignments:
            readid = alignment.split(',')[0]  # .rsplit('|',1)[0]
            l_locusreads.append(readid)
            transcriptid = alignment.split(',')[1]
            d_readlocus_transcripts[n_locus][readid].append(transcriptid)
        l_locireads.append(set(l_locusreads))
    fh_merged_bed.close()

    # for n_locus, line in enumerate(fh_merged_bed):
    #     # chr14     64814786    64814804    -   tag_1308593|1|r,ENSMUST00000176386;tag_1308594|2|r,ENSMUST00000176386
    #     f = line.rstrip('\n').split('\t')
    #     locuslen = float(f[2])-float(f[1])+1
    #     l_locilen.append(locuslen)
    #     pos = ':'.join([f[0], f[1], f[2], f[3]])
    #     l_locipos.append(pos)
    #     alignments = f[4].split(';')  # in the description column of the bed file,alignments are seperated by ';'
    #     l_locusreads = []
    #     for alignment in alignments:
    #         readid = alignment.split(',')[0]
    #         l_locusreads.append(readid)
    #         transcriptid = alignment.split(',')[1]
    #         d_readlocus_transcripts[n_locus][readid].append(transcriptid)
    #     l_locireads.append(set(l_locusreads))
    # fh_merged_bed.close()

    print("size of d_readlocus_transcripts:", sys.getsizeof(d_readlocus_transcripts))

    print(len(l_locireads), len(d_readlocus_transcripts))

    l_mmgloci = []
    n_group = 0
    print("There are a total of " + str(len(l_locireads)) + " loci")
    print("size of l_locireads: " + str(sys.getsizeof(l_locireads)))
    # creating MMGs is an iterative process.
    # for each locus check how many reads of it are shared among all other loci
    # with significant overlap, merge the loci into MMG
    # a single locus can be part of multiple MMGs
    for n_locus, l_locusreads in enumerate(l_locireads):
        isnewmmglocus = False
        c_mmg_member = 0
        for mmgid, l_loci in enumerate(l_mmgloci):
            for locusid in l_loci:
                n_common_reads = len(set(l_locireads[locusid]) & set(l_locusreads))
                # skip the whole mmg even if one of the locus doesn't overlap enough
                if n_common_reads / float(len(set(l_locireads[locusid]))) < crg_share\
                        or n_common_reads / float(len(set(l_locusreads))) < crg_share:
                    isnewmmglocus = True
                    break
                else:
                    isnewmmglocus = False
            if not isnewmmglocus:
                # print(str(n_locus) + " belongs to " + str(mmgid))
                l_mmgloci[mmgid].append(n_locus)
                c_mmg_member += 1
        # n_locus is not a member of any mmg, hence create a new mmg
        if c_mmg_member == 0:
            # print(str(n_locus) + " makes up a new  " + str(len(l_mmgloci)))
            l_mmgloci.append([n_locus])
            n_group += 1

    print("There are a total of " + str(len(l_mmgloci)) + " mmgs")
    isolated_loci = []  # loci that are separated from MMG because of not enough overall share
    integrated_loci = []  # loci already belong to any MMG
    d_mmg_locus_overlap = defaultdict(lambda: defaultdict(float))
    d_locus_mmg_overlap = defaultdict(lambda: defaultdict(float))

    d_mmg_loci = defaultdict(list)
    d_locus_mmgs = defaultdict(list)
    for mmgid, l_loci in enumerate(l_mmgloci):
        mmg_reads = []
        for locusid in l_loci:
            mmg_reads.extend(l_locireads[locusid])
        for locusid in l_loci:
            locus_share = len(set(l_locireads[locusid]) & set(mmg_reads)) / float(len(set(mmg_reads)))
            # in the end there might be multiple groups with only same single locus because of this eimination
            if locus_share >= crg_share:
                # for each locus choose only one mmg that shares most reads with
                # if there are more than 1 best keep them all
                best_share = True
                for mmg in d_locus_mmg_overlap[locusid]:
                    if d_locus_mmg_overlap[locusid][mmg] > locus_share:
                        best_share = False
                if not best_share:
                    continue
                d_mmg_locus_overlap[mmgid][locusid] = locus_share
                d_locus_mmg_overlap[locusid][mmgid] = locus_share
                d_mmg_loci[mmgid].append(locusid)
                d_locus_mmgs[locusid].append(mmgid)
                # add to integrated list to cross check with isolated list before creating new mmgs
                integrated_loci.append(locusid)
            else:
                # print("locus_"+str(locusid), "group_"+str(mmgid), locus_share)
                # add to isloated list to create new mmgs out of it later
                isolated_loci.append(locusid)

    mmg_index = len(l_mmgloci)
    # every locus that is in isolated_loci makes its own mmg
    for locusid in sorted(set(isolated_loci)):
        if locusid in set(integrated_loci):
            # already member of some mmg, safely ignore it
            continue
        d_mmg_locus_overlap[mmg_index][locusid] = 1.0
        d_locus_mmg_overlap[locusid][mmg_index] = 1.0
        d_mmg_loci[mmg_index].append(locusid)
        d_locus_mmgs[locusid].append(mmg_index)
        mmg_index += 1
        print("Creating a new mmg " + str(mmg_index))

    print("size of d_mmg_locus_overlap: ", sys.getsizeof(d_mmg_locus_overlap))

    mmg_reads = defaultdict(list)
    reads_mmgs = defaultdict(list)
    for mmgid in d_mmg_locus_overlap.keys():
        l_mmgreads = []
        for locusid in d_mmg_locus_overlap[mmgid]:
            l_mmgreads.extend(l_locireads[locusid])
            mmg_reads[mmgid].extend(l_locireads[locusid])
        mmg_reads[mmgid] = set(sorted(l_mmgreads))
        reads_mmgs[tuple(set(sorted(l_mmgreads)))].append(mmgid)

    d_new = defaultdict(lambda: defaultdict(float))
    c = 0
    for l_mmgreads in reads_mmgs.keys():
        for mmgid in reads_mmgs[l_mmgreads]:
            for locusid in d_mmg_loci[mmgid]:
                locus_share = len(set(l_locireads[locusid]) & set(l_mmgreads)) / float(len(set(l_mmgreads)))
                d_new[c][locusid] = locus_share
        c += 1

    #
    #     z = x.copy()
    #     z.update(y)

    # # Finally, merge any two mmgs that have significantly overlapping reads
    # l_group_pairs = list(itertools.combinations(d_mmg_locus_overlap, r=2))
    # d_mmg_locus_overlap_new = defaultdict(lambda: defaultdict(float))
    # for group1, group2 in l_group_pairs:
    #     mmg_reads1 = []
    #     for locusid in d_mmg_locus_overlap[group1]:
    #         mmg_reads1.extend(l_locireads[locusid])
    #     mmg_reads2 = []
    #     for locusid in d_mmg_locus_overlap[group2]:
    #         mmg_reads2.extend(l_locireads[locusid])
    #     mmg_reads = []
    #     if len(set(mmg_reads1) & set(mmg_reads2))/float(len(set(mmg_reads1) | set(mmg_reads2))) >= crg_share:
    #         mmg_reads = set(mmg_reads1) | set(mmg_reads2)
    #
    print("Writing MMGs")
    d_duplicate_entries = {}
    fh_groups_file = open(groups_file, "w")
    for n_group in sorted(d_new.keys()):
        for n_locus in d_new[n_group]:
            # if n_locus not in d_readlocus_transcripts:
            #     continue
            for readid in l_locireads[n_locus]:
                # if there are duplicate entries, they are written to output
                entry = "\t".join([readid,
                                   "locus_" + str(n_locus),
                                   "group_" + str(n_group),
                                   l_locipos[n_locus],
                                   str(d_new[n_group][n_locus]),
                                   ";".join(set(d_readlocus_transcripts[n_locus][readid]))])
                if entry in d_duplicate_entries:
                    print("Duplicate entry found: " + entry)
                    continue
                fh_groups_file.write(entry + "\n")
                d_duplicate_entries[entry] = 1
    fh_groups_file.close()
    print("MMGs written")


def tpm(d_group_expression, d_group_locipos):
    total_rpk = 0
    d_group_tpm = defaultdict(float)
    for mmgid in sorted(d_group_expression.keys()):
        mmg_expression = d_group_expression[mmgid]
        mmg_len = utilities.median(sorted(d_group_locipos[mmgid].values()))/1000.0  # length in kbs
        # print(mmgid, d_group_expression[mmgid], sorted(d_group_locipos[mmgid].values()), mmg_len)
        rpk = mmg_expression/mmg_len
        d_group_tpm[mmgid] = rpk
        total_rpk += rpk
    millions_of_rpk = total_rpk/1000000.0
    for mmgid in sorted(d_group_expression.keys()):
        group_tpm = d_group_tpm[mmgid]/millions_of_rpk
        d_group_tpm[mmgid] = group_tpm
    return d_group_tpm


def quantify_mmgs(loci_groups_file, em_threshold):
    d_group_total_counts = defaultdict(float)
    d_group_reads = defaultdict(list)
    d_read_groups = defaultdict(list)
    d_group_locipos = defaultdict(lambda: defaultdict(str))
    d_read_group_share = defaultdict(float)
    d_group_expression = defaultdict(float)
    fh_loci_groups_file = open(loci_groups_file, "r")
    for line in fh_loci_groups_file:
        f = line.rstrip("\n").split("\t")
        readid = f[0]
        locusid = f[1]
        groupid = f[2]
        [chrom, start, end, strand] = f[3].split(':')
        locuslength = int(end) - int(start) + 1
        # a single locus can belong to multiple mmgs
        # one read can be part of multiple mmgs
        d_group_reads[groupid].append(readid)
        d_read_groups[readid].append(groupid)
        d_group_locipos[groupid][locusid] = locuslength
    fh_loci_groups_file.close()
    print("no. of keys in d_read_groups: " + str(len(d_read_groups)))

    for readid in sorted(d_read_groups.keys()):
        for groupid in set(d_read_groups[readid]):
            d_read_group_share[readid + "\t" + groupid] = 1 / float(len(set(d_read_groups[readid])))

    for groupid in sorted(d_group_reads.keys()):
        count = 0.0
        for readid in sorted(set(d_group_reads[groupid])):
            count += d_read_group_share[readid + "\t" + groupid]
        d_group_total_counts[groupid] = count

    sys.setrecursionlimit(100)
    d_res = maximize(d_read_group_share, d_group_reads, d_read_groups, d_group_total_counts, em_threshold)

    for readid in sorted(d_read_groups.keys()):
        # now d_read_groups should contain only multi mapped reads because because
        # uniquely mapped reads were already removed
        # set() because a read can occur 2 times at a locus in file with 2 tx ids
        # one read can be part of multiple mmgs
        for multimapped_group in set(d_read_groups[readid]):
            d_group_expression[multimapped_group] += d_res[readid+"\t"+multimapped_group]

    d_group_tpm = tpm(d_group_expression, d_group_locipos)

    return d_group_expression, d_group_tpm