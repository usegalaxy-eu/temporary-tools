import os
import sys
import re
from collections import defaultdict
import pysam
import pkg_resources
# pkg_resources.require("pysam==0.8.3")


def file_prefixes(sampledir, firstpart, secondpart):
    if firstpart == "miRNA":
        parentdir = os.path.join(sampledir, 'miRNA_vs_transcriptome')
        first_prefix = os.path.join(parentdir, 'vs_miRNA')
        pairs_prefix = os.path.join(parentdir, 'miRNA')
        if secondpart == "miRNA":
            second_prefix = os.path.join(parentdir, 'clipped_vs_mirna')
            pairs_prefix += '_vs_miRNA'
        elif secondpart == "transcriptome":
            second_prefix = os.path.join(parentdir, 'clipped_vs_transcriptome')
            pairs_prefix += '_vs_transcriptome'
        else:
            sys.exit("Unknown hybrid first part")
    elif firstpart == "transcriptome":
        parentdir = os.path.join(sampledir, 'transcriptome_vs_transcriptome')
        first_prefix = os.path.join(parentdir, 'vs_transcriptome')
        pairs_prefix = os.path.join(parentdir, 'transcriptome')
        if secondpart == "transcriptome":
            second_prefix = os.path.join(parentdir, 'clipped_vs_transcriptome')
            pairs_prefix += '_vs_transcriptome'
        else:
            sys.exit("Unknown hybrid second part")
    else:
        sys.exit("Unknown hybrid part")

    return first_prefix, second_prefix, pairs_prefix


def transcript_to_genomic_pos(transcripts_bed, f_geneexonbed, f_txexonbed, p_intersectbed):
        id2pol = {}
        id2chr = {}
        exid2start = {}
        exid2end = {}
        print("Read in genomic exon features .. ")
        fh_gene_exon_bed = open(f_geneexonbed, "r")
        for line in fh_gene_exon_bed:
            # chrName, s, e, exid, pol = line.rstrip('\n').split('\t')[0,1,2,3,5]
            f = line.rstrip('\n').split('\t')
            chrname = f[0]
            s = int(f[1])
            e = int(f[2])
            exid = f[3]
            pol = f[5]
            match = re.match("(.+?)_e", exid)
            transcriptid = match.group(1)
            id2pol[transcriptid] = pol
            id2chr[transcriptid] = chrname
            exid2start[exid] = s
            exid2end[exid] = e
        fh_gene_exon_bed.close()
        print("done")
        print("Transcripts with exons:   " + str(len(id2pol)))

        overlapout = transcripts_bed.replace(".bed", ".overlap.txt")
        wojunctions = transcripts_bed.replace(".bed", ".genomic.bed")

        print("Calculating bed overlap .. ")
        intersect_command = (p_intersectbed + " -a " + transcripts_bed + " -b " + f_txexonbed + " -wb > " + overlapout)
        print(intersect_command)
        os.system(intersect_command)

        fh_overlapout = open(overlapout, "r")
        n_out = 0
        n_plus = 0
        n_minus = 0

        d_idsseen = {}
        d_idsdouble = {}
        d_withjunctions = defaultdict(list)

        # ENSMUST00000023614      5301    5319    2175|tag_1004650|1,ENSMUST00000023614,ENSMUSG00000022897,5301,5319 \
        #       0       +       ENSMUST00000023614      1817    5754    ENSMUST00000023614_e012 0       +
        for line in fh_overlapout:
            # transcriptid, s, e, readid, exonstart, exonend, exonid, pol
            f = line.rstrip('\n').split('\t')
            transcriptid = f[0]
            txstart = int(f[1])
            txend = int(f[2])
            readid = f[3]
            exonstart = int(f[7])
            exonid = f[9]
            pol = f[11]
            # Calculate genomic hit positions.
            # Plus strand case.
            if pol == "+":
                genomicstart = exid2start[exonid] + (txstart - exonstart)
                genomicend = exid2start[exonid] + (txend - exonstart)
                n_plus += 1
            # Minus strand case.
            elif pol == "-":
                genomicstart = exid2end[exonid] - (txend - exonstart)
                genomicend = exid2end[exonid] - (txstart - exonstart)
                n_minus += 1
            else:
                sys.exit("ERROR: Check polarity " + pol + "in line: " + line)

            # Store ID.
            if readid in d_idsseen:
                d_idsdouble[readid] = d_idsdouble.get(readid, 0) + 1
            else:
                d_idsseen[readid] = d_idsseen.get(readid, 0) + 1

            n_out += 1
            d_withjunctions[readid].append(id2chr[transcriptid] + "\t" +
                                           str(genomicstart) + "\t" +
                                           str(genomicend) + "\t" +
                                           readid +
                                           "\t0\t" +
                                           pol + '\n')
            if str(id2pol[transcriptid]) != pol:
                sys.exit("Different strands!!")
        fh_overlapout.close()

        n_double = len(d_idsdouble)

        print("Total alignments:        " + str(n_out))
        print("Hits on plus strand:     " + str(n_plus))
        print("Hits on minus strand:    " + str(n_minus))
        print("Hits over exon borders:  " + str(n_double))

        # Write alignments spanning exon borders to a file
        fh_wojunctions = open(wojunctions, "w")
        n_final = 0
        print("Filtering out exon border reads .. ")
        for readid, lines in d_withjunctions.items():
            # Alignments over exon borders have more than 2 lines in hash
            # This way we can identify them and filter them out.
            # TODO: Considering them as 2 seperate alignments at genomic level doesn't harm. Just need to make sure that
            # TODO they should not be treated as multi-mapped read. Give them a special ID while quantification
            if len(lines) > 1:
                # print readid
                continue
            n_final += 1
            for line in lines:
                fh_wojunctions.write(line)
        fh_wojunctions.close()

        print("done")
        print("Alignments after filtering:  " + str(n_final))
        os.system("rm " + overlapout)

        return


def get_numberof_matches(alignment):
    # sum up all aligned bases and subtract the number of mismatches
    # insertions and deletions are not considered here
    alignedbases = sum(x[1] for x in alignment.cigar if x[0] == 0)
    mismatches = alignment.opt("NM")
    n_matches = alignedbases - mismatches
    return n_matches


def get_mapped_parts(l_matchedlengths):
    l_sorted_indices = sorted(range(len(l_matchedlengths)), key=lambda x: l_matchedlengths[x])
    c = 1
    prev = 0
    max1_indices = []
    max2_indices = []
    for i in reversed(l_sorted_indices):
        if prev == 0:
            max1_indices.append(i)
        elif l_matchedlengths[i] == prev and c == 1:
            max1_indices.append(i)
        elif l_matchedlengths[i] == prev and c == 2:
            max2_indices.append(i)
        elif l_matchedlengths[i] != prev:
            if c == 2:
                break
            max2_indices.append(i)
            c += 1
        prev = l_matchedlengths[i]
    if len(max1_indices) > 1:
        max2_indices = max1_indices

    dist = 0
    left = ""
    right = ""
    for i in max1_indices:
        for j in max2_indices:
            if abs(i - j) > dist:
                dist = abs(i - j)
                left = i
                right = j
                if i > j:
                    left = j
                    right = i
    return left, right


def overlap(f, s):
    return max(0, min(f[1], s[1]) - max(f[0], s[0]))


def overlap_on_same_reference(first_read_pos, second_read_pos, min_looplen):
    overlapping = False
    f = first_read_pos.split(':')
    s = second_read_pos.split(':')
    reference1 = f.pop(0)
    reference2 = s.pop(0)
    f = list(map(int, f))
    s = list(map(int, s))
    if reference1 == reference2:  # same chromosomes
        # 3 is the minimum loop length
        if min(f[1], s[1]) - max(f[0], s[0]) >= -min_looplen:
            overlapping = True
        # else:
        #     print("overlapping on same reference", first_read_pos, second_read_pos)
    return overlapping


def get_read_pos(bam):
    d_left_readpositions = defaultdict(list)
    d_right_readpositions = defaultdict(list)
    fh_bam = pysam.Samfile(bam, "rb")
    for alignment in fh_bam.fetch(until_eof=True):
        readid = alignment.qname
        tagid = readid.split('|', 1)[1]
        # these are the positions on the aligned part of the read, not reference positions!
        read_start = alignment.query_alignment_start
        read_end = alignment.query_alignment_end
        referenceid = fh_bam.getrname(alignment.tid)
        # NOTE: one of the following should work
        if readid.split('|')[-1] == "l":
            d_left_readpositions[tagid].append(referenceid
                                               + ":" + str(read_start)
                                               + ":" + str(read_end))
        elif readid.split('|')[-1] == "r":
            d_right_readpositions[tagid].append(referenceid
                                                + ":" + str(read_start)
                                                + ":" + str(read_end))
        # TODO: beware, this might be a single RNA duplex, just don't ignore
        else:
            if readid == "4715|tag_6262|19":
                print("here: ", readid, referenceid, str(read_start),  str(read_end))
    fh_bam.close()
    return d_left_readpositions, d_right_readpositions


def median(x):
    n = len(x)
    mid = int(n/2)
    if not n%2:
        return (x[mid-1] + x[mid]) / 2.0
    return x[mid]


def chunks(ids, n):
    return [ids[i:i+n] for i in range(0, len(ids), n)]


def parse_bed(bed):
    d_read_genomic_pos = defaultdict(lambda: defaultdict(str))
    d_read_tx_pos = defaultdict(lambda: defaultdict(str))
    fh_bed = open(bed, "r")
    for line in fh_bed:
        f = line.rstrip('\n').split('\t')
        pos = ':'.join([f[0], f[1], f[2], f[5]])
        desc = f[3].split(',')  # in the description column of the bed file,alignments are seperated by ';'
        readid = desc[0]
        transcriptid = desc[1]
        # at this level reads have unique ids preceeded by a serialnumber 
        d_read_genomic_pos[readid][transcriptid] = pos
        d_read_tx_pos[readid][transcriptid] = desc[3] + "\t" + desc[4]
    return d_read_genomic_pos, d_read_tx_pos


def parse_files(counts_file, groups_file):
    d_group_expression = {}
    d_group_tpm = {}
    print(counts_file)
    fh_counts = open(counts_file, "r")
    for line in fh_counts:
        f = line.rstrip('\n').split('\t')
        groupid = f[0]
        count = f[1]
        group_tpm = f[2]
        d_group_expression[groupid] = count
        d_group_tpm[groupid] = float(group_tpm)
    fh_counts.close()

    d_readgroups = defaultdict(list)
    d_grouploci = defaultdict(list)
    d_locipos = {}
    d_locusgroup_share = defaultdict(lambda: defaultdict(float))
    d_locus_transcripts = defaultdict(lambda: defaultdict(list))
    fh_groups = open(groups_file, "r")
    for line in fh_groups:
        f = line.rstrip('\n').split('\t')
        readid = f[0]
        locusid = f[1]
        groupid = f[2]
        locuspos = f[3]
        locusshare = f[4]
        # a single locus can be part of 2 different mmgs with different expression levels
        d_readgroups[readid].append(groupid)
        d_grouploci[groupid].append(locusid)
        d_locipos[locusid] = locuspos
        d_locusgroup_share[locusid][groupid] = locusshare
        d_locus_transcripts[readid][locusid].extend(f[5].split(";"))
    fh_counts.close()

    return d_group_expression, d_group_tpm, d_readgroups, d_grouploci, d_locipos, d_locusgroup_share, d_locus_transcripts


def pairspos_to_bed(pairs_prefix, bed1, bed2):
    fh_readpairs = open(pairs_prefix + ".loci.pairs", "r")
    fh_bed1 = open(bed1, "w")
    fh_bed2 = open(bed2, "w")
    for line in fh_readpairs:
        f = line.rstrip('\n').split('\t')
        chimeraid = f[0]
        pos1 = f[16]
        f1 = pos1.split(':')
        # TODO : add 20nt flanking
        f1[1] = str(int(f1[1]) - 0)
        f1[2] = str(int(f1[2]) + 0)
        f1.insert(3, chimeraid)
        f1.insert(4, "0")
        fh_bed1.write('\t'.join(f1)+"\n")

        pos2 = f[17]
        f2 = pos2.split(':')
        # TODO : add 20nt flanking
        f2[1] = str(int(f2[1]) - 0)
        f2[2] = str(int(f2[2]) + 0)
        f2.insert(3, chimeraid)
        f2.insert(4, "0")
        fh_bed2.write('\t'.join(f2)+"\n")

    fh_readpairs.close()
    fh_bed1.close()
    fh_bed2.close()
    return


def run_intarna(d_seq_records1, d_seq_records2, ids, matrix_file):
    os.system("[ -e " + matrix_file + "] && rm " + matrix_file + ";  touch " + matrix_file)
    for seqId in ids:
        record1 = d_seq_records1[seqId]
        record2 = d_seq_records2[seqId]
        os.system("python " +
                  os.path.dirname(os.path.realpath(__file__)) + "/intarna_predict_hybrid_return_base_pairs.py " +
                  " -m " + str(record1.seq) +
                  " -t " + str(record2.seq) +
                  " >> " + matrix_file)
    return


def run_hybridmin(d_seq_records1, d_seq_records2, ids, hybrid_file, pairs_prefix, p_hybmin):
    fh_hybrids = open(hybrid_file, "w")
    for seqId in ids:
        record1 = d_seq_records1[seqId]
        record2 = d_seq_records2[seqId]
        hybrid_prefix = pairs_prefix + "." + seqId
        with open(hybrid_prefix + ".s1", "w") as f1:
            f1.write(str(record1.seq))
        with open(hybrid_prefix + ".s2", "w") as f2:
            f2.write(str(record2.seq))
        os.system(p_hybmin + " -o " + hybrid_prefix + " " + hybrid_prefix + ".s1 " + hybrid_prefix + ".s2 > /dev/null") # + hybrid_prefix + ".log"

        mfe = None
        with open(hybrid_prefix + ".dG") as fh_dg:
            next(fh_dg)
            for line in fh_dg:
                mfe = line.split("\t")[1]
        i = []
        j = []
        with open(hybrid_prefix + ".37.plot") as fh_plot:
            next(fh_plot)
            for line in fh_plot:
                f = line.split("\t")
                i.append(int(f[0]))
                j.append(int(f[1]))
        dotbracket = ""
        for index in range(1, len(str(record1.seq))+1):
            if index in i:
                dotbracket += "("
            else:
                dotbracket += "."
        dotbracket += "&"
        for index in range(1, len(str(record2.seq))+1):
            if index in j:
                dotbracket += ")"
            else:
                dotbracket += "."
        if not i:
            i.append(-1)
            j.append(-1)
        fh_hybrids.write("\t".join([seqId,
                                   str(record1.seq) + "&" + str(record2.seq),
                                   dotbracket,
                                   str(i[0]) + "&" + str(j[-1]),
                                   mfe,
                                   "\n"]))
        os.system("rm " + hybrid_prefix + ".*")
    fh_hybrids.close()
    return


def populate_db_dicts(d_txid_geneid, d_geneid_genetype, d_geneid_genename,
                      d_txid_len,
                      d_txid_goterms,
                      d_mirna_family,
                      d_mature_mir_pos,
                      f_txgenemap,
                      f_txexonbed,
                      f_goterms,
                      f_mirfammap,
                      f_mirgtf):

    fh_txgenemap = open(f_txgenemap, "r")
    for line in fh_txgenemap:
        fields = line.rstrip('\n').split('\t')
        txid = fields[0]
        geneid = fields[1]
        genetype = fields[2]
        genename = fields[3]
        d_txid_geneid[txid] = geneid
        d_geneid_genetype[geneid] = genetype
        d_geneid_genename[geneid] = genename
    fh_txgenemap.close()

    fh_txexonbed = open(f_txexonbed, "r")
    for line in fh_txexonbed:
        fields = line.rstrip('\n').split('\t')
        txid = fields[0]
        txlen = fields[2]
        d_txid_len[txid] = txlen
    fh_txexonbed.close()

    fh_go_terms = open(f_goterms, "r")
    for line in fh_go_terms:
        fields = line.rstrip('\n').split('\t')
        geneid = fields[0]
        txid = fields[1]
        goterm = fields[2]
        d_txid_goterms[txid+"\t"+geneid].append(goterm)
    fh_go_terms.close()

    fh_mirna_family = open(f_mirfammap, "r")
    for l in fh_mirna_family:
        fields = l.rstrip('\n').split('\t')
        mirnaid = fields[0] + "_" + fields[1]
        familyid = fields[3]
        # if familyid == "NONE":
        #     d_mirna_family[mirnaid] = mirnaid
        # else:
        d_mirna_family[mirnaid] = familyid
    fh_mirna_family.close()

    fh_mirbasegtf = open(f_mirgtf, "r")
    for line in fh_mirbasegtf:
        if line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        if fields[2] == "miRNA_primary_transcript":
            continue
        chrname = fields[0]
        startpos = fields[3]
        endpos = fields[4]
        strand = fields[6]
        description = fields[8]
        m = re.match(".+Name=(.+);Derives_from=(.+)", description)
        maturename = m.group(1)
        precursorid = m.group(2)
        if precursorid not in d_mature_mir_pos:
            d_mature_mir_pos[precursorid] = {}
        d_mature_mir_pos[precursorid][maturename.lower()] = chrname + '\t' + startpos + '\t' + endpos + '\t' + strand
    fh_mirbasegtf.close()


