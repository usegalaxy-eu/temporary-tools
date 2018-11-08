# coding=utf-8
import utilities
import quantifier
# import multi_sample

import argparse
import os
import sys
from collections import defaultdict
import pysam
# from Bio import SeqIO
import multiprocessing
import logging
import pkg_resources
import itertools
# pkg_resources.require("pysam==0.8.3")

p_bowtie2 = "/usr/local/tools/bowtie2/2.2.6/iuc/package_bowtie_2_2_6/6d6cca69a34a/bowtie2"
p_intersectbed = "/usr/local/tools/bedtools/default/bin/intersectBed"
p_fastafrombed = "/usr/local/tools/bedtools/default/bin/fastaFromBed"
p_hybmin = "/home/videmp/install/oligoarrayaux-3.8/build/bin/hybrid-min"

# Mouse db
f_miridx = "/data/5/galaxy_import/galaxy_user_data/pospisilik/miconomou/eclash/db/miRNA/miRNA"
# transcriptome index not intersected with miRNAs
f_txidx = "/data/5/galaxy_import/galaxy_user_data/pospisilik/miconomou/eclash/db/ensembl/transcriptome"
f_txgenemap = "/data/5/galaxy_import/galaxy_user_data/pospisilik/miconomou/eclash/db/whole_transcriptome.map"
# f_txgenemap = "/data/5/galaxy_import/galaxy_user_data/pospisilik/miconomou/eclash/db/whole_transcriptome_w_tx_info.map"
f_mirgtf = "/data/5/galaxy_import/galaxy_user_data/pospisilik/miconomou/eclash/db/miRNA/mmu.gff3"
f_genomicfasta = "/data/db/reference_genomes/mm10/seq/mm10.fa"
f_geneexonbed = "/data/5/galaxy_import/galaxy_user_data/pospisilik/miconomou/eclash/db/genomic_exons.filtered.bed"
f_mirfammap = "/data/5/galaxy_import/galaxy_user_data/pospisilik/miconomou/eclash/db/miRNA/miRNA_miFAM.map.new"
f_txexonbed = "/data/5/galaxy_import/galaxy_user_data/pospisilik/miconomou/eclash/db/transcriptomic_exons.filtered.bed"
f_goterms = "/data/5/galaxy_import/galaxy_user_data/pospisilik/miconomou/eclash/db/mm10_biomart_goterms.txt"

# Human db
# f_miridx = "/scratch/bi02/videmp/ProjectPospisilik/db/hg38/miRNA/miRNA"
# f_txidx = "/scratch/bi02/videmp/ProjectPospisilik/db/hg38/transcriptome/transcriptome"
# f_txgenemap = "/scratch/bi02/videmp/ProjectPospisilik/db/hg38/whole_transcriptome.map"
# f_mirgtf = "/scratch/bi02/videmp/ProjectPospisilik/db/hg38/miRNA/hsa.gff3"
# f_genomicfasta = "/data/db/reference_genomes/hg38/seq/hg38.fa"
# f_geneexonbed = "/scratch/bi02/videmp/ProjectPospisilik/db/hg38/genomic_exons.bed"
# f_mirfammap = "/scratch/bi02/videmp/ProjectPospisilik/db/hg38/miRNA/miRNA_miFAM.map"
# f_txexonbed = "/scratch/bi02/videmp/ProjectPospisilik/db/hg38/transcriptomic_exons.bed"
# f_goterms = "/scratch/bi02/videmp/ProjectPospisilik/db/hg38/hg38_biomart_goterms.txt"

d_txid_geneid = {}  # maps [transcript id]-->[gene id]
d_geneid_genetype = {}  # maps [gene id]-->[gene biotype]
d_geneid_genename = {}  # maps [gene id]-->[gene name]
d_low_support_transcripts = {}
d_mirna_family = {}  # maps [miRNA id]-->[miRNA family id]
d_txid_len = {}
d_txid_goterms = defaultdict(list)
d_mature_mir_pos = {}

utilities.populate_db_dicts(d_txid_geneid, d_geneid_genetype, d_geneid_genename, d_txid_len,
                            d_txid_goterms,d_mirna_family, d_mature_mir_pos,
                            f_txgenemap, f_txexonbed, f_goterms, f_mirfammap, f_mirgtf)


class ChiRA:

    def __init__(self, insample=None, outputdir=None, threads=None, stranded=None,
                 min_looplen=None, tpm_cutoff=None,score_cutoff=None, merge_overlap=None,
                 crg_share=None, em_threshold=None, so_len_threshold=None, hybridize=None):
        self.insample = insample
        self.outputdir = outputdir
        self.threads = threads
        self.stranded = stranded
        self.min_looplen = min_looplen
        self.tpm_cutoff = tpm_cutoff
        self.score_cutoff = score_cutoff
        self.merge_overlap = merge_overlap
        self.crg_share = crg_share
        self.em_threshold = em_threshold
        self.so_len_threshold = so_len_threshold
        self.hybridize = hybridize

    def map_reads(self, mappedto, mapagainst):
        """
            Funtion that maps the reads to the miRNAs or transcriptome. Different parameters
            are used for miRNA and transcriptome mapping.
        """

        if mappedto:
            logging.info("Mapping " + mappedto + " mapped softclipped sequences to " + mapagainst)
        else:
            logging.info("Mapping reads to " + mapagainst)

        alignment_fileprefix = unmapped_fasta = log = refindex = None
        # f(x) = 2 + 0.2 * sqrt(x) where x is the length of the aligned read part.
        # Sets a function governing the interval between seed substrings to use during multiseed alignment
        i = "S,1,1.5"
        rdg = "5,30"  # read gap costs, huge penalties for extension leads to softclipped alignments instead of gaps
        rfg = "5,3"  # reference gap open and extension penalties
        n = "1"
        if not mappedto:
            if mapagainst == "miRNA":
                analysisdir = os.path.join(self.outputdir, 'miRNA_vs_transcriptome')
                query_fasta = self.insample
                alignment_fileprefix = os.path.join(analysisdir, 'vs_miRNA')
                unmapped_fasta = os.path.join(analysisdir, 'vs_miRNA.unmapped.fasta')
                log = os.path.join(analysisdir, 'vs_miRNA.bowtie2.log')
                refindex = f_miridx
                # i = "S,0,0"  # f(x) = 1 + 0.2 * sqrt(x) where x is the length of the aligned read part
                # rdg = "5,1"  # read gap open and extension penalties
                # rfg = "5,1"  # reference gap open and extension penalties
            elif mapagainst == "transcriptome":
                analysisdir = os.path.join(self.outputdir, 'transcriptome_vs_transcriptome')
                query_fasta = os.path.join(self.outputdir, 'miRNA_vs_transcriptome', 'vs_miRNA.unmapped.fasta')
                alignment_fileprefix = os.path.join(analysisdir, 'vs_transcriptome')
                unmapped_fasta = os.path.join(analysisdir, 'vs_transcriptome.unmapped.fasta')
                log = os.path.join(analysisdir, 'vs_transcriptome.bowtie2.log')
                refindex = f_txidx
            else:
                system.exit("ERROR: map against unknown!")

        elif mappedto == "miRNA":
            analysisdir = os.path.join(self.outputdir, 'miRNA_vs_transcriptome')
            if mapagainst == "transcriptome":
                # TODO: change this back or does it make sense?
                # query_fasta = os.path.join(analysisdir, "vs_miRNA.miRNA.nonChimeras.fa")
                query_fasta = os.path.join(analysisdir, "vs_miRNA.clipped.fa")
                # map the soft clipped reads which were not mapped to miRNAs.
                # i.e., first preference for the miRNA targets then transcriptome targets
    #            query_fasta = os.path.join(analysisdir, "clipped_vs_mirna.unmapped.fasta")
                alignment_fileprefix = os.path.join(analysisdir, 'clipped_vs_transcriptome')
                unmapped_fasta = os.path.join(analysisdir, 'clipped_vs_transcriptome.unmapped.fasta')
                log = os.path.join(analysisdir, 'clipped_vs_transcriptome.bowtie2.log')
                refindex = f_txidx
                n = "0"
                i = "L,0,0.3"
                rdg = "5,3"
                rfg = "5,3"
            else:
                system.exit("ERROR: map against unknown!")

        opt_stranded = ""
        if self.stranded == "fw":
            opt_stranded = "--norc"
        elif self.stranded == "rc":
            opt_stranded = "--nofw"

        bowtie2command = \
            (p_bowtie2 +
                ' --local'
                ' --ignore-quals ' +
                opt_stranded +
                ' -k 100'
                # TODO: report all alignments otherwise -k 100 to report top 100 alignments
                # ' -a'
                ' -D 20'
                ' -R 3'
                ' -N ' + n +  # number of mismatches to allowed in a seed alignment
                ' -L 18'  # seed substring length. For example set to 10 to increase sensityvity, but takes longer
                ' --ma 1'  # match bonus
                ' --np 0'
                ' --mp 1,1'  # mismatch penalty, same as match bonus
                ' --rdg ' + rdg +
                ' --rfg ' + rfg +
                ' -i ' + i +
                ' --score-min L,18,0' +  # simply set to 18 matches at least
                ' -p ' + str(self.threads) +
                ' -x ' + refindex +
                ' -f -U ' + query_fasta +
                ' -S ' + alignment_fileprefix + ".sam"
                ' --un ' + unmapped_fasta +
                '  2>&1 | tee ' + log)
        print(bowtie2command)
        os.system(bowtie2command)

        # SAM to BAM --> sort --> index, remove SAM, remove unsorted BAM
        # Call it .temp.bam, because we process it further
        print(pysam.__version__)
        pysam.view("-b", "-S", "-F", "4", "-h",
                   "-o", alignment_fileprefix + ".unsorted.bam",
                   alignment_fileprefix + ".sam",
                   catch_stdout=False)
        os.system("rm " + alignment_fileprefix + ".sam")
        pysam.sort("-m", "1G", "-@", str(self.threads),
                   alignment_fileprefix + ".unsorted.bam",
                   "-T", alignment_fileprefix + ".temp",
                   "-o", alignment_fileprefix + ".temp.bam")  # TODO: 1G memory used!
        pysam.index(alignment_fileprefix + ".temp.bam")
        os.system("rm " + alignment_fileprefix + ".unsorted.bam")

        logging.info("Done")
        return

    def write_softclipped(self, mappedto):
        logging.info("Writing soft clipped parts of " + mappedto + " mapped reads")

        n_total_alignments = 0
        n_softclipped_alignments = 0
        n_discarded_alignments = 0
        n_singleton_alignments = 0

        outputdir = ""
        file_name = ""
        if mappedto == "miRNA":
            outputdir = os.path.join(self.outputdir, 'miRNA_vs_transcriptome')
            file_name = "vs_miRNA"
        elif mappedto == "transcriptome":
            outputdir = os.path.join(self.outputdir, 'transcriptome_vs_transcriptome')
            file_name = "vs_transcriptome"
        else:
            print("Mapped to unknown!")
            exit()

        # TODO: check whether a read gives both intra and inter molecular hybrids
        # Contains [readid"\t"referenceid]-->[longesetClippedSeq]. Used to write the fasta file  with clipped sequences.
        # readid is always prefixed with a counter, inorder to give unique id for multiply mapped reads
        d_clipped_alignment_seq = defaultdict(str)
        # To consider best alignment for same read mapped multiple times on a single reference
        d_read_reference_score_inter = defaultdict(float)
        d_read_reference_score_intra = defaultdict(float)
        d_read_reference_mismatches_intra = defaultdict(float)
        # Stores [readid"\t"referenceid]--> [alignmentrefstart+"\t"+alignedcigar], which have soft clipped squences
        d_filtered_alignments = defaultdict(str)
        l_aligned_tags = []
        l_softclipped_tags = []
        l_singleton_tags = []

        first_prefix = os.path.join(outputdir, file_name)
        first_bam = first_prefix + '.temp.bam'
        filtered_bam = first_prefix + '.filtered.bam'
        unmapped_fasta = first_prefix + '.unmapped.fasta'
        clipped_fasta = first_prefix + '.clipped.fa'
        clip_log = first_prefix + '.clip.log'

        # make a copy of low support transcripts which is used to avoid printing of
        # warning message about same transcript multiple times.
        # NOTE: Let the quantification deal with the low support transcripts. Do not read them to any dictionary
    #    lowSupportTranscriptsCopy = d_low_support_transcripts

        fh_first_bam = pysam.AlignmentFile(first_bam, "rb")

        for alignment in fh_first_bam.fetch():  # until_eof=True also fetches the unmapped reads
            readid = alignment.qname
            alignedcigar = alignment.cigarstring
            cigartuples = alignment.cigar
            readseq = alignment.seq
            referenceid = fh_first_bam.getrname(alignment.tid)
            alignscore = alignment.opt("AS")
            alignmentrefstart = alignment.reference_start

            n_total_alignments += 1
            l_aligned_tags.append(readid)

            # cigar types: M=0, I=1, D=2, S=4
            # skip to next alignment if cigar has no softclipping
            # soft clipping can be seen in ends of cigar
            # if cigartuples[0][0] != 4 and cigartuples[-1][0] != 4:
            #     continue
            issoftclipped = False
            leftclippedseq = ""
            rightclippedseq = ""
            # consider alignment only if there are an overall at least 18 matches to reference
            if utilities.get_numberof_matches(alignment) < 18:
                continue
            pos = 0
            for (cigartype, cigarlength) in cigartuples:
                # trim out the "S" part if it is longer than 18nt
                # set issoftclipped to True only if softclipped and have length > 18
                # In some cases there are more than one softclipped parts, in that case find and choose the longest
                # eg: 18S30M23S
                if int(cigartype) == 4 and int(cigarlength) >= 18:
                    if pos > 0:
                        rightclippedseq = readseq[pos:pos+cigarlength]
                    else:
                        leftclippedseq = readseq[pos:pos+cigarlength]
                    issoftclipped = True
                if int(cigartype) != 2:
                    pos += int(cigarlength)

            # cigar and alignment start uniquely defines each alignment
            alignment_key = "\t".join([readid, referenceid, str(alignmentrefstart), alignedcigar])

            if issoftclipped:
                n_softclipped_alignments += 1
                # an alignment can have 2 valid (>18nt) soft clipped parts (20S33M22S), bu these cases are rare
                # in this case choose the longest soft clipped part
                # TODO: instead of choosing longest soft clipped, consider both possibilities
                longestclippedseq = rightclippedseq
                mappedpart = "l"
                if len(leftclippedseq) > len(rightclippedseq):
                    longestclippedseq = leftclippedseq
                    mappedpart = "r"
                d_clipped_alignment_seq[alignment_key] = longestclippedseq
                d_filtered_alignments[alignment_key] = mappedpart
                l_softclipped_tags.append(readid)
            else:
                # if there is no soft clipping, it could have intra molecular hybrid
                # count number of matched parts with length > 18
                # consider it as intra-molecuar interaction if there are at least 2 significant matched parts
                l_matchedlengths = []
                for (cigartype, cigarlength) in cigartuples:
                    if int(cigartype) == 0 and int(cigarlength) >= 18:
                        l_matchedlengths.append(cigarlength)
                if len(l_matchedlengths) < 2:
                    n_discarded_alignments += 1
                    continue
                left_mapped_part, right_mapped_part = utilities.get_mapped_parts(l_matchedlengths)
                # sort and take best 2 matched lengths and sum them up.
                # TODO: consider NM tag, but with more than 3 matched parts it is hard to distribute NM
                total_matched_length = left_mapped_part + right_mapped_part
                mismatches = alignment.opt("NM")
                if readid + "\t" + referenceid in d_read_reference_score_intra:
                    if int(total_matched_length) < int(d_read_reference_score_intra[readid + "\t" + referenceid]):
                        n_discarded_alignments += 1
                        continue
                    elif int(total_matched_length) == int(d_read_reference_score_intra[readid + "\t" + referenceid]):
                        if int(mismatches) < int(d_read_reference_mismatches_intra[readid + "\t" + referenceid]):
                            n_discarded_alignments += 1
                            continue
                d_read_reference_score_intra[readid + "\t" + referenceid] = total_matched_length
                d_read_reference_mismatches_intra[readid + "\t" + referenceid] = mismatches
                # unlike inter-molecular hybrids, these have more than 2 mapped parts
                d_filtered_alignments[alignment_key] = None
                # d_filtered_alignments_intra[readid + "\t" +referenceid] = ["\t".join([str(left_mapped_part),
                #                                                                   alignedcigar,
                #                                                                    "l"]),
                #                                                            "\t".join([str(right_mapped_part),
                #                                                                       alignedcigar,
                #                                                                       "r"])]
        fh_first_bam.close()

        d_read_reference_score_inter.clear()
        d_read_reference_score_intra.clear()

        fh_first_bam = pysam.AlignmentFile(first_bam, "rb")
        fh_filtered_bam = pysam.AlignmentFile(filtered_bam, "wb", template=fh_first_bam)
        fh_clipped_fasta = open(clipped_fasta, "w")
        counter = 1
        for alignment in fh_first_bam.fetch():
            readid = alignment.qname
            referenceid = fh_first_bam.getrname(alignment.tid)
            alignedcigar = alignment.cigarstring
            alignmentrefstart = alignment.reference_start
            alignment_key = "\t".join([readid, referenceid, str(alignmentrefstart), alignedcigar])
            # consider alignments that are filtered
            if alignment_key in d_filtered_alignments:
                mappedpart = d_filtered_alignments[alignment_key]
                if mappedpart:
                    unmappedpart = "r"
                    if mappedpart == "r":
                        unmappedpart = "l"
                    alignment.qname = str(counter) + "|" + readid + '|' + mappedpart
                    fh_clipped_fasta.write('>' + str(counter) + "|" + readid + '|' + unmappedpart + "\t" + referenceid
                                           + '\n' + d_clipped_alignment_seq[alignment_key] + '\n')
                # if there is no softclip there can be multiple aligned parts
                else:
                    alignment.qname = str(counter) + '|' + readid
                fh_filtered_bam.write(alignment)
                counter += 1
            else:
                n_singleton_alignments += 1
                l_singleton_tags.append(readid)

        fh_first_bam.close()
        fh_filtered_bam.close()
        fh_clipped_fasta.close()
        d_filtered_alignments.clear()

        n_aligned_tags = float(len(set(l_aligned_tags)))
        n_total_alignments = float(n_total_alignments)
        n_discarded_tags = len(set(l_aligned_tags).difference(set(l_softclipped_tags)))

        # convert tags to reads for reporting. sum up the read counts in the read ids
        n_aligned_reads = 0.0
        for tagid in set(l_aligned_tags):
            n_aligned_reads += int(tagid.split('|')[1])
        n_softclipped_reads = 0.0
        for tagid in set(l_softclipped_tags):
            n_softclipped_reads += int(tagid.split('|')[1])
        n_singleton_reads = 0.0
        for tagid in set(l_singleton_tags):
            n_singleton_reads += int(tagid.split('|')[1])
        n_discarded_reads = 0.0
        for tagid in set(l_aligned_tags).difference(set(l_softclipped_tags)):
            n_discarded_reads += int(tagid.split('|')[1])

        fh_unmapped_fasta = open(unmapped_fasta, "r")
        n_unaligned_tags = 0.0
        n_unaligned_reads = 0.0
        for line in fh_unmapped_fasta:
            if not line.startswith('>'):
                continue
            else:
                n_unaligned_reads += int(line.split('|')[1])
                n_unaligned_tags += 1
        n_total_tags = float(n_aligned_tags + n_unaligned_tags)
        n_total_reads = float(n_aligned_reads + n_unaligned_reads)

        fh_log = open(clip_log, "w")
        fh_log.write("total input reads:            %d\n" % n_total_reads)
        fh_log.write("total unaligned reads:        %d (%.2f%%)\n" % (n_unaligned_reads,
                                                                      n_unaligned_reads / n_total_reads * 100))
        fh_log.write("total aligned reads:          %d (%.2f%%)\n" % (n_aligned_reads,
                                                                      n_aligned_reads/n_total_reads*100))
        fh_log.write("possible chimeric reads:      %d (%.2f%%)\n" % (n_softclipped_reads,
                                                                      n_softclipped_reads/n_total_reads*100))
        fh_log.write("Singleton reads:              %d (%.2f%%)\n" % (n_singleton_reads,
                                                                      n_singleton_reads/n_total_reads*100))
        fh_log.write("discarded reads:              %d (%.2f%%)\n" % (n_discarded_reads,
                                                                      n_discarded_reads/n_total_reads*100))
        fh_log.write("\n")
        fh_log.write("total input tags:             %d\n" % n_total_tags)
        fh_log.write("total unaligned tags:         %d (%.2f%%)\n" % (n_unaligned_tags,
                                                                      n_unaligned_tags / n_total_tags * 100))
        fh_log.write("total aligned tags:           %d (%.2f%%)\n" % (n_aligned_tags,
                                                                      n_aligned_tags/n_total_tags*100))
        fh_log.write("possible chimeric tags:       %d (%.2f%%)\n" % (len(set(l_softclipped_tags)),
                                                                      len(set(l_softclipped_tags))/n_total_tags*100))
        fh_log.write("Singleton tags:               %d (%.2f%%)\n" % (len(set(l_singleton_tags)),
                                                                      len(set(l_singleton_tags))/n_total_tags*100))
        fh_log.write("discarded tags:               %d (%.2f%%)\n" % (n_discarded_tags,
                                                                      n_discarded_tags/n_total_tags*100))
        fh_log.write("\n")
        fh_log.write("total alignments:             %d\n" % n_total_alignments)
        fh_log.write("softclipped alignments:       %d (%.2f%%)\n" % (n_softclipped_alignments,
                                                                      n_softclipped_alignments/n_total_alignments*100))
        fh_log.write("singleton alignments:         %d (%.2f%%)\n" % (n_singleton_alignments,
                                                                      n_singleton_alignments/n_total_alignments*100))
        fh_log.write("discarded alignments:         %d (%.2f%%)\n" % (n_discarded_alignments,
                                                                      n_discarded_alignments/n_total_alignments*100))
        fh_log.close()

        # TODO: 1G memory used!
        pysam.sort("-m", "1G", "-@", str(self.threads), filtered_bam, "-T", first_prefix, "-o", first_prefix + ".bam")
        pysam.index(first_prefix + ".bam")
        os.system("rm " + filtered_bam)

        logging.info("Done")
        return

    def get_chimeric_alignments(self, d_left_readpositions, d_right_readpositions):
        d_chimeric_alignments = {}
        # go through one of the dicts, because each tag should be present in both the dicts
        for tagid_left in d_left_readpositions.keys():
            tagid_right = tagid_left[:-1] + "r"
            if tagid_right not in d_right_readpositions:
                continue
            l_transcript_pairs = list(itertools.product(d_left_readpositions[tagid_left],
                                                        d_right_readpositions[tagid_right]))
            # TODO: hard cutoff?
            # TODO: Not all the multi mappings are considered!!
            if len(l_transcript_pairs) > 1000:
                continue
            for pair in l_transcript_pairs:
                first_readpos = pair[0]
                second_readpos = pair[1]
                [transcriptid1, start1, end1] = first_readpos.split(':')
                [transcriptid2, start2, end2] = second_readpos.split(':')
                # In case of mir-mir interactions,soft clipped parts are not mapped separately.
                # Hence, make sure that the both mapped parts of a read excluding the overlapping are
                # at least 18nt long. If it is on same transcript, the parts shold not overlap
                # MMMMMMMMMMSSSSSSS
                # SSSSSSSMMMMMMMMMM
                first_only_matches = int(end1) - int(start1) + 1 - utilities.overlap([int(start1), int(end1)],
                                                                                     [int(start2), int(end2)])
                second_only_matches = int(end2) - int(start2) + 1 - utilities.overlap([int(start1), int(end1)],
                                                                                      [int(start2), int(end2)])
                if first_only_matches >= 18 and second_only_matches >= 18 and\
                        not utilities.overlap_on_same_reference(first_readpos, second_readpos, self.min_looplen):
                    d_chimeric_alignments[tagid_left + "\t" + transcriptid1] = 1
                    d_chimeric_alignments[tagid_right + "\t" + transcriptid2] = 1
        return d_chimeric_alignments

    def separate_alignments(self, firstpart, secondpart):
        logging.info("Separating " + firstpart + "-" + secondpart + " alignments from others")
        first_prefix, second_prefix, pairs_prefix = utilities.file_prefixes(self.outputdir, firstpart, secondpart)
        d_left_readpositions, d_right_readpositions = utilities.get_read_pos(first_prefix + '.bam')

        print("Done extracting read positions:: "
              + str(len(d_left_readpositions))
              + "," + str(len(d_right_readpositions)))

        d_chimeric_alignments = self.get_chimeric_alignments(d_left_readpositions, d_right_readpositions)
        d_left_readpositions.clear()
        d_right_readpositions.clear()
        print("Done extracting chimeric reads:: " + str(len(d_chimeric_alignments)))

        fh_first_bam = pysam.Samfile(first_prefix + '.bam', "rb")
        fh_first_chimeric_bam = pysam.Samfile(first_prefix + '.' + firstpart + '.bam', "wb", template=fh_first_bam)
        fh_first_chimeric_bed = open(first_prefix + '.' + firstpart + '.bed', "w")
        for alignment in fh_first_bam.fetch(until_eof=True):
            alignment.qname = alignment.qname.split('|',1)[1]
            readid = alignment.qname
            referenceid = fh_first_bam.getrname(alignment.tid)
            alignmentrefstart = alignment.reference_start
            alignmentrefend = alignment.reference_end
            bedentry = '\t'.join([referenceid,
                                  str(alignmentrefstart),
                                  str(alignmentrefend),
                                  ','.join([readid,
                                            referenceid,
                                            d_txid_geneid[referenceid],
                                            str(alignmentrefstart),
                                            str(alignmentrefend)]),
                                  "0",
                                  "+"])
            if readid+"\t"+referenceid in d_chimeric_alignments:
                fh_first_chimeric_bam.write(alignment)
                fh_first_chimeric_bed.write(bedentry + "\n")
        fh_first_bam.close()
        fh_first_chimeric_bam.close()
        fh_first_chimeric_bed.close()
        d_chimeric_alignments.clear()
        print("running on " + first_prefix + '.' + firstpart + '.bed')
        utilities.transcript_to_genomic_pos(first_prefix + '.' + firstpart + '.bed',
                                            f_geneexonbed,
                                            f_txexonbed,
                                            p_intersectbed)

        logging.info("Done")
        return

    def merge_loci_and_quantify(self, firstpart, secondpart):
        logging.info("Merging to genomic loci and quantifiying of " + firstpart + "-" + secondpart + " alignments")
        first_prefix, second_prefix, pairs_prefix = utilities.file_prefixes(self.outputdir, firstpart, secondpart)

        if firstpart == "miRNA":
            if secondpart == "miRNA":
                first_prefix += ".miRNA"
            elif secondpart == "transcriptome":
                first_prefix += ".transcriptome"
        elif firstpart == "transcriptome":
            first_prefix += ".transcriptome"

        first_bed_merged = first_prefix + ".genomic.merged.bed"
        first_loci_groups_file = first_prefix + ".loci.groups.txt"
        first_counts = first_prefix + ".group.counts"

        second_bed_merged = second_prefix + ".genomic.merged.bed"
        second_loci_groups_file = second_prefix + ".loci.groups.txt"
        second_counts = second_prefix + ".group.counts"

        if firstpart == secondpart:
            first_bed = first_prefix + ".genomic.bed"
            second_bed = second_prefix + ".genomic.bed"
        else:
            first_bed = first_prefix + ".genomic.filtered.bed"
            second_bed = second_prefix + ".genomic.filtered.bed"

        # Merge the loci in the first part
        print("Merging first part loci...")
        quantifier.merge_bed(first_bed, first_bed_merged, self.merge_overlap)
        print("done")
        quantifier.create_mmgs(first_bed_merged, first_loci_groups_file, self.crg_share)
        d_first_group_counts, d_first_group_tpms = quantifier.quantify_mmgs(first_loci_groups_file, self.em_threshold)
        fh_first_counts = open(first_counts, "w")
        for locusmmg in sorted(d_first_group_counts.keys()):
            fh_first_counts.write(locusmmg + "\t"
                                  + str(d_first_group_counts[locusmmg]) + "\t"
                                  + str(d_first_group_tpms[locusmmg]) + "\n")
        fh_first_counts.close()
        # Only merge the loci in the second part if it is mir vs transcriptome
        if firstpart != secondpart:
            print("Merging second part loci...")
            quantifier.merge_bed(second_bed, second_bed_merged, self.merge_overlap)
            print("done")
            quantifier.create_mmgs(second_bed_merged, second_loci_groups_file, self.crg_share)
            d_second_group_expression, d_second_group_tpms = quantifier.quantify_mmgs(second_loci_groups_file, self.em_threshold)
            fh_second_counts = open(second_counts, "w")
            for locusmmg in sorted(d_second_group_expression.keys()):
                fh_second_counts.write(locusmmg + "\t"
                                       + str(d_second_group_expression[locusmmg])
                                       + "\t" + str(d_second_group_tpms[locusmmg]) + "\n")
            fh_second_counts.close()

        logging.info("Done")
        return

    def write_second_part(self, firstpart, secondpart):
        logging.info("Separating " + firstpart + "-" + secondpart + " alignments")
        first_prefix, second_prefix, pairs_prefix = utilities.file_prefixes(self.outputdir, firstpart, secondpart)

        second_bam = second_prefix + '.temp.bam'
        print(second_bam)
        fh_second_bam = pysam.Samfile(second_bam, "rb")
        filtered_bam = second_prefix + '.bam'
        fh_secondbed = open(second_prefix + '.bed', "w")

        # contains [readid"\t"referenceid]-->[score].
        # For a single read with multiple alignments on same reference, best alignment is chosen by its score.
        # Hence, there is only one alignment per a read-reference pair
        d_read_reference_score_inter = {}
        d_filtered_alignments = {}
        secondpart_readids = {}  # stores readIds of all qualified second parts, used to match with first part read ids

        print("secondbam: " + second_bam)
        # This block of code populates the above given hashes
        for alignment in fh_second_bam.fetch(until_eof=True):
            if alignment.is_unmapped:
                continue
            readid = alignment.qname
            referenceid = fh_second_bam.getrname(alignment.tid)
            alignscore = alignment.opt("AS")
            alignmentrefstart = alignment.reference_start
            alignedcigar = alignment.cigarstring

            if readid+"\t"+referenceid in d_read_reference_score_inter:
                if int(alignscore) <= int(d_read_reference_score_inter[readid+"\t"+referenceid]):
                    continue
            d_read_reference_score_inter[readid+"\t"+referenceid] = alignscore
            d_filtered_alignments[readid+"\t"+referenceid] = str(alignmentrefstart)+"\t"+alignedcigar
            if readid.split("|")[-1] == "r":
                readid = readid[:-1] + "l"
            elif readid.split("|")[-1] == "l":
                readid = readid[:-1] + "r"
            secondpart_readids[readid] = 1
        fh_second_bam.close()
        d_read_reference_score_inter.clear()

        print("filtered_bam: " + filtered_bam)
        fh_filtered_bam = pysam.Samfile(filtered_bam, "wb", template=fh_second_bam)
        fh_second_bam = pysam.Samfile(second_bam, "rb")
        counter = 1
        for alignment in fh_second_bam.fetch(until_eof=True):
            readid = alignment.qname
            referenceid = fh_second_bam.getrname(alignment.tid)
            alignmentrefstart = alignment.reference_start
            alignmentrefend = alignment.reference_end
            alignedcigar = alignment.cigarstring
            if readid+"\t"+referenceid in d_filtered_alignments:
                [filteredrefstart, filteredcigar] = d_filtered_alignments[readid+"\t"+referenceid].split('\t')
                if int(filteredrefstart) == int(alignmentrefstart) and filteredcigar == alignedcigar:
                    bedentry = '\t'.join([referenceid,
                                          str(alignmentrefstart),
                                          str(alignmentrefend),
                                          ','.join([readid,
                                                    referenceid,
                                                    d_txid_geneid[referenceid],
                                                    str(alignmentrefstart),
                                                    str(alignmentrefend)]),
                                          "0",
                                          "+"])
                    fh_secondbed.write(bedentry+"\n")
                    # alignment.qname = str(counter) + '|' + alignment.qname
                    fh_filtered_bam.write(alignment)
                    counter += 1
        fh_second_bam.close()
        fh_filtered_bam.close()
        fh_secondbed.close()
        d_filtered_alignments.clear()

        print("first_prefix.bam: " + first_prefix + '.bam')
        print("first_prefix.secondpart.bam: " + first_prefix + '.' + secondpart + '.bam')
        fh_first_bam = pysam.Samfile(first_prefix + '.bam', "rb")
        fh_first_bamout = pysam.Samfile(first_prefix + '.' + secondpart + '.bam', "wb", template=fh_first_bam)
        fh_first_bedout = open(first_prefix + '.' + secondpart + '.bed', "w")
        for alignment in fh_first_bam.fetch(until_eof=True):
            readid = alignment.qname
            if readid in secondpart_readids:
                referenceid = fh_first_bam.getrname(alignment.tid)
                alignmentrefstart = alignment.reference_start
                alignmentrefend = alignment.reference_end
                bedentry = '\t'.join([referenceid,
                                      str(alignmentrefstart),
                                      str(alignmentrefend),
                                      ','.join([readid,
                                                referenceid,
                                                d_txid_geneid[referenceid],
                                                str(alignmentrefstart),
                                                str(alignmentrefend)]),
                                      "0",
                                      "+"])
                fh_first_bamout.write(alignment)
                fh_first_bedout.write(bedentry+"\n")
        fh_first_bam.close()
        fh_first_bamout.close()
        fh_first_bedout.close()
        secondpart_readids.clear()

        print("running on " + first_prefix + '.' + secondpart + '.bed')
        utilities.transcript_to_genomic_pos(first_prefix + '.' + secondpart + '.bed', f_geneexonbed, f_txexonbed, p_intersectbed)
        print("running on " + second_prefix + '.bed')
        utilities.transcript_to_genomic_pos(second_prefix + '.bed', f_geneexonbed, f_txexonbed, p_intersectbed)

        logging.info("Done")
        return

    def filter_multimappings(self, firstpart, secondpart):
        """
            Filter multi-mappings based on the length of the alignments. Extract the longest combined alignment for
            each read. Then keep the alignments which are significantly longer compared to the combined longest alignment.
            Write them to filtered BED files.
        """
        logging.info("Filtering " + firstpart + "-" + secondpart + " multimappings")
        first_prefix, second_prefix, pairs_prefix = utilities.file_prefixes(self.outputdir, firstpart, secondpart)
        if firstpart == "miRNA":
            if secondpart == "miRNA":
                first_prefix += ".miRNA"
            elif secondpart == "transcriptome":
                first_prefix += ".transcriptome"
        elif firstpart == "transcriptome":
            first_prefix += ".transcriptome"

        fh_first_bed = open(first_prefix + ".genomic.bed")
        d_first_readlen = defaultdict(lambda: defaultdict(list))
        for line in fh_first_bed:
            f = line.split("\t")
            readid = f[3].split(",")[0]
            # consider tagid instead of readid because tagid represents actual intial read
            # for tagid consider "|l" and "|r" part separately
            tagid = "|".join(readid.split('|')[1:])
            # referenceid contains transcriptid,geneid,txstart,txend. Even if there are multiple such referenceids for
            # a single tagid, the align length should be the same
            referenceid = ",".join(f[3].split(",")[1:])
            align_len = int(f[2]) - int(f[1]) + 1
            d_first_readlen[tagid][referenceid] = align_len
        fh_first_bed.close()

        fh_second_bed = open(second_prefix + ".genomic.bed")
        d_second_readlen = defaultdict(lambda: defaultdict(list))
        for line in fh_second_bed:
            f = line.split("\t")
            readid = f[3].split(",")[0]
            tagid = "|".join(readid.split('|')[1:])
            referenceid = ",".join(f[3].split(",")[1:])
            align_len = int(f[2]) - int(f[1]) + 1
            d_second_readlen[tagid][referenceid] = align_len
        fh_second_bed.close()

        d_combined_long_alignments1 = defaultdict(list)
        d_combined_long_alignments2 = defaultdict(list)

        for tagid1 in d_first_readlen.keys():
            # consider reads that are present in both parts (hybrids)
            # if it is "l" in first it should be "r" in second and viceversa
            if tagid1.split("|")[-1] == "r":
                tagid2 = tagid1[:-1] + "l"
            elif tagid1.split("|")[-1] == "l":
                tagid2 = tagid1[:-1] + "r"
            if tagid2 not in d_second_readlen:
                continue
            longest1 = 0
            longest2 = 0
            longest_combined = 0
            # get the longest combined alignment for each read
            for referenceid in d_first_readlen[tagid1]:
                if d_first_readlen[tagid1][referenceid] > longest1:
                    longest1 = d_first_readlen[tagid1][referenceid]
            for referenceid in d_second_readlen[tagid2]:
                if d_second_readlen[tagid2][referenceid] > longest2:
                    longest2 = d_second_readlen[tagid2][referenceid]

            if longest1 + longest2 > longest_combined:
                longest_combined = longest1 + longest2
            # TODO make this cut-off as an option
            longest_combined = longest_combined * self.so_len_threshold

            # keep the alignments which are significantly longer compared to the combined longest alignment
            for referenceid1 in d_first_readlen[tagid1]:
                for referenceid2 in d_second_readlen[tagid2]:
                    if d_first_readlen[tagid1][referenceid1] + d_second_readlen[tagid2][referenceid2] >= longest_combined:
                        d_combined_long_alignments1[tagid1].append(referenceid1)
                        d_combined_long_alignments2[tagid2].append(referenceid2)

        n_discarded_alignments = 0
        fh_first_bed = open(first_prefix + ".genomic.bed")
        fh_first_bed_filtered = open(first_prefix + ".genomic.filtered.bed", "w")
        for line in fh_first_bed:
            f = line.split("\t")
            readid = f[3].split(",")[0]
            tagid = "|".join(readid.split('|')[1:])
            referenceid = ",".join(f[3].split(",")[1:])
            if tagid in d_combined_long_alignments1 and referenceid in d_combined_long_alignments1[tagid]:
                fh_first_bed_filtered.write(line)
            else:
                n_discarded_alignments += 1
        fh_first_bed.close()
        fh_first_bed_filtered.close()
        print("discarded alignments from first part:  " + str(n_discarded_alignments))

        fh_second_bed = open(second_prefix + ".genomic.bed")
        fh_second_bed_filtered = open(second_prefix + ".genomic.filtered.bed", "w")
        for line in fh_second_bed:
            f = line.split("\t")
            readid = f[3].split(",")[0]
            tagid = "|".join(readid.split('|')[1:])
            referenceid = ",".join(f[3].split(",")[1:])
            if tagid in d_combined_long_alignments2 and referenceid in d_combined_long_alignments2[tagid]:
                fh_second_bed_filtered.write(line)
            else:
                n_discarded_alignments += 1
        fh_second_bed.close()
        fh_second_bed_filtered.close()
        print("discarded alignments from second part: " + str(n_discarded_alignments))

        logging.info("Done")
        return

    def write_loci_pairs(self, firstpart, secondpart):
        logging.info("Writing " + firstpart + "-" + secondpart + " interactions")
        first_prefix, second_prefix, pairs_prefix = utilities.file_prefixes(self.outputdir, firstpart, secondpart)

        first_prefix += "." + secondpart
        if firstpart == secondpart:
            second_prefix = first_prefix

        first_bed = first_prefix + ".genomic.bed"
        first_counts_file = first_prefix + ".group.counts"
        first_groups_file = first_prefix + ".loci.groups.txt"
        second_bed = second_prefix + ".genomic.bed"
        second_counts_file = second_prefix + ".group.counts"
        second_groups_file = second_prefix + ".loci.groups.txt"
        locipairs = pairs_prefix + ".loci.pairs"

        d_first_group_expression, d_first_group_tpm, d_first_readgroups, d_first_grouploci, d_first_locipos, d_first_locusgroup_share, d_first_locus_transcripts = \
            utilities.parse_files(first_counts_file, first_groups_file)
        d_second_group_expression, d_second_group_tpm, d_second_readgroups, d_second_grouploci, d_second_locipos, d_second_locusgroup_share, d_second_locus_transcripts = \
            utilities.parse_files(second_counts_file, second_groups_file)
        d_first_read_genomic_pos, d_first_read_tx_pos = utilities.parse_bed(first_bed)
        d_second_read_genomic_pos, d_second_read_tx_pos = utilities.parse_bed(second_bed)

        first_uniq_tpms = sorted(list(set(d_first_group_tpm.values())))
        first_tpm_threshold = first_uniq_tpms[int(self.tpm_cutoff*len(first_uniq_tpms))]

        second_uniq_tpms = sorted(list(set(d_second_group_tpm.values())))
        second_tpm_threshold = second_uniq_tpms[int(self.tpm_cutoff*len(second_uniq_tpms))]

        d_first_readgroups_expression = defaultdict(lambda: defaultdict(float))
        d_second_readgroups_expression = defaultdict(lambda: defaultdict(float))

        print("groups less than " + str(first_tpm_threshold) + " and " + str(second_tpm_threshold) + " are discarded")

        # discard groups that don't meet tpm threshold
        # compute the normalized expression over all multimapped groups
        for readid in set(d_first_readgroups.keys()):
            first_groups_expression = 0.0
            for first_group in set(d_first_readgroups[readid]):
                if d_first_group_tpm[first_group] < first_tpm_threshold:
                    continue
                first_groups_expression += float(d_first_group_expression[first_group])
            for first_group in set(d_first_readgroups[readid]):
                if first_groups_expression != 0.0:
                    first_group_expression = float(d_first_group_expression[first_group])
                    first_group_expression_normalized = first_group_expression / first_groups_expression
                    d_first_readgroups_expression[readid][first_group] = first_group_expression_normalized
                else:
                    print(readid, first_group)

        for readid in set(d_second_readgroups.keys()):
            second_groups_expression = 0.0
            for second_group in set(d_second_readgroups[readid]):
                if d_second_group_tpm[second_group] < second_tpm_threshold:
                    continue
                second_groups_expression += float(d_second_group_expression[second_group])
            for second_group in set(d_second_readgroups[readid]):
                if second_groups_expression != 0.0:
                    second_group_expression = float(d_second_group_expression[second_group])
                    second_group_expression_normalized = second_group_expression / second_groups_expression
                    d_second_readgroups_expression[readid][second_group] = second_group_expression_normalized
                else:
                    print(readid, second_group)

        n_chimeras = 0
        fh_locipairs = open(locipairs, "w")
        d_duplicate_pairs = {}
        readid2 = ""
        for readid1 in d_first_readgroups_expression.keys():
            readid = readid1
            if readid1.split("|")[-1] == "r":
                readid2 = readid1[:-1] + "l"
            elif readid1.split("|")[-1] == "l":
                readid2 = readid1[:-1] + "r"
            if readid2 not in d_second_readgroups_expression.keys():
                continue
            if firstpart != secondpart:
                # this will overwrite the values of miRNAs from type to family
                # remove serial number infront of the tagid, this id will be used for writing final table
                readid = readid1.split('|',1)[1]

            d_geneid_genetype.update(d_mirna_family)

            for first_group in set(d_first_readgroups_expression[readid1]):
                first_group_expression_normalized = d_first_readgroups_expression[readid1][first_group]
                for second_group in set(d_second_readgroups_expression[readid2]):
                    second_group_expression_normalized = d_second_readgroups_expression[readid2][second_group]
                    for first_locus in set(d_first_grouploci[first_group]):
                        # TODO: distributed final score for each locus based on it's share?
                        first_locus_score = float(first_group_expression_normalized) * \
                                            float(d_first_locusgroup_share[first_locus][first_group])
                        if first_locus_score < self.score_cutoff:
                            continue
                        for second_locus in set(d_second_grouploci[second_group]):
                            second_locus_score = float(second_group_expression_normalized) * \
                                                 float(d_second_locusgroup_share[second_locus][second_group])
                            if second_locus_score < self.score_cutoff:
                                continue
                            combined_score_normalized = first_locus_score * second_locus_score
                            for first_transcriptid in set(d_first_locus_transcripts[readid1][first_locus]):
                                for second_transcriptid in set(d_second_locus_transcripts[readid2][second_locus]):
                                    if d_first_locipos[first_locus] == d_second_locipos[second_locus]:
                                        continue

                                    if ",".join([readid1,
                                                 first_locus,
                                                 second_locus,
                                                 first_transcriptid,
                                                 second_transcriptid]) in d_duplicate_pairs \
                                            or ",".join([readid2,
                                                         second_locus,
                                                         first_locus,
                                                         second_transcriptid,
                                                         first_transcriptid]) in d_duplicate_pairs:
                                        continue
                                    d_duplicate_pairs[",".join([readid1,
                                                                first_locus,
                                                                second_locus,
                                                                first_transcriptid,
                                                                second_transcriptid])] = 1
                                    n_chimeras += 1
                                    chimeraid = str(n_chimeras)
                                    interaction = "\t".join([chimeraid,
                                                             readid,
                                                             first_transcriptid,
                                                             second_transcriptid,
                                                             d_txid_geneid[first_transcriptid],
                                                             d_txid_geneid[second_transcriptid],
                                                             d_geneid_genename[d_txid_geneid[first_transcriptid]],
                                                             d_geneid_genename[d_txid_geneid[second_transcriptid]],
                                                             d_geneid_genetype[d_txid_geneid[first_transcriptid]],
                                                             d_geneid_genetype[d_txid_geneid[second_transcriptid]],
                                                             d_first_read_tx_pos[readid1][first_transcriptid],
                                                             d_txid_len[first_transcriptid],
                                                             d_second_read_tx_pos[readid2][second_transcriptid],
                                                             d_txid_len[second_transcriptid],
                                                             d_first_read_genomic_pos[readid1][first_transcriptid],
                                                             d_second_read_genomic_pos[readid2][second_transcriptid],
                                                             first_locus,
                                                             second_locus,
                                                             d_first_locipos[first_locus],
                                                             d_second_locipos[second_locus],
                                                             first_group,
                                                             second_group,
                                                             str(d_first_group_tpm[first_group]),
                                                             str(d_second_group_tpm[second_group]),
                                                             str(first_locus_score),
                                                             str(second_locus_score),
                                                             str(combined_score_normalized)])
                                    fh_locipairs.write(interaction + "\n")
        fh_locipairs.close()

        logging.info("Done")
        return

    def validate_groups(self, firstpart, secondpart):
        first_prefix, second_prefix, pairs_prefix = utilities.file_prefixes(self.outputdir, firstpart, secondpart)

        first_prefix += "." + secondpart
        if firstpart == secondpart:
            second_prefix = first_prefix

        first_groups_file = first_prefix + ".loci.groups.txt"
        second_groups_file = second_prefix + ".loci.groups.txt"

        fh_first_groups_file = open(first_groups_file)
        fh_second_groups_file = open(second_groups_file)

        d_first_group_transcripts = defaultdict(list)
        d_second_group_transcripts = defaultdict(list)

        for line in fh_first_groups_file:
            f = line.rstrip("\n").split("\t")
            d_first_group_transcripts[f[2]].extend(set(f[5].split(";")))

        for line in fh_second_groups_file:
            f = line.rstrip("\n").split("\t")
            d_second_group_transcripts[f[2]].extend(set(f[5].split(";")))

        for group in sorted(set(d_first_group_transcripts.keys())):
            l_all_goterms = []
            n_all_transcripts = 0.0
            n_goterm_transcripts = 0.0
            for transcriptid in set(d_first_group_transcripts[group]):
                n_all_transcripts += 1
                l_current_goterms = []
                # print(transcriptid+"\t"+d_txid_geneid[transcriptid])
                # print(d_txid_goterms[transcriptid + "\t" + d_txid_geneid[transcriptid]])
                l_current_goterms.extend(d_mirna_family[transcriptid])
                if len(set(list(set(l_all_goterms).intersection(l_current_goterms)))) > 0:
                    n_goterm_transcripts += 1
                l_all_goterms.extend(d_mirna_family[transcriptid])
            if n_all_transcripts > 2:
                print(group, n_goterm_transcripts , n_all_transcripts, n_goterm_transcripts / n_all_transcripts)


        # TODO: currently, putting all transcripts into a single list. But multiple transcripts may belong to the same
        # genomic locus.Check if at least one of the locus transcripts has a goterm in common instead of individual transcpx
        for group in sorted(set(d_second_group_transcripts.keys())):
            l_all_goterms = []
            n_all_transcripts = 0.0
            n_goterm_transcripts = 0.0
            for transcriptid in set(d_second_group_transcripts[group]):
                n_all_transcripts += 1
                l_current_goterms = []
                # print(transcriptid+"\t"+d_txid_geneid[transcriptid])
                # print(d_txid_goterms[transcriptid + "\t" + d_txid_geneid[transcriptid]])
                l_current_goterms.extend(d_txid_goterms[transcriptid + "\t" + d_txid_geneid[transcriptid]])
                if len(set(list(set(l_all_goterms).intersection(l_current_goterms)))) > 0:
                    n_goterm_transcripts += 1
                # print(set(list(set(l_all_goterms).intersection(l_current_goterms))))
                l_all_goterms.extend(d_txid_goterms[transcriptid + "\t" + d_txid_geneid[transcriptid]])
                # # if transcriptid+"\t"+d_txid_geneid[transcriptid] not in d_txid_geneid:
                # #     continue
                # print(d_txid_goterms[transcriptid + "\t" + d_txid_geneid[transcriptid]])
                # d_goterms[d_txid_goterms[transcriptid+"\t"+d_txid_geneid[transcriptid]]] += 1
            #print(len(d_goterms))
            # if n_all_transcripts > 2:
            #     print(group, n_goterm_transcripts , n_all_transcripts, n_goterm_transcripts / n_all_transcripts)

            # print(group)
            # l_all_gotersms = []
            # l_current_goterms = []
            # d_goterms = defaultdict(int)
            # for transcriptid in d_second_group_transcripts[group]:
            #     if transcriptid+"\t"+d_txid_geneid[transcriptid] not in d_txid_geneid:
            #         continue
            #     print(d_txid_goterms[transcriptid + "\t" + d_txid_geneid[transcriptid]])
            #     l_current_goterms.append(d_txid_goterms[transcriptid+"\t"+d_txid_geneid[transcriptid]])
            #     print(list(set(d_all_gotersms.keys()).intersection(l_current_goterms)))
            #     l_all_gotersms.append(d_txid_goterms[transcriptid+"\t"+d_txid_geneid[transcriptid]])
            #     print(l_all_gotersms)
            # # if len(d_goterms) > 2:
            # #     print(group, max(d_goterms.values())/float(sum(d_goterms.values())))

    def fasta_to_dict(fasta):
        seqid = None
        seq = ""
        fasta_dict = {}
        with open (fasta) as fh_fasta:
            for line in fh_fasta:
                if line.startswith(">"):
                    fasta_dict[seqid] = seq
                    seqid = line
                    seq = ""
                else:
                    seq += trim(line)
        return fasta_dict


    def hybridize_pairs(self, firstpart, secondpart):
        first_prefix, second_prefix, pairs_prefix = utilities.file_prefixes(self.outputdir, firstpart, secondpart)
        bed1 = pairs_prefix + "." + firstpart + ".bed"
        bed2 = pairs_prefix + "." + secondpart + ".bed"
        fasta1 = pairs_prefix + "." + firstpart + ".fasta"
        fasta2 = pairs_prefix + "." + secondpart + ".fasta"
        utilities.pairspos_to_bed(pairs_prefix, bed1, bed2)
        os.system(p_fastafrombed + " -s -name -fi " + f_genomicfasta + " -bed " + bed1 + " -fo " + fasta1)
        os.system(p_fastafrombed + " -s -name -fi " + f_genomicfasta + " -bed " + bed2 + " -fo " + fasta2)

        d_seq_records1 = fasta_to_dict(fasta1)
        d_seq_records2 = fasta_to_dict(fasta2)
        l_seqids = d_seq_records1.keys()
        chunk_size = len(l_seqids)/(self.threads-1)
        slices = utilities.chunks(l_seqids, chunk_size)
        print(len(l_seqids), str(chunk_size), len(slices))
        manager = multiprocessing.Manager()
        jobs = []
        for i, s in enumerate(slices):
            hybrid_file = pairs_prefix + ".hybrid."+str(i)
            print(len(l_seqids), len(s), hybrid_file)
            j = multiprocessing.Process(target=utilities.run_hybridmin,
                                        args=(d_seq_records1, d_seq_records2, s, hybrid_file, pairs_prefix, p_hybmin))
            jobs.append(j)
        for j in jobs:
            j.start()
        for j in jobs:
            j.join()
        os.system("cat " + pairs_prefix + ".hybrid.* > " + pairs_prefix + ".hybrids")
        os.system("rm " + pairs_prefix + ".hybrid.*")
        return

    def write_hybrids(self, firstpart, secondpart):
        first_prefix, second_prefix, pairs_prefix = utilities.file_prefixes(self.outputdir, firstpart, secondpart)
        locipairs = pairs_prefix + ".loci.pairs"
        hybrids = pairs_prefix + ".hybrids"
        interactome = pairs_prefix + ".interactome"

        d_hybrids = {}
        with open(hybrids) as fh_hybrids:
            for line in fh_hybrids:
                f = line.split("\t")
                d_hybrids[f[0]] = "\t".join(f[1:])

        fh_interactome = open(interactome, "w")
        with open(locipairs) as fh_locipairs:
            for line in fh_locipairs:
                fh_interactome.write(line.rstrip("\n") + "\t" + d_hybrids[line.split("\t")[0]])
        fh_interactome.close()

    def parse_arguments(self, argv=None):
        def score_float(x):
            x = float(x)
            if x < 0.0 or x > 1.0:
                raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
            return x

        parser = argparse.ArgumentParser(description='Chimeric read annotator for interactome data',
                                         usage='%(prog)s [-h] [-v,--version]',
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        parser.add_argument('-i', '--in', action='store', dest='insample', required=True,
                            metavar='', help='Input sample path')

        parser.add_argument('-o', '--out', action='store', dest='outputdir', required=True, metavar='',
                            help='Output directory path for the whole analysis')

        parser.add_argument('-p', '--threads', action='store', type=int, default=1, metavar='',
                            dest='threads',
                            help='Number of threads to use')

        parser.add_argument("-s", '--stranded', type=str, choices=["fw", "rc", "both"], default='fw', metavar='',
                            dest='stranded',
                            help="Strand-specificity of input samples. "
                                 "fw = map to transcript strand"
                                 "rc = map to reverse compliment of transcript strand"
                                 "both = try to map on both strnads")

        parser.add_argument('-m', '--min_looplen', action='store', type=int, default=3, metavar='',
                            dest='min_looplen',
                            help='Set minimum loop length in case of single RNA duplex')

        parser.add_argument('-tc', '--tpm_cutoff', action='store', type=score_float, default=0.1, metavar='',
                            dest='tpm_cutoff',
                            help='Transcripts with less than this percentile TPMs will be discarded in '
                                 'the final output. [0-1.0]')

        parser.add_argument('-sc', '--score_cutoff', action='store', type=score_float, default=0.0, metavar='',
                            dest='score_cutoff',
                            help='Hybrids with less than this score will be discarded in the final output. [0-1.0]')

        parser.add_argument('-mo', '--merge_overlap', action='store', type=score_float, default=0.7, metavar='',
                            dest='merge_overlap',
                            help='Minimum percentage overlap among BED entries inorder to merge. [0-1.0]')

        parser.add_argument('-cs', '--crg_share', action='store', type=score_float, default=0.7, metavar='',
                            dest='crg_share',
                            help='Minimum percentage overlap of a locus on a CRG inorder to merge it '
                                 'into a CRG. [0-1.0]')

        parser.add_argument('-e', '--em_threshold', action='store', type=float, default=0.01, metavar='',
                            dest='em_threshold',
                            help='The maximum difference of transcripts expression between two consecutive iterations '
                                 'of EM algorithm to converge.')

        parser.add_argument('-so', '--suboptimal', action='store', type=float, default=0.8, metavar='',
                            dest='so_len_threshold',
                            help='In case of multimappings, hybrids with total alignment length greater than this % of'
                                 'the longest hybrid are considered.')

        parser.add_argument("-f", '--hybridize', action='store_true', help="Hybridize the predicted chimeras")

        parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')

        args = parser.parse_args()
        self.insample = args.insample
        self.outputdir = args.outputdir
        self.threads = args.threads
        self.stranded = args.stranded
        self.min_looplen = args.min_looplen
        self.tpm_cutoff = args.tpm_cutoff  # we need percentile later
        self.score_cutoff = args.score_cutoff
        self.merge_overlap = args.merge_overlap
        self.crg_share = args.crg_share
        self.em_threshold = args.em_threshold
        self.so_len_threshold = args.so_len_threshold
        self.hybridize = args.hybridize

        print('Sample                           : ' + args.insample)
        print('Output directory                 : ' + args.outputdir)
        print('Number of threads                : ' + str(args.threads))
        print('Stranded                         : ' + args.stranded)
        print('Minimum loop length              : ' + str(args.min_looplen))
        print('TPM cutoff                       : ' + str(args.tpm_cutoff))
        print('Score cutoff                     : ' + str(args.score_cutoff))
        print('Merge overlap                    : ' + str(args.merge_overlap))
        print('CRG share                        : ' + str(args.crg_share))
        print('EM threshold                     : ' + str(args.em_threshold))
        print('Suboptimal hybrid threshold      : ' + str(args.so_len_threshold))


if __name__ == "__main__":

    chira = ChiRA()
    chira.parse_arguments()

    print("===================================================================")
    if not os.path.exists(os.path.join(chira.outputdir, 'miRNA_vs_transcriptome')):
        os.makedirs(os.path.join(chira.outputdir, 'miRNA_vs_transcriptome'))
    if not os.path.exists(os.path.join(chira.outputdir, 'transcriptome_vs_transcriptome')):
        os.makedirs(os.path.join(chira.outputdir, 'transcriptome_vs_transcriptome'))

    logging.basicConfig(level=logging.INFO,
                        filename=os.path.join(chira.outputdir, 'run.log'),  # log to this file
                        filemode='w',
                        format='%(asctime)s %(message)s')

    chira.map_reads(None, "miRNA")
    chira.write_softclipped("miRNA")
    chira.separate_alignments("miRNA", "miRNA")
    chira.merge_loci_and_quantify("miRNA", "miRNA")
    chira.write_loci_pairs("miRNA", "miRNA")

    chira.map_reads("miRNA", "transcriptome")
    chira.write_second_part("miRNA", "transcriptome")
    chira.filter_multimappings("miRNA", "transcriptome")
    chira.merge_loci_and_quantify("miRNA", "transcriptome")
    chira.write_loci_pairs("miRNA", "transcriptome")

    chira.map_reads(None, "transcriptome")
    chira.write_softclipped("transcriptome")
    chira.separate_alignments("transcriptome", "transcriptome")
    chira.merge_loci_and_quantify("transcriptome", "transcriptome")
    chira.write_loci_pairs("transcriptome", "transcriptome")

    if chira.hybridize:
        chira.hybridize_pairs("miRNA", "miRNA")
        chira.write_hybrids("miRNA", "miRNA")

        chira.hybridize_pairs("miRNA", "transcriptome")
        chira.write_hybrids("miRNA", "transcriptome")

        chira.hybridize_pairs("transcriptome", "transcriptome")
        chira.write_hybrids("transcriptome", "transcriptome")

    #
    # # multi_sample.reproducible_interactions_at_loci_level(chira.outputdir, replicates, "miRNA", "miRNA")
    # #     reproducible_interactions(outputdir, replicates, "miRNA", "transcriptome")
    # #    reproducible_interactions(outputdir, replicates, "transcriptome", "transcriptome")
    #
    #     # /usr/local/tools/p_bowtie2/2.2.6/iuc/package_bowtie_2_2_6/6d6cca69a34a/p_bowtie2 --local --ignore-quals --norc -k 100 -N 1 --rdg 10,10 --rfg 10,10 --score-min L,36,0 -p 3 -x /data/5/galaxy_import/galaxy_user_data/pospisilik/miconomou/eclash/db/miRNA/miRNA -f -U /scratch/bi02/videmp/ProjectPospisilik/fastq/clash5/trimmed/Bcells2.trimmed.uniq.fasta -S vs_miRNA.sam --un vs_miRNA.unmapped.fasta; cat vs_miRNA.sam | grep -o -P 'NM:i:\d+' | sort | uniq -c | sort -n
    #
    # # . /usr/local/tools/pysam/0.8.3/iuc/package_python_2_7_pysam_0_8_3/7defbb2032a1/env.sh
    # # . /usr/local/tools/biopython/1.66/env.sh
    # # . /usr/local/tools/bedtools/2.24/env.sh
    # # . /usr/local/tools/sailfish/0.7.6/iuc/package_sailfish_0_7_6/04996d8b8587/env.sh
