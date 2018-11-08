import argparse
import os
import re

def gtfToBED(gtf, bed):
    GTF = open(gtf, "r")
    BED = open(bed, "w")
    for line in  GTF:
        f = line.rstrip('\n').split('\t')
        if f[0].startswith('#'): continue
        if f[2] != "exon": continue
        pos = '\t'.join([f[0],f[3],f[4]])
        strand = f[6]
        descDict = dict(item.replace("\";","").split(" \"") for item in filter(None,f[8].rstrip().split("\"; ")))
        if not "gene_id" in descDict:
            sys.stderr.write("attribute \"gene_id\" not found for position :" + pos)
            continue
        elif not "transcript_id" in descDict:
            sys.stderr.write("attribute \"transcript_id\" not found for gene :" + descDict["gene_id"])
            continue
        elif not "exon_number" in descDict:
            sys.stderr.write("attribute \"exon_number\" not found for gene :" + descDict["gene_id"] + " and transcript:" + descDict["transcript_id"])
            continue
        desc = descDict["gene_id"]+":"+descDict["transcript_id"]+":"+descDict["exon_number"]
        BED.write("\t".join([pos, desc, "0", strand])+"\n")
    GTF.close()
    BED.close()

def main():
    parser = argparse.ArgumentParser(description='Maps the DEXSeq couting bins back to the transcript and exon ids')
    parser.add_argument('-i','--input', help='DEXSeq output',required=True)
    parser.add_argument('-g','--gtf',help='Original annotation gtf file used for DEXSeq prepare tool', required=True)
    parser.add_argument('-o','--output', help='Output table',required=True)
    args = parser.parse_args()
    print ("DEXSeq output file: %s" % args.input)
    print ("Annotation file: %s" % args.gtf)
    print ("Mapped output file: %s" % args.output)

    os.system("awk '{print $12\"\t\"$13\"\t\"$14\"\t\"$1\"\t0\t\"$16}' " + args.input + " | sort -k1,1 -k2n,2 > output.bed")
    gtfToBED(args.gtf, "annotation.bed")
    os.system("intersectBed -wo -s -a output.bed -b annotation.bed > overlap.txt")

    binExonMap = {}
    OVERLAP = open("overlap.txt", "r")
    for line in OVERLAP:
        binId = line.split('\t')[3]
        exonId = line.split('\t')[9]
        binExonMap.setdefault(binId, []).append(exonId)
    OVERLAP.close()

    IN = open(args.input, "r")
    OUT = open(args.output, "w")
    for line in IN:
        binId = line.split('\t')[0]
        OUT.write(line.rstrip('\n')+'\t'+','.join(binExonMap[binId])+'\n')
    IN.close()
    OUT.close()

if __name__ == "__main__":
    main()

