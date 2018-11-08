import argparse
import os
import re

def gtfToDict(gtf, itype):
    GTF = open(gtf, "r")
    geneStartPosMap = {}
    geneEndPosMap = {}
    annotationMap = {}
    for line in  GTF:
        f = line.rstrip('\n').split('\t')
        if f[2] != "exon": continue
        pos = '\t'.join([f[0],f[3],f[4]])
        strand = f[6]
        descDict = dict(item.replace("\";","").split(" \"") for item in filter(None,f[8].rstrip().split("\"; ")))
        if not "gene_id" in descDict:
            sys.stderr.write("attribute \"gene_id\" not found for position :" + pos)
            continue
        if descDict["gene_id"] in geneStartPosMap:
            if int(f[3]) < int(geneStartPosMap[descDict["gene_id"]]):
                geneStartPosMap[descDict["gene_id"]] = int(f[3])
            if int(f[4]) > int(geneEndPosMap[descDict["gene_id"]]):
                geneEndPosMap[descDict["gene_id"]] = int(f[4])
        else:
            geneStartPosMap[descDict["gene_id"]] = int(f[3])
            geneEndPosMap[descDict["gene_id"]] = int(f[4])
        if "gene_biotype" not in descDict:
            descDict["gene_biotype"] = "NA"
        if "gene_name" not in descDict:
            descDict["gene_name"] = "NA"
        if itype == "deseq":
            desc = descDict["gene_biotype"]+"\t"+descDict["gene_name"]+"\t"+f[0]+"\t"+str(geneStartPosMap[descDict["gene_id"]])+"\t"+str(geneEndPosMap[descDict["gene_id"]])+"\t"+strand
        elif itype == "dexseq":
            desc = descDict["gene_biotype"]+"\t"+descDict["gene_name"]
        annotationMap[descDict["gene_id"]] = desc
    GTF.close()
    return annotationMap

def main():
    parser = argparse.ArgumentParser(description='Maps the DEXSeq couting bins back to the transcript and exon ids')
    parser.add_argument('-i','--input', help='DESeq2/DEXSeq output',required=True)
    parser.add_argument('-t','--itype', help='Type intput',required=True)
    parser.add_argument('-g','--gtf',help='Original annotation gtf file used', required=True)
    parser.add_argument('-o','--output', help='Output table',required=True)
    args = parser.parse_args()
    print ("DE(X)Seq output file: %s" % args.input)
    print ("Type of output file: %s" % args.input)
    print ("Annotation file: %s" % args.gtf)
    print ("Annotated output file: %s" % args.output)

    annotationMap = gtfToDict(args.gtf, args.itype)

    IN = open(args.input, "r")
    OUT = open(args.output, "w")
    for line in IN:
        if args.itype == "deseq":
            geneId = line.split('\t')[0]
            OUT.write(line.rstrip('\n')+'\t'+annotationMap[geneId]+'\n')
        elif args.itype == "dexseq":
            geneIds = line.split('\t')[1].split('+')
            geneNames = []
            geneTypes = []
            for geneId in geneIds:
                [geneName, geneType] = annotationMap[geneId].split('\t')
                geneNames.append(geneName)
                geneTypes.append(geneType)
            OUT.write(line.rstrip('\n')+'\t'+'+'.join(geneTypes)+'\t'+'+'.join(geneNames)+'\n')
    IN.close()
    OUT.close()

if __name__ == "__main__":
    main()

