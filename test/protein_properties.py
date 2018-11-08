#!/usr/bin/env python

import sys
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

sys.stdout.write("ID\tMW\tIP\tgravy\tlength\tinstability\tmonoisotpoic\tSequence\n")

def calc_property( prop ):
    try:
        return prop()
    except:
        return 'Error: undefined AA'

for record in SeqIO.parse(sys.stdin, "fasta"):
    a = ProteinAnalysis(str(record.seq))

    properties = list()
    properties.append(record.id)
    properties.append( calc_property(a.molecular_weight) )
    properties.append( calc_property(a.isoelectric_point) )
    properties.append( calc_property(a.gravy) )
    properties.append(a.length)
    properties.append( calc_property(a.instability_index))
    properties.append( calc_property(a.aromaticity))
    # always last column to make the output more readable
    properties.append(a.sequence)
    sys.stdout.write( '\t'.join(map(str, properties))+"\n" )

