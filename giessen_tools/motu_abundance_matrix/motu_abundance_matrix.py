#!/usr/bin/env python3
"""
program to create an abundance matrix for each mOTU over each sample.
Needs the cluster file from the -uc setting from vsearch as input and
a tabular blast file in the output format "-outfmt '6 std staxids sscinames'.
"""


__author__ = "sgriep"
__projekt__ = "motu_abundance_matrix"
__date__ = "2023-01-11"
__version__ = "1.0"


import sys
import argparse
from pathlib import Path

def main():
    """
    this is where the whole script takes place. parsing arguments,
    reading in clusters, writing clusternumbers into blast output
    and merge cluster counts with best hit species per sample
    """

    args = parse_args()

    (clusters, read2cluster, matrix_header) = parse_clusters(args.input_cluster)

    clusters = parse_and_rewrite_blast(args.input_blast, args.output_blast, clusters, read2cluster)

    write_otu_matrix(clusters, matrix_header, args.output_otutable)

    return 0

def parse_args():
    """
    parse arguments passed through CLI call
    """

    parser = argparse.ArgumentParser()

    parser.add_argument("-ic", "--input_cluster",
                        help="cluster file in VSEARCH cluster format from -uc setting")
    parser.add_argument("-ib", "--input_blast",
                        help="tab separated blast file with the following entrys:\n"
                        "-outfmt '6 std staxids sscinames'")
    parser.add_argument("-ot", "--output_otutable",
                        help="path to store the calculated mOTU table")
    parser.add_argument("-ob", "--output_blast",
                        help="path to store the modified blast table file")

    return parser.parse_args()

def get_read_name(seq_id):
    """
    extract read name from sequence id with information separated by ';'
    """

    seq_id = seq_id.split(';')[0]
    has_centroid = seq_id.split('=')

    if len(has_centroid) > 1:
        return has_centroid[1]

    return seq_id

def get_sample_name(seq_id):
    """
    extract sample name from sequence id
    """

    seq_id = get_read_name(seq_id)
    sample_name = seq_id.split('___')[0]

    return sample_name

def get_read_size(read_id):
    """
    extract read abundance from sequence id
    """

    read_parts = read_id.split(';')[1:]

    for part in read_parts:
        (key, val) = part.split('=')

        if key == 'size':
            return int(val)

    return 1

def parse_clusters(input_cluster):
    """
    parses vsearch clusters into to dict
    """

    clusters = {}
    read2cluster = {}
    header_samples = set()
    fh_input_cluster = open(str(Path(input_cluster).absolute()))

    for line in fh_input_cluster:
        cols = line.strip().split("\t")

        if cols[0] == "S":
            cluster_id = cols[1]
            cluster_read = cols[8]
            cluster_name = get_read_name(cluster_read)
            member_sample_name = get_sample_name(cluster_read)
            member_read_size = get_read_size(cluster_read)

            if cluster_name not in clusters:
                header_samples.add(member_sample_name)
                clusters[cluster_name] = {
                    'id': cluster_id,
                    'name': cluster_name,
                    'matrix': {
                        member_sample_name: member_read_size,
                    },
                    'total': member_read_size,
                }
                read2cluster[cluster_name] = cluster_name

        elif cols[0] == "H":
            cluster_id = cols[1]
            member_read = cols[8]
            cluster_read = cols[9]
            member_name = get_read_name(member_read)
            cluster_name = get_read_name(cluster_read)
            member_sample_name = get_sample_name(member_read)
            member_read_size = get_read_size(member_read)

            if cluster_name in clusters:

                if member_sample_name in clusters[cluster_name]['matrix']:
                    clusters[cluster_name]['matrix'][member_sample_name] += member_read_size

                else:
                    header_samples.add(member_sample_name)
                    clusters[cluster_name]['matrix'][member_sample_name] = member_read_size

                clusters[cluster_name]['total'] += member_read_size
                read2cluster[member_name] = cluster_name

            else:
                header_samples.add(member_sample_name)
                clusters[cluster_name] = {
                    'id': cluster_id,
                    'name': cluster_name,
                    'matrix': {
                        member_sample_name: member_read_size,
                    },
                    'total': member_read_size,
                }
                read2cluster[cluster_name] = cluster_name

                read2cluster[member_name] = cluster_name

        elif cols[0] == "C":
            cluster_id = cols[1]
            cluster_size = cols[2]
            cluster_read = cols[8]
            cluster_name = get_read_name(cluster_read)

            if cluster_name in clusters:

                if cluster_id == clusters[cluster_name]['id'] and cluster_size == clusters[cluster_name]['total']:
                    continue

            else:
                sys.stderr.write(f'cluster {cluster_id} missed through parsing\n')

        else:
            continue

    fh_input_cluster.close()

    header = sorted(header_samples, key=sort_by_sample_name)

    return (clusters, read2cluster, header)

def handle_s_lines():
    pass

def handle_h_lines():
    pass

def parse_and_rewrite_blast(input_blast, output_blast, clusters, read2cluster):
    """
    parses the blast output into to dict containing only the best hit(s)
    and writes out blast table with substituting the query name with cluster number
    """

    taxas = {}
    input_blast_path = Path(input_blast)
    output_blast_path = Path(output_blast)
    fh_input_blast = open(str(input_blast_path.absolute()))
    fh_output_blast = open(str(output_blast_path.absolute()), 'w')
    
    header_blast_output = [
        'qseqid (OUT) (Clusternumber)',
        'sseqid (Genbank)',
        'pident (identity)',
        'length',
        'mismatch',
        'gapopen',
        'qstart',
        'qend',
        'sstart',
        'send',
        'evalue',
        'bitscore',
        'staxids',
        'ssciname',
    ]
    fh_output_blast.write('{header}\n'.format(header='\t'.join(header_blast_output)))

    for line in fh_input_blast:
        cols = line.strip().split('\t')

        # ignore header lines
        if cols[0].startswith("#"):
            continue

        # strips size info from read name
        read_name = get_read_name(cols[0])
        cluster_name = read2cluster[read_name]
        pident = cols[2]
        bitscore = cols[11].strip()
        taxa = cols[13].strip()

        output_cols = cols[1:]
        fh_output_blast.write('{clusternumber}\t'.format(clusternumber=clusters[cluster_name]['id']))
        fh_output_blast.write('{cols}\n'.format(cols='\t'.join(output_cols)))

        # skip if current taxa already stored
        ##if read_name in taxas:
        ##    continue

        ##else:
        ##    taxas[read_name] = taxa

        # current cluster got already a blast hit
        if 'best_hit' in clusters[cluster_name]:

            # the current bitscore equals the previously stored best hit
            if float(bitscore) == float(clusters[cluster_name]['best_hit']['bitscore']):

                # the current taxa is not stored in sscinames of the cluster yet
                if taxa not in clusters[cluster_name]['best_hit']['sscinames']:
                    clusters[cluster_name]['best_hit']['sscinames'].add(taxa)

            # the current bitscore is greater than the previously stored best hit
            elif float(bitscore) > float(clusters[cluster_name]['best_hit']['bitscore']):
                # replace the best hit for this cluster
                clusters[cluster_name]['best_hit'] = {
                    'sscinames': set(),
                    'pident': pident,
                    'bitscore': bitscore,
                }
                clusters[cluster_name]['best_hit']['sscinames'].add(taxa)

        # a blast hit for the current cluster appears for the very first time
        else:
            # add blast hit as new best hit for this cluster
            clusters[cluster_name]['best_hit'] = {
                'sscinames': set(),
                'pident': pident,
                'bitscore': bitscore,
            }
            clusters[cluster_name]['best_hit']['sscinames'].add(taxa)

    fh_input_blast.close()
    fh_output_blast.close()

    return clusters

def sort_by_cluster_size(item):
    """
    sort by cluster size
    """

    return int(item[1]['total'])

def sort_by_sample_name(item):
    """
    sort by sample name
    """

    #sample = int(item.split('_')[0].replace('S', ''))
    sample = item.split('_')[0]
    taxa = item.split('_')[-1]

    return (sample, taxa)

def write_otu_matrix(clusters, header, output_otutable):
    """
    write otu matrix into file
    """

    prefix_header = ['#OTU_ID', 'ssciname', 'pident', 'bitscore']
    header_line = []
    header_line.extend(prefix_header)
    header_line.extend(header)
    fh_output_otu = open(output_otutable, 'w')

    fh_output_otu.write('{header}\n'.format(header='\t'.join(header_line)))

    for name, cluster in sorted(clusters.items(), key=sort_by_cluster_size, reverse=True):

        if 'best_hit' not in cluster:
            continue

        fh_output_otu.write('{id}'.format(id=cluster['id']))
        fh_output_otu.write('\t{ssciname}'.format(ssciname=','.join(sorted(cluster['best_hit']['sscinames']))))
        fh_output_otu.write('\t{pident}'.format(pident=cluster['best_hit']['pident']))
        fh_output_otu.write('\t{bitscore}'.format(bitscore=cluster['best_hit']['bitscore']))

        for matrix_header in header:

            if matrix_header not in cluster['matrix']:
                fh_output_otu.write('\t0')
                continue

            fh_output_otu.write('\t{count}'.format(count=cluster['matrix'][matrix_header]))

        fh_output_otu.write('\n')

    fh_output_otu.close()

if __name__ == "__main__":
    sys.exit(main())
