#!/usr/bin/env python3

"""
Table
=====

Generates a CSV table from a given GenBank file and a BioSyc entry file.

Usage:
    $ ./table.py sequence.gb pathway-genes -o table.csv
"""

import argparse
import csv
import sys
from Bio import SeqIO
from biocyc import read_entries

BIOSYQ = 'data/pathway-genes'
SEQUENCE = 'data/sequence.gb'

def get_cds_features(record, entries):
    accessions = [e.gene_accession for e in entries]
    return filter(lambda f: f.type == 'CDS' and \
                            f.qualifiers['locus_tag'][0] in accessions,
                  record.features)

def write_table(csvfile, features, entries):
    table = csv.writer(csvfile, delimiter=',', quotechar='"',
                       quoting=csv.QUOTE_MINIMAL)
    table.writerow([
        'Gene ID',
        'Locus tag',
        'Gene name',
        'Strand',
        'Protein ID',
        'Protein name',
        'Aminoacid number',
        'Location',
        'EC',
        'Function'
    ])

    for f in features:
        entry = list(filter(lambda e: e.gene_accession == f.qualifiers['locus_tag'][0],
                       entries))[0]
        table.writerow([
            entry.gene_id,
            f.qualifiers['locus_tag'][0],
            entry.gene_name,
            f.location.strand,
            f.qualifiers['protein_id'][0],
            f.qualifiers['product'][0],
            len(f.qualifiers['translation'][0]),
            f.location,
            entry.reaction_ec,
            f.qualifiers['function'][0] if 'function' in f.qualifiers else ''
        ])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("sequence", type=argparse.FileType('r'),
                        default=SEQUENCE, nargs='?',
                        help="file containing the sequence in Genbank format")
    parser.add_argument("biosyc", type=argparse.FileType('r'),
                        default=BIOSYQ, nargs='?',
                        help="file downloaded from BioSyq")
    parser.add_argument("-o", "--outfile", type=argparse.FileType('w'),
                        default=sys.stdout,
                        help="output file to write the table")
    args = parser.parse_args()
    record = SeqIO.read(args.sequence, 'genbank')
    entries = read_entries(args.biosyc)
    features = get_cds_features(record, entries)
    write_table(args.outfile, features, entries)

if __name__ == "__main__":
    main()
