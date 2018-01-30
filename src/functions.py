#!/usr/bin/env python3

"""
Functions
=========

Extracts the functions of a given set of genes.
Exports the results as a CVS matrix where the rows correspond to the Gene
and each column to each function. If that gene contains that function then the
matrix value will be 1 or else otherwise.

Usage:
    $ ./functions.py sequence.gb pathway-genes -o functions.csv
"""

import argparse
import csv
import re
import sys
from Bio import SeqIO
from biocyc import read_entries
from table import get_cds_features

BIOSYQ = 'data/pathway-genes'
SEQUENCE = 'data/sequence.gb'

def get_functions(features, entries):
    functions = []

    for f in features:
        gene_id = list(filter(lambda e:
                              e.gene_accession == f.qualifiers['locus_tag'][0],
                              entries))[0].gene_id
        fs = []

        if 'function' in f.qualifiers:
            fs = re.split(r'\s*[,/]\s*', f.qualifiers['function'][0])

        functions.append((gene_id, f.qualifiers['locus_tag'][0], fs))

    return functions


def write_table(csvfile, functions):
    header = ["GeneID", "Accesion"] + list(set([f for _, _, fs in functions for f in fs]))

    table = csv.writer(csvfile, delimiter=',', quotechar='"',
                       quoting=csv.QUOTE_MINIMAL)
    table.writerow(header)

    for gene_id, accession, fs in functions:
        row = [gene_id, accession] + ['X' if h in fs else '' for h in header[2:]]
        table.writerow(row)

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
    functions = get_functions(features, entries)

    write_table(args.outfile, functions)

if __name__ == "__main__":
    main()
