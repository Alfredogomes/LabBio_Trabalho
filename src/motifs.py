#!/usr/bin/env python3

"""
Motifs
======
"""

import argparse
import sys
from Bio import SeqIO
from Bio.ExPASy import get_sprot_raw
from Bio.ExPASy import ScanProsite

SEQUENCES = 'data/translations.fasta'

def get_ids(fasta):
    ids = []

    for record in SeqIO.parse(fasta, 'fasta'):
        print(record.id)
        ids.append(record.id)

    return ids


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("sequences", type=argparse.FileType('r'),
                        default=SEQUENCES, nargs='?',
                        help="Fasta file containing the protein sequences")
    parser.add_argument("-o", "--outfile", type=argparse.FileType('w'),
                        default=sys.stdout,
                        help="output file to write the motifs")
    args = parser.parse_args()
    
    for protein in get_ids(args.sequences):
        with get_sprot_raw(protein) as handle:
            record = ScanProsite.read(handle)
            print(record)

if __name__ == "__main__":
    main()
