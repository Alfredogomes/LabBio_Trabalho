#!/usr/bin/env python3

"""
Motifs
======
"""

import argparse
import sys
from Bio import SeqIO
from Bio.ExPASy import ScanProsite

SEQUENCES = 'data/translations.fasta'

def get_motifs(records, outfile):
    for record in SeqIO.parse(records, 'fasta'):
        seq = record.seq
        keywords = {'CC':'/SKIP-FLAG=FALSE;'}
        response = ScanProsite.scan(seq, 'http://www.expasy.org', 'xml', **keywords)
        obj = ScanProsite.read(response)
        outfile.write(str(obj) + "\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("sequences", type=argparse.FileType('r'),
                        default=SEQUENCES, nargs='?',
                        help="Fasta file containing the protein sequences")
    parser.add_argument("-o", "--outfile", type=argparse.FileType('w'),
                        default=sys.stdout,
                        help="output file to write the motifs")
    args = parser.parse_args()

    get_motifs(args.sequences, args.outfile)
   
if __name__ == "__main__":
    main()
