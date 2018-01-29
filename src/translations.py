#!/usr/bin/env python3

"""
Translations
============

Creates a FASTA output, containing all the selected genes translations.

Usage:
    $ ./translations sequence.gb pathway-genes  -o translations.fasta
"""

import argparse
import sys
from Bio import SeqIO
from biocyc import read_entries
from table import get_cds_features

BIOSYQ = 'data/pathway-genes'
SEQUENCE = 'data/sequence.gb'

def get_translations(features):
    translations = []

    for feature in features:
        translations.append((feature.qualifiers['locus_tag'][0],
                             feature.qualifiers['translation'][0]))

    return translations


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("sequence", type=argparse.FileType('r'),
                        default=SEQUENCE, nargs='?',
                        help="file containing the sequence in Genbank format")
    parser.add_argument("biosyc", default=BIOSYQ, nargs='?',
                        type=argparse.FileType('r'),
                        help="file downloaded from BioSyq")
    parser.add_argument("-o", "--outfile", type=argparse.FileType('w'),
                        default=sys.stdout,
                        help="output file to write the fasta sequences")
    args = parser.parse_args()

    record = SeqIO.read(args.sequence, 'genbank')
    entries = read_entries(args.biosyc)

    translations = get_translations(get_cds_features(record, entries))

    for locus_tag, translation in translations:
        args.outfile.write("> {}\n{}\n".format(locus_tag, translation))

if __name__ == "__main__":
    main()
