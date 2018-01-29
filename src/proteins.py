#!/usr/bin/env python3

"""
Proteins
========

Extracts all the protein codes from a BLAST result XML file.

Usage:
    $ ./proteins resultsBlast.xml -o proteins.txt

"""

import argparse
import re
import sys
from Bio.Blast import NCBIXML

RESULT = 'data/resultsBlast.xml'

def get_proteins(blast):
    records = NCBIXML.parse(blast)
    proteins = []

    for record in records:
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < 0.01:
                    title = re.search(r'^(.*?)\|(\d+)\|(.*?)\|(.*?)\|.*', alignment.title).group(4)
                    proteins.append(title)
    return proteins

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("result", type=argparse.FileType('r'),
                        default=RESULT, nargs='?',
                        help="XML blast result file")
    parser.add_argument("-o", "--outfile", type=argparse.FileType('w'),
                        default=sys.stdout,
                        help="output file to write the protein list")
    args = parser.parse_args()

    for protein in get_proteins(args.result):
        print(protein)

if __name__ == "__main__":
    main()

