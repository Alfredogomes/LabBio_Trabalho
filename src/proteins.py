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

    for r in records:
        alignments = sorted(r.alignments, key=lambda a: a.hsps[0].expect)[0:10]

        for a in alignments:
            title = re.search(r'^(.*?)\|(\d+)\|(.*?)\|(.*?)\|.*',
                              a.title).group(4)
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
        args.outfile.write(protein + "\n")

if __name__ == "__main__":
    main()

