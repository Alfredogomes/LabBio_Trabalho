#!/usr/bin/env python3

"""
Uniprot
=======

For a given protein file, fetches from Uniprot the accession, keywords and go
of each protein.

Usage:
    $ ./uniprot proteins -o uniprot.csv
"""

import argparse
import csv
import sys
import re
import requests
from Bio import ExPASy
from Bio import SeqIO
from Bio import SwissProt

PROTEINS = 'data/proteins'

def get_info(proteins):
    info = []

    for protein in proteins:
        record = requests.get("http://www.uniprot.org/uniprot/" +
                              protein.rstrip() + ".xml")
        try:
            ac = re.findall(r'<accession>(.*?)</accession>', record.text)[0]
            kw = ','.join(re.findall(r'<keyword id=".*">(.*)</.*>',
                                     record.text))
            go = ','.join(re.findall(r'<property type="term" value="(.*)"/>',
                                     record.text))
            info.append((ac, kw, go))
        except Exception as e:
            print(e)

    return info

def write_table(csvfile, info):
    table = csv.writer(csvfile, delimiter=',', quotechar='"',
                       quoting=csv.QUOTE_MINIMAL)
    table.writerow(['Accession', 'Keywords', 'Go'])

    for (ac, kw, go) in info:
        table.writerow([ac, kw, go])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("proteins", type=argparse.FileType('r'),
                        default=PROTEINS, nargs='?',
                        help="File containing a list of proteins")
    parser.add_argument("-o", "--outfile", type=argparse.FileType('w'),
                        default=sys.stdout,
                        help="output file to write the table")
    args = parser.parse_args()
    info = get_info(args.proteins)
    write_table(args.outfile, info)

if __name__ == "__main__":
    main()
