#!/usr/bin/env python3

"""
Genome
======

Fecthes the genome with a given ID from NCBI either in FASTA or Genbank format.

Examples
--------

We can use this utility to quickly fecth a genome from NCBI:

    $ ./genome.py NZ_CP013742 --email 'your@email.com'
    $ ./genome.py NZ_CP013742 -f fasta -o ../data/sequence.fasta
"""

import argparse
import sys
from Bio import Entrez
from Bio import SeqIO

# This is the default id to fetch
ID='NZ_CP013742'

def fetch(record_id, db, outformat, outfile, email):
    Entrez.email = email
    handle = Entrez.efetch(db=db, id=record_id, rettype=outformat,
                           retmode='text')
    outfile.write(handle.read())
    handle.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("id", default=ID, nargs='?', help="Record id")
    parser.add_argument("-o", "--outfile", default=sys.stdout,
                        type=argparse.FileType('w'),
                        help="File to output the results")
    parser.add_argument("-d", "--db", choices=['nucleotide', 'protein'],
                        default='nucleotide',
                        help="database to fetch the results")
    parser.add_argument("-f", "--format", choices=['fasta', 'gb'],
                        default='gb', help="Format of the output file")
    parser.add_argument("--email", default='',
                        help="Email address to use with NCBI")
    args = parser.parse_args()
    fetch(args.id, args.db, args.format, args.outfile, args.email)

if __name__ == "__main__":
    main()
