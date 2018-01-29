#!/usr/bin/env python3

"""
Domains
=======

Tries to find more info about the genes whose functions are unknown.
"""

import argparse
import csv
import sys
from Bio import Entrez

FUNCTIONS = 'data/functions.csv'

def without_function(csvfile):
    """ Returns the list of genes without any known function """
    genes = []
    reader = csv.reader(csvfile, delimiter=',')

    for row in reader:
        if ''.join(row[2:]) == '':
            genes.append(row[0])

    return genes

def fetch(gene_id, email):
    Entrez.email = email
    handle = Entrez.efetch(db='cdd', id=gene_id, rettype='gb', retmode='text')
    print(handle.read())
    handle.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("functions", type=argparse.FileType('r'),
                        default=FUNCTIONS, nargs='?',
                        help="csv file containing the gene functions")
    parser.add_argument("-o", "--outfile", type=argparse.FileType('w'),
                        default=sys.stdout,
                        help="output file to write the new gene functions")
    parser.add_argument("--email", default='',
                        help="Email address to use with NCBI")
    args = parser.parse_args()

    for gene in without_function(args.functions):
        fetch(gene, args.email)

if __name__ == "__main__":
    main()
