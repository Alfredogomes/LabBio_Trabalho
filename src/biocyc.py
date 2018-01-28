#!/usr/bin/env python3

"""
BioSyc
======

Parses a BioSyc entry file, and allows to apply filters to its result.
The output is a set of entries in JSON format.

Examples
--------

Suppose we have a file downloaded from Biosyc, containing all gene entries
named 'pathway-genes', like this one:
<https://biocyc.org/GCF_001941585/pathway-genes?object=PWY-841>

We can use this utility to parse the file and filter the results:

 - show all the genes:
    $ ./biocyc.py pathway-genes

 - find all the genes with a given ID:
    $ ./biocyc.py --gene G1FLP-2393 pathway-genes

 - you can combine the filters, to get more control over the result:
    $ ./biocyc.py --reaction ATPSYN-RXN  --reaction-ec 3.6.3.14 pathway-genes
"""

import argparse
import csv

# This is the default file that this script will parse, if no one is specified
# Obtained from:
# https://biocyc.org/GCF_001941585/pathway-genes?object=PWY-841
FILE = 'pathway-genes'

class Entry:
    def __init__(self, gene_id, gene_accession, gene_name, reaction_id,
                 reaction_ec, enzymatic_activity, evidence):
        self.gene_id = gene_id
        self.gene_accession = gene_accession
        self.gene_name = gene_name
        self.reaction_id = reaction_id
        self.reaction_ec = reaction_ec
        self.enzymatic_activity = enzymatic_activity
        self.evidence = evidence

    def __str__(self):
        return '{\n' + \
               '  "Gene ID": "{}",\n'.format(self.gene_id) + \
               '  "Gene Accession": "{}",\n'.format(self.gene_accession) + \
               '  "Gene name": "{}",\n'.format(self.gene_name) + \
               '  "Reaction id": "{}",\n'.format(self.reaction_id) + \
               '  "Reaction EC": "{}",\n'.format(self.reaction_ec) + \
               '  "Enzymatic activity": "{}",\n'.format(self.enzymatic_activity) + \
               '  "Evidence": "{}",\n'.format(self.evidence) + \
               '}'

def read_entries(csvfile):
    entries = []

    fields = ['Gene ID', 'Gene Accession', 'Gene name', 'Reaction id',
              'Reaction EC', 'Enzymatic activity', 'Evidence']
    reader = csv.DictReader(csvfile, fieldnames=fields, delimiter='\t')

    # Skip first two lines
    next(reader)
    next(reader)

    for row in reader:
        entry = Entry(gene_id=row['Gene ID'],
                      gene_accession=row['Gene Accession'],
                      gene_name=row['Gene name'],
                      reaction_id=row['Reaction id'],
                      reaction_ec=row['Reaction EC'],
                      enzymatic_activity=row['Enzymatic activity'],
                      evidence=row['Evidence'])
        entries.append(entry)

    return entries

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", type=argparse.FileType('r'), default=FILE, nargs='?',
                        help="biocyc genes file that contains the gene entries")
    parser.add_argument("-g", "--gene", help="filter by gene id")
    parser.add_argument("-a", "--accession", help="filter by gene accession")
    parser.add_argument("-n", "--name", help="filter by gene name")
    parser.add_argument("-r", "--reaction", help="filter by reaction id")
    parser.add_argument("--reaction-ec", help="filter by reaction ec")
    parser.add_argument("--enzymatic-activity", help="filter by enzymatic activity")
    parser.add_argument("-e", "--evidence", help="filter by evidence")
    args = parser.parse_args()
    entries = read_entries(args.file)

    if args.gene is not None:
        entries = filter(lambda e: e.gene_id == args.gene, entries)
    if args.accession is not None:
        entries = filter(lambda e: e.gene_accession == args.accession, entries)
    if args.name is not None:
        entries = filter(lambda e: e.gene_name == args.name, entries)
    if args.reaction is not None:
        entries = filter(lambda e: e.reaction_id == args.reaction, entries)
    if args.reaction_ec is not None:
        entries = filter(lambda e: e.reaction_ec == args.reaction_ec, entries)
    if args.enzymatic_activity is not None:
        entries = filter(lambda e: e.enzymatic_activity == args.enzymatic_activity, entries)

    for entry in entries:
        print(entry)

if __name__ == "__main__":
    main()
