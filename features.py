#!/usr/bin/env python3

"""
Reads the features from a sequence, and exports them to plain text or HTML.
"""

from Bio import SeqIO

GENOME = 'genome.gb'
SEQUENCE = 'sequence.txt'

record = SeqIO.read(GENOME, 'genbank')

for feature in record.features:
    print(feature)
