#!/usr/bin/env python3

"""
Fetches the genome from NCBI nucleotide database, and writes it to Genbank
format. It also extracts the desired sequence and exports to TXT.
"""

from Bio import Entrez
from Bio import SeqIO

EMAIL = 'biolg7@gmail.com'
ACCESSION_ID = 'NC_002942.5'
START = '3386197'
STOP = '3392758'
GENBANK = 'genome.gb'
SEQUENCE = 'sequence.txt'

# Fetch the genome from NCBI nucleotide database
Entrez.email = EMAIL
handle = Entrez.efetch(db='nucleotide', rettype='gb', retmode='text',
              seq_start=START, seq_stop=STOP, id=ACCESSION_ID)
genome = SeqIO.read(handle, 'gb')

# Write the genome to Genbank format
SeqIO.write(genome, GENBANK, 'gb')

# Write the sequence to a TXT file
with open(SEQUENCE, 'w') as f:
    f.write(str(genome.seq))

handle.close()
