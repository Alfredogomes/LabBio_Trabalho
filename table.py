"""
Generates a table with the description of the gene/protein annotations.
It contains the following columns:

 - Gene identification (geneID NCBI, accession number NCBI, locus tag, gene
   name, strand)
 - Protein identification (Uniprot ID, revision, accession number NCBI,
   protein name)
 - Protein properties (aminoacid number, celular location)
 - Associated term list of GeneOntology
 - Associated EC numbers, or TC numbers
 - Description (text about protein's role)
 - Comments

"""

#!/usr/bin/env python3

