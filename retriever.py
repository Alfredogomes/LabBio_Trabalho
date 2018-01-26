#!/usr/bin/env python3

from Bio import SeqIO
import csv

# Parte 1:
# =======

# Ficheiro obtido em:
# <https://biocyc.org/GCF_001941585/pathway-genes?object=PWY-841>
# Este ficheiro foi convertido para CSV para facilitar a sua leitura
GENES = 'pathway-genes.csv'
ACCESSIONS = []

with open(GENES, 'r') as csvfile:
    gene_reader = csv.DictReader(csvfile)
    for row in gene_reader:
        ACCESSIONS.append(row['Gene Accession'])


# Parte 2: Leitura das features dos genes selecionados
# ====================================================

# Sequência GenBank obtida a partir de:
# <https://www.ncbi.nlm.nih.gov/nuccore/NZ_CP013742>
record = SeqIO.read('sequence.gb', 'genbank')
features = []

for feature in record.features:
    if feature.type == 'CDS' and feature.qualifiers['locus_tag'][0] in ACCESSIONS:
        features.append(feature)

# Parte 3: Escrita dos resultados numa tabela CSV
# ===============================================

with open('table.csv', 'w') as csvfile:
    table = csv.writer(csvfile, delimiter=',', quotechar='"',
                       quoting=csv.QUOTE_MINIMAL)
    table.writerow([
        'ID do gene',
        'Locus tag',
        'Nome do gene',
        'Strand',
        'ID da Proteina',
        'Nome da proteína',
        'Nº aminoácidos',
        'Localização',
        'EC',
        'Função'
    ])

    for feature in features:
        table.writerow([
            'GeneID',
            feature.qualifiers['locus_tag'][0],
            'Nome do gene',
            feature.location.strand,
            feature.qualifiers['protein_id'][0],
            feature.qualifiers['product'][0],
            len(feature.qualifiers['translation'][0]),
            feature.location,
            'EC',
            feature.qualifiers['function'][0] if 'function' in feature.qualifiers else ''
        ])
