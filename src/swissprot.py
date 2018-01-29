#!/usr/bin/env python3

"""
Swissprot
=========
"""

import argparse
import csv
import sys
from Bio import ExPASy
from Bio import SwissProt

def write_table(csvfile, proteins):
    table = csv.writer(csvfile, delimiter=',', quotechar='"',
                       quoting=csv.QUOTE_ALL)
    table.writerow([
        'Proteína',
        'Nome',
        'Classe',
        'Tipo de molecula',
        'Tamanho da sequencia',
        'Codigo de acesso',
        'Criado',
        'Descrição',
        'Nome do gene',
        'Organismo',
        'Classificacao do organismo',
        'ID da taxonomia',
        'Comentários',
        'Palavras-chave',
        'Características',
        'Sequência',
        'Informação da sequência',
        'Referências'
    ])

    for protein in proteins:
        handle = ExPASy.get_sprot_raw(protein.strip())

        try:
            record = SwissProt.read(handle)
            table.writerow([
                protein.strip(),
                record.entry_name,
                record.data_class,
                record.molecule_type,
                record.sequence_length,
                record.accessions[0],
                record.created[0],
                record.description,
                record.gene_name,
                record.organism,
                ','.join(record.organism_classification),
                record.taxonomy_id[0],
                ','.join(record.comments),
                ','.join(record.keywords),
                record.features,
                record.sequence,
                record.seqinfo,
                ','.join([str(r.references) for r in record.references]),
            ])
        except Exception as e:
            print(e)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("proteins", type=argparse.FileType('r'),
                        default=sys.stdin, nargs='?',
                        help="file containing the list of proteins")
    parser.add_argument("-o", "--outfile", type=argparse.FileType('w'),
                        default=sys.stdout,
                        help="output file to write the proteins info")
    args = parser.parse_args()

    write_table(args.outfile, args.proteins)

if __name__ == "__main__":
    main()
