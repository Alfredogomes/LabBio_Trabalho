#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio import SeqIO
from Bio.ExPASy import ScanProsite
import os

def scanFiles():
    
    results = ''
    # percorre os ficheiros fasta na directoria Prot_Biocyc e, para cada um,
    # retira a sequencia e faz uma pesquisa no prosite sobre esta.
    for file in os.listdir("Prot_Biocyc"):
        #obtem a sequencia presente no ficheiro
        record = SeqIO.read( "Prot_Biocyc/"+file, "fasta")
        seq = record.seq
        #prepara as keywords usadas no pedido de scan
        keywords = {'CC':'/SKIP-FLAG=FALSE;'}
        #efectua o pedido de pesquisa
        response = ScanProsite.scan(seq, 'http://www.expasy.org', 'xml', **keywords );
        #obtem-se a resposta
        obj = ScanProsite.read(response)
        #guarda-se os resultados numa string para adicionar ao ficheiro de resultados
        results += file+"\n"+str(obj)+"\n\n"
    
    #abre um ficheiro e grava os resultados obtidos na pesquisa
    fileRecord = open("results.txt", "w")
    fileRecord.write(results)
    fileRecord.close()
     
if __name__ == "__main__":
    scanFiles()
