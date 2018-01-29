# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 14:51:25 2018

@author: Joao
"""

# -*- coding: utf-8 -*-

#Procura os artigos em PubMed relacionados  com o organismo 
#Legionella pneumophila Philadelphia-1.
#Retorna a quantidade de artigos
import sys
from Bio import Entrez


Entrez.email = "biolg07@gmail.com" 
handle = Entrez.egquery(term="Legionella pneumophila Philadelphia-1")
record = Entrez.read(handle)
for row in record["eGQueryResult"]:
   if row["DbName"]=="pubmed":
       var = row["Count"]

print(var)

#IDs dos artigos da PubMed (IdList)
handle = Entrez.esearch(db="pubmed", term="Legionella pneumophila Philadelphia-1", retmax=var)
record = Entrez.read(handle)
handle.close()
idlist = record["IdList"]
print(idlist)


#Download dos registos da Medline
from Bio import Medline
handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="text")
registos = Medline.parse(handle)


valorOriginal = sys.stdout
f = open('artigos.txt', 'w')
sys.stdout = f

registos = list(registos)
for registo in registos:
    print("Titulo:", registo.get("TI", "?"))
    print("Autores:", registo.get("AU", "?"))
    print("Fontes:", registo.get("SO", "?"))
    print("")

sys.stdout = valorOriginal
f.close()

