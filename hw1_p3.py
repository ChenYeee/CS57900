#!/usr/bin/python

with open("dna.txt", "r") as infile:
    rnaArray = []
    for line in infile:
        rna = line.replace("T", "U")
        rnaArray.append(rna)
infile.close()

codemap = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"s", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"-", "UAG":"-",
    "UGU":"C", "UGC":"C", "UGA":"-", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}

f = open('output.txt', 'w')
proteinArray = []

for rna in rnaArray:
    i = 0
    protein = ""
    while i + 3 <= len(rna):
        coden = rna[i: i+3]
        protein += codemap[coden]
        i += 3
    print >> f, protein
    proteinArray.append(protein)

i = 0
for aa in proteinArray[0]:
    if aa != proteinArray[1][i]:
        print >> f, i, aa, proteinArray[1][i]
    i += 1
f.close()
