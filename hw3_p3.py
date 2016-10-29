#!usr/bin/python
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio import SeqIO
import csv

VALUE_THRESH = 80

def blast(record, mode):
    if mode == 1:
        result_handle = NCBIWWW.qblast("blastn", "nr", record.format("fasta"), entrez_query='all [filter] NOT(environmental samples[organism] OR metagenomes[orgn]) AND bacteria[organism] NOT Gram-positive bacteria[organism] NOT uncultured bacterium[organism] NOT predicted[title]', hitlist_size=250, megablast=True)
        blast_record = NCBIXML.read(result_handle)
        parseResult(blast_record, 'ouput_hs.csv', record)
    else:
        result_handle = NCBIWWW.qblast("blastn", "nr", record.format("fasta"), entrez_query='all [filter] NOT(environmental samples[organism] OR metagenomes[orgn]) AND bacteria[organism] NOT Gram-positive bacteria[organism] NOT uncultured bacterium[organism] NOT predicted[title]', hitlist_size=250, megablast=False, nucl_penalty=-1, nucl_reward=1)
        blast_record = NCBIXML.read(result_handle)
        parseResult(blast_record, 'ouput_ss.csv', record)

def parseResult(blast_record, filename, record):
    specieses = {}
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            # calculate the %coverage
            coverage = float((hsp.align_length - hsp.gaps)/(len(record)*0.01))
            # calculate the %identity
            identity = float(hsp.identities/(len(hsp.match)*0.01))
            if identity >= VALUE_THRESH and coverage >= VALUE_THRESH:
                species = alignment.title.split()[1]+' '+alignment.title.split()[2]
                specieses[species] = [species, alignment.title.split('|')[4], coverage, identity]
    write2file(filename, specieses)

def write2file(filename, specieses):
            with open(filename, 'w') as csvfile:
                fieldnames = ['Species', 'Gene description', 'Query Coverage', 'Percent Identity']
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                for i, j in specieses.items():
                    writer.writerow({'Species':j[0], 'Gene description':j[1], 'Query Coverage':j[2], 'Percent Identity':j[3]})

def main():
    record = SeqIO.read("ermb.fasta", format="fasta")
    blast(record, 1);
    blast(record, 0);

if __name__ == "__main__":
    main()
