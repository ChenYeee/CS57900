import numpy as np
from Bio import pairwise2
from Bio import SeqIO

def calc_align(ori_seq1, ori_seq2, matrix):

    gap_open = [-6, -10]
    gap_ext = [-5, -4]

    for i in gap_open:
        for j in gap_ext:
            print "Gap open:", i, ", Gap extension:", j
            alignments = pairwise2.align.globalds(ori_seq1, ori_seq2, matrix, i, j)
            align_str_array = np.array([list(rec) for rec in alignments], np.character)
            seq1 = align_str_array[0][0].tostring()
            seq2 = align_str_array[0][1].tostring()

            score = align_str_array[0][2]

            match = 0
            gap = 0
            transition = 0
            transversion = 0

            for k in range(len(seq1)):
                if seq1[k] == '-':
                    gap += 1
                if seq2[k] == '-':
                    gap += 1

                if seq1[k] == seq2[k]:
                    match += 1
                elif (seq1[k] == 'A' and seq2[k] == 'G') or (seq1[k] == 'G' and seq2[k] == 'A') or (seq1[k] == 'C' and seq2[k] == 'T') or (seq1[k] == 'T' and seq2[k] == 'C'):
                    transition += 1
                elif seq1[k] != '-' and seq2[k] != '-':
                    transversion += 1
                else:
                    pass

            print "Score: ", score
            print "#match: ", match
            print "#gap: ", gap
            print "#transition: ", transition
            print "#transversion: ", transversion
            print "\n"

def main():
    seq1918 = SeqIO.read("SC1918.txt", "fasta")
    seq2007 = SeqIO.read("Brisbane2007.txt", "fasta")
    seq2009 = SeqIO.read("California2009.txt", "fasta")

    matrix1 = {('A', 'A'): 5, ('A', 'T'): -4, ('A', 'G'): -4, ('A', 'C'): -4,
    ('T', 'A'): -4, ('T', 'T'): 5, ('T', 'G'): -4, ('T', 'C'): -4,
    ('G', 'A'): -4, ('G', 'T'): -4, ('G', 'G'): 5, ('G', 'C'): -4,
    ('C', 'A'): -4, ('C', 'T'): -4, ('C', 'G'): -4, ('C', 'C'): 5}

    matrix2 = {('A', 'A'): 5, ('A', 'T'): -5, ('A', 'G'): -1, ('A', 'C'): -5,
    ('T', 'A'): -5, ('T', 'T'): 5, ('T', 'G'): -5, ('T', 'C'): -1,
    ('G', 'A'): -1, ('G', 'T'): -5, ('G', 'G'): 5, ('G', 'C'): -5,
    ('C', 'A'): -5, ('C', 'T'): -1, ('C', 'G'): -5, ('C', 'C'): 5}

    matrix  = [matrix1, matrix2]

    for i in range(2):
        print "Under Matrix ", i + 1, ": \n"
        print "Sequences: SC1918 & Brisbane2007"
        calc_align(seq1918.seq, seq2007.seq, matrix[i])
        print "\n\n\nSequences: SC1918 & California2009"
        calc_align(seq1918.seq, seq2009.seq, matrix[i])
        print "\n\n\nSequences: Brisbane2007 & California2009"
        calc_align(seq2007.seq, seq2009.seq, matrix[i])



if __name__ == "__main__":
    main()
