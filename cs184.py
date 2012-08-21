#-------------------------------------------------------------------------------
# Name:        CS184
# Purpose:
#
# Author:      Michael Vu
#
# Created:     28/09/2011
# Copyright:   (c) Vum2 2011
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python
# xhx@ics.uci.edu

import sys
from string import maketrans

def readSeqsFromFile(filepath):
    f = open(filepath, 'r')
    sqs = [];
    for line in f:
        if line[:1] == '>':
            sqs.append('')
        elif(len(sqs) == 0):
            continue
        else:
            sqs[-1] += line[:len(line)-1]
    f.close()
    #print sqs
    return sqs


# combines every possible character from DNA molecules to create a codon table
def buildCodonTable():
    codons = []
    bases = 'T', 'C', 'A', 'G'
    amino_acids = 'FFLLSSSSYY!!CC!WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    for a in bases:
        for b in bases:
            for c in bases:
                codons.append(a+b+c)
    codon_table = dict(zip(codons,amino_acids))
    return codon_table

codon_table = buildCodonTable()

def codonTable():
    return codon_table

# assumes input length is multiple of 3
def translateSeq(sq):
    codon_table = codonTable()
    s = ''
    # traverse the sequence 3 characters at a time
    for x in range(0,len(sq)/3):
        s += codon_table[sq[x*3:x*3+3]]
    #print s
    return s

def findFrames(sq):
    start_codon = 'ATG'
    stop_codons = 'TAA', 'TAG', 'TGA'
    genes = []
    last = ''
    for x in range(0, len(sq)):
        # traverse 3 characters at a time
        if (sq[x:x+3] == start_codon):
            for y in range (1, (len(sq)-x)/3):
                if (sq[x+y*3:x+y*3+3] in stop_codons):
                    # make sure it's at least 200 amino acids long
                    if (y*3 >= 600):
                        genes.append((x, y*3))
                    else:
                        last = (x, y*3)
                    break
            else: # end of sequence
                if (len(sq)-x >= 600):
                    genes.append((x, len(sq)-x))
    #print genes
    return genes

def reverseComplement(sq):
    into = "TCAG"
    outto = "AGTC"
    trans = maketrans(into, outto)
    return sq.translate(trans)[::-1]

def readFindTrans(filepath):
    sequence_count = 0
    longest_gene_across = 0
    shortest_gene_across = 0
    # read file
    for idx, sq in enumerate(readSeqsFromFile(filepath)):
        sequence_count += 1
        gene_count = 0
        longest_gene = 0
        shortest_gene = 0
        print "SEQ #" + repr(idx+1)
        # forward string
        print "=== FORWARD STRING ==="
        for j, gene in enumerate(findFrames(sq)):
            gene_count += 1
            trans = translateSeq(sq[gene[0]:gene[0]+gene[1]])
            l = len(trans)
            if ( l > longest_gene):
                longest_gene = l
            if (l < shortest_gene or shortest_gene == 0):
                shortest_gene = l
            print 'Location: ' + repr(gene[0]) + '\nGene # ' + repr(j +1) + ': ' + trans
        # reverse string
        print "=== REVERSE COMPLEMENT ==="
        reverse = reverseComplement(sq)
        #print reverse
        for j, gene in enumerate(findFrames(reverse)):
            gene_count += 1
            trans = translateSeq(reverse[gene[0]:gene[0]+gene[1]])
            l = len(trans)
            if ( l > longest_gene):
                longest_gene = l
            if (l < shortest_gene):
                shortest_gene = l
            print 'Location: ' + repr(gene[0]) + '\nGene # ' + repr(j + 1) + ': ' + trans
        # Begin sequence statistics
        print '=== SEQ #' + repr(idx+1) + ' STATISTICS ==='
        print 'Number of genes: ' + repr(gene_count)
        print 'Length of shortest gene (in amino acids): ' + repr(shortest_gene)
        print 'Length of longest gene (in amino acids): ' + repr(longest_gene)
        if (longest_gene > longest_gene_across):
            longest_gene_across = longest_gene
        if (shortest_gene < shortest_gene_across or shortest_gene_across == 0):
            shortest_gene_across = shortest_gene
    #Begin overall statistics
    print '=== OVERALL STATISTICS ==='
    print 'Number of sequences: ' + repr(sequence_count)
    print 'Length of shortest gene across all sequences (in amino acids): ' + repr(shortest_gene_across)
    print 'Length of longest gene across all sequences (in amino acids): ' + repr(longest_gene_across)



def main():
    if (len(sys.argv) > 1):
        readFindTrans(sys.argv[1])
    else:
        readFindTrans(raw_input("Specify input filepath:"))
    pass

if __name__ == '__main__':
    main()