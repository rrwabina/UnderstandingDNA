import numpy as np
import json
import Bio
from Bio.Seq import Seq
from Bio import SeqIO
import fastaparser
import reverse_translator as biotools

'''
        The FASTA file format is a standard text-based format for representing nucleotide and 
        aminoacid sequences (usual file extensions include: .fasta, .fna, .ffn, .faa and .frn). 
        FastaParser is able to parse such files and extract the biological sequences within 
        into Python objects.
'''
filename = '/root/NLU/Project/data/test_cases.fa'
def parser_fasta(filename):
        name, seq = None, []
        for line in filename:
            line = line.rstrip()
            if line.startswith(">"):
                if name: yield (name, ''.join(seq))
                name, seq = line, []
            else:
                seq.append(line)
        if name: yield (name, ''.join(seq))

def read_fasta(filename):
    sequence = []
    with open(filename) as fp:
        for name, seq in parser_fasta(fp):
            sequence.append(seq)
    return sequence[1]

def prior_nucleotide(sequence, molecule = 'dna'):
    '''
    Measures the priors of oligonucleotide counts of DNA/RNA sequence
    '''
    total_seq = len(sequence)
    total_nucleotides = {}
    nucleotides = ['adenine', 'cytosine', 'guanine', 'thymine']
    if molecule == 'dna':
        bases = ['A', 'C', 'G', 'T']
        for base, nuclei in zip(bases, nucleotides):
            total_nucleotides[nuclei] = sequence.count(base)
    elif molecule == 'rna':
        bases = ['A', 'C', 'G', 'U']
        for base, nuclei in zip(bases, nucleotides):
            total_nucleotides[nuclei] = sequence.count(base)
    else:
        return 'Wrong input!'
    prior = {}
    for nuclei in nucleotides:
        prior[nuclei] = np.round(total_nucleotides[nuclei]/total_seq, 4)
    return prior

def prior_codon(sequence, ngram = 3, dictionary = True):
    sequence = Seq(sequence).transcribe()
    sequence = str(sequence)
    codon = [sequence[i: i + ngram] for i in range(0, len(sequence), ngram)]
    total_codon = len(codon)
    prior = {}
    codon_list = ['UUU', 'UUC', 'UUA', 'UUG', 'UCU', 'UAU', 'UGU',
                  'CUU', 'CUC', 'CUA', 'CUG', 'UCC', 'UAC', 'UGC',
                  'AUU', 'AUC', 'AUA', 'AUG', 'UCA', 'UAA', 'UGA',
                  'GUU', 'GUC', 'GUA', 'GUG', 'UCG', 'UAG', 'UGG',
                  'CCU', 'CAU', 'CGU', 'CCC', 'CAC', 'CGC', 'CCA',
                  'CAA', 'CGA', 'CCG', 'CAG', 'CGG', 'ACU', 'AAU',
                  'AGU', 'ACC', 'AAC', 'AGC', 'ACA', 'AAA', 'AGA',
                  'ACG', 'AAG', 'AGG', 'GCU', 'GAU', 'GGU', 'GCC',
                  'GAC', 'GGC', 'GCA', 'GAA', 'GGA', 'GCG', 'GAG',
                  'GGG']
    for codon_ in codon_list:
        prior[codon_] = np.round(codon.count(codon_)/total_codon, 4)
    return json.dumps(prior, indent = 2)
    
# sequence = read_fasta(filename) 
# total_nucleotides = prior_nucleotide(sequence)
# #print((prior_codon(sequence)))
# print(total_nucleotides)       
    



