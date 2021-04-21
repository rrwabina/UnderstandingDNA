import numpy as np
import itertools
from Bio.Seq import Seq
from Bio.Data import CodonTable
from dnachisel import biotools

def generate_dna(length = 1000, gc_share = None, probas = None, seed = None):
    '''
    Generate random DNA sequence with 1000 nucleotide bases
    '''
    if seed is not None:
        np.random.seed(seed)
    if gc_share is not None:
        g_or_c = gc_share / 2
        not_g_or_c = (1 - gc_share) / 2.0
        probas = {'G': g_or_c, 'C': g_or_c, 'A': not_g_or_c, 'T': not_g_or_c}
    if probas is None:
        sequence = np.random.choice(list('ATCG'), length)
    else:
        bases, probas = zip(*probas.items())
        sequence = np.random.choice(bases, length, p = probas)
    return ''.join(sequence)

def generate_proteins(length = 100, seed = None):
    '''
    Generate random amino acid sequence with default length of 100
    Uses 1-letter amino acid code from IUPAC
    '''
    if seed is not None:
        np.random.seed(seed)
    amino_acids = ['ACEDGFIHKLNQPSRTWVY']
    amino_acids_choice = np.random.choice(amino_acids, length - 2)
    return ''.join(amino_acids_choice)

def dna_complement(dna_sequence):
    '''
    Input: DNA sequence in normal format
    Output: Complementary sequence of DNA during Transcription process
    '''
    if len(dna_sequence) <= 100:
        return "".join([COMPLEMENTS[nuc] for nuc in dna_sequence])
    return str(Seq(dna_sequence).complement())

def reverse_dna_complement(sequence):
    '''
    Return the reverse-complement of DNA sequence, reverse-transcription process
    '''
    return complement(sequence[::-1])

def dna_translate(dna_sequence, start_codon = False, delete_stop_codon = True):
    '''
    Translate a DNA sequence into amino acids
    '''
    dna_sequence = dna_sequence[:len(dna_sequence)-(len(dna_sequence) % 3)]  
    aa = biotools.translate(dna_sequence, assume_start_codon = start_codon)
    if delete_stop_codon is True:
        aa = aa.replace('*', '')     
    return aa

def reverse_proteins(protein_sequence, randomize_codons = False):
    '''
    Input: Protein sequence in normal format
    Output: DNA sequence using stochastic codon-optimization
    '''
    return biotools.reverse_translate(protein_sequence, randomize_codons)


