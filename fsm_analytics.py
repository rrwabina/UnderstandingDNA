import numpy as np
import pandas as pandas
import get_data
import hmm_viterbi as hmm
import reverse_translator
import finite_state_machine as fsm

filename = '/root/NLU/Project/data/test_cases.fa'
def read_dna(filename, pos = 0):
    '''
    First state of the finite state machine
    '''
    dna_sequence = get_data.read_fasta(filename, pos)
    if dna_sequence == dna_sequence:
        newState = 'Extract Bases' 
    return (newState, dna_sequence)


def extract_base(sequence):
    'Second state of FSM'
    bases = []
    for base in sequence:
        bases.append(base)
    newState = 'Convert Bases'
    return (newState, bases)

def convert_base(bases):
    'Third state of FSM'
    'Transcribes DNA to mRNA'
    cbases = []
    for base in bases:
        if base == 'A':
            cbase = 'U'
        elif base == 'C':
            cbase = 'G'
        elif base == 'G':
            cbase = 'C'
        else:
            cbase = 'A'
        cbases.append(cbase)
    if cbases == cbases:
        newState = 'mRNA'
    return (newState, cbases)

def rna_sequence(cbases):
    'Fourth state in FSM'
    if cbases[0] == 'AUG':
        sequence = cbases
        newState = 'Stop Codon'

    elif cbases[0] == 'A' or cbases[0] == 'U' or cbases[0] == 'C' or cbases[0] == 'G':  
        sequence = ''.join([str(cbase) for cbase in cbases])
        newState = 'Start Codon'
    else:
        newState = 'Final State'
        sequence = cbases
    return (newState, sequence)

def start_codon(sequence):
    ngram = 3
    codons = [sequence[i: i + ngram] for i in range(0, len(sequence), ngram)]
    start = 'AUG'
    start_seq = []
    if codons[0] == start:
        for codon in codons:
            if codon == 'AUG':
                start_seq.append(codon)
        newState = 'Final State'
    else:
        for codon in codons:
            if codon == start:
                start_seq.append(codon)
                codons.remove(codon)
                codons = ['AUG'] + codons
        newState = 'mRNA'
    return (newState, codons)

def stop_codon_state(sequence):
    stop_codons = ['UGA', 'UAA', 'UAG']
    stop_codons_list = []
    stop_index = []
    for seq in sequence:
        if seq == 'UGA' or seq == 'UAA' or seq == 'UAG':
            stop_codons_list.append(seq)
            stop_index.append(sequence.index(seq))
    for i in stop_index:
        sequence.pop(i)        
    newState = 'Convert Codon'
    return (newState, sequence)

def convert_codon(codons):
    '''
    State 8 translates the extracted codon to its corresponding amino acid 
    based on the genetic code. This amino acid is transferred to the final state
    '''
    codon_list = {"UUU" : "F", "CUU" : "L", "AUU" : "I", "GUU" : "V",
           "UUC" : "F", "CUC" : "L", "AUC" : "I", "GUC" : "V",
           "UUA" : "L", "CUA" : "L", "AUA" : "I", "GUA" : "V",
           "UUG" : "L", "CUG" : "L", "AUG" : "M", "GUG" : "V",
           "UCU" : "S", "CCU" : "P", "ACU" : "T", "GCU" : "A",
           "UCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
           "UCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
           "UCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
           "UAU" : "Y", "CAU" : "H", "AAU" : "N", "GAU" : "D",
           "UAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
           "UAA" : "STOP", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
           "UAG" : "STOP", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
           "UGU" : "C", "CGU" : "R", "AGU" : "S", "GGU" : "G",
           "UGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
           "UGA" : "STOP", "CGA" : "R", "AGA" : "R", "GGA" : "G",
           "UGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G" 
           }
    protein_string = ''
    for codon in codons:
        if codon_list[codon] == 'STOP':
            # print('STOP')
            break
        protein_string += codon_list[codon]
    print('Proteins: ', protein_string)
    newState = 'Final State'
    return (newState, protein_string)
        

m = fsm.FSM()
m.add_state('Read DNA', read_dna, end_state = 0)
m.add_state('Extract Bases', extract_base)
m.add_state('Convert Bases', convert_base)
m.add_state('mRNA', rna_sequence)
m.add_state('Start Codon', start_codon)
m.add_state('Stop Codon', stop_codon_state)
m.add_state('Convert Codon', convert_codon)
m.add_state('Final State', None, end_state = 1)
m.set_start('Read DNA')
m.run(filename)